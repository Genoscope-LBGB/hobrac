"""Offline tests for version discovery, CLI wiring and Snakemake reruns."""

import io
import os
import shlex
import subprocess
import sys
from email.message import Message
from pathlib import Path
from types import SimpleNamespace
from urllib.error import HTTPError, URLError

import pytest

import hobrac.main as hobrac_main
import hobrac.mash_database as database
from hobrac.command_line import get_args


def analysis_args(monkeypatch, *extra):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "hobrac",
            "-a",
            "assembly.fa",
            "-n",
            "Test species",
            "-t",
            "1",
            "-e",
            "local",
            *extra,
        ],
    )
    return get_args()


@pytest.mark.parametrize("version", ["latest", "20260629", "20270101"])
def test_version_url(version):
    assert database.mash_db_url(version) == database.MASH_DB_URL + version + "/"


@pytest.mark.parametrize(
    "version",
    [
        "",
        "../20260629",
        "2026-06-29",
        "20260629/",
        "latest;echo bad",
        "https://example.org",
        "２０２６０６２９",
    ],
)
def test_invalid_versions_rejected(version, monkeypatch):
    with pytest.raises(ValueError):
        database.mash_db_url(version)
    with pytest.raises(SystemExit) as exc:
        analysis_args(monkeypatch, "--mash-db-version", version)
    assert exc.value.code == 2


def test_resolve_latest_uses_version_file(monkeypatch):
    calls = []

    def fetch(url, timeout):
        calls.append((url, timeout))
        return io.BytesIO(b"20260629\n")

    monkeypatch.setattr(database, "urlopen", fetch)
    assert database.resolve_mash_db_version("latest") == "20260629"
    assert calls == [(database.MASH_DB_URL + "latest/VERSION", 30)]


def test_explicit_version_needs_no_network(monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail("Pinned releases must not look up latest")

    monkeypatch.setattr(database, "urlopen", unexpected)
    assert database.resolve_mash_db_version("20260629") == "20260629"


@pytest.mark.parametrize(
    "payload",
    [
        b"",
        b"latest\n",
        b"2026-06-29",
        b"../20260629",
        b"<html>Error</html>",
        b"20260629\n20270101",
        b"20260629" + b" " * 100,
        b"\xff",
    ],
)
def test_bad_version_file_fails_without_fallback(monkeypatch, payload):
    monkeypatch.setattr(database, "urlopen", lambda *a, **kw: io.BytesIO(payload))
    with pytest.raises(RuntimeError, match="MASH database version"):
        database.resolve_mash_db_version("latest")


def test_listing_filters_and_sorts_directory_links(monkeypatch):
    html = b"""
    <a href="?C=N;O=D">Name</a>
    <a href="../">Parent Directory</a>
    <a href="20260629/">First release</a>
    <a href="latest/">latest</a>
    <a href="20270101/">New release</a>
    <a href="./20260629/">Duplicate</a>
    <a href="Chordata.msh">Old unversioned sketch</a>
    <a href="20270101">Not a directory</a>
    <a href="20270101.txt">A file</a>
    <a href="nested/20280101/">Not a release directory</a>
    <a href="https://example.org/20290101/">External directory</a>
    <a>No href</a>
    """
    calls = []

    def fetch(url, timeout):
        calls.append((url, timeout))
        return io.BytesIO(html)

    monkeypatch.setattr(database, "urlopen", fetch)
    assert database.list_mash_db_versions() == ["latest", "20270101", "20260629"]
    assert calls == [(database.MASH_DB_URL, 30)]


def test_listing_does_not_invent_latest(monkeypatch):
    monkeypatch.setattr(
        database, "urlopen", lambda *a, **kw: io.BytesIO(b'<a href="20260629/">x</a>')
    )
    assert database.list_mash_db_versions() == ["20260629"]


@pytest.mark.parametrize("html", [b"", b'<a href="Chordata.msh">Chordata</a>'])
def test_empty_or_legacy_index_is_an_error(monkeypatch, html):
    monkeypatch.setattr(database, "urlopen", lambda *a, **kw: io.BytesIO(html))
    with pytest.raises(RuntimeError, match="No MASH database versions found"):
        database.list_mash_db_versions()


@pytest.mark.parametrize(
    "error",
    [
        URLError("offline"),
        TimeoutError("timed out"),
        HTTPError(database.MASH_DB_URL, 404, "Not Found", Message(), None),
    ],
)
def test_listing_network_errors_are_actionable(monkeypatch, error):
    def fetch(*args, **kwargs):
        raise error

    monkeypatch.setattr(database, "urlopen", fetch)
    with pytest.raises(RuntimeError, match="Could not list MASH database versions"):
        database.list_mash_db_versions()
    with pytest.raises(RuntimeError, match="Could not resolve latest"):
        database.resolve_mash_db_version("latest")


def test_standalone_listing_needs_no_inputs_or_tools(monkeypatch, tmp_path, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["hobrac", "--list-mash-db-versions"])
    monkeypatch.setattr(
        hobrac_main, "list_mash_db_versions", lambda: ["latest", "20260629"]
    )

    def unexpected(*args, **kwargs):
        pytest.fail("Listing must not prepare or launch an analysis")

    monkeypatch.setattr(hobrac_main, "resolve_mash_db_version", unexpected)
    monkeypatch.setattr(hobrac_main, "check_dependencies", unexpected)
    monkeypatch.setattr(hobrac_main, "create_dir", unexpected)
    monkeypatch.setattr(hobrac_main, "generate_snakemake_command", unexpected)
    hobrac_main.main()
    assert capsys.readouterr().out == "latest\n20260629\n"
    assert list(tmp_path.iterdir()) == []


def test_listing_failure_exits_cleanly(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["hobrac", "--list-mash-db-versions"])

    def fail():
        raise RuntimeError("Server unavailable")

    monkeypatch.setattr(hobrac_main, "list_mash_db_versions", fail)
    with pytest.raises(SystemExit) as exc:
        hobrac_main.main()
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert not captured.out
    assert captured.err == "Error: Server unavailable\n"


@pytest.mark.parametrize("argv", [[], ["-a", "assembly.fa"], ["-t", "1"]])
def test_analysis_still_requires_inputs(monkeypatch, argv, capsys):
    monkeypatch.setattr(sys, "argv", ["hobrac", *argv])
    with pytest.raises(SystemExit) as exc:
        get_args()
    assert exc.value.code == 2
    assert "arguments are required" in capsys.readouterr().err


@pytest.mark.parametrize("version", [None, "20260629"])
def test_cli_passes_database_version_to_snakemake(monkeypatch, version):
    extra = ["--mash-db-version", version] if version else []
    args = analysis_args(monkeypatch, *extra)
    args.busco_reference_override_paths = {}
    assert args.mash_db_version == (version or "latest")
    command = shlex.split(hobrac_main.generate_snakemake_command(args))
    assert f"mash_db_version={version or 'latest'}" in command


@pytest.mark.parametrize("manual,stop", [(False, False), (True, False), (True, True)])
def test_main_resolves_only_when_mash_is_needed(tmp_path, monkeypatch, manual, stop):
    monkeypatch.chdir(tmp_path)
    Path("assembly.fa").write_text(">query\nACGT\n")
    Path("reference.fa").write_text(">reference\nACGT\n")
    extra = ["-r", "reference.fa"] if manual else []
    if stop:
        extra.append("--stop-after-mash")
    args = analysis_args(monkeypatch, *extra)
    monkeypatch.setattr(hobrac_main, "get_args", lambda: args)
    monkeypatch.setattr(hobrac_main, "check_dependencies", lambda **kw: None)
    lookups, commands = [], []

    def resolve(version):
        lookups.append(version)
        return "20260629"

    def launch(command, **kwargs):
        commands.append(shlex.split(command))
        return SimpleNamespace(returncode=0, wait=lambda: 0)

    monkeypatch.setattr(hobrac_main, "resolve_mash_db_version", resolve)
    monkeypatch.setattr(hobrac_main, "subprocess", SimpleNamespace(Popen=launch))
    with pytest.raises(SystemExit) as exc:
        hobrac_main.main()
    assert exc.value.code == 0
    expected = "20260629" if not manual or stop else "latest"
    assert lookups == (["latest"] if not manual or stop else [])
    assert f"mash_db_version={expected}" in commands[0]


def test_resolution_failure_stops_before_creating_output(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    analysis_args(monkeypatch)

    def fail(*args, **kwargs):
        raise URLError("offline")

    monkeypatch.setattr(database, "urlopen", fail)
    with pytest.raises(SystemExit) as exc:
        hobrac_main.main()
    assert exc.value.code == 1
    assert "Could not resolve latest" in capsys.readouterr().err
    assert not Path("hobrac_analysis").exists()


def test_workflow_downloads_version_and_reruns_on_change(tmp_path, monkeypatch):
    """Execute the real workflow with fake tools; never contact the database server."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "assembly.fa").write_text(">query\nACGT\n")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    tools = {
        "taxonkit": "#!/bin/sh\ncat >/dev/null\nprintf '1\\tChordata\\n'\n",
        "wget": (
            "#!/bin/sh\nset -eu\n"
            "printf '%s\\n' \"$1\" >> downloads.txt\n"
            'printf \'%s\\n\' "$1" > "$3"\n'
        ),
        "mash": (
            '#!/bin/sh\ncase "$1" in\n'
            'info) cat "$2" ;;\n'
            'dist) cat "$4" ;;\n'
            "*) exit 1 ;;\nesac\n"
        ),
    }
    for name, content in tools.items():
        path = bindir / name
        path.write_text(content)
        path.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ["PATH"])
    args = analysis_args(monkeypatch)
    args.busco_reference_override_paths = {}
    command = shlex.split(hobrac_main.generate_snakemake_command(args))
    config = command[command.index("--config") + 1 :]
    config = [item for item in config if not item.startswith("mash_db_version=")]

    # Mock HTTP inside the Snakemake subprocess, including changing the alias.
    wrapper = tmp_path / "Snakefile"
    wrapper.write_text(
        "import io\n"
        "import hobrac.mash_database as database\n"
        "def fetch_version(url, timeout):\n"
        "    assert url == database.MASH_DB_URL + 'latest/VERSION'\n"
        "    with open('resolutions.txt', 'a') as log:\n"
        "        log.write(url + '\\n')\n"
        "    with open('latest_version.txt', 'rb') as version:\n"
        "        return io.BytesIO(version.read())\n"
        "database.urlopen = fetch_version\n"
        f"include: {hobrac_main.snakefile_path!r}\n"
    )
    Path("latest_version.txt").write_text("20260629\n")

    def run(version=None):
        # Omitting the config key also tests direct-Snakemake backward compatibility.
        selected = [] if version is None else [f"mash_db_version={version}"]
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "snakemake",
                "--snakefile",
                str(wrapper),
                "--cores",
                "1",
                "--scheduler",
                "greedy",
                "--printshellcmds",
                "mash/mash.dist",
                "--config",
                *config,
                *selected,
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    run()
    version_file = Path("mash/mash_db.msh.version")
    assert version_file.read_text() == "20260629\n"
    assert Path("mash/mash_db.msh.phylum").read_text() == "Chordata\n"
    first_url = database.MASH_DB_URL + "20260629/Chordata.msh\n"
    assert Path("mash/mash.dist").read_text() == first_url
    lookup = database.MASH_DB_URL + "latest/VERSION\n"
    assert Path("resolutions.txt").read_text() == lookup

    # Moving latest must trigger a download from the new immutable release.
    Path("latest_version.txt").write_text("20270101\n")
    run()
    next_url = database.MASH_DB_URL + "20270101/Chordata.msh\n"
    assert version_file.read_text() == "20270101\n"
    assert Path("mash/mash.dist").read_text() == next_url
    assert Path("downloads.txt").read_text() == first_url + next_url
    assert Path("resolutions.txt").read_text() == lookup * 2

    # Selecting the same resolved date reuses the sketch without a lookup.
    run("20270101")
    assert Path("downloads.txt").read_text() == first_url + next_url
    assert Path("resolutions.txt").read_text() == lookup * 2

    # Explicitly pinning an older release still triggers a download.
    run("20260629")
    assert version_file.read_text() == "20260629\n"
    assert Path("mash/mash.dist").read_text() == first_url
    assert Path("downloads.txt").read_text() == first_url + next_url + first_url
    assert Path("resolutions.txt").read_text() == lookup * 2
