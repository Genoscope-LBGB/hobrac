"""Version selection and discovery for the public MASH database."""

import re
from html.parser import HTMLParser
from urllib.error import URLError
from urllib.request import urlopen

MASH_DB_URL = "https://www.genoscope.cns.fr/lbgb/mash_versioned/"


def validate_mash_db_version(value: str) -> str:
    """Accept the moving alias or a dated release, never an arbitrary URL/path."""
    if value != "latest" and not re.fullmatch(r"[0-9]{8}", value):
        raise ValueError("MASH database version must be 'latest' or YYYYMMDD")
    return value


def mash_db_url(version: str) -> str:
    return f"{MASH_DB_URL}{validate_mash_db_version(version)}/"


def resolve_mash_db_version(version: str) -> str:
    """Resolve latest to an immutable release; explicit dates need no network access."""
    validate_mash_db_version(version)
    if version != "latest":
        return version

    url = mash_db_url("latest") + "VERSION"
    try:
        with urlopen(url, timeout=30) as response:
            payload = response.read(65)
        resolved = payload.decode("utf-8").strip()
    except (OSError, URLError, UnicodeError) as exc:
        raise RuntimeError(
            f"Could not resolve latest MASH database version from {url}: {exc}"
        ) from exc
    if len(payload) > 64 or not re.fullmatch(r"[0-9]{8}", resolved):
        raise RuntimeError(
            f"Invalid MASH database version in {url}; expected a YYYYMMDD release."
        )
    return resolved


class _VersionIndexParser(HTMLParser):
    """Read release-directory links from the server's Apache directory index."""

    def __init__(self):
        super().__init__()
        self.versions = set()

    def handle_starttag(self, tag, attrs):
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                match = re.fullmatch(r"(?:\./)?(latest|[0-9]{8})/", value)
                if match:
                    self.versions.add(match.group(1))


def list_mash_db_versions() -> list[str]:
    """List published releases (newest first), with latest first if advertised."""
    try:
        with urlopen(MASH_DB_URL, timeout=30) as response:
            index = response.read().decode("utf-8")
    except (OSError, URLError, UnicodeError) as exc:
        raise RuntimeError(
            f"Could not list MASH database versions at {MASH_DB_URL}: {exc}"
        ) from exc

    parser = _VersionIndexParser()
    parser.feed(index)
    if not parser.versions:
        raise RuntimeError(
            f"No MASH database versions found at {MASH_DB_URL}; "
            "expected YYYYMMDD/ release directories and a latest/ alias."
        )
    return sorted(parser.versions, reverse=True)
