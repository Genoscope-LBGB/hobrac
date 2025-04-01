import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="HoBRAC",
    version="0.1.0",
    author=["Jean-Marc Aury", "Benjamin Istace"],
    author_email=["jmaury@genoscope.cns.fr", "bistace@genoscope.cns.fr"],
    description="HoBRAC is a rdbioseq tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data = { "hobrac": ["workflow/*"] },
    install_requires=[
        "snakemake", 
        "snakemake-executor-plugin-slurm",
        "find_reference_genomes",
    ],
    dependency_links=['http://github.com/user/repo/tarball/master#egg=package-1.0'],
    entry_points={
        "console_scripts": [
            "hobrac=hobrac.main:main",
        ],  
    },
    
)
