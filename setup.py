import os
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.dist import Distribution


class PlatformDistribution(Distribution):
    """The bundled FAST/metacount executables are native Linux binaries."""

    def has_ext_modules(self):
        return True

PACKAGE_ROOT = Path(os.path.realpath(__file__)).parent
NAME = "metapathways".lower()
ENTRY_POINTS =  [
    'MetaPathways=metapathways.pipeline:main',
    'metapathways=metapathways.pipeline:main',
]
with open(os.path.join(PACKAGE_ROOT, "metapathways", "_version.py")) as fp:
    dat = fp.readlines()
    ver_dict = {x.strip().split(' = ')[0].strip('__'):x.strip().split(' = ')[1]
                for x in dat
                }
VERSION = ver_dict['version'].strip('"')
CLASSIFIERS = [
    "Development Status :: 2 - Pre-Alpha",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    f"License :: {ver_dict['license']}",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if __name__ == "__main__":
    setup(
        distclass=PlatformDistribution,
        name=NAME,
        version=VERSION,
        author=ver_dict['author'],
        author_email=ver_dict['contact'],
        description=(
            "MetaPathways is a modular pipeline to build PGDBs"
            " from Metagenomic sequences."
        ),
        license=ver_dict['license'],
        keywords="metagenomics pipeline",
        url="https://github.com/hallamlab/MetaPathways/",
        packages=find_packages(include=["metapathways", "metapathways.*"]),
        package_data={'': ['metapathways/regtests/**/*'], 'metapathways': ['report_assets/*', 'resources/*.json']},
        scripts=["bin/metapathways-install-deps.sh",
                "bin/metapathways-data-install.sh",
                "dev/metacount",
                "dev/fastal",
                "dev/fastdb",
                "dev/pgdb_build_wf.py",
                "dev/pgdb_build_single.py",
                "dev/run-pathway-tools-and-copy-pgdb.sh",
                "dev/run-pathway-tools-and-copy-pgdb_taxprune.sh",
                "dev/run-pathway-tools-and-copy-pgdb_local.sh",
                "dev/run-pathway-tools-and-copy-pgdb_taxprune_local.sh",
                "dev/abund_calc.py",
                "dev/gff2gtf.py",
                "dev/ptRNAscan.py",
                "dev/metacyc_add_ont2pwy.py",
                "dev/pwy_abundance.py",
                "dev/install-pathway-tools.sh",
                "dev/run-pathway-tools.sh",
                "dev/ptools-init.dat"],
        entry_points={"console_scripts": ENTRY_POINTS},
        long_description=read("README.md"),
        include_package_data=True,
        classifiers=CLASSIFIERS,
        extras_require={
            "test": ["pytest", "pytest-cov", "tox"],
        },
        python_requires=">=3.10",
        install_requires=[line.strip() for line in read("requirements-workflow.txt").splitlines()
                          if line.strip() and not line.startswith("#")],
    )
