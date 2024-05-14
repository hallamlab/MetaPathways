import os
from pathlib import Path
from setuptools import setup, find_packages

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
        url="https://bitbucket.org/BCB2/metapathways/",
        packages=find_packages(),
        scripts=["bin/metapathways-install-deps.sh",
                "bin/metapathways-data-install.sh",
                "bin/metacount",
                "bin/fastal",
                "bin/fastdb",
                "bin/pgdb_build_wf.py",
                "bin/pgdb_build_single.py",
                "bin/run-pathway-tools-and-copy-pgdb.sh",
                "bin/run-pathway-tools-and-copy-pgdb_taxprune.sh",
                "bin/run-pathway-tools-and-copy-pgdb_local.sh",
                "bin/run-pathway-tools-and-copy-pgdb_taxprune_local.sh",
                "dev/abund_calc.py",
                "dev/gff2gtf.py",
                "bin/ptRNAscan.py",
                "bin/metacyc_add_ont2pwy.py",
                "bin/pwy_abundance.py",
                "bin/install-pathway-tools.sh",
                "bin/run-pathway-tools.sh",
                "bin/ptools-init.dat"],
        entry_points={"console_scripts": ENTRY_POINTS},
        long_description=read("README.md"),
        include_package_data=True,
        data_files=[('data_file_test', ['README.md', 'Makefile'])],
        classifiers=CLASSIFIERS,
        extras_require={
            "test": ["pytest", "pytest-cov", "tox"],
        },
        python_requires=">=3.10",
        install_requires=[],
    )
