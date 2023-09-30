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
    _, v = fp.read().strip().split(" = ")
VERSION = v.strip('"')

CLASSIFIERS = [
    "Development Status :: 2 - Pre-Alpha",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: MIT License",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if __name__ == "__main__":
    setup(
        name=NAME,
        version=VERSION,
        author="Kishori Mohan Konwar",
        author_email="kishori82@gmail.com",
        description=(
            "MetaPathways is a modular pipeline to build PGDBs"
            " from Metagenomic sequences."
        ),
        license="MIT",
        keywords="metagenomics pipeline",
        url="https://bitbucket.org/BCB2/metapathways/",
        packages=find_packages(),
        scripts=["bin/metapathways-install-deps.sh",
                "bin/metapathways-data-install.sh",
                "bin/metacount",
                "bin/fastal",
                "bin/fastdb",
                "dev/pgdb_build_wf.py",
                "dev/run-pathway-tools-and-copy-pgdb-singularity.sh",
                "dev/run-pathway-tools-and-copy-pgdb-singularity_taxprune.sh",
                "dev/abund_calc.py",
                "dev/gff2gtf.py"],
        entry_points={"console_scripts": ENTRY_POINTS},
        long_description=read("README.md"),
        #package_data={'resources': ['Dsignal', 'TPCsignal', 'template_param.txt']},
        include_package_data=True,
        data_files=[('data_file_test', ['README.md', 'Makefile'])],
        classifiers=CLASSIFIERS,
        extras_require={
            "test": ["pytest", "pytest-cov", "tox"],
        },
        python_requires=">=3.10",
        install_requires=[
            "camelot_frs @ git+https://bitbucket.org/tomeraltman/camelot-frs@dev#egg=camelot-frs"
        ],
    )
