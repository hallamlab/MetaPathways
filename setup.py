import os
from setuptools import setup, find_packages


package_root = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(package_root, "metapathways", "_version.py")) as fp:
    _, v = fp.read().strip().split(" = ")
version = v.strip('"')

CLASSIFIERS = [
    "Development Status :: 2 - Pre-Alpha",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: MIT License",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.5",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="MetaPathways",
    version=version,
    author="Kishori Mohan Konwar",
    author_email="kishori82@gmail.com",
    description=(
        "MetaPathways is a modular pipeline to build PGDBs"
        " from Metagenomic sequences."
    ),
    license="MIT",
    keywords="metagenomics pipeline",
    url="http://packages.python.org/",
    download_url="https://github.com/kishori82/MetaPathways_Python.3.0/archive/kmk-develop.zip",
    packages=find_packages(),
    scripts=["bin/metapathways-install-deps.sh",
             "bin/metapathways-data-install.sh",
             "bin/metacount",
             "bin/fastal",
             "bin/fastdb",
             "dev/pgdb_build_wf.py",
             "dev/run-pathway-tools-and-copy-pgdb-singularity.sh"],
    install_requires=["pyfastx",
                      "docopt",
                      'camelot_frs @ git+https://bitbucket.org/tomeraltman/camelot-frs@dev#egg=camelot-frs',
                      'pandas',
                      'tqdm',
                      'pprodigal',
                      'sexpdata'
                      ],
    entry_points={"console_scripts": ["MetaPathways=metapathways.pipeline:main"]},
    long_description=read("README.md"),
    #package_data={'resources': ['Dsignal', 'TPCsignal', 'template_param.txt']},
    include_package_data=True,
    data_files=[('data_file_test', ['README.md', 'Makefile'])],
    classifiers=CLASSIFIERS,
    extras_require={
        "test": ["pytest", "pytest-cov", "tox"],
    },
    python_requires=">3.5.2",
)
