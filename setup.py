import os
from setuptools import setup

CLASSIFIERS = [
    "Development Status :: 3 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3.0+",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "MetaPathways",
    version = "3.5.0",
    author = "Kishori Mohan Konwar",
    author_email = "kishori82@gmail.com",
    description = ("MetaPathways is a modular pipeline to build PGDBs"\
                   " from Metagenomic sequences."),
    license = "MIT",
    keywords = "metagenomics pipeline",
    url = "http://packages.python.org/",
    package_dir = {"": "src"},
    packages = ["metapathways", 
                "metapathways/pipeline",
                "metapathways/test", 
                "metapathways/utils",
                "metapathways/taxonomy",
                "metapathways/parsers",
                "metapathways/scripts",
                "metapathways/annotate",
                "metapathways/diagnostics"
    ],
    py_modules = ['metapathways/modules'],
    scripts = ['bin/compress_by_ec'],
    entry_points={
        "console_scripts": [
            "MetaPathways = metapathways." + "MetaPathways:main"
         ]
    },
    long_description=read('README.md'),
    include_package_data=True,
    classifiers = CLASSIFIERS
)
