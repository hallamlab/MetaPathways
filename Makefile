#### Assumptions:
## * AWS CLI installed & configured
## * sudo apt-get install make python2.7

#### Make configuration:

## Use Bash as default shell, and in strict mode:
SHELL := /bin/bash
.SHELLFLAGS = -ec

## If the parent env doesn't ste TMPDIR, do it ourselves:
TMPDIR ?= /tmp

## Users can override this variable from the command line,
## to install MP binaries somewhere other than /usr/local,
## if they lack root privileges:
DESTDIR ?= /usr/local

## This makes all recipe lines execute within a shared shell process:
## https://www.gnu.org/software/make/manual/html_node/One-Shell.html#One-Shell
.ONESHELL:

## If a recipe contains an error, delete the target:
## https://www.gnu.org/software/make/manual/html_node/Special-Targets.html#Special-Targets
.DELETE_ON_ERROR:

## This is necessary to make sure that these intermediate files aren't clobbered:
.SECONDARY:


### Local Definitions:
PYTHON ?= python3

### Install Singularity

OS := linux
ARCH := amd64
GO-VERSION := 1.17.6
SY-VERSION := 3.9.3
singularity-install:
	sudo apt-get update
	sudo apt-get install -y \
	   build-essential \
	   libseccomp-dev \
	   pkg-config \
	   squashfs-tools \
	   cryptsetup
	cd $(TMPDIR)
	wget https://dl.google.com/go/go$(GO-VERSION).$(OS)-$(ARCH).tar.gz # Downloads the required Go package
	sudo tar -C /usr/local -xzvf go$(GO-VERSION).$(OS)-$(ARCH).tar.gz  # Extracts the archive
	rm go$(GO-VERSION).$(OS)-$(ARCH).tar.gz                            # Deletes the ``tar`` file
	export PATH=$$PATH:/usr/local/go/bin
	wget https://github.com/sylabs/singularity/releases/download/v${SY-VERSION}/singularity-ce-${SY-VERSION}.tar.gz
	tar -xzf singularity-ce-${SY-VERSION}.tar.gz
	rm singularity-ce-${SY-VERSION}.tar.gz
	cd singularity-ce-${SY-VERSION}
	./mconfig
	make -C builddir
	sudo make -C builddir install


### Container Automation
docker-start:
	sudo systemctl start docker

## If git_branch is empty, it probably means that we are building off of a tagged version.
## So, we then grab the tag string:
docker-build: #pre-docker-builds
	git_branch=$$(git symbolic-ref --short -q HEAD) \
		|| git_branch=$$(git describe --tags)
	cd docker
	sudo docker build --network=host \
			--build-arg git_branch=$$git_branch \
			-t quay.io/hallamlab/metapathways:$$git_branch .

docker-run:
	git_branch=$$(git symbolic-ref --short -q HEAD)
	sudo docker run -it --network=host --rm \
		-v $(CURDIR):/input \
		-v $(CURDIR)/out:/output \
		quay.io/hallamlab/metapathways:$$git_branch /bin/bash

docker-test:
	cp $(CURDIR)/regtests/input/A1.fasta /tmp
	git_branch=$$(git symbolic-ref --short -q HEAD)
	sudo docker run -it -v /tmp:/input \
		quay.io/hallamlab/metapathways:$$git_branch \
		MetaPathways -v \
			-i /opt/mp_repo/tests/data/input/lagoon-sample2.fasta \
			-o /input/A1_MP_out/ \
			-p /opt/mp_repo/resources/template_param.txt \
			-d /opt/mp_repo/tests/data/ref_data


docker-deploy:
	git_branch=$$(git symbolic-ref --short -q HEAD)
	sudo docker push quay.io/hallamlab/metapathways:$$git_branch

docker-fetch:
	sudo docker pull quay.io/hallamlab/metapathways

singularity-local-build:
	git_branch=$$(git symbolic-ref --short -q HEAD) \
		|| git_branch=$$(git describe --tags)
	sudo singularity build metapathways-$$git_branch.sif docker-daemon://quay.io/hallamlab/metapathways:$$git_branch

singularity-local-shell:
	singularity shell metapathways-dev.sif

singularity-docker-build:
	sudo singularity build metapathways-dev.sif docker://quay.io/hallamlab/metapathways:dev

singularity-docker-shell:
	singularity shell docker://quay.io/hallamlab/metapathways:dev


### Conda Packaging
##

## Conda Installation:
conda-install: conda-install-deps #extensions-install 

## Install conda build tools:
conda-build-init:
	conda init bash
	conda install --yes conda-build

conda-install-deps:
	bin/metapathways-install-deps.sh


### Python PyPI Packaging:
##
##

### Python installs:

create-package: clean-package
	$(PYTHON) -m pip install --user --upgrade setuptools wheel twine
	$(PYTHON) setup.py sdist bdist_wheel --universal
	cp -r dist docker/dist

clean-package:
	rm -rf dist MetaPathways.egg-info build docker/dist

install-package:
	$(PYTHON) -m pip install --user .

install-dev-package:
	$(PYTHON) -m pip install --user --upgrade -e .

install-dist-package:
	$(PYTHON) -m pip install --user dist/MetaPathways-*-py*-none-any.whl

deploy-package-to-pypi:
	twine upload dist/*

### Build & Install Extensions
##
##

extensions-build:
	$(MAKE) -C extensions clean
	$(MAKE) -C extensions

## Deprecated, see Python package.
extensions-build-install: extensions-build
	mkdir -p $(DESTDIR)/bin
	cp extensions/FAST/fast*            $(DESTDIR)/bin
	cp extensions/metacount/metacount   $(DESTDIR)/bin

## Deprecated, these binaries need to get revisited.
extensions-install:
	mkdir -p $(DESTDIR)/bin
	$(MAKE) -C extensions/metacount clean
	$(MAKE) -C extensions/metacount 
	cp extensions/FAST/fastal $(DESTDIR)/bin/fastal
	cp extensions/FAST/fastdb $(DESTDIR)/bin/fastdb
	cp extensions/metacount/metacount $(DESTDIR)/bin/metacount
	chmod 755 $(DESTDIR)/bin/fastal
	chmod 755 $(DESTDIR)/bin/fastdb
	chmod 755 $(DESTDIR)/bin/metacount

# The location of the expat directory
CC=gcc  
LEX=lex  
LEXFLAGS=-lfl
CFLAGS=-C


#example: 
#     export  METAPATHWAYS_DB=../fogdogdatabases
#     export  PTOOLS_DIR=../ptools/
#
#     a) make install-without-ptools  METAPATHWAYS_DB=../fogdogdatabases  
#     this will get the files uploaded by wholebiome into the path in METAPATHWAYS_DB but NOT the ptools
#
#     b) make mp-regression-tests:
#
#     c) make install-with-ptools
#     this will get the files uploaded by koonkie into the path in METAPATHWAYS_DB and the ptools.tar.gz into the PTOOLS_DIR
#

OS_PLATFORM=linux
#should be the same as the EXECUTABLES_DIR in the template_config.txt file

NCBI_BLAST=ncbi-blast-2.10.1+-x64-linux.tar.gz
NCBI_BLAST_VER=ncbi-blast-2.10.1+
BINARY_FOLDER=executables/$(OS_PLATFORM)


BLASTP=$(BINARY_FOLDER)/blastp
LASTAL=$(BINARY_FOLDER)/lastal+
RPKM=$(BINARY_FOLDER)/rpkm
BWA=$(BINARY_FOLDER)/bwa
FAST=$(BINARY_FOLDER)/fastal
PRODIGAL=$(BINARY_FOLDER)/prodigal

GIT_SUBMODULE_UPDATE=gitupdate
# Alias for target 'all', for compliance with FogDog deliverables standard:

#all: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA)  $(RPKM)
all: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA) $(RPKM) $(BLASTP) METAPATHWAYS_DB_FETCH
#pre-docker-builds: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA) $(RPKM) $(BLASTP) 


install-without-ptools: all METAPATHWAYS_DB_FETCH

install-with-ptools: all METAPATHWAYS_DB_FETCH 

.PHONY: METAPATHWAYS_DB_FETCH
METAPATHWAYS_DB_FETCH:
	@if [ -z $(METAPATHWAYS_DB) ]; then  echo "Variable METAPATHWAYS_DB not set. Set it as export METPATHWAYS_DB=<path>" ;  exit 1; fi
	@if [ ! -d $(METAPATHWAYS_DB) ]; then  echo "Fetching the database from S3 to $(METAPATHWAYS_DB)....";  mkdir $(METAPATHWAYS_DB); fi
	@if [ ! -d $(METAPATHWAYS_DB)/functional ]
	then
		cd $(METAPATHWAYS_DB)
		wget $(METAPATHWAYS_DB_URL)
		unzip MetaPathways_DBs.zip
		rm MetaPathways_DBs.zip
	fi

NOT_USED:
	@if [ ! -d $(METAPATHWAYS_DB) ]; then \
		mkdir $(METAPATHWAYS_DB); \
		echo  "Fetching the databases...."  \
		aws s3 cp s3://wbfogdog/a2ac7fc4db0bfae6c05ca12a5818792d/Metapathways_DBs_2016-04.tar.xz ${METAPATHWAYS_DB}/; \
		echo  "Unzipping the database...." 
		tar -xvJf ${METAPATHWAYS_DB}/Metapathways_DBs_2016-04.tar.xz  --directory $(METAPATHWAYS_DB);  \
		mv  ${METAPATHWAYS_DB}/MetaPathways_DBs/* $(METAPATHWAYS_DB)/;  \
	fi


.PHONY: $(GIT_SUBMODULE_UPDATE) 
$(GIT_SUBMODULE_UPDATE):
	@echo git submodule update  rpkm
	git submodule update  --init executables/source/rpkm 
	@echo git submodule update  bwa
	git submodule update  --init executables/source/bwa 
	@echo git submodule update  FAST
	git submodule update  --init executables/source/FAST 
	@echo git submodule update  prodigal
	git submodule update  --init executables/source/prodigal 

$(RPKM):  
	$(MAKE) $(CFLAGS) executables/source/rpkm 
	mv executables/source/rpkm/rpkm $(BINARY_FOLDER)/

$(BWA):  
	$(MAKE) $(CFLAGS) executables/source/bwa 
	mv executables/source/bwa/bwa $(BINARY_FOLDER)/

$(PRODIGAL):  
	$(MAKE) $(CFLAGS) executables/source/prodigal 
	mv executables/source/prodigal/prodigal $(BINARY_FOLDER)/

$(FAST):  
	$(MAKE) $(CFLAGS) executables/source/FAST
	mv executables/source/FAST/fastal $(BINARY_FOLDER)/
	mv executables/source/FAST/fastdb $(BINARY_FOLDER)/

$(BLASTP): $(NCBI_BLAST) 
	@echo -n "Extracting the binaries for BLAST...." 
	tar --extract --file=$(NCBI_BLAST)  $(NCBI_BLAST_VER)/bin
	mv $(NCBI_BLAST_VER)/bin/blastx  executables/$(OS_PLATFORM)/
	mv $(NCBI_BLAST_VER)/bin/blastp  executables/$(OS_PLATFORM)/
	mv $(NCBI_BLAST_VER)/bin/blastn  executables/$(OS_PLATFORM)/
	mv $(NCBI_BLAST_VER)/bin/makeblastdb  executables/$(OS_PLATFORM)/
	rm -rf  $(NCBI_BLAST_VER)
	rm -rf  $(NCBI_BLAST)
	@echo "done" 

$(NCBI_BLAST):
	@echo -n "Downloading BLAST from NCBI website...." 
	wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$(NCBI_BLAST)
	@echo "done" 


$(METAPATHWAYS_DB_TAR):
	@echo  "Fetching the databases...." 
	aws s3 cp s3://wbfogdog/a2ac7fc4db0bfae6c05ca12a5818792d/Metapathways_DBs_2016-04.tar.xz .

$(METAPATHWAYS_DB): $(METAPATHWAYS_DB_TAR)
	@echo  "Unzipping the database...." 
	tar -xvJf Metapathways_DBs_2016-04.tar.xz


$(BINARY_FOLDER): 
	@if [ ! -d $(BINARY_FOLDER) ]; then mkdir $(BINARY_FOLDER); fi


### Utilities:
clean:
	$(MAKE) $(CFLAGS) executables/source/rpkm clean
	$(MAKE) $(CFLAGS) executables/source/prodigal.v2_00 clean
	$(MAKE) $(CFLAGS) executables/source/bwa clean

remove:
	rm -rf ../$(OS_PLATFORM)/bwa  
	rm -rf ../$(OS_PLATFORM)/prodigal
	rm -rf ../$(OS_PLATFORM)/rpkm 

