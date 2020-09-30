#### Assumptions:
## * AWS CLI installed & configured
## * sudo apt-get install make python2.7

#### Make configuration:

## Use Bash as default shell, and in strict mode:
SHELL := /bin/bash
.SHELLFLAGS = -ec

## If the parent env doesn't ste TMPDIR, do it ourselves:
TMPDIR ?= /tmp


## This makes all recipe lines execute within a shared shell process:
## https://www.gnu.org/software/make/manual/html_node/One-Shell.html#One-Shell
.ONESHELL:

## If a recipe contains an error, delete the target:
## https://www.gnu.org/software/make/manual/html_node/Special-Targets.html#Special-Targets
.DELETE_ON_ERROR:

## This is necessary to make sure that these intermediate files aren't clobbered:
.SECONDARY:



### Docker Automation
docker-start:
	sudo systemctl start docker

docker-build: pre-docker-builds
	sudo docker build --network=host -t taltman/metapathways:taltman_dev .

docker-run:
	sudo docker run -it --rm -v $(CURDIR):/input -v $(CURDIR)/out:/output taltman/darth:maul bash 

docker-deploy:
	sudo docker login
	sudo docker push taltman/darth:maul

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
TRNASCAN=$(BINARY_FOLDER)/trnascan-1.4
FAST=$(BINARY_FOLDER)/fastal
PRODIGAL=$(BINARY_FOLDER)/prodigal

METAPATHWAYS_DB_DEFAULT=../fogdogdatabases
METAPATHWAYS_DB_TAR=Metapathways_DBs_2016-04.tar.xz
METAPATHWAYS_DB_URL=https://www.dropbox.com/s/ye3kpve041e0r39/MetaPathways_DBs.zip


GIT_SUBMODULE_UPDATE=gitupdate
# Alias for target 'all', for compliance with FogDog deliverables standard:

#all: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA) $(TRNASCAN)  $(RPKM)
all: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA) $(TRNASCAN)  $(RPKM) $(BLASTP) METAPATHWAYS_DB_FETCH
pre-docker-builds: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA) $(TRNASCAN)  $(RPKM) $(BLASTP) 


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
	@echo git submodule update  trnascan
	git submodule update  --init executables/source/trnascan 
	@echo git submodule update  rpkm
	git submodule update  --init executables/source/rpkm 
	@echo git submodule update  bwa
	git submodule update  --init executables/source/bwa 
	@echo git submodule update  FAST
	git submodule update  --init executables/source/FAST 
	@echo git submodule update  prodigal
	git submodule update  --init executables/source/prodigal 

$(TRNASCAN):  
	$(MAKE) $(CFLAGS) executables/source/trnascan 
	mv executables/source/trnascan/trnascan-1.4 $(BINARY_FOLDER)/

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
	$(MAKE) $(CFLAGS) executables/source/trnascan clean
	$(MAKE) $(CFLAGS) executables/source/rpkm clean
	$(MAKE) $(CFLAGS) executables/source/prodigal.v2_00 clean
	$(MAKE) $(CFLAGS) executables/source/FAST clean
	$(MAKE) $(CFLAGS) executables/source/bwa clean

remove:
	rm -rf  ../$(OS_PLATFORM)/trnascan-1.4 
	rm -rf ../$(OS_PLATFORM)/fastal  
	rm -rf ../$(OS_PLATFORM)/fastdb  
	rm -rf ../$(OS_PLATFORM)/bwa  
	rm -rf ../$(OS_PLATFORM)/prodigal
	rm -rf ../$(OS_PLATFORM)/rpkm 

### Testing:

mp-regression-tests:
	./run_regtests.sh
	@exit $$?

## Top-level test target
test: test-mp-regression-tests

no-ptools-unit-test:
	mkdir -p test
	source MetaPathwaysrc
	touch executables/linux/FGS+
	touch executables/linux/ptools
	touch /tmp/mp_db_dir/MetaPathways_DBs/ncbi_tree/RefSeq-release80.catalog
	time ./MetaPathways.py -i regtests/input/B1.fasta -o test/B1_MPout/ -p test/mp_param.txt -c test/mp_config.txt

docker-test:
	cp $(CURDIR)/regtests/input/A1.fasta /tmp
	mkdir -p /tmp/mp_db_dir/MetaPathways_DBs/functional/formatted
	mkdir -p /tmp/mp_db_dir/MetaPathways_DBs/taxonomic/formatted
	mkdir -p /tmp/mp_db_dir/MetaPathways_DBs/ncbi_tree/formatted
	touch /tmp/mp_db_dir/MetaPathways_DBs/ncbi_tree/RefSeq-release80.catalog
	sudo docker run -it -v /tmp:/input taltman/metapathways:taltman_dev /root/mp_repo/MetaPathways.py -i /input/A1.fasta -o /input/A1_MP_out/ -p /root/mp_repo/resources/docker_param.txt -c /root/mp_repo/resources/docker_config.txt
