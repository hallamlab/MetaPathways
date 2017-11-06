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

NCBI_BLAST=ncbi-blast-2.6.0+-x64-linux.tar.gz
NCBI_BLAST_VER=ncbi-blast-2.6.0+
BINARY_FOLDER=executables/$(OS_PLATFORM)


BLASTP=$(BINARY_FOLDER)/blastp
LASTAL=$(BINARY_FOLDER)/lastal+
RPKM=$(BINARY_FOLDER)/rpkm
BWA=$(BINARY_FOLDER)/bwa
TRNASCAN=$(BINARY_FOLDER)/trnascan-1.4
FAST=$(BINARY_FOLDER)/fastal
PRODIGAL=$(BINARY_FOLDER)/prodigal

METAPATHWAYS_DB_DEFAULT=../fogdogdatabases
PTOOLS_DIR=../
METAPATHWAYS_DB_TAR=Metapathways_DBs_2016-04.tar.xz



GIT_SUBMODULE_UPDATE=gitupdate

## Alias for target 'all', for compliance with FogDog deliverables standard:
all: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA) $(TRNASCAN)  $(RPKM)
#all: $(GIT_SUBMODULE_UPDATE) $(BINARY_FOLDER) $(PRODIGAL)  $(FAST)  $(BWA) $(TRNASCAN)  $(RPKM) $(MICROBE_CENSUS) METAPATHWAYS_DB_FETCH PTOOLS_FETCH PTOOLS_INSTALL
#all: PTOOLS_FETCH

install-without-ptools: all METAPATHWAYS_DB_FETCH

install-with-ptools: all METAPATHWAYS_DB_FETCH 

.PHONY: METAPATHWAYS_DB_FETCH
METAPATHWAYS_DB_FETCH:
	@if [ -z $(METAPATHWAYS_DB) ]; then  echo "Variable METAPATHWAYS_DB not set. Set it as export METPATHWAYS_DB=<path>" ;  exit 1; fi
	@if [ ! -d $(METAPATHWAYS_DB) ]; then  echo "Fetching the database from S3 to $(METAPATHWAYS_DB)....";  mkdir $(METAPATHWAYS_DB); fi
	@if [ ! -d $(METAPATHWAYS_DB)/functional ]; then  aws s3 sync s3://fogdogdatabases  $(METAPATHWAYS_DB)/; fi

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


mp-regression-tests:
	./run_regtests.sh
	@exit $$?


## Top-level test target
test: test-mp-regression-tests


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

