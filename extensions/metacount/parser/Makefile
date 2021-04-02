VERSION ?= 1.0.14

.PHONY: package

SUBDIRS=general bam samtools fastq

include Makefiles/Makefile.base


clean:$(SUBDIRS)
	rm -f $(STAT_GEN_LIB_OPT)
	rm -f $(STAT_GEN_LIB_DEBUG)
	rm -f $(STAT_GEN_LIB_PROFILE)

# general depends on samtools
general: samtools

# other subdirectories depend on general
bam: general

RELEASE_FILE?=libbamparser.$(VERSION).tgz

# Package the library.
package : 
# the touch gets rid of a tar warning
	touch $(RELEASE_FILE)
	tar chvz --exclude="*~" --exclude=$(RELEASE_FILE) --exclude='obj/*' --exclude='*.a'  --exclude='include/*' --exclude='bin/*' --exclude='test/results/*' --exclude-vcs -f $(RELEASE_FILE) --transform 's,^,libbamparser_$(VERSION)/,' * --show-transformed-names 
