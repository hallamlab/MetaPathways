ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

all:
	make -C trnascan
	make -C metacount
	make -C rpkm

# Removing this command from the "all" target for now:
#make -C FAST


.PHOHY: clean
clean:
	make -C FAST clean
	make -C trnascan clean
	make -C metacount clean
	make -C rpkm clean

.PHONY: install
install:
	mkdir -p $(DESTDIR)$(PREFIX)/bin
#
#	install FAST/fastal $(DESTDIR)$(PREFIX)/bin/
#	install FAST/fastdb $(DESTDIR)$(PREFIX)/bin/
#
	install trnascan/trnascan-1.4 $(DESTDIR)$(PREFIX)/bin/
#
	install metacount/metacount $(DESTDIR)$(PREFIX)/bin/
#
	install rpkm/rpkm $(DESTDIR)$(PREFIX)/bin/

.PHONY: uninstall
uninstall:
#	rm -f $(DESTDIR)$(PREFIX)/bin/fastal  $(DESTDIR)$(PREFIX)/bin/fastdb
	rm -f $(DESTDIR)$(PREFIX)/bin/trnascan-1.4
	rm -f $(DESTDIR)$(PREFIX)/bin/metacount
	rm -f $(DESTDIR)$(PREFIX)/bin/rpkm
