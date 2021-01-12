ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

all:
	make -C FAST
	make -C trnascan
	make -C metacount

.PHOHY: clean
clean:
	make -C FAST clean
	make -C trnascan clean
	make -C metacount clean

.PHONY: install
install: bwa/bwa
	mkdir -p $(DESTDIR)$(PREFIX)/bin
#
	install FAST/fastal $(DESTDIR)$(PREFIX)/bin/
	install FAST/fastdb $(DESTDIR)$(PREFIX)/bin/
#
	install trnascan/trnascan-1.4 $(DESTDIR)$(PREFIX)/bin/
#
	install metacount/metacount $(DESTDIR)$(PREFIX)/bin/

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/fastal  $(DESTDIR)$(PREFIX)/bin/fastdb
	rm -f $(DESTDIR)$(PREFIX)/bin/trnascan-1.4
	rm -f $(DESTDIR)$(PREFIX)/bin/metacount
