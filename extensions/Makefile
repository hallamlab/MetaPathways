ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

all:
	make -C metacount
	make -C rpkm

.PHOHY: clean
clean:
	make -C metacount clean
	make -C rpkm clean

.PHONY: install
install:
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	install metacount/metacount $(DESTDIR)$(PREFIX)/bin/
	install rpkm/rpkm $(DESTDIR)$(PREFIX)/bin/

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/metacount
	rm -f $(DESTDIR)$(PREFIX)/bin/rpkm
