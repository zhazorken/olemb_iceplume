sun:
	@cp make.rules.sun make.rules
	-@mv Makefile Makefile.bak
	@cp ic.inc.in ic.inc
	@/usr/ccs/lib/cpp -DsunMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v "^#" | sed -e "s/^%/#/" > Makefile

win32:
	@cp make.rules.win32 startup.mk
	-@mv -f Makefile Makefile.bak
	@sed -f fortpowr.sed ic.inc.in > ic.inc
	@cl /EP /D"win32MAKE" Makefile.in \
		| grep -v "^#" | sed -f make.sed > Makefile

sgi:
	-@mv Makefile Makefile.bak
	@cp ic.inc.in ic.inc
	@cc -E -DsgiMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v '^#' | sed -e "s/^%/#/" > Makefile
