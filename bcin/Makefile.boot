sun:
	@cp make.rules.sun make.rules
	-@mv Makefile Makefile.bak
	@/usr/ccs/lib/cpp -DsunMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v "^#" | sed -e "s/^%/#/" > Makefile

win32:
	@cp make.rules.win32 startup.mk
	-@mv -f Makefile Makefile.bak
	@cl /EP /D"win32MAKE" Makefile.in \
		| grep -v "^#" | sed -f make.sed > Makefile

sgi:
	-@mv Makefile Makefile.bak
	@cc -E -DsgiMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v '^#' | sed -e "s/^%/#/" > Makefile
