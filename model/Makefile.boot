#
# Makefile.boot - used to make system dependent files
#
sun sun-f77:
	@cp make.rules.sun make.rules
	-@mv Makefile Makefile.bak
	@cp olemf.inc.in olemf.inc
	@cp moddef.inc.in moddef.inc
	@cp varinit.inc.in varinit.inc
	@/usr/ccs/lib/cpp -DsunMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v '^#' | sed -e "s/^%/#/" > Makefile

sun-mpi:
	@cp make.rules.sun make.rules
	-@mv Makefile Makefile.bak
	@cp olemf.inc.in olemf.inc
	@cp moddef.inc.in moddef.inc
	@cp varinit.inc.in varinit.inc
	@/usr/ccs/lib/cpp -DsunmpiMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v '^#' | sed -e "s/^%/#/" > Makefile

exemplar:
	-@mv Makefile Makefile.bak
	@cp olemf.inc.in olemf.inc
	@cp moddef.inc.in moddef.inc
	@cp varinit.inc.in varinit.inc
	@/opt/langtools/lbin/cpp -DexemplarMAKE -Usun -Ucray \
		-Uunix Makefile.in \
		| grep -v '^#' | sed -e "s/^%/#/" > Makefile

win32:
	@cp make.rules.win32 startup.mk
	-@mv -f Makefile Makefile.bak
	@cp olemf.inc.in olemf.inc
	@cp moddef.inc.in moddef.inc
	@cp varinit.inc.in varinit.inc
	@cl /EP /D"win32MAKE" Makefile.in \
		| grep -v "^#" | sed -f make.sed > Makefile

sgi:
	-@mv Makefile Makefile.bak
	@cp olemf.inc.in olemf.inc
	@cp moddef.inc.in moddef.inc
	@cp varinit.inc.in varinit.inc
	@cc -E -DsgiMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v '^#' | sed -e "s/^%/#/" > Makefile

ibm:
	@cp make.rules.ibm rules.mk
	-@mv Makefile Makefile.bak
	@cp olemf.inc.in olemf.inc
	@cp moddef.inc.in moddef.inc
	@cp varinit.inc.in varinit.inc
	@/usr/ccs/lib/cpp -DibmMAKE -Usun -Ucray -Uunix Makefile.in \
		| grep -v '^#' | sed -e "s/^%/#/" > Makefile
