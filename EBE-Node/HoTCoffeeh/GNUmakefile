all	: svwr.e cfwr.e

svwr.e  :
		(cd svwr; make -f GNUmakefile_svwr distclean; make -f GNUmakefile_svwr install)
cfwr.e	:
		(cd cfwr; make -f GNUmakefile && make -f GNUmakefile install)

distclean	:
	(cd svwr; make -f GNUmakefile_svwr distclean)
	(cd cfwr; make -f GNUmakefile distclean)
	rm *.e
