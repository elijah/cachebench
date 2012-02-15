all:
	cd blasbench; $(MAKE)
	cd cachebench; $(MAKE)
	cd mpbench; $(MAKE)

clean:
	cd blasbench; $(MAKE) clean
	cd cachebench; $(MAKE) clean
	cd mpbench; $(MAKE) clean
	-rm -f *~

clobber: clean
	cd blasbench; $(MAKE) clobber
	cd cachebench; $(MAKE) clobber
	cd mpbench; $(MAKE) clobber

