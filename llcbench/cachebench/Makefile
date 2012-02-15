SHELL	= /bin/sh
HOSTNAME = uname -n
FULLNAME = uname -smn

include ../user.def
include ../sys.def

CBOPTS  = -m $(CB_Memsize) -e $(CB_RepeatCount) -x $(CB_Resolution) -d $(CB_SecondsPerIteration)
DATA	= $(CB_Datatype)
TAG	= ""

cachebench: cachebench.c
	$(CB_CC) $(CB_CFLAGS) -DUSE_$(DATA) cachebench.c $(CB_LIBS) -o $@

run: read.dat write.dat rmw.dat handread.dat handwrite.dat handrmw.dat memset.dat memcpy.dat
	@sed -e "s|TITLE|Memory Hierarchy Performance|g" < cachegraph.gp > results/${Result_Prefix}_cache.gp; echo "plot \"${Result_Prefix}_cache_read.dat\" title 'read' with linespoints \\" >> results/${Result_Prefix}_cache.gp; for i in write rmw handread handwrite handrmw memset memcpy; do echo ", \"${Result_Prefix}_cache_$$i.dat\" title '$$i' with linespoints \\" >> results/${Result_Prefix}_cache.gp; done;
	@cd results; tar cf ${Result_Prefix}-cachebench-datafiles.tar ${Result_Prefix}_cache*.dat ${Result_Prefix}_cache.gp; echo ""; echo "Datafiles and GNUplot scripts are located in the results directory."; echo "";

graphs graph: 
	@cd results; gnuplot <  ${Result_Prefix}_cache.gp > /dev/null
	cd results; gnuplot < ${Result_Prefix}_cache.gp > ${Result_Prefix}_cache.ps
	@echo ""; echo "Memory heirarchy graph is in results/${Result_Prefix}_cache.ps"

	@cd results; tar cf ${Result_Prefix}-cachebench-graphs.tar *.ps; echo ""
clean:
	-rm -f *.o *~ core */*~

clobber: clean
	-rm -f cachebench

doc:
	cd doc; $(MAKE)

read.dat: cachebench 
	@echo "Measuring Read..."
	./cachebench  $(CBOPTS) -r > results/${Result_Prefix}_cache_$@
write.dat: cachebench 
	@echo "Measuring Write..."
	./cachebench  $(CBOPTS) -w > results/${Result_Prefix}_cache_$@
rmw.dat: cachebench 
	@echo "Measuring RMW..."
	./cachebench  $(CBOPTS) -b > results/${Result_Prefix}_cache_$@
handread.dat: cachebench 
	@echo "Measuring Tuned Read..."
	./cachebench  $(CBOPTS) -tr > results/${Result_Prefix}_cache_$@
handwrite.dat: cachebench 
	@echo "Measuring Tuned Write..."
	./cachebench  $(CBOPTS) -tw > results/${Result_Prefix}_cache_$@
handrmw.dat: cachebench 
	@echo "Measuring Tuned RMW..."
	./cachebench  $(CBOPTS) -tb > results/${Result_Prefix}_cache_$@
memset.dat: cachebench 
	@echo "Measuring memset()..."
	./cachebench  $(CBOPTS) -s > results/${Result_Prefix}_cache_$@
memcpy.dat: cachebench 
	@echo "Measuring memcpy()..."
	./cachebench  $(CBOPTS) -p > results/${Result_Prefix}_cache_$@

