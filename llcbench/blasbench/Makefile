SHELL	= /bin/sh

include ../user.def
include ../sys.def

BBOPTS	= -i $(BB_Iterations) -e $(BB_RepeatCount) -m $(BB_Memsize) -x $(BB_Resolution) -c -d $(BB_ReportMflops) -o

all: vblasbench

vblasbench: bb.o 
	$(BB_CC) $(BB_LDFLAGS) -o $@ bb.o $(VBLASLIB) $(BB_LIBS)

blasbench: bb.o 
	-if [ -f "$(BLASLIB)" ]; then $(BB_F77) $(BB_LDFLAGS) -o $@ bb.o $(BLASLIB) $(BB_LIBS); fi;

bb.o:
	$(BB_CC) $(BB_CFLAGS) -c bb.c

blasperf: 
	@cd results; gnuplot < ${Result_Prefix}_blasref.gp > /dev/null
	-cd results; gnuplot < ${Result_Prefix}_blasref.gp > ${Result_Prefix}_blasref.ps
	@echo ""; @echo "UTK BLAS graph is in results/${Result_Prefix}_blasref.ps"; echo

vblasperf: 
	@cd results; gnuplot < ${Result_Prefix}_blas.gp > /dev/null
	-cd results; gnuplot < ${Result_Prefix}_blas.gp > ${Result_Prefix}_blas.ps
	@echo ""; echo "Vendor BLAS graph is in results/${Result_Prefix}_blas.ps"; echo ""

compare: 
	@cd results; gnuplot < ${Result_Prefix}_blascomp.gp > /dev/null
	-cd results; gnuplot < ${Result_Prefix}_blascomp.gp > ${Result_Prefix}_blascomp.ps
	echo ""; @echo "BLAS comparison graph is in results/${Result_Prefix}_blascomp.ps"; echo

graph graphs: vblasperf
	-@cd results; tar cf ${Result_Prefix}-blasbench-graphs.tar ${Result_Prefix}_blas*.ps; 

compare.gp: 
	-@sed -e "s|TITLE|Vendor vs. UTK BLAS performance|g" < blasgraph.gp > results/$@; echo "plot \"${Result_Prefix}_blas_ix}.dat\" title 'daxpy' with linespoints \\" >> results/$@; for i in vdaxpy dgemv vdgemv dgemm vdgemm; do echo ", \"$$i-${Result_Prefix}.dat\" title '$$i' with linespoints \\" >> results/$@; done;

custom.gp: 
	-@sed -e "s|TITLE|UTK BLAS performance|g" < blasgraph.gp > results/$@; echo "plot \"daxpy-${Result_Prefix}.dat\" title 'saxpy' with linespoints \\" >> results/$@; for i in dgemv dgemm; do echo ", \"$$i-${Result_Prefix}.dat\" title '$$i' with linespoints \\" >> results/$@; done;

run: vdaxpy vdgemv vdgemm
	-@sed -e "s|TITLE|Vendor BLAS performance|g" < blasgraph.gp > results/${Result_Prefix}_blas.gp; echo "plot \"${Result_Prefix}_blas_vdaxpy.dat\" title 'vdaxpy' with linespoints \\" >> results/${Result_Prefix}_blas.gp; for i in vdgemv vdgemm; do echo ", \"${Result_Prefix}_blas_$$i.dat\" title '$$i' with linespoints \\" >> results/${Result_Prefix}_blas.gp; done;
	-@cd results; tar cf ${Result_Prefix}-blasbench-datafiles.tar ${Result_Prefix}_blas*.dat ${Result_Prefix}_blas.gp; 
	@echo ""; echo "Datafiles and GNUPlot scripts are located in the results directory."; echo "";

saxpy: blasbench
	-if [ -x blasbench ]; then blasbench $(BBOPTS) -vs > results/${Result_Prefix}_blas_$@.dat; fi
daxpy: blasbench
	-if [ -x blasbench ]; then blasbench  $(BBOPTS) -v > results/${Result_Prefix}_blas_$@.dat; fi
sgemv: blasbench
	-if [ -x blasbench ]; then blasbench  $(BBOPTS) -as > results/${Result_Prefix}_blas_$@.dat; fi
dgemv: blasbench
	-if [ -x blasbench ]; then blasbench  $(BBOPTS) -a > results/${Result_Prefix}_blas_$@.dat; fi
sgemm: blasbench
	-if [ -x blasbench ]; then blasbench  $(BBOPTS) -ts > results/${Result_Prefix}_blas_$@.dat; fi
dgemm: blasbench
	-if [ -x blasbench ]; then blasbench  $(BBOPTS) -t > results/${Result_Prefix}_blas_$@.dat; fi
vsaxpy: vblasbench
	@echo "Measuring SAXPY..."
	vblasbench  $(BBOPTS) -vs > results/${Result_Prefix}_blas_$@.dat
vdaxpy: vblasbench
	@echo "Measuring DAXPY..."
	vblasbench  $(BBOPTS) -v > results/${Result_Prefix}_blas_$@.dat
vsgemv: vblasbench
	@echo "Measuring SGEMV..."
	vblasbench  $(BBOPTS) -as > results/${Result_Prefix}_blas_$@.dat
vdgemv: vblasbench
	@echo "Measuring DGEMV..."
	vblasbench  $(BBOPTS) -a > results/${Result_Prefix}_blas_$@.dat
vsgemm: vblasbench
	@echo "Measuring SGEMM..."
	vblasbench  $(BBOPTS) -ts > results/${Result_Prefix}_blas_$@.dat
vdgemm: vblasbench
	@echo "Measuring DGEMM..."
	vblasbench  $(BBOPTS) -t > results/${Result_Prefix}_blas_$@.dat

clean:
	-rm -f *.o *~ core */*~
clobber: clean
	-rm -f blasbench vblasbench
doc:
	cd doc; $(MAKE)
