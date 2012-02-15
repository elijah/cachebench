SHELL = /bin/sh

include sys.def

dum: 
	@echo "Please use one of the following targets:"
	@echo
	@echo "For all three : compile, run, graph, clean, clobber, reconfig"
	@echo "For BlasBench : blas-bench, blas-run, blas-graph"
	@echo "For CacheBench: cache-bench, cache-run, cache-graph"
	@echo "For MPBench   : mp-bench, mp-run, mp-graph"

compile bench: blas-bench mp-bench cache-bench

run: blas-run cache-run mp-run

graphs graph: blas-graph cache-graph mp-graph

blas-bench:
	cd blasbench; make  
blas-run: 
	cd blasbench; make run
blas-graph:
	cd blasbench; make graph
cache-bench:
	cd cachebench; make 
cache-run:
	cd cachebench; make run
cache-graph:
	cd cachebench; make graph
mp-bench:
	cd mpbench; make 
mp-run:
	cd mpbench; make run
mp-graph:
	cd mpbench; make graph

sys.def reconfig:
	-rm -f sys.def
	ln -s conf/sys.default sys.def

clean:
	cd mpbench; make clean
	cd blasbench; make clean
	cd cachebench; make clean
	cd doc; make clean
	rm -f *~ */*~

clobber: reconfig
	cd mpbench; make clobber
	cd blasbench; make clobber
	cd cachebench; make clobber
	cd doc; make clobber
	rm -f *~ */*~ results/*

dist: reconfig clean
	cd ..; tar cvf llcbench.tar llcbench; gzip llcbench.tar

install: dist
	cd 
	cp ../llcbench.tar.gz $(HOME)/www-home/llcbench
	cp www/*.html $(HOME)/www-home/llcbench
