all: mpbench.ps clean

mpbench.ps: mpbench.tex
	latex mpbench.tex
	latex mpbench.tex
	dvips -o mpbench.ps mpbench.dvi
	ps2pdf mpbench.ps mpbench.pdf

clean:
	rm -f *.log *.aux *~ *.dvi

clobber: clean
	rm -f *.ps *.pdf
