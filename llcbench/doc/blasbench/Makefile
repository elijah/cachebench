all: blasbench.ps clean

blasbench.ps: blasbench.tex
	latex blasbench.tex
	latex blasbench.tex
	dvips -o blasbench.ps blasbench.dvi
	ps2pdf blasbench.ps blasbench.pdf

clean:
	rm -f *.log *.aux *~ *.dvi

clobber: clean
	rm -f *.ps *.pdf
