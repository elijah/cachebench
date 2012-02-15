all: cachebench.ps clean

cachebench.ps: cachebench.tex
	latex cachebench.tex
	latex cachebench.tex
	dvips -o cachebench.ps cachebench.dvi
	ps2pdf cachebench.ps cachebench.pdf

clean:
	rm -f *.log *.aux *~ *.dvi

clobber:
	rm -f *.ps *.pdf