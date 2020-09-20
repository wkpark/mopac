all: mopac.pdf

mopac.dvi: mopac.tex
	latex mopac.tex
	makeindex mopac.idx
	latex mopac.tex

mopac.pdf: mopac.tex
	pdflatex mopac.tex
	makeindex mopac.idx
	pdflatex mopac.tex

clean:
	rm -f mopac.aux mopac.bbl mopac.dvi mopac.pdf mopac.toc mopac.log mopac.ind mopac.log mopac.toc mopac.idx mopac.out mopac.ilg
