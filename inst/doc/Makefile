
all: fabia.pdf clean

clean:
	rm -f *.log *.aux *.bak *~ *.out *.tex 
	rm -rf .Rd2dvi*

fabia.pdf: fabia.tex Makefile bioinf.cls bioinf-article.txi
	pdflatex fabia 
	bibtex fabia
	pdflatex fabia

fabia.tex: fabia.Rnw Makefile
	R CMD Sweave fabia.Rnw