
all: manual_7.1.pdf 

manual_7.1.pdf:manual_7.1.tex
	pdflatex  manual_7.1.tex
	bibtex manual_7.1
	pdflatex   manual_7.1.tex
	pdflatex   manual_7.1.tex
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg .DS_Store