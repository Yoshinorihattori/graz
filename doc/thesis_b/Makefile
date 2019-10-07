TARGET = thesis.pdf

all : $(TARGET)

thesis.pdf : grad_thesis.sty thesis.tex 1_introduction.tex 2_appararus.tex
	platex thesis.tex
	platex thesis.tex
	dvipdfmx -p a4 thesis.dvi
	rm -f *.dvi *.aux *.log *.toc *~

clean :
	rm -f *.dvi *.aux *.log *.toc *~
	rm -f $(TARGET)
