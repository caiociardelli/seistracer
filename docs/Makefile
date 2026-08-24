all:
	pdflatex Manual.tex
	biber Manual
	pdflatex Manual.tex   # second pass for TOC / references / lastpage
	pdflatex Manual.tex   # third pass (just in case) for TOC / references / lastpage

clean:
	rm -f *.aux *.log *.out *.toc *.pdf