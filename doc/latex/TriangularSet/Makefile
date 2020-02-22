filename=triangularset-doc

all : pdf

pdf: 
	pdflatex ${filename}
	pdflatex ${filename}

clean:
	bash -c 'rm -f ${filename}.{ps,pdf,log,aux,out,dvi,bbl,blg,snm,bcf,toc,nav,run.xml}'
