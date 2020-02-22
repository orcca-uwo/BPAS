
all : regularchain-doc zerodimensionalregularchain-doc

regularchain-doc: 
	pdflatex $@
	pdflatex $@

zerodimensionalregularchain-doc:
	pdflatex $@
	pdflatex $@

clean: clean-rc clean-zdrc

clean-rc:
	bash -c 'rm -f regularchain-doc.{ps,pdf,log,aux,out,dvi,bbl,blg,snm,bcf,toc,nav,run.xml}'

clean-zdrc:
	bash -c 'rm -f zerodimensionalregularchain-doc.{ps,pdf,log,aux,out,dvi,bbl,blg,snm,bcf,toc,nav,run.xml}'

