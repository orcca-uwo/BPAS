## BPAS/Main/src/Makefile

DOCS = RegularChain TriangularSet 
DOCS_DIR = $(PWD)

ALLTARGS = $(DOCS:%=%.all)
CLEANTARGS = $(DOCS:%=%.clean)

all: $(ALLTARGS)

clean: $(CLEANTARGS)

%.all:
	(cd $(DOCS_DIR)/$* ; make)
	
%.clean:
	(cd $(DOCS_DIR)/$* ; make clean)

	
