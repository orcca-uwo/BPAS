## BPAS/Main/include/Makefile
SHELL = bash
BPAS_HOME = $(PWD)/..
BPAS_INCLUDE =  $(BPAS_HOME)/include/
BPAS_DISTRIB = $(BPAS_HOME)/BPAS/
EXCLUDE_FILE = $(PWD)/Makefile.distrib.exclude

distrib_32bit:
	mkdir -p $(BPAS_DISTRIB)/include
	(cp -pf bpas.h global.h ring.h polynomial.h $(BPAS_DISTRIB)/include/)
	mkdir -p $(BPAS_DISTRIB)/include/IntegerPolynomial/
	(cd $(BPAS_INCLUDE)/IntegerPolynomial ; cp -pf *.h $(BPAS_DISTRIB)/include/IntegerPolynomial/ ; cp -pfr Multiplication $(BPAS_DISTRIB)/include/IntegerPolynomial/ ; rm -f $(BPAS_DISTRIB)/include/IntegerPolynomial/Multiplication/*64*bit* )
	mkdir -p $(BPAS_DISTRIB)/include/FFT/src/ 
	(cd $(BPAS_INCLUDE)/FFT/ ; cp -pf *.h $(BPAS_DISTRIB)/include/FFT/ )
	(cd $(BPAS_INCLUDE)/FFT/src/ ; cp -pfr * $(BPAS_DISTRIB)/include/FFT/src/ )
	mkdir -p $(BPAS_DISTRIB)/include/RationalNumberPolynomial/
	(cd $(BPAS_INCLUDE)/RationalNumberPolynomial/ ; cp -pfr * $(BPAS_DISTRIB)/include/RationalNumberPolynomial/ )
	mkdir -p $(BPAS_DISTRIB)/include/Parser/
	(cd $(BPAS_INCLUDE)/Parser/ ; cp -pfr * $(BPAS_DISTRIB)/include/Parser/ )
	mkdir -p $(BPAS_DISTRIB)/include/Polynomial/
	(cd $(BPAS_INCLUDE)/Polynomial/ ; cp -pfr * $(BPAS_DISTRIB)/include/Polynomial/ )
	mkdir -p $(BPAS_DISTRIB)/include/TriangularSet/
	(cd $(BPAS_INCLUDE)/TriangularSet/ ; cp -pfr * $(BPAS_DISTRIB)/include/TriangularSet/ )
	mkdir -p $(BPAS_DISTRIB)/include/RationalFunction/
	(cd $(BPAS_INCLUDE)/RationalFunction/ ; cp -pfr * $(BPAS_DISTRIB)/include/RationalFunction/ )

distrib_64bit:
	mkdir -p $(BPAS_DISTRIB)/include
	(rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/)

	#DyadicRationalNumber
	mkdir -p $(BPAS_DISTRIB)/include/DyadicRationalNumber/Multiplication
	(cd $(BPAS_INCLUDE)/DyadicRationalNumber ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/DyadicRationalNumber)
	(cd $(BPAS_INCLUDE)/DyadicRationalNumber/Multiplication ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/DyadicRationalNumber/Multiplication)

	#ExpressionTree
	mkdir -p $(BPAS_DISTRIB)/include/ExpressionTree
	(cd $(BPAS_INCLUDE)/ExpressionTree ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/ExpressionTree)

	#FFT
	mkdir -p $(BPAS_DISTRIB)/include/FFT/src/
	(cd $(BPAS_INCLUDE)/FFT ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/FFT)
	(cd $(BPAS_INCLUDE)/FFT/src ; rsync -avhWr --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *  $(BPAS_DISTRIB)/include/FFT/src)

	#FiniteFields
	mkdir -p $(BPAS_DISTRIB)/include/FiniteFields
	(cd $(BPAS_INCLUDE)/FiniteFields ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/FiniteFields)

	#IntegerPolynomial
	mkdir -p $(BPAS_DISTRIB)/include/IntegerPolynomial/Multiplication
	(cd $(BPAS_INCLUDE)/IntegerPolynomial ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/IntegerPolynomial)
	(cd $(BPAS_INCLUDE)/IntegerPolynomial/Multiplication ; rsync -avhWr --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *  $(BPAS_DISTRIB)/include/IntegerPolynomial/Multiplication)

	#Interval
	mkdir -p $(BPAS_DISTRIB)/include/Interval
	(cd $(BPAS_INCLUDE)/Interval ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/Interval)

	#ModularPolynomial
	mkdir -p $(BPAS_DISTRIB)/include/ModularPolynomial
	(cd $(BPAS_INCLUDE)/ModularPolynomial ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/ModularPolynomial)
	
	#RationalFunction
	mkdir -p $(BPAS_DISTRIB)/include/RationalFunction
	(cd $(BPAS_INCLUDE)/RationalFunction ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/RationalFunction)
	
	#RationalNumberPolynomial
	mkdir -p $(BPAS_DISTRIB)/include/RationalNumberPolynomial
	(cd $(BPAS_INCLUDE)/RationalNumberPolynomial ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/RationalNumberPolynomial)
	
	#Parser
	mkdir -p $(BPAS_DISTRIB)/include/Parser
	(cd $(BPAS_INCLUDE)/Parser ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/Parser)

	#RegularChain
	mkdir -p $(BPAS_DISTRIB)/include/RegularChain
	(cd $(BPAS_INCLUDE)/RegularChain ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/RegularChain)
	
	#Ring
	mkdir -p $(BPAS_DISTRIB)/include/Ring
	(cd $(BPAS_INCLUDE)/Ring ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/Ring)
	
	#RingPolynomial
	mkdir -p $(BPAS_DISTRIB)/include/RingPolynomial
	(cd $(BPAS_INCLUDE)/RingPolynomial ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/RingPolynomial)

	#Symbol
	mkdir -p $(BPAS_DISTRIB)/include/Symbol
	(cd $(BPAS_INCLUDE)/Symbol ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/Symbol)
	
	#TriangularSet
	mkdir -p $(BPAS_DISTRIB)/include/TriangularSet
	(cd $(BPAS_INCLUDE)/TriangularSet ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/TriangularSet)
	
	#Utils
	mkdir -p $(BPAS_DISTRIB)/include/Utils
	(cd $(BPAS_INCLUDE)/Utils ; rsync -avhW --no-compress --ignore-missing-args --progress --exclude-from=$(EXCLUDE_FILE) *.h *.hpp  $(BPAS_DISTRIB)/include/Utils)
	
	#Finally, clean any empty dirs. This maybe that an entire directory had files excluded
	(cd $(BPAS_DISTRIB)/include ; find . -mindepth 1 -type d -print | xargs rmdir --ignore-fail-on-non-empty 2>/dev/null)

