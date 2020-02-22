#Makefile for mul

CILKPP = g++ 
COMPILEARG = -std=c++11 -c -O2 -fcilkplus -fpic -DLINUXINTEL64=1 -march=native
#-Wall
ifeq ($(DEBUG),1)
	COMPILEARG+= -DTDEBUG=1
endif
INCLUDE_DIR = ../../include/RingPolynomial



#############################################################################################
###################### Build, Compile ### ###################################################
#############################################################################################

all: dupolynomial.o 

%.o: %.cpp 
	$(CILKPP) $(COMPILEARG) $<

serial: COMPILEARG += -DSERIAL=1
serial: dupolynomial.o 

test:

clean:
	rm -rf *.o
