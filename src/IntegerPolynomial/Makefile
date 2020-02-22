#Makefile for mul

CILKPP = g++ 
CC = gcc

INCLUDES = ../../include
COMPILEARG = -std=c++11 -I$(INCLUDES) -c -O2 -fcilkplus -fpic -DLINUXINTEL64=1 -march=native
CCOMPILEARG = -std=gnu11 -I$(INCLUDES) -c -O2 -Wall -Wno-unused-function -fpic -DLINUXINTEL64=1 -march=native

#-Wall
ifeq ($(DEBUG),1)
	COMPILEARG+= -DTDEBUG=1
endif
INCLUDE_DIR = ../../include/IntegerPolynomial

vpath %.cpp Multiplication/ ./
vpath %.h $(INCLUDE_DIR)/ $(INCLUDE_DIR)/Multiplication/

#############################################################################################
###################### Build, Compile ### ###################################################
#############################################################################################

all: DUZP_Support.o SMZP_Support.o SMZP_Support_Unpacked.o SMZP_Support_Test.o SMZP_Support_Recursive.o SMZP_Support_Recursive_Unpacked.o mzpolynomial.o uzpolynomial.o Poly.o BivariatePoly.o Mul.o MulNaive.o MulKS.o MulDnC.o MulToom4.o MulToom8.o MulSSA.o modpoly.o

%.o: %.cpp
	$(CILKPP) $(COMPILEARG) $<

%.o: %.c
	$(CC) $(CCOMPILEARG) $<

serial: COMPILEARG += -DSERIAL=1
serial: all #uzpolynomial.o Poly.o BivariatePoly.o Mul.o MulNaive.o MulKS.o MulDnC.o MulToom4.o MulToom8.o MulSSA.o modpoly.o

BLAD: COMPILEARG += -DWITH_BLAD=1
BLAD: CCOMPILEARG += -DWITH_BLAD=1
BLAD: all

test:

clean:
	rm -rf *.o
