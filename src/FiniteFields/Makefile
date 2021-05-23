CILKPP  = g++
CC = g++ #-pg -g3 -g
INCLUDES = ../../include
COMPILEARG = -I$(INCLUDES) -lgmp -c -O3 -fpic -w -fcilkplus -DLINUXINTEL64=1 -march=native

CCOMPILEARG= -I$(INCLUDES) -lgmp -c -O3 -fpic -w -DLINUXINTEL64 -march=native
								 #-DPROFILING_ENABLED -DPROFILING_GMP_ENABLED=1
#-Wall
vpath %.cpp ./ 
vpath %.c ./ 
vpath %.h $(INCLUDES)

#######################################
SRC_LIST:=\
	SmallPrimeField\
	BigPrimeField\
	GeneralizedFermatPrimeField\
	BigPrimeField_Support\
	SmallPrimeField_Support\
	small_prime_field_fft\

#######################################
OBJS:=$(SRC_LIST:=.o)
#######################################
all: $(OBJS)
	make -f Makefile.new
#######################################
serial-gmp: $(OBJS) 
	make -f Makefile.new serial-gmp
#######################################
%.o: %.cpp
	$(CILKPP) $(COMPILEARG) $<

#GFPF_Support.o: GFPF_Support.c
#	$(CC) $(CCOMPILEARG) $<

BigPrimeField_Support.o: BigPrimeField_Support.cpp
	$(CC) $(CCOMPILEARG) $<

SmallPrimeField_Support.o: SmallPrimeField_Support.c
	$(CC) $(CCOMPILEARG) $<

test:


clean:
	rm -rf *.o
	make -f Makefile.new clean

serial: serial-gmp
	make -f Makefile.new serial
serial: COMPILEARG += -DSERIAL=1
serial: SmallPrimeField.o BigPrimeField.o GeneralizedFermatPrimeField.o BigPrimeField_Support.o SmallPrimeField_Support.o

