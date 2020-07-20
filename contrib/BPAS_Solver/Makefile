CC = g++
INCLUDES = ../../include
CPPCOMPILEARG = -std=c++11 -O2 -I$(INCLUDES) -fcilkplus -DLINUXINTEL64=1 -march=x86-64 -DWITH_MAPLE
LDARGS = -L../../build/ -lmodpnLINUXINTEL64 -lmaplec -static-libstdc++ -static-libgcc -lcilkrts -lgmpxx -lgmp -lm -lpthread 

vpath %.cpp ./

BPAS_LIB = ../../build/libbpas.a
NTL_LIB = /usr/lib/libntl.a
TARGET = BPAS_Solver.bin

all: ${TARGET}

BPAS_Solver.bin: BPAS_Solver.o
	$(CC) -o $@ $< $(BPAS_LIB)  $(LDARGS)

%.o: %.cpp
	$(CC) -c $(CPPCOMPILEARG) $<

clean:
	rm -f *.o ${TARGET}


