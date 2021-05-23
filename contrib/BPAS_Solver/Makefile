CC = g++
INCLUDES = ../../include
CPPCOMPILEARG = -std=c++11 -O2 -I$(INCLUDES) -fcilkplus -DLINUXINTEL64=1 -march=x86-64 
LDARGS = -L../../build/ -lmodpnLINUXINTEL64 -lntl -static-libstdc++ -static-libgcc -lcilkrts -lgmpxx -lgmp -lm -lpthread 

vpath %.cpp ./

BPAS_LIB = ../../build/libbpas.a
TARGET = BPAS_Solver.bin

all: ${TARGET}

with_maple: clean
with_maple: CPPCOMPILEARG += -DWITH_MAPLE=1
with_maple: LDARGS += -lmaplec
with_maple: ${TARGET}

with_maple_serial: clean
with_maple_serial: CPPCOMPILEARG += -DWITH_MAPLE=1 -DSERIAL=1
with_maple_serial: LDARGS += -lmaplec
with_maple_serial: ${TARGET}


BPAS_Solver.bin: BPAS_Solver.o
	$(CC) -o $@ $< $(BPAS_LIB)  $(LDARGS)

%.o: %.cpp
	$(CC) -c $(CPPCOMPILEARG) $<

clean:
	rm -f *.o ${TARGET}


