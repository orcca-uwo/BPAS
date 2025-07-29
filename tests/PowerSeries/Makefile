

CPPC = g++-7
CC = gcc-7
INCLUDES = ../../include
CCOMPILEARG = -I$(INCLUDES) -std=gnu99 -c -g -Wall -Wno-unused-function -march=native -DLINUXINTEL64=1
CPPCOMPILEARG = -I$(INCLUDES) -std=c++11 -c -g -Wall -Wno-unused-function -march=native -DLINUXINTEL64=1 -fcilkplus 
CLIBARG = -lgmpxx -lgmp -lm -lpthread
CPPLIBARG =  -lntl -lgmpxx -lgmp -lm -lpthread -lcilkrts

BPAS_LIB = ../../build/libbpas.a
MODPN_LIB = ../../build/libmodpnLINUXINTEL64.so

LODEoUPS_DEMO_src = LinearODEOverUnivariatePowerSeries_Demo_Tests

LODEoUPS_DEMO_BIN = LODEoUPS_demo_tests.bin

all: clean dobuild

dobuild: test

test: $(LODEoUPS_DEMO_BIN)

$(LODEoUPS_DEMO_src).o : $(LODEoUPS_DEMO_src).cpp
	$(CPPC) $(CPPCOMPILEARG) LinearODEOverUnivariatePowerSeries_Demo_Tests.cpp	

$(LODEoUPS_DEMO_BIN) : $(LODEoUPS_DEMO_src).o $(BPAS_LIB) $(MODPN_LIB)
	$(CPPC) -o $@ $^ $(CPPLIBARG)			

ODE_solving_demo_tests: $(LODEoUPS_DEMO_BIN)
	./$(LODEoUPS_DEMO_BIN)

clean :
	rm -rf *.o *.bin
