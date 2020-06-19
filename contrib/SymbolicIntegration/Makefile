# Use: ./intTest

CILKPP	= g++ 
LIBARG  = -lbpas -lgmpxx -lgmp -lcilkrts -lmodpnLINUXINTEL64 -lmps -lpthread
#COMPILEARG = -c -O2 -g -fcilkplus -DLINUXINTEL64=1   
COMPILEARG = -c -O2 -I../../include -fpic -fcilkplus -DLINUXINTEL64=1 -march=native
TARGET  = intTest
#TARGET2 = issacDemo

all: clean $(TARGET) run_test
	#$(TARGET2)

$(TARGET): test.o realSymbolicNumericIntegratePFD.o floatIntegrate.o approximateIntegration.o polyParse.o
	$(CILKPP) -o $@ $^ $(LIBARG)

#$(TARGET2): issacDemo.o floatIntegrate.o approximateIntegration.o polyParse.o
#	$(CILKPP) -o $@ $^ $(LIBARG)

#../../include/Polynomial/upolynomial.h
%.o: %.cpp 
	$(CILKPP) $(COMPILEARG) $<

serial: COMPILEARG += -DSERIAL=1
serial: $(TARGET) $(TARGET2)

test: clean_test $(TARGET)
	./$(TARGET)

run_test:
	./$(TARGET)
	
#demo:
#	./$(TARGET2)

clean_test:
	rm -rf $(TARGET) test.o

clean:
	rm -rf $(TARGET) *.out *~ *.log *.o *.dat data.txt #$(TARGET2)
