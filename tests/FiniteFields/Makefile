
#######################################
SRC_LIST:=\
	gfpf_six_step_fft_test\
	small_prime_field_six_step_fft\

#######################################
INCLUDES = ../../include/FiniteFields
CC  = g++ -L$(LIB_PATH) #-pg -g3 -g
LIB_PATH:=../../src/FiniteFields
LIBS= -lm -lgmp -lcilkrts -lbpas
#######################################
COMPILER_OPTIONS=-fpic -fcilkplus -flto\
								 -DLINUXINTEL64\
								 -march=native -Ofast\
								 -I$(INCLUDES)\
								 -DVERBOSE=0\
								 -DSMALL_FFT_STEPS_TIMING_ENABLED=0\
								 -DVERIFICATION_ENABLED=0\
#######################################
#######################################
COMPILE_ARGS=$(CC) $(COMPILER_OPTIONS) $(LIBS) #-Wl,--start-group
#######################################
OBJS:=$(SRC_LIST:=.o)

#######################################
all : clean test

#######################################
check: parallel

#######################################
test: parallel
	./gfpf_six_step_fft_test.bin 32 2	
	./gfpf_six_step_fft_test.bin 32 3

test-validate: test ;

#######################################
parallel: clean
parallel: parallel_src 
parallel_src:
	(cd $(BPAS_HOME)/src/FiniteFields;\
		make clean;\
		cd ..;\
		make;\
		)
parallel: $(OBJS)
#######################################
#######################################
serial-gmp: clean
serial-gmp: parallel_src_serial_gmp 
parallel_src_serial_gmp:
	(cd $(BPAS_HOME)/src/FiniteFields;\
		make clean;\
		make serial-gmp;\
		cd ..;\
		make;\
		)
serial-gmp: $(OBJS)
#######################################
#######################################
%.o: %.c
	@echo gcc [$<] ...	
	$(COMPILE_ARGS) $< $(LIBS) -o $(@:%.o=%.bin) 

#######################################
SERIAL_OBJS:=$(SRC_LIST:=.serial)
#######################################
serial: clean
serial: serial_src
serial_src: 
	(cd $(BPAS_HOME)/src/FiniteFields;\
		make clean;\
		make serial;\
		cd ..;\
		make;\
		)

#######################################
serial: $(SERIAL_OBJS)
#######################################
%.serial: %.c
	@echo gcc serial [$<] ...
	@$(COMPILE_ARGS) -DENABLE_CILK=0 $< -c -o $@
	@cp $@ $(@:%.serial=%.o)
	@$(COMPILE_ARGS) -DENABLE_CILK=0 $< $(LIBS) -o $(@:%.serial=%.bin) 

#######################################
benchmark-serial: serial
	./benchmark.sh serial

#######################################
benchmark-serial-gmp: serial-gmp
	./benchmark.sh serial-gmp

#######################################
benchmark-parallel: parallel
	./benchmark.sh parallel

#######################################
RM_LIST= bin log o data out tmp serial
RM_STR=$(RM_LIST:%=*.%) callgrind.out.*

clean:
	echo $(RM_STR)
	rm -rf $(RM_STR)
