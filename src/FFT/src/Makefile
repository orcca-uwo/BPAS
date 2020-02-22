## BPAS/Main/src/FFT/src.64/Makefile
CILKPP = g++
COMPILEARG = -std=c++11 -c -O2 -fcilkplus -fpic -DLINUXINTEL64=1 -march=native
THRESHOLD=1024
ifeq ($(XEON5600),1)
	COMPILEARG+= -DXEON5600=1
else
	COMPILEARG+= -DAMD=1
endif
ifeq ($(ASM),1)
	COMPILEARG+= -DASM=1
endif
ifeq ($(DEBUG),1)
	COMPILEARG+= -DDEBUG=1
	COMPILEARG+= -DMODDEBUG=1
endif
ifdef HTHRESHOLD
	THRESHOLD=$(HTHRESHOLD)
endif
ifeq ($(COMPOSESHUFFLE),1)
	COMPILEARG+= -DCOMPOSESHUFFLE=1
endif
ifeq ($(FURERCOMPLEXITY),1)
	COMPILEARG+= -DFURER=1
endif
ifeq ($(SPEMONT),1)
	COMPILEARG+= -DSPECIALMONT=1
endif
## COMPILEARG = -g -O2 -fPIC -fcilkplus -DLINUXINTEL64=1

INCLUDE_DIR = ../../../include/FFT/src

vpath %.cpp  ./
vpath %.h $(INCLUDE_DIR)/ $(INCLUDE_DIR)/modpn_hfiles

#############################################################################################
###################### Build, Compile ### ###################################################
#############################################################################################

ifeq ($(wildcard fft_iter*),)
all: GENERATE_FFT
endif

ifeq ($(wildcard fft_furer*),)
all: GENERATE_FURER
endif

ifeq ($(wildcard tft_tree*),)
all: GENERATE_TFT_RELAX
endif

#ifeq ($(wildcard custom_DFT*),)
#all: GENERATE_DFT
#endif

#all : dft8.o dft16.o dft32.o dft64.o dft_general.o dft_fermat_fwd_16.o dft_fermat_inverse_16.o dft_general_fermat.o custom_DFT.o 
all: twocon_general.o twocon_basic.o fft_iter1.o fft_iter2.o fft_furer1.o fft_furer2.o arraybitreversal.o modpn_export.o ser_general_routine.o ser_basic_routine.o PrimeField.o general_routine.o basic_routine.o pbpas_basic.o tft_tree1.o tft_tree2.o transpose.o

%.o: %.cpp
	$(CILKPP) $(COMPILEARG) $<
#this part is for fft####################
fft_furer1.o : fft_furer1.cpp
	$(CILKPP)  $(COMPILEARG) fft_furer1.cpp
fft_furer2.o : fft_furer2.cpp
	$(CILKPP)  $(COMPILEARG) fft_furer2.cpp
fft_iter1.o : fft_iter1.cpp
	$(CILKPP)  $(COMPILEARG) fft_iter1.cpp
fft_iter2.o : fft_iter2.cpp
	$(CILKPP)  $(COMPILEARG) fft_iter2.cpp
arraybitreversal.o : arraybitreversal.cpp
	$(CILKPP) $(COMPILEARG) arraybitreversal.cpp

GENERATE_FFT:generate_fft.py
	python generate_fft.py $(THRESHOLD)
GENERATE_FURER:generate_fft_furer.py
	python generate_fft_furer.py $(THRESHOLD)

### for six step dft #######
#custom_DFT.o : custom_DFT.cpp
#	$(CILKPP) $(COMPILEARG) custom_DFT.cpp
GENERATE_DFT:generate_dft_wc.py
	python generate_dft_wc.py

###for tft##################
tft_tree1.o : tft_tree1.cpp
	$(CILKPP)  $(COMPILEARG) tft_tree1.cpp
tft_tree2.o : tft_tree2.cpp
	$(CILKPP)  $(COMPILEARG) tft_tree2.cpp
GENERATE_TFT_RELAX: generate_tft_relax.py
	python generate_tft_relax.py $(THRESHOLD)


serial: COMPILEARG += -DSERIAL=1
serial: GENERATE_FFT GENERATE_FURER twocon_general.o twocon_basic.o fft_iter1.o fft_iter2.o fft_furer1.o fft_furer2.o arraybitreversal.o modpn_export.o ser_general_routine.o ser_basic_routine.o PrimeField.o general_routine.o basic_routine.o pbpas_basic.o GENERATE_TFT_RELAX tft_tree1.o tft_tree2.o transpose.o

test:

clean:
	rm -rf *.o fft_furer* fft_iter* tft_tree* custom_DFT*
