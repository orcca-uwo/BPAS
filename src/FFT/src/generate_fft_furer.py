import math

def generate_roots():
	roots = open("roots","r")
	c_roots = open("ROOTS.c","w")
	k=0
	for line in roots:
		c_roots.write(("#define W%i "%k)+line)
		k=k+1
	return

def logof(n):
	l = 0
	while n>1:
		l=l+1
		n>>=1
	return l

def RevBitIncr(a,b,pow):
	l = logof(pow)+1
	u = logof(b)
	v = u%l
	c = b-(b>>v)
	b=b>>v
	a=a+b
	if((a&c) != 0 or b==1):
		return a
	c = b - (b>>l)
	a = a-(b<<v)
	while True:
		b = b>>l
		a = a+b
		if( (a&c) != 0 or b ==1):
			break
		c = c>>l
		a = a-(b<<l)
	return a

def RevBidMap(THRESHOLD,pow):
	RevBitMap = [k for k in range(0,THRESHOLD)]
	if THRESHOLD==1:
		return RevBitMap
	i = 0
	j = 0
	while i != THRESHOLD:
		RevBitMap[i] = j
		i = i+1
		j = RevBitIncr(j,THRESHOLD,pow)
	return RevBitMap

def generate_Array_Bit_Reversal(code,header,Array,name=""):
	s = "void ArrayBitReversalSpe" + name +"(sfixn *A)"
	header.write(s+";\n")
	s +="{\n"
	s += "\tsfixn t,u;\n"
	code.write(s)
	n = len(Array)
	explored = [1 for k in range(0,n)]
	for i in range(0,n):
		if(explored[i] == 1 and Array[i] != i):
			j = Array[i]
			explored[i] = 0
			s=""
			s+="\tt = A[%i];\n"%(i)
			while(explored[j] == 1):
				s+="\tu = A[%i];\n"%(j)
				s+="\tA[%i] = t;\n"%(j)
				s+="\tt = u;\n"
				explored[j] = 0
				j = Array[j]
			s+="\tA[%i] = t;\n"%j
			code.write(s);
	s = "}\n"
	code.write(s)

def generate_reversal_files(H,pow):
	reversalheader = open("../../../include/FFT/src/arraybitreversalfurer.h","w")
	reversalcode = open("arraybitreversalfurer.cpp","w")
	reversalheader.write("#include \"modpn.h\"\n");
	reversalheader.write("namespace PBPAS{\n")
	reversalcode.write("#include \"../../../include/FFT/src/arraybitreversalfurer.h\"\n")
	reversalcode.write("namespace PBPAS{\n")
	l = logof(pow)+1
	generate_Array_Bit_Reversal(reversalcode,reversalheader,RevBidMap(H,1))
	H=H>>1
	while H>=8:
		generate_Array_Bit_Reversal(reversalcode,reversalheader,RevBidMap(H,1),"%i"%H)
		H=H>>1
	reversalheader.write("}")
	reversalcode.write("}")
	reversalheader.close()
	return

def generate_16_point_fft(H,code):
	s=""
	s+="\tfor(int k=0;k<%i;k+=16){\n"%H
	s+="\t\t//FFT_16POINT(A+k,Wp);\n"
	s+="\t\tsfixn u = A[k+0];\n"
	s+="\t\tsfixn t = A[k+1];\n"
	s+="\t\tA[k+0] = AddModSpe(u,t);\n"
	s+="\t\tA[k+1] = SubModSpe(u,t);\n"
	s+="\t\tu = A[k+2];\n"
	s+="\t\tt = A[k+3];\n"
	s+="\t\tA[k+2] = AddModSpe(u,t);\n"
	s+="\t\tA[k+3] = SubModSpe(u,t);\n"
	s+="\t\tu = A[k+4];\n"
	s+="\t\tt = A[k+5];\n"
	s+="\t\tA[k+4] = AddModSpe(u,t);\n"
	s+="\t\tA[k+5] = SubModSpe(u,t);\n"
	s+="\t\tu = A[k+6];\n"
	s+="\t\tt = A[k+7];\n"
	s+="\t\tA[k+6] = AddModSpe(u,t);\n"
	s+="\t\tA[k+7] = SubModSpe(u,t);\n"
	s+="\t\tu = A[k+8];\n"
	s+="\t\tt = A[k+9];\n"
	s+="\t\tA[k+8] = AddModSpe(u,t);\n"
	s+="\t\tA[k+9] = SubModSpe(u,t);\n"
	s+="\t\tu = A[k+10];\n"
	s+="\t\tt = A[k+11];\n"
	s+="\t\tA[k+10] = AddModSpe(u,t);\n"
	s+="\t\tA[k+11] = SubModSpe(u,t);\n"
	s+="\t\tu = A[k+12];\n"
	s+="\t\tt = A[k+13];\n"
	s+="\t\tA[k+12] = AddModSpe(u,t);\n"
	s+="\t\tA[k+13] = SubModSpe(u,t);\n"
	s+="\t\tu = A[k+14];\n"
	s+="\t\tt = A[k+15];\n"
	s+="\t\tA[k+14] = AddModSpe(u,t);\n"
	s+="\t\tA[k+15] = SubModSpe(u,t);\n"
	s+="\t\tA[k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+3],*(Wp-3));\n"
	s+="\t\tA[k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+7],*(Wp-3));\n"
	s+="\t\tA[k+11] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+11],*(Wp-3));\n"
	s+="\t\tA[k+15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+15],*(Wp-3));\n"
	s+="\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k,A+k+2);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+4,A+k+6);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+8,A+k+10);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+12,A+k+14);\n"
	s+="\t\tA[k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+5],*(Wp-11));\n"
	s+="\t\tA[k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+6],*(Wp-10));\n"
	s+="\t\tA[k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+7],*(Wp-9));\n"
	s+="\t\tA[k+13] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+13],*(Wp-11));\n"
	s+="\t\tA[k+14] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+14],*(Wp-10));\n"
	s+="\t\tA[k+15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+15],*(Wp-9));\n"
	s+="\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k,A+k+4);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+2,A+k+6);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+8,A+k+12);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+10,A+k+14);\n"
	s+="\n"
	s+="\t\tA[k+9] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+9],*(Wp-27));\n"
	s+="\t\tA[k+10] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+10],*(Wp-26));\n"
	s+="\t\tA[k+11] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+11],*(Wp-25));\n"
	s+="\t\tA[k+12] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+12],*(Wp-24));\n"
	s+="\t\tA[k+13] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+13],*(Wp-23));\n"
	s+="\t\tA[k+14] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+14],*(Wp-22));\n"
	s+="\t\tA[k+15] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+15],*(Wp-21));\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k,A+k+8);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+2,A+k+10);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+4,A+k+12);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+6,A+k+14);\n"
	s+="\t}\n"
	code.write(s)


#def writeVerySmallFFTR(code,power,num):
#	log = logof(power)
#	code.write("inline void very_small_butterflies%i(sfixn* A, int next, int k){\n"%(num))
#	if num==0:
#		num=log
#	root = 2**(log-2)
#	if 2**log<8:
#		for k in range(0,2**(num-1)):
#			code.write("\t\tfor(int j=%i*next+k;j<%i*next+k;j+=%i){\n"%(2*k,2*k+1,2**log))
#			code.write("\t\t\tsfixn t,u;\n")
#			for j in range(0,2**log):
#				code.write("\t\t\tu = AddModSpe(*(A+j+%i),*(A+j+%i+next));\n"%(j,j))
#				code.write("\t\t\tt = SubModSpe(*(A+j+%i),*(A+j+%i+next));\n"%(j,j))
#				code.write("\t\t\t*(A+j+%i) = u;\n\t\t\t*(A+j+%i+next) = t;\n"%(j,j))
#			code.write("\t\t}\n")
#		for i in range(2,num+1):
#			for k in range(0,(2**num)/(2**i)):
#				for l in range(2**(i-1),2**(i)):
#					t = l-2**(i-1)
#					code.write("\t\tfor(int j=%i*next+k;j<%i*next+k;j+=%i){\n"%(k*2**i+t,k*2**i+t+1,2**log))
#					code.write("\t\t\tsfixn t,u;\n")
#					if(l-2**(i-1) != 0):
#						for j in range(0,2**log):
#							code.write("\t\t\t*(A+j+%i+%i*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+%i+%i*next),MY_ROOT%i);\n"%(j,2**(i-1),j,2**(i-1),(l-2**(i-1))*root))
#					for j in range(0,2**log):
#						code.write("\t\t\tu = AddModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
#						code.write("\t\t\tt = SubModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
#						code.write("\t\t\t*(A+j+%i) = u;\n\t\t\t*(A+j+%i+%i*next) = t;\n"%(j,j,2**(i-1)))
#					code.write("\t\t}\n")
#			root = root>>1;
#		code.write("}\n")
#		return
#
#	for k in range(0,2**(num-1)):
#		code.write("\t\tfor(int j=%i*next+k;j<%i*next+k;j+=8){\n"%(2*k,2*k+1))
#		code.write("\t\t\tsfixn t,u;\n")
#		for j in range(0,8):
#			code.write("\t\t\tu = AddModSpe(*(A+j+%i),*(A+j+%i+next));\n"%(j,j))
#			code.write("\t\t\tt = SubModSpe(*(A+j+%i),*(A+j+%i+next));\n"%(j,j))
#			code.write("\t\t\t*(A+j+%i) = u;\n\t\t\t*(A+j+%i+next) = t;\n"%(j,j))
#		code.write("\t\t}\n")
#	for i in range(2,num+1):
#		for k in range(0,(2**num)/(2**i)):
#			for l in range(2**(i-1),2**(i)):
#				t = l-2**(i-1)
#				code.write("\t\tfor(int j=%i*next+k;j<%i*next+k;j+=8){\n"%(k*2**i+t,k*2**i+t+1))
#				code.write("\t\t\tsfixn t,u;\n")
#				if(l-2**(i-1) != 0):
#					for j in range(0,8):
#						code.write("\t\t\t*(A+j+%i+%i*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+%i+%i*next),MY_ROOT%i);\n"%(j,2**(i-1),j,2**(i-1),(l-2**(i-1))*root))
#				for j in range(0,8):
#					code.write("\t\t\tu = AddModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
#					code.write("\t\t\tt = SubModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
#					code.write("\t\t\t*(A+j+%i) = u;\n\t\t\t*(A+j+%i+%i*next) = t;\n"%(j,j,2**(i-1)))
#				code.write("\t\t}\n")
#		root = root>>1;
#	code.write("}\n")

def writeSmallFFTR(code,power,num):
	log = logof(power)
	code.write("inline void small_butterflies%i(sfixn* A, int next, int k){\n"%(num))
	if num==0:
		num=log
	root = 2**(log-2)


#--------------------IF RADIX LESS THAN 8, unrolling differ---------------
	if 2**log<8:
		for i in range(2,num+1):
			for k in range(0,(2**num)/(2**i)):
				for l in range(2**(i-1),2**(i)):
					t = l-2**(i-1)
					code.write("\t\tfor(int j=%i*next+k;j<%i*next+k;j+=%i){\n"%(k*2**i+t,k*2**i+t+1,2**log))
					code.write("\t\t\tsfixn t,u;\n")
					if(l-2**(i-1) != 0):
						for j in range(0,2**log):
							code.write("\t\t\t*(A+j+%i+%i*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+%i+%i*next),MY_ROOT%i);\n"%(j,2**(i-1),j,2**(i-1),(l-2**(i-1))*root))
					for j in range(0,2**log):
						code.write("\t\t\tu = AddModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
						code.write("\t\t\tt = SubModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
						code.write("\t\t\t*(A+j+%i) = u;\n\t\t\t*(A+j+%i+%i*next) = t;\n"%(j,j,2**(i-1)))
					code.write("\t\t}\n")
			root = root>>1;
		code.write("}\n")
		return

#----------------------OTHERWISE CLASSICAL unrolling---------------------
	for i in range(2,num+1):
		for k in range(0,int((2**num)/(2**i))):
			for l in range(2**(i-1),2**(i)):
				t = l-2**(i-1)
				code.write("\tfor(int j=%i*next+k;j<%i*next+k;j+=8){\n"%(k*2**i+t,k*2**i+t+1))
				code.write("\t\tsfixn t,u;\n")
				if(l-2**(i-1) != 0):
					for j in range(0,8):
						code.write("\t\t*(A+j+%i+%i*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+%i+%i*next),MY_ROOT%i);\n"%(j,2**(i-1),j,2**(i-1),(l-2**(i-1))*root))
				for j in range(0,8):
					code.write("\t\tu = AddModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
					code.write("\t\tt = SubModSpe(*(A+j+%i),*(A+j+%i+%i*next));\n"%(j,j,2**(i-1)))
					code.write("\t\t*(A+j+%i) = u;\n\t\t*(A+j+%i+%i*next) = t;\n"%(j,j,2**(i-1)))
				code.write("\t}\n")
		root = root>>1;
	#h = 8
	#code.write("\tsfixn t,u,j;\n")
	#code.write("\tfor(int i=k;i<next+k;i+=%i){\n"%h)
	#for i in range(2,num+1):
	#	for k in range(0,(2**num)/(2**i)):
	#		for l in range(2**(i-1),2**(i)):
	#			t = l-2**(i-1)
	#			#code.write("\tfor(int j=%i*next+k;j<%i*next+k;j+=8){\n"%(k*2**i+t,k*2**i+t+1))
	#			code.write("\n\n\t\tj=%i*next;\n"%(k*2**i+t))
	#			if(l-2**(i-1) != 0):
	#				for j in range(0,h):
	#					code.write("\t\t*(A+j+i+%i+%i*next) = MontMulModSpe_OPT3_AS_GENE_INLINE(*(A+j+i+%i+%i*next),MY_ROOT%i);\n"%(j,2**(i-1),j,2**(i-1),(l-2**(i-1))*root))
	#			for j in range(0,h):
	#				code.write("\t\tu = AddModSpe(*(A+j+i+%i),*(A+j+i+%i+%i*next));\n"%(j,j,2**(i-1)))
	#				code.write("\t\tt = SubModSpe(*(A+j+i+%i),*(A+j+i+%i+%i*next));\n"%(j,j,2**(i-1)))
	#				code.write("\t\t*(A+j+i+%i) = u;\n\t\t*(A+j+i+%i+%i*next) = t;\n"%(j,j,2**(i-1)))
	#			#code.write("\t}\n")
	#	root = root>>1;
	#code.write("\t}\n")



	code.write("}\n")


def generate_small_fft(code,num,log):
	code.write("inline void small_fft%i(sfixn* u0"%(num%log))
	if num==0:
		num=log
	for k in range(1,2**num):
		code.write(", sfixn* u%i"%k)
	Array = RevBidMap(2**num,1)
	code.write("){\n")
	code.write("\tsfixn t,u;\n")
	root = 2**(log-2)
	for k in range(0,2**(num-1)):
		code.write("\tu = AddModSpe(*u%i,*u%i);\n"%(2*k,2*k+1))
		code.write("\tt = SubModSpe(*u%i,*u%i);\n"%(2*k,2*k+1))
		code.write("\t*u%i = u;\n\t*u%i = t;\n"%(2*k,2*k+1))
	for i in range(2,num+1):
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**num)/(2**i))):
				code.write("\t*u%i = MontMulModSpe_OPT3_AS_GENE_INLINE(*u%i,MY_ROOT%i);\n"%(((2**i)*k+l),(k*2**i+l),(l-2**(i-1))*root))
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**num)/(2**i))):
				code.write("\tu = AddModSpe(*u%i,*u%i);\n"%(k*2**i+l,k*2**i+l+2**(i-1)))
				code.write("\tt = SubModSpe(*u%i,*u%i);\n"%(k*2**i+l,k*2**i+l+2**(i-1)))
				code.write("\t*u%i = u;\n\t*u%i = t;\n"%(k*2**i+l,k*2**i+l+2**(i-1)))
		root = root>>1;
	code.write("}\n")



def generate_dft_iterative(size,pow,code,name=""):
	code.write("void DFT_iter"+ name+"(sfixn *A,sfixn *W){\n")
	decal = 0
	add = 2*pow
	l = logof(pow)+1
	while add<(size/(2*pow)):
		add = add<<l
		decal = decal + add
	decal = decal + size
	code.write("\tsfixn *Wp = W + %i;\n"%(decal))
	code.write("\tsfixn* Wt;\n")
	if size>=2*pow:
		code.write("\tfor(int k=0;k<%i;k+=MY_POW){\n"%size)
		code.write("\tsmall_fft0((A+k)")
		for j in range(1,2*pow):
			code.write(",(A+%i+k)"%(j))
		code.write(");\n")
		code.write("\t}\n")
	else:
		l = logof(size)
		code.write("\t\tsmall_fft%i((A+%i)"%(l,k*2*pow))
		for j in range(1,size):
			code.write(",(A+%i)"%j)
		code.write(");\n")
		return
	l = logof(pow)+1
	m=size>>l
	roots = (2*pow)
	Array = RevBidMap(2*pow,1)
	tab ="{%i"%Array[0]
	for k in range(1,2*pow):
		tab = tab+",%i"%Array[k]
	tab = tab+"}"
	code.write(("\tunsigned long my_array[%i]="%(2*pow))+tab+";\n")

	while m>=2*pow:
		jump = size/m
		roots = roots * 2*pow
		code.write("\tWp = Wp- %i;\n"%roots)
		code.write("\tfor(int k=0;k<%i;k+=%i){\n"%(size,jump*2*pow))
		code.write("\t\tunsigned long pos = 1;\n")



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#--------------------MULTIPLICATION BY TWIDDLE FACTORS-------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


		if jump<8:
			code.write("\t\t\tunsigned long add = my_array[pos];\n")
			code.write("\t\t\tunsigned long start = 0;\n")
			for t in range(0,jump):
				code.write("\t\t\t\tA[%i+%i+k] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+%i+k],*(Wp+start));\n"%(t,jump,t,jump))
				code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tsfixn t,u;\n")
			for t in range(0,jump):
				code.write("\t\t\t\tu=AddModSpe(A[%i+k],A[%i+k+%i]);\n"%(t,t,jump))
				code.write("\t\t\t\tt=SubModSpe(A[%i+k],A[%i+k+%i]);\n"%(t,t,jump))
				code.write("\t\t\t\tA[%i+k] = u;\n"%t)
				code.write("\t\t\t\tA[%i+k+%i] = t;\n"%(t,jump))
		else :
			code.write("\t\t\tunsigned long add = my_array[pos];\n")
			code.write("\t\t\tunsigned long start = add;\n")
			code.write("\t\t\t\tA[%i+k+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+k+1],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+k+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+k+2],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+k+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+k+3],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+k+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+k+4],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+k+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+k+5],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+k+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+k+6],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+k+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+k+7],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")

#--------------------ADDSUB----------------------


			code.write("\t\t\t\tsfixn t,u;\n")
			for k in range(0,8):
				code.write("\t\t\t\tu=AddModSpe(A[%i+k],A[%i+%i+k]);\n"%(k,k,jump))
				code.write("\t\t\t\tt=SubModSpe(A[%i+k],A[%i+%i+k]);\n"%(k,k,jump))
				code.write("\t\t\t\tA[%i+k]=u;\n"%k)
				code.write("\t\t\t\tA[%i+%i+k]=t;\n"%(k,jump))
#--------------------MAIN LOOP------------------------


			code.write("\t\t\tfor(int i=8+k;i<%i+k;i+=8){\n"%(jump))
			code.write("\t\t\t\tA[%i+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+1],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+2],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+3],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+4],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+5],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+6],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")
			code.write("\t\t\t\tA[%i+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+7],*(Wp+start));\n"%(jump,jump))
			code.write("\t\t\t\tstart+=add;\n")

#--------------------ADDSUB------------------------

			for k in range(0,8):
				code.write("\t\t\t\tu=AddModSpe(A[%i+i],A[%i+i+%i]);\n"%(k,k,jump))
				code.write("\t\t\t\tt=SubModSpe(A[%i+i],A[%i+i+%i]);\n"%(k,k,jump))
				code.write("\t\t\t\tA[%i+i]=u;\n"%k)
				code.write("\t\t\t\tA[%i+i+%i]=t;\n"%(k,jump))

			code.write("\t\t\t}\n")
		code.write("\t\tpos = 2;\n")
		code.write("\t\tfor(int j=k+2*%i;j<k+%i;j+=2*%i,pos+=2){\n"%(jump,jump*2*pow,jump))
		if jump<8:
			code.write("\t\t\tunsigned long add1 = my_array[pos];\n")
			code.write("\t\t\tunsigned long start1 = 0;\n")
			for t in range(0,jump):
				code.write("\t\t\t\tA[%i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+j],*(Wp+start1));\n"%(t,t))
				code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\tunsigned long add2 = my_array[pos+1];\n")
			code.write("\t\t\tunsigned long start2 = 0;\n")
			for t in range(0,jump):
				code.write("\t\t\t\tA[%i+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+%i+j],*(Wp+start2));\n"%(t,jump,t,jump))
				code.write("\t\t\t\tstart2+=add2;\n")

#----------------ADDSUB-----------------------------
			for t in range(0,jump):
				code.write("\t\t\t\tu=AddModSpe(A[%i+j],A[%i+j+%i]);\n"%(t,t,jump))
				code.write("\t\t\t\tt=SubModSpe(A[%i+j],A[%i+j+%i]);\n"%(t,t,jump))
				code.write("\t\t\t\tA[%i+j] = u;\n"%t)
				code.write("\t\t\t\tA[%i+j+%i] = t;\n"%(t,jump))


		else :
			code.write("\t\t\tunsigned long add1 = my_array[pos];\n")
			code.write("\t\t\tunsigned long start1 = add1;\n")
			code.write("\t\t\t\tA[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")



			code.write("\t\t\tunsigned long add2 = my_array[pos+1];\n")
			code.write("\t\t\tunsigned long start2 = add2;\n")
			code.write("\t\t\t\tA[1+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[2+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[3+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[4+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[5+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[6+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[7+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")


#-----------------------ADDSUB-------------------

			for k in range(0,8):
				code.write("\t\t\t\tu=AddModSpe(A[%i+j],A[%i+j+%i]);\n"%(k,k,jump))
				code.write("\t\t\t\tt=SubModSpe(A[%i+j],A[%i+j+%i]);\n"%(k,k,jump))
				code.write("\t\t\t\tA[%i+j]=u;\n"%k)
				code.write("\t\t\t\tA[%i+j+%i]=t;\n"%(k,jump))


			code.write("\t\t\tfor(int i=8;i<%i;i+=8){\n"%(jump))
			code.write("\t\t\t\tA[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")
			code.write("\t\t\t\tA[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));\n")
			code.write("\t\t\t\tstart1+=add1;\n")

#---------------------------------------------

			code.write("\t\t\t\tA[i+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[i+1+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[i+2+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[i+3+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[i+4+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[i+5+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[i+6+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")
			code.write("\t\t\t\tA[i+7+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\t\tstart2+=add2;\n")

#-----------------------ADDSUB-------------------

			for k in range(0,8):
				code.write("\t\t\t\tu=AddModSpe(A[%i+i+j],A[%i+j+i+%i]);\n"%(k,k,jump))
				code.write("\t\t\t\tt=SubModSpe(A[%i+i+j],A[%i+j+i+%i]);\n"%(k,k,jump))
				code.write("\t\t\t\tA[%i+j+i]=u;\n"%k)
				code.write("\t\t\t\tA[%i+j+i+%i]=t;\n"%(k,jump))
			code.write("\t\t\t}\n")


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

		code.write("\t\t}\n")
		code.write("\t\tunsigned long next = %i;\n"%jump)
		code.write("\t\tsmall_butterflies0(A,next,k);\n")
		code.write("\t}\n")
		m=m>>l



#-------------LAST LEVEL IN THE BUTTERFLY----------------------------------
	if m<2*pow and m>1:
		jump = size/m
		roots = roots*m
		code.write("\tWp = Wp- %i;\n"%roots)
		code.write("\tunsigned long pos = 1;\n")



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#--------------------MULTIPLICATION BY TWIDDLE FACTORS-------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

		if jump<8:
			code.write("\tunsigned long add = my_array[pos]>>%i;\n"%(l-logof(m)))
			code.write("\tunsigned long start = 0;\n")
			for t in range(0,jump):
				code.write("\tA[%i+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+%i],*(Wp+start));\n"%(t,jump,t,jump))
				code.write("\tstart+=add;\n")
			code.write("\tsfixn t,u;\n")
			for t in range(0,jump):
				code.write("\tu=AddModSpe(A[%i+j],A[%i+j+%i]);\n"%(t,t,jump))
				code.write("\tt=SubModSpe(A[%i+j],A[%i+j+%i]);\n"%(t,t,jump))
				code.write("\tA[%i+j] = u;\n"%t)
				code.write("\tA[%i+j+%i] = t;\n"%(t,jump))
		else :
			code.write("\tunsigned long add = my_array[pos]>>%i;\n"%(l-logof(m)))
			code.write("\tunsigned long start = add;\n")
			code.write("\tA[%i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+1],*(Wp+start));\n"%(jump,jump))
			code.write("\tstart+=add;\n")
			code.write("\tA[%i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+2],*(Wp+start));\n"%(jump,jump))
			code.write("\tstart+=add;\n")
			code.write("\tA[%i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+3],*(Wp+start));\n"%(jump,jump))
			code.write("\tstart+=add;\n")
			code.write("\tA[%i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+4],*(Wp+start));\n"%(jump,jump))
			code.write("\tstart+=add;\n")
			code.write("\tA[%i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+5],*(Wp+start));\n"%(jump,jump))
			code.write("\tstart+=add;\n")
			code.write("\tA[%i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+6],*(Wp+start));\n"%(jump,jump))
			code.write("\tstart+=add;\n")
			code.write("\tA[%i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+7],*(Wp+start));\n"%(jump,jump))
			code.write("\tstart+=add;\n")

#--------------------ADDSUB----------------------


			code.write("\tsfixn t,u;\n")
			for k in range(0,8):
				code.write("\tu=AddModSpe(A[%i],A[%i+%i]);\n"%(k,k,jump))
				code.write("\tt=SubModSpe(A[%i],A[%i+%i]);\n"%(k,k,jump))
				code.write("\tA[%i]=u;\n"%k)
				code.write("\tA[%i+%i]=t;\n"%(k,jump))
#--------------------MAIN LOOP------------------------


			code.write("\tfor(int i=8;i<%i;i+=8){\n"%(jump))
			code.write("\t\tA[%i+i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")
			code.write("\t\tA[%i+i+1] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+1],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")
			code.write("\t\tA[%i+i+2] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+2],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")
			code.write("\t\tA[%i+i+3] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+3],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")
			code.write("\t\tA[%i+i+4] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+4],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")
			code.write("\t\tA[%i+i+5] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+5],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")
			code.write("\t\tA[%i+i+6] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+6],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")
			code.write("\t\tA[%i+i+7] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+i+7],*(Wp+start));\n"%(jump,jump))
			code.write("\t\tstart+=add;\n")

#--------------------ADDSUB------------------------

			for k in range(0,8):
				code.write("\t\tu=AddModSpe(A[%i+i],A[%i+i+%i]);\n"%(k,k,jump))
				code.write("\t\tt=SubModSpe(A[%i+i],A[%i+i+%i]);\n"%(k,k,jump))
				code.write("\t\tA[%i+i]=u;\n"%k)
				code.write("\t\tA[%i+i+%i]=t;\n"%(k,jump))

			code.write("\t}\n")
		code.write("\tpos = 2;\n")
		code.write("\tfor(int j=2*%i;j<%i;j+=2*%i,pos+=2){\n"%(jump,size,jump))
		if jump<8:
			code.write("\t\tunsigned long add1 = my_array[pos]>>%i;\n"%(l-logof(m)))
			code.write("\t\tunsigned long start1 = 0;\n")
			for t in range(0,jump):
				code.write("\t\tA[%i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+j],*(Wp+start1));\n"%(t,t))
				code.write("\t\tstart1+=add1;\n")
			code.write("\t\tunsigned long add2 = my_array[pos+1]>>%i;\n"%(l-logof(m)))
			code.write("\t\tunsigned long start2 = 0;\n")
			for t in range(0,jump):
				code.write("\t\tA[%i+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i+%i+j],*(Wp+start2));\n"%(t,jump,t,jump))
				code.write("\t\tstart2+=add2;\n")

#----------------ADDSUB-----------------------------
			for t in range(0,jump):
				code.write("\t\tu=AddModSpe(A[%i+j],A[%i+j+%i]);\n"%(t,t,jump))
				code.write("\t\tt=SubModSpe(A[%i+j],A[%i+j+%i]);\n"%(t,t,jump))
				code.write("\t\tA[%i+j] = u;\n"%t)
				code.write("\t\tA[%i+j+%i] = t;\n"%(t,jump))


		else :
			code.write("\t\tunsigned long add1 = my_array[pos]>>%i;\n"%(l-logof(m)))
			code.write("\t\tunsigned long start1 = add1;\n")
			code.write("\t\tA[1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j],*(Wp+start1));\n")
			code.write("\t\tstart1+=add1;\n")
			code.write("\t\tA[2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j],*(Wp+start1));\n")
			code.write("\t\tstart1+=add1;\n")
			code.write("\t\tA[3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j],*(Wp+start1));\n")
			code.write("\t\tstart1+=add1;\n")
			code.write("\t\tA[4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j],*(Wp+start1));\n")
			code.write("\t\tstart1+=add1;\n")
			code.write("\t\tA[5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j],*(Wp+start1));\n")
			code.write("\t\tstart1+=add1;\n")
			code.write("\t\tA[6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j],*(Wp+start1));\n")
			code.write("\t\tstart1+=add1;\n")
			code.write("\t\tA[7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j],*(Wp+start1));\n")
			code.write("\t\tstart1+=add1;\n")



			code.write("\t\tunsigned long add2 = my_array[pos+1]>>%i;\n"%(l-logof(m)))
			code.write("\t\tunsigned long start2 = add2;\n")
			code.write("\t\tA[1+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[1+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\tstart2+=add2;\n")
			code.write("\t\tA[2+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[2+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\tstart2+=add2;\n")
			code.write("\t\tA[3+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[3+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\tstart2+=add2;\n")
			code.write("\t\tA[4+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[4+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\tstart2+=add2;\n")
			code.write("\t\tA[5+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[5+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\tstart2+=add2;\n")
			code.write("\t\tA[6+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[6+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\tstart2+=add2;\n")
			code.write("\t\tA[7+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[7+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\tstart2+=add2;\n")


#-----------------------ADDSUB-------------------

			for k in range(0,8):
				code.write("\t\tu=AddModSpe(A[%i+j],A[%i+j+%i]);\n"%(k,k,jump))
				code.write("\t\tt=SubModSpe(A[%i+j],A[%i+j+%i]);\n"%(k,k,jump))
				code.write("\t\tA[%i+j]=u;\n"%k)
				code.write("\t\tA[%i+j+%i]=t;\n"%(k,jump))


			code.write("\t\tfor(int i=8;i<%i;i+=8){\n"%(jump))
			code.write("\t\t\tA[i+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")
			code.write("\t\t\tA[i+1+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")
			code.write("\t\t\tA[i+2+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")
			code.write("\t\t\tA[i+3+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")
			code.write("\t\t\tA[i+4+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")
			code.write("\t\t\tA[i+5+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")
			code.write("\t\t\tA[i+6+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")
			code.write("\t\t\tA[i+7+j] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j],*(Wp+start1));\n")
			code.write("\t\t\tstart1+=add1;\n")

#---------------------------------------------

			code.write("\t\t\tA[i+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")
			code.write("\t\t\tA[i+1+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+1+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")
			code.write("\t\t\tA[i+2+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+2+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")
			code.write("\t\t\tA[i+3+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+3+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")
			code.write("\t\t\tA[i+4+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+4+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")
			code.write("\t\t\tA[i+5+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+5+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")
			code.write("\t\t\tA[i+6+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+6+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")
			code.write("\t\t\tA[i+7+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[i+7+j+%i],*(Wp+start2));\n"%(jump,jump))
			code.write("\t\t\tstart2+=add2;\n")

#-----------------------ADDSUB-------------------

			for k in range(0,8):
				code.write("\t\t\tu=AddModSpe(A[%i+i+j],A[%i+j+i+%i]);\n"%(k,k,jump))
				code.write("\t\t\tt=SubModSpe(A[%i+i+j],A[%i+j+i+%i]);\n"%(k,k,jump))
				code.write("\t\t\tA[%i+j+i]=u;\n"%k)
				code.write("\t\t\tA[%i+j+i+%i]=t;\n"%(k,jump))
			code.write("\t\t}\n")
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


		code.write("\t}\n")
		code.write("\tunsigned long next = %i;\n"%jump)
		if m==1:
			m = 2*pow
		code.write("\tsmall_butterflies%i(A,next,0);\n"%(logof(m)%l))
	code.write("}\n")

def extended_gcd(a,b):
	r = a
	u = 1
	v = 0
	r_ = b
	u_ = 0
	v_ = 1
	while r_!=0:
		q = r/r_
		(r,u,v,r_,u_,v_) = (r_,u_,v_,r-q*r_,u-q*u_,v-q*v_)
	return (r,u,v)

#-----------Writes the code for the original Montgomery trick-----------------
def writeRegularMontMul(code,num):
#	code.write("inline sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){\n")
	code.write("sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){\n")
	code.write("\tasm(\"mulq %2\\n\\t\"\n")
	code.write("\t\"movq %%rax,%%rsi\\n\\t\"\n")
	code.write("\t\"movq %%rdx,%%rdi\\n\\t\"\n")
	code.write("\t\"imulq %3,%%rax\\n\\t\"\n")
	code.write("\t\"mulq %4\\n\\t\"\n")
	code.write("\t\"add %%rsi,%%rax\\n\\t\"\n")
	code.write("\t\"adc %%rdi,%%rdx\\n\\t\"\n")
	code.write("\t\"subq %4,%%rdx\\n\\t\"\n")
	code.write("\t\"mov %%rdx,%%rax\\n\\t\"\n")
	code.write("\t\"sar $63,%%rax\\n\\t\"\n")
	code.write("\t\"andq %4,%%rax\\n\\t\"\n")
	code.write("\t\"addq %%rax,%%rdx\\n\\t\"\n")
	code.write("\t: \"=d\" (a)\n")
	code.write("\t: \"a\"(a),\"rm\"(b),\"b\"((sfixn) FURER_INV_PRIME%i),\"c\"((sfixn) FURER_MY_PRIME%i)\n"%(num,num))
	code.write("\t:\"rsi\",\"rdi\");\n")
	code.write("\treturn a;\n")
	code.write("}\n")

#----------Writes the code of the Montgomery Mult with two short multiplications-----------------
def writeSpecialMontMul(code,c,Npow,Rpow,num):
	lpow = 64-Rpow + Npow
#	code.write("inline sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){\n")
	code.write("sfixn MontMulModSpe_OPT3_AS_GENE_INLINE(sfixn a,sfixn b){\n")
	code.write("\tasm(\"mulq %2\\n\\t\"\n")
	code.write("\t\"movq %%rax,%%rsi\\n\\t\"\n")
	code.write("\t\"movq %%rdx,%%rdi\\n\\t\"\n")
	code.write("\t\"imulq %3,%%rax\\n\\t\"\n")

#---------------3RD MULTIPLICATION-------------------------------------
	code.write("\t\"add %%rax,%%rsi\\n\\t\"\n")
	code.write("\t\"adc $0,%%rdi\\n\\t\"\n")
	code.write("\t\"shr $%i"%(64-Rpow) +",%%rax\\n\\t\"\n")
	code.write("\t\"imulq $%i"%c +",%%rax\\n\\t\"\n")
	code.write("\t\"rolq $%i"%lpow+",%%rax\\n\\t\"\n")
	code.write("\t\"movq %%rax,%%rdx\\n\\t\"\n")
	code.write("\t\"andq %5,%%rax\\n\\t\"\n")
	code.write("\t\"andq %6,%%rdx\\n\\t\"\n")
#---------------------------------------------------------------------
	code.write("\t\"add %%rsi,%%rax\\n\\t\"\n")
	code.write("\t\"adc %%rdi,%%rdx\\n\\t\"\n")
	code.write("\t\"subq %4,%%rdx\\n\\t\"\n")
	code.write("\t\"mov %%rdx,%%rax\\n\\t\"\n")
	code.write("\t\"sar $63,%%rax\\n\\t\"\n")
	code.write("\t\"andq %4,%%rax\\n\\t\"\n")
	code.write("\t\"addq %%rax,%%rdx\\n\\t\"\n")
	code.write("\t: \"=d\" (a)\n")
	code.write("\t: \"a\"(a),\"d\"(b),\"b\"((sfixn) FURER_INV_PRIME%i),\"c\"((sfixn) FURER_MY_PRIME%i),\"rm\"(%iu),\"rm\"(%i)\n"%(num,num,2**64-2**lpow,2**lpow-1))
	code.write("\t:\"rsi\",\"rdi\");\n")
	code.write("\treturn a;\n")
	code.write("}\n")



def generate_everything(p,b,gen,power,H,num):
	file = "fft_furer%i"%num
	p_1 = p-1
	Npow = 0
	while (p_1%(1<<Npow)==0):
		Npow=Npow+1
	Npow = Npow - 1
	c = p_1>>Npow
	Rpow = int(math.log(p,2)+1)
	c_sft = c<<Npow
	(r,v,u)= extended_gcd(2**Rpow,p)
	y = u%(2**Rpow)
	x = ((u-y)/(2**Rpow))*p
	if y>0:
		y = y-2**Rpow
		x = x+p
	R = (2**Rpow) %p
	R_2 = (2**(2*Rpow))%p
	base_rpow = 64-Rpow
	root = pow(int(gen),int((p-1)/(2*power)),p)
	header = open("../../../include/FFT/src/"+file+".h","w")
	header.write("#include \"modpn.h\"\n")
	header.write("#ifndef FURER_FFTSPE%i\n"%num)
	header.write("#define FURER_FFTSPE%i\n"%num)

#-----------------Constants used in FFT-----------------------------------------

	header.write("#define FURER_MY_PRIME%i %iu\n"%(num,p))
	header.write("#define FURER_INV_PRIME%i %iu\n"%(num,-y))
	header.write("#define FURER_C_SFT%i %i\n"%(num,c_sft))
	header.write("#define FURER_SEE%i %i\n"%(num,c))
	header.write("#define FURER_RINV%i %iu\n"%(num,v%p))
	header.write("#define FURER_RSFT%i %i\n"%(num,64-Rpow))
	header.write("#define FURER_NPOW%i %i\n"%(num,Npow))
	header.write("namespace FURERPBPAS%i{\n"%num)
	code = open(file+".cpp","w")
	code.write("#include \"../../../include/FFT/src/"+file+".h\"\n")
	code.write("#include \"../../../include/FFT/src/arraybitreversal.h\"\n")
	code.write("#include \"../../../include/FFT/src/modpn.h\"\n")
	code.write("#include <iostream>\n")
	code.write("#include <string.h>\n")
	code.write("#define FFT_THRESHOLD %i\n"%(H))
	code.write("#define FFT_THRESHOLD_LOG %i\n"%(len(bin(H))-3))
	code.write("#define COMPLETE %iu\n"%(2**64-p))
	mask = H
	bound = 2**64
	while mask < bound:
		mask = (mask*power*2+H)
	mask = mask% 2**64
	masks=[(mask<<i)%2**64 for i in range(0,logof(power)+1)]

#---------------Construction of the masks to determine the first Shuffle------------------
	for i in range(0,logof(power)+1):
		code.write("#define MY_MASK%i %iu\n"%(i,masks[i]))
#---------------Roots of unity of small order are encoded as constants-------------------
	for i in range(0,2*power):
		code.write("#define MY_ROOT%i %iu\n"%(i,(((2**Rpow)*root**i)%p)<<(64-Rpow)))

	code.write("#define LOG_OF_POW %i\n"%(logof(power)+1))
	code.write("#define MY_POW %i\n"%(2*power))
	code.write("#define CONST_R2 %iu\n"%R_2)
	code.write("#define CONST_R %iu\n"%R)
	code.write("#define BASE_RPOW %i\n"%base_rpow)
	code.write("#define GENERATOR %i\n"%gen)
	template = open("generate_fft_template_furer.cpp","r")
	header.write("void Shuffle2(int n, sfixn* A,sfixn* B);\n")
	header.write("void Shuffle(int n, sfixn* A,sfixn* B);\n")
	header.write("void RootsTableFurer(int n, int r,sfixn *T);\n")
	header.write("sfixn testDFT(int n,int index,sfixn* A,sfixn *W);\n")


#--------------We replace some expressions met in the template-----------------------
	for line in template:
		if "void DFT_eff" in line:
			code.write(line.replace("DFT_eff","DFT_eff_p%i"%num))
			header.write(line.replace("{",";\n").replace("DFT_eff","DFT_eff_p%i"%num))
		elif "void InvDFT_eff" in line:
			code.write(line.replace("InvDFT_eff","InvDFT_eff_p%i"%num))
			header.write(line.replace("{",";\n").replace("InvDFT_eff","InvDFT_eff_p%i"%num))
		elif "void InvDFTKeepMont_eff" in line:
			code.write(line.replace("InvDFTKeepMont_eff","InvDFTKeepMont_eff_p%i"%num))
			header.write(line.replace("{",";\n").replace("InvDFTKeepMont_eff","InvDFTKeepMont_eff_p%i"%num))
		elif "DFT_rec" in line:
			code.write(line.replace("DFT_rec","DFT_rec_p%i"%num))
		elif "DFT_eff" in line:
			code.write(line.replace("DFT_eff","DFT_eff_p%i"%num))
		elif "namespace FURERPBPAS" in line:
			code.write(line.replace("FURERPBPAS","FURERPBPAS%i"%num))
		elif line == "DFT_ITERATIVE\n":
			generate_dft_iterative(H,power,code)
			tmpH = H>>1
#---------- We generate here all small iterative DFT ---------------------------
			while(tmpH>16):
				generate_dft_iterative(tmpH,power,code,"%i"%tmpH)
				tmpH=tmpH>>1
		elif "FIRST_LEVEL_SHUFFLE" in line:
			for i in range(0,logof(power)+1):
				code.write("\t\tif(m & MY_MASK%i){\n"%i)
				if i!=0:
					code.write("\t\t\tnext=m>>%i;\n"%i)
					code.write("\t\t\tadjust=%i;\n"%(logof(power)+1-i))
				else:
					code.write("\t\t\tnext=m>>%i;\n"%(logof(power)+1))
				code.write("\t\t}\n")
		elif "COMPUTE_POS_IN_SHUFFLES" in line:
			for i in range(1,logof(power)+2):
				code.write("\t\t\tif(m & MY_MASK%i){\n"%(i%(logof(power)+1)))
				code.write("\t\t\t\tm>>=%i;\n}\n"%i)
		elif "SMALL_FFTs" in line:
			for k in range(1,logof(power)+2):
				code.write("\t\tif(n& MY_MASK%i){\n"%(k%(logof(power)+1)))
				code.write("\t\tsmall_butterflies%i(A,next,0);\n"%(k%(logof(power)+1)))
				code.write("\t\t}\n")
		elif "SMALL_FFTR" in line:
			code.write("\t\tsmall_butterflies0(A,next,0);\n")
		elif "GENERATE_SMALL_FFT" in line:
			for k in range(0,logof(power)+1):
				generate_small_fft(code,k,logof(power)+1)
				writeSmallFFTR(code,2*power,k)
		elif "FFT_CASES" in line:
			tmpH = 32
			while tmpH<H:
				code.write("\t\t\tcase %i:\n"%tmpH)
				code.write("\t\t\t\tPBPAS::ArrayBitReversalSpe%i(A);\n"%tmpH)
				code.write("\t\t\t\tDFT_iter%i(A,W);\n"%tmpH)
				code.write("\t\t\t\tbreak;\n")
				tmpH=tmpH<<1

		elif "MONTGOMERYMUL" in line:
			if c*2**Rpow >=2**64:
				writeRegularMontMul(code,num)
			else:
				code.write("#if SPECIALMONT\n")
				writeSpecialMontMul(code,c,Npow,Rpow,num)
				code.write("#else\n")
				writeRegularMontMul(code,num)
				code.write("#endif\n")
		else:
			s = line
			if "MY_PRIME" in s:
				s=(s.replace("MY_PRIME","FURER_MY_PRIME%i"%num))
			if "INV_PRIME" in s:
			      s=(s.replace("INV_PRIME","FURER_INV_PRIME%i"%num))
			if "C_SFT" in s:
			      s=(s.replace("C_SFT","FURER_C_SFT%i"%num))
			if "SEE" in s:
			      s=(s.replace("SEE","FURER_SEE%i"%num))
			if "RSFT" in s:
			      s=(s.replace("RSFT","FURER_RSFT%i"%num))
			if "RINV" in s:
				s=(s.replace("RINV","FURER_RINV%i"%num))
			if "NPOW" in s:
				s=(s.replace("NPOW","FURER_NPOW%i"%num))
			if "ARRAY_FOR_MULS" in s :
				Array = RevBidMap(2*power,1)
				tab ="{%i"%Array[0]
				for k in range(1,2*power):
					tab = tab+",%i"%Array[k]
				tab = tab+"}"
				s=(s.replace("ARRAY_FOR_MULS",tab))
			code.write(s)
	header.write("}\n")
	header.write("#endif\n")
	header.close()
	template.close()
	code.close();

import sys
H=int(sys.argv[1])
config = open(".config_furer","r")
g = []
b = []
p = []
power = int(config.readline())
files = []
k=0
for line in config:
	if(k%2==0):
		p.append(int(line))
	else:
		g.append(int(line))
	k=k+1
k=0
for x in g:
	b.append(pow(int(x),int((p[k]-1)/(2*power)),p[k]))
	k=k+1
for i in range(0,len(b)):
	h = logof(H)
	generate_everything(p[i],b[i],g[i],power,H,i+1)
generate_reversal_files(H,power)
config.close()
