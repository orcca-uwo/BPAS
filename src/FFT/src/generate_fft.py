import math

def generate_roots():
	roots = open("roots","r")
	c_roots = open("ROOTS.c","w")
	k=0
	for line in roots:
		c_roots.write(("#define W%i "%k)+line)
		k=k+1
	return

def RevBitIncr(a,b):
	while True:
		b = b>>1
		a = a^b
		if( (a&b) != 0 or b ==1):
			break
	return a

def RevBidMap(THRESHOLD):
	RevBitMap = [k for k in range(0,THRESHOLD)]
	if THRESHOLD==1:
		return RevBitMap
	i = 0
	j = 0
	while i != THRESHOLD:
		if i<j:
			t = RevBitMap[i]
			RevBitMap[i] = j
			RevBitMap[j] = t
		i = i+1
		j = RevBitIncr(j,THRESHOLD)
	return RevBitMap

def generate_Array_Bit_Reversal(code,header,Array,name=""):
	s = "void ArrayBitReversalSpe" + name +"(sfixn *A)"
	header.write(s+";\n")
	s +="{\n"
	s += "\tsfixn t0,t1,t2,t3,t4,t5,t6,t7;\n"
	code.write(s)
	n = len(Array)
	for i in range(0,n,8):
		s=""
		for k in range(0,8):
			j = Array[i+k]
			if i+k<j:
				s+="\tt%i = A[%i];\n"%(k,i+k)
				s+="\tA[%i] = A[%i];\n"%(i+k,j)
				s+="\tA[%i] = t%i;\n"%(j,k)
		code.write(s);
	s = "}\n"
	code.write(s)
def generate_reversal_files(H):
	reversalheader = open("../../../include/FFT/src/arraybitreversal.h","w")
	reversalcode = open("arraybitreversal.cpp","w")
	reversalheader.write("#include \"modpn.h\"\n");
	reversalheader.write("namespace PBPAS{\n")
	reversalcode.write("#include \"../../../include/FFT/src/arraybitreversal.h\"\n")
	reversalcode.write("namespace PBPAS{\n")
	generate_Array_Bit_Reversal(reversalcode,reversalheader,RevBidMap(H))
	H=H>>1
	while H>=8:
		generate_Array_Bit_Reversal(reversalcode,reversalheader,RevBidMap(H),"%i"%H)
		H=H>>1
	reversalheader.write("}")
	reversalcode.write("}")
	reversalheader.close()
	return
def generate_first_block_fft(log,H,code):
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
	s+="\t\tA[k+3] = MontMulModSpe_OPT3_AS_GENE(A[k+3],*(Wp-3));\n"
	s+="\t\tA[k+7] = MontMulModSpe_OPT3_AS_GENE(A[k+7],*(Wp-3));\n"
	s+="\t\tA[k+11] = MontMulModSpe_OPT3_AS_GENE(A[k+11],*(Wp-3));\n"
	s+="\t\tA[k+15] = MontMulModSpe_OPT3_AS_GENE(A[k+15],*(Wp-3));\n"
	s+="\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k,A+k+2);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+4,A+k+6);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+8,A+k+10);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+12,A+k+14);\n"
	s+="\t\tA[k+5] = MontMulModSpe_OPT3_AS_GENE(A[k+5],*(Wp-11));\n"
	s+="\t\tA[k+6] = MontMulModSpe_OPT3_AS_GENE(A[k+6],*(Wp-10));\n"
	s+="\t\tA[k+7] = MontMulModSpe_OPT3_AS_GENE(A[k+7],*(Wp-9));\n"
	s+="\t\tA[k+13] = MontMulModSpe_OPT3_AS_GENE(A[k+13],*(Wp-11));\n"
	s+="\t\tA[k+14] = MontMulModSpe_OPT3_AS_GENE(A[k+14],*(Wp-10));\n"
	s+="\t\tA[k+15] = MontMulModSpe_OPT3_AS_GENE(A[k+15],*(Wp-9));\n"
	s+="\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k,A+k+4);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+2,A+k+6);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+8,A+k+12);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+10,A+k+14);\n"
	s+="\n"
	s+="\t\tA[k+9] = MontMulModSpe_OPT3_AS_GENE(A[k+9],*(Wp-27));\n"
	s+="\t\tA[k+10] = MontMulModSpe_OPT3_AS_GENE(A[k+10],*(Wp-26));\n"
	s+="\t\tA[k+11] = MontMulModSpe_OPT3_AS_GENE(A[k+11],*(Wp-25));\n"
	s+="\t\tA[k+12] = MontMulModSpe_OPT3_AS_GENE(A[k+12],*(Wp-24));\n"
	s+="\t\tA[k+13] = MontMulModSpe_OPT3_AS_GENE(A[k+13],*(Wp-23));\n"
	s+="\t\tA[k+14] = MontMulModSpe_OPT3_AS_GENE(A[k+14],*(Wp-22));\n"
	s+="\t\tA[k+15] = MontMulModSpe_OPT3_AS_GENE(A[k+15],*(Wp-21));\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k,A+k+8);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+2,A+k+10);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+4,A+k+12);\n"
	s+="\t\tAddSubSpeSSEModInplace(A+k+6,A+k+14);\n"
	s+="\t}\n"
	code.write(s)
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

def generate_last_block(log,H,depth,code):
	jump = (2**(depth-1))
	block_size = 2**log
	s="\tWp = Wp-(1L<<%i);\n"%depth
	s+="\tWt = Wp;\n"
	s+="\tWp = Wt;\n"
	#s+="asm(\"\":::\"memory\");\n"
	for i in range(1,log+1):
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				if ((l-2**(i-1))*jump!=0):
					s+="\tA[%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i],*(Wp+%i));\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\tA[%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[%i],*(Wp+%i));\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\tAddSubSpeSSEModInplace(A+%i,A+%i);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
		if i!=log:
			s+="\tWp = Wp-(1L<<%i);\n"%(depth+i)
	#s+="asm(\"\":::\"memory\");\n"
	s+="\tfor(int j=2;j<%i;j+=2){\n"%jump
	s+="\t\tWp = Wt;\n"
	for i in range(1,log+1):
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tA[j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+%i],*(Wp+j+%i));\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tA[j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[j+%i],*(Wp+j+%i));\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tAddSubSpeSSEModInplace(A+j+%i,A+j+%i);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
		if i!=log:
			s+="\t\tWp = Wp-(1L<<%i);\n"%(depth+i)
	s+="\t}\n"
	code.write(s)

def generate_block(log,H,depth,code):
#-------------Blocks of size 2**log
#-------------H : size of input vector
#-------------depth : number of the row where we are traversing datas
#--------------code : file where to write
	jump = (2**(depth-1))
	block_size = 2**log
	s="\tWp = Wp-(1L<<%i);\n"%depth
	s+="\tWt = Wp;\n"
	s+="\tfor(int k=0;k<%i;k+=%i){\n"%(H,jump<<log)
	s+="\t\tWp = Wt;\n"
	for i in range(1,log+1):
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				if ((l-2**(i-1))*jump!=0):
					s+="\t\tA[k+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+%i],*(Wp+%i));\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tA[k+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+%i],*(Wp+%i));\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tAddSubSpeSSEModInplace(A+k+%i,A+k+%i);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
		if i!=log:
			s+="\t\tWp = Wp-(1L<<%i);\n"%(depth+i)
	s+="\t\tfor(int j=2;j<%i;j+=2){\n"%jump
	s+="\t\t\tWp = Wt;\n"
	for i in range(1,log+1):
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\t\tA[k+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+%i],*(Wp+j+%i));\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\t\tA[k+j+%i] = MontMulModSpe_OPT3_AS_GENE_INLINE(A[k+j+%i],*(Wp+j+%i));\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\t\tAddSubSpeSSEModInplace(A+k+j+%i,A+k+j+%i);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
		if i!=log:
			s+="\t\t\tWp = Wp-(1L<<%i);\n"%(depth+i)
	s+="\t\t}\n"
	s+="\t}\n"
	code.write(s)

def generate_dft_iterative(size,code,name=""):
	code.write("void DFT_iter"+ name+"(sfixn *A,sfixn *W){\n")
	code.write("\tsfixn *Wp = W + (%i<<1)-4;\n"%size)
	code.write("\tsfixn* Wt;\n")
	generate_16_point_fft(size,code)
	code.write("\tWp = Wp - 28;\n")
	pos = 16
	depth = 5
	block = 3
	while(pos*(2**block)<size):
		generate_block(block,size,depth,code)
		decal = 0
		for k in range(1,block):
			decal = decal + 2**(depth+k)
		code.write("\tWp = Wt -%i;\n"%(decal))
		pos = pos*2**block
		depth = depth+block
	rem = size/pos
	log = 0
	while rem>1:
		log += 1
		rem = rem//2
	generate_last_block(log,size,depth,code)
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

def generate_everything(file,p,H,num):
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

	header = open("../../../include/FFT/src/"+file+".h","w")
	header.write("#include \"modpn.h\"\n")
	header.write("#ifndef FFTSPE%i\n"%num)
	header.write("#define FFTSPE%i\n"%num)
	header.write("#define MY_PRIME%i %i\n"%(num,p))
	header.write("#define INV_PRIME%i %i\n"%(num,-y))
	header.write("#define C_SFT%i %i\n"%(num,c_sft))
	header.write("#define SEE%i %i\n"%(num,c))
	header.write("#define RINV%i %i\n"%(num,v%p))
	header.write("#define RSFT%i %i\n"%(num,64-Rpow))
	header.write("#define NPOW%i %i\n"%(num,Npow))
	header.write("namespace PBPAS%i{\n"%num)
	code = open(file+".cpp","w")
	code.write("#include \"../../../include/FFT/src/"+file+".h\"\n")
	code.write("#include \"../../../include/FFT/src/arraybitreversal.h\"\n")
	code.write("#include \"../../../include/FFT/src/modpn.h\"\n")
	code.write("#include <iostream>\n")
	code.write("#include <string.h>\n")
	code.write("#define FFT_THRESHOLD %i\n"%(H))
	code.write("#define FFT_THRESHOLD_LOG %i\n"%(len(bin(H))-2))
#---------------NOT USED ANYMORE------------------------
	mask = 1317624576693539401 #12297829382473034410
	mask1 = mask
	mask2 = mask<<1
	mask3 = mask<<2
	if (H&mask1)==0:
		mask1 = mask<<1
		mask2 = mask<<2
		mask3 = mask

	if (H&mask1)==0:
		mask1  = mask<<2
		mask2  = mask
		mask3  = mask<<1
	code.write("#define MY_MASK1 %i\n"%mask1)
	code.write("#define MY_MASK2 %i\n"%mask2)
	code.write("#define MY_MASK3 %i\n"%mask3)
#---------------------------------------------------------

	template = open("generate_fft_template.cpp","r")
	header.write("void Shuffle2(int n, sfixn* A,sfixn* B);\n")
	header.write("void Shuffle(int n, sfixn* A,sfixn* B);\n")
	header.write("sfixn testDFT(int n,int index,sfixn* A,sfixn *W);\n")
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
		elif "namespace PBPAS" in line:
			code.write(line.replace("PBPAS","PBPAS%i"%num))
		elif line == "DFT_ITERATIVE\n":
			generate_dft_iterative(H,code)
			tmpH = H>>1
			while(tmpH>16):
				generate_dft_iterative(tmpH,code,"%i"%tmpH)
				tmpH=tmpH>>1
		elif "FFT_CASES" in line:
			tmpH = 32
			while tmpH<H:
				code.write("\t\t\tcase %i:\n"%tmpH)
				code.write("\t\t\t\tPBPAS::ArrayBitReversalSpe%i(A);\n"%tmpH)
				code.write("\t\t\t\tDFT_iter%i(A,W);\n"%tmpH)
				code.write("\t\t\t\tbreak;\n")
				tmpH=tmpH<<1


		else:
			s = line
			if "MY_PRIME" in s:
				s=(s.replace("MY_PRIME","MY_PRIME%i"%num))
			if "INV_PRIME" in s:
			      s=(s.replace("INV_PRIME","INV_PRIME%i"%num))
			if "C_SFT" in s:
			      s=(s.replace("C_SFT","C_SFT%i"%num))
			if "SEE" in s:
			      s=(s.replace("SEE","SEE%i"%num))
			if "RSFT" in s:
			      s=(s.replace("RSFT","RSFT%i"%num))
			if "RINV" in s:
				s=(s.replace("RINV","RINV%i"%num))
			if "NPOW" in s:
				s=(s.replace("NPOW","NPOW%i"%num))
			code.write(s)
	header.write("}\n")
	header.write("#endif\n")
	header.close()
	template.close()
	code.close();

import sys
H=int(sys.argv[1])
config = open(".config","r")
p = []
files = []
k=0

for line in config:
	if(k%2==0):
		files.append(line.replace("\n",""))
	else:
		p.append(int(line))
	k=k+1
for i in range(0,len(files)):
	generate_everything(files[i],p[i],H,i+1)
generate_reversal_files(H)
config.close()
