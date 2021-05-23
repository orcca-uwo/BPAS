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
	s = "void ArrayBitReversalSpec" + name +"(long long int *A)"
	header.write(s+";\n")
	s +="{\n"
	s += "\tlong long int t0,t1,t2,t3,t4,t5,t6,t7;\n"
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
def generate_reversal_files(arrayfile,H):
	reversalheader = open("../../../include/FFT/src/"+arrayfile+".h","w")
	reversalcode = open(arrayfile+".c","w")
	reversalheader.write("//This header genreated for FFT_THRESHOLD=%i\n\n"%H);
	# reversalheader.write("namespace PBPAS{\n")
	reversalcode.write("#include \"../../../include/FFT/src/arraybitreversal_%i.h\"\n"%H);
	# reversalcode.write("namespace PBPAS{\n")
	generate_Array_Bit_Reversal(reversalcode,reversalheader,RevBidMap(H),"_H%i"%H)
	tmpH=H>>1
	while tmpH>=8:
		generate_Array_Bit_Reversal(reversalcode,reversalheader,RevBidMap(tmpH),"_H%i_%i"%(H,tmpH))
		tmpH=tmpH>>1
	# reversalheader.write("}")
	# reversalcode.write("}")
	reversalheader.close()
	return

def generate_16_point_fft(H,code):
	s=""
	s+="\tfor(int k=0;k<%i;k+=16){\n"%H
	s+="\t\t//FFT_16POINT(A+k,Wp);\n"
	s+="\t\tlong long int u = A[k+0];\n"
	s+="\t\tlong long int t = A[k+1];\n"
	s+="\t\tA[k+0] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+1] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tu = A[k+2];\n"
	s+="\t\tt = A[k+3];\n"
	s+="\t\tA[k+2] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+3] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tu = A[k+4];\n"
	s+="\t\tt = A[k+5];\n"
	s+="\t\tA[k+4] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+5] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tu = A[k+6];\n"
	s+="\t\tt = A[k+7];\n"
	s+="\t\tA[k+6] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+7] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tu = A[k+8];\n"
	s+="\t\tt = A[k+9];\n"
	s+="\t\tA[k+8] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+9] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tu = A[k+10];\n"
	s+="\t\tt = A[k+11];\n"
	s+="\t\tA[k+10] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+11] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tu = A[k+12];\n"
	s+="\t\tt = A[k+13];\n"
	s+="\t\tA[k+12] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+13] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tu = A[k+14];\n"
	s+="\t\tt = A[k+15];\n"
	s+="\t\tA[k+14] = smallprimefield_AddMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+15] = smallprimefield_SubMod(u,t,MY_PRIME);\n"
	s+="\t\tA[k+3] = smallprimefield_Mul_INLINE(A[k+3],*(Wp-3),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+7] = smallprimefield_Mul_INLINE(A[k+7],*(Wp-3),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+11] = smallprimefield_Mul_INLINE(A[k+11],*(Wp-3),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+15] = smallprimefield_Mul_INLINE(A[k+15],*(Wp-3),INV_PRIME,MY_PRIME);\n"
	s+="\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k,A+k+2,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+4,A+k+6,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+8,A+k+10,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+12,A+k+14,MY_PRIME);\n"
	s+="\t\tA[k+5] = smallprimefield_Mul_INLINE(A[k+5],*(Wp-11),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+6] = smallprimefield_Mul_INLINE(A[k+6],*(Wp-10),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+7] = smallprimefield_Mul_INLINE(A[k+7],*(Wp-9),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+13] = smallprimefield_Mul_INLINE(A[k+13],*(Wp-11),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+14] = smallprimefield_Mul_INLINE(A[k+14],*(Wp-10),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+15] = smallprimefield_Mul_INLINE(A[k+15],*(Wp-9),INV_PRIME,MY_PRIME);\n"
	s+="\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k,A+k+4,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+2,A+k+6,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+8,A+k+12,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+10,A+k+14,MY_PRIME);\n"
	s+="\n"
	s+="\t\tA[k+9] = smallprimefield_Mul_INLINE(A[k+9],*(Wp-27),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+10] = smallprimefield_Mul_INLINE(A[k+10],*(Wp-26),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+11] = smallprimefield_Mul_INLINE(A[k+11],*(Wp-25),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+12] = smallprimefield_Mul_INLINE(A[k+12],*(Wp-24),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+13] = smallprimefield_Mul_INLINE(A[k+13],*(Wp-23),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+14] = smallprimefield_Mul_INLINE(A[k+14],*(Wp-22),INV_PRIME,MY_PRIME);\n"
	s+="\t\tA[k+15] = smallprimefield_Mul_INLINE(A[k+15],*(Wp-21),INV_PRIME,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k,A+k+8,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+2,A+k+10,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+4,A+k+12,MY_PRIME);\n"
	s+="\t\tAddSubSpecSSEModInplace(A+k+6,A+k+14,MY_PRIME);\n"
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
					s+="\tA[%i] = smallprimefield_Mul_INLINE(A[%i],*(Wp+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\tA[%i] = smallprimefield_Mul_INLINE(A[%i],*(Wp+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\tAddSubSpecSSEModInplace(A+%i,A+%i,MY_PRIME);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
		if i!=log:
			s+="\tWp = Wp-(1L<<%i);\n"%(depth+i)
	#s+="asm(\"\":::\"memory\");\n"
	s+="\tfor(int j=2;j<%i;j+=2){\n"%jump
	s+="\t\tWp = Wt;\n"
	for i in range(1,log+1):
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tA[j+%i] = smallprimefield_Mul_INLINE(A[j+%i],*(Wp+j+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tA[j+%i] = smallprimefield_Mul_INLINE(A[j+%i],*(Wp+j+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tAddSubSpecSSEModInplace(A+j+%i,A+j+%i,MY_PRIME);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
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
					s+="\t\tA[k+%i] = smallprimefield_Mul_INLINE(A[k+%i],*(Wp+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tA[k+%i] = smallprimefield_Mul_INLINE(A[k+%i],*(Wp+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\tAddSubSpecSSEModInplace(A+k+%i,A+k+%i,MY_PRIME);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
		if i!=log:
			s+="\t\tWp = Wp-(1L<<%i);\n"%(depth+i)
	s+="\t\tfor(int j=2;j<%i;j+=2){\n"%jump
	s+="\t\t\tWp = Wt;\n"
	for i in range(1,log+1):
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\t\tA[k+j+%i] = smallprimefield_Mul_INLINE(A[k+j+%i],*(Wp+j+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump,(k*2**i+l)*jump,(l-2**(i-1))*jump)
		for l in range(2**(i-1),2**(i)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\t\tA[k+j+%i] = smallprimefield_Mul_INLINE(A[k+j+%i],*(Wp+j+%i),INV_PRIME,MY_PRIME);\n"%(((2**i)*k+l)*jump+1,(k*2**i+l)*jump+1,(l-2**(i-1))*jump+1)
		for l in range(0,2**(i-1)):
			for k in range(0,int((2**log)/(2**i))):
				s+="\t\t\tAddSubSpecSSEModInplace(A+k+j+%i,A+k+j+%i,MY_PRIME);\n"%((k*2**i+l)*jump,(k*2**i+l+2**(i-1))*jump)
		if i!=log:
			s+="\t\t\tWp = Wp-(1L<<%i);\n"%(depth+i)
	s+="\t\t}\n"
	s+="\t}\n"
	code.write(s)

def generate_dft_iterative(size,code,name=""):
	code.write("void DFT_iter"+ name+"(long long int *A,long long int *W,long long int INV_PRIME, long long int MY_PRIME){\n")
	code.write("\tlong long int *Wp = W + (%i<<1)-4;\n"%size)
	code.write("\tlong long int* Wt;\n")
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
		log=log+1
		# rem=rem>>1
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

def generate_everything(file,arrayfile,H):
	# p_1 = p-1
	# Npow = 0
	# while (p_1%(1<<Npow)==0):
	# 	Npow=Npow+1
	# Npow = Npow - 1
	# c = p_1>>Npow
	# Rpow = int(math.log(p,2)+1)
	# c_sft = c<<Npow
	# (r,v,u)= extended_gcd(2**Rpow,p)
	# y = u%(2**Rpow)
	# x = ((u-y)/(2**Rpow))*p
	# if y>0:
	# 	y = y-2**Rpow
	# 	x = x+p

	header = open("../../../include/FFT/src/"+file+".h","w")
	header.write("#ifndef _FFT_SPEC_%i_\n"%H);
	header.write("#define _FFT_SPEC_%i_\n\n"%H);
	header.write("#include \"../../FiniteFields/SmallPrimeField_Support.h\"\n")

	code = open(file+".c","w")
	code.write("#include \"../../../include/FFT/src/"+file+".h\"\n")
	code.write("#include \"../../../include/FFT/src/"+arrayfile+".h\"\n")
	code.write("#define FFT_THRESHOLD %i\n"%(H))
	code.write("#define FFT_THRESHOLD_LOG %i\n"%(len(bin(H))-2))

	template = open("generate_fft_iter_gen_template.c","r")
	for line in template:
		if "void DFT_eff" in line:
			code.write(line.replace("DFT_eff","DFT_eff_H%i"%H))
			header.write(line.replace("{",";\n").replace("DFT_eff","DFT_eff_H%i"%H))
		elif "void InvDFT_eff" in line:
			code.write(line.replace("InvDFT_eff","InvDFT_eff_H%i"%H))
			header.write(line.replace("{",";\n").replace("InvDFT_eff","InvDFT_eff_H%i"%H))
		elif "void InvDFTKeepMont_eff" in line:
			code.write(line.replace("InvDFTKeepMont_eff","InvDFTKeepMont_eff_H%i"%H))
			header.write(line.replace("{",";\n").replace("InvDFTKeepMont_eff","InvDFTKeepMont_eff_H%i"%H))
		elif "DFT_iter" in line:
			code.write(line.replace("DFT_iter","DFT_iter_H%i"%H))
		elif "DFT_rec" in line:
			code.write(line.replace("DFT_rec","DFT_rec_H%i"%H))
		elif "DFT_eff" in line:
			code.write(line.replace("DFT_eff","DFT_eff_H%i"%H))
		elif "ArrayBitReversalSpec8" in line:
			code.write(line.replace("ArrayBitReversalSpec8","ArrayBitReversalSpec_H%i_8"%H))
		elif "ArrayBitReversalSpec16" in line:
			code.write(line.replace("ArrayBitReversalSpec16","ArrayBitReversalSpec_H%i_16"%H))
		elif "ArrayBitReversalSpec" in line:
			code.write(line.replace("ArrayBitReversalSpec","ArrayBitReversalSpec_H%i"%H))
		elif line == "DFT_ITERATIVE\n":
			generate_dft_iterative(H,code,"_H%i"%H);
			tmpH = H>>1
			while(tmpH>16):
				generate_dft_iterative(tmpH,code,"_H%i_%i"%(H,tmpH));
				tmpH=tmpH>>1
		elif "FFT_CASES" in line:
			tmpH = 32
			while tmpH<H:
				code.write("\t\t\tcase %i:\n"%tmpH)
				code.write("\t\t\t\tArrayBitReversalSpec_H%i_%i(A);\n"%(H,tmpH))
				code.write("\t\t\t\tDFT_iter_H%i_%i(A,W,INV_PRIME,MY_PRIME);\n"%(H,tmpH));
				code.write("\t\t\t\tbreak;\n")
				tmpH=tmpH<<1
		else :
			code.write(line);

	header.write("#endif\n")
	header.close()
	template.close()
	code.close();

import sys
H=int(sys.argv[1])

outfile = "fft_iter_%i"%H;
arrayfile = "arraybitreversal_%i"%H;
generate_everything(outfile, arrayfile,H)
generate_reversal_files(arrayfile,H)
