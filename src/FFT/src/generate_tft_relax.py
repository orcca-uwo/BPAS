import math

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
	header.write("#include <math.h>\n")
	header.write("#include <algorithm>\n")
	header.write("#include <iostream>\n")
	header.write("#include \"modpn.h\"\n")
	header.write("#include \"transpose.h\"\n")
	header.write("#include \"fft_iter%i.h\"\n"%num)
	header.write("#ifndef TFTSPE%i\n"%num)
	header.write("#define TFTSPE%i\n"%num)
	header.write("#define MY_PRIME%i %i\n"%(num,p))
	header.write("#define INV_PRIME%i %i\n"%(num,-y))
	header.write("#define C_SFT%i %i\n"%(num,c_sft))
	header.write("#define SEE%i %i\n"%(num,c))
	header.write("#define RINV%i %i\n"%(num,v%p))
	header.write("#define RSFT%i %i\n"%(num,64-Rpow))
	header.write("#define R2 2559286960657440491\n")
	header.write("#define Mont_two 3458764513820540920\n")
	#header.write("#define Mont_invtwo -9223372036854775808\n")

	header.write("#define TFT_grainsize 2048\n")


	header.write("#define NPOW%i %i\n"%(num,Npow))
	header.write("using namespace PBPAS%i;\n"%num)
	header.write("namespace TFT_tree%i{\n"%num)
	code = open(file+".cpp","w")
	code.write("#include \"../../../include/FFT/src/"+file+".h\"\n")
	code.write("#include \"../../../include/FFT/src/arraybitreversal.h\"\n")
	code.write("#include \"../../../include/FFT/src/modpn.h\"\n")
	code.write("#include \"../../../include/FFT/src/general_routine.h\"\n")
	#code.write("#include \"../../../src/FFT/src/fft_iter%i.cpp\"\n"%num)
	code.write("#include <cilk/cilk.h>\n")
	code.write("#include <cilk/cilk_api.h>\n")
	code.write("#include <string.h>\n")
	#code.write("#include <cilktools/cilkview.h>\n")




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

	template = open("generate_tft_tree_template.cpp","r")
	header.write("void Shuffle2(int n, sfixn* A,sfixn* B);\n")
	header.write("void Shuffle(int n, sfixn* A,sfixn* B);\n")
	header.write("sfixn testDFT(int n,int index,sfixn* A,sfixn *W);\n")
	for line in template:
		if "void DFT_eff" in line:
			code.write(line.replace("DFT_eff","DFT_eff_p%i"%num))
			header.write(line.replace("{",";\n").replace("DFT_eff","DFT_eff_p%i"%num))


		elif "void InvDFTKeepMont_eff" in line:
			code.write(line.replace("InvDFTKeepMont_eff","InvDFTKeepMont_eff_p%i"%num))
			header.write(line.replace("{",";\n").replace("InvDFTKeepMont_eff","InvDFTKeepMont_eff_p%i"%num))

##########         for tft functions         #########################
		elif "void Mont_ITFT_Core" in line:
			code.write(line.replace("Mont_ITFT_Core","Mont_ITFT_Core_p%i"%num))
			header.write(line.replace("{",";\n").replace("Mont_ITFT_Core","Mont_ITFT_Core_p%i"%num))
		elif "Mont_ITFT_Core" in line:
			code.write(line.replace("Mont_ITFT_Core","Mont_ITFT_Core_p%i"%num))

		elif "static inline void  TFT_AddSubSpeSSEModInplace" in line:
			code.write(line.replace("TFT_AddSubSpeSSEModInplace","TFT_AddSubSpeSSEModInplace"))
			header.write(line.replace("{",";\n").replace("TFT_AddSubSpeSSEModInplace","TFT_AddSubSpeSSEModInplace"))
		elif "TFT_AddSubSpeSSEModInplace" in line:
			code.write(line.replace("TFT_AddSubSpeSSEModInplace","TFT_AddSubSpeSSEModInplace"))

		elif "void inline TFT_4POINT" in line:
			code.write(line.replace("TFT_4POINT","TFT_4POINT_p%i"%num))
			header.write(line.replace("{",";\n").replace("TFT_4POINT","TFT_4POINT_p%i"%num))
		elif "TFT_4POINT" in line:
			code.write(line.replace("TFT_4POINT","TFT_4POINT_p%i"%num))

		elif "void inline TFT_8POINT" in line:
			code.write(line.replace("TFT_8POINT","TFT_8POINT_p%i"%num))
			header.write(line.replace("{",";\n").replace("TFT_8POINT","TFT_8POINT_p%i"%num))
		elif "TFT_8POINT" in line:
			code.write(line.replace("TFT_8POINT","TFT_8POINT_p%i"%num))

		elif "void inline TFT_16POINT" in line:
			code.write(line.replace("TFT_16POINT","TFT_16POINT_p%i"%num))
			header.write(line.replace("{",";\n").replace("TFT_16POINT","TFT_16POINT_p%i"%num))
		elif "TFT_16POINT" in line:
			code.write(line.replace("TFT_16POINT","TFT_16POINT_p%i"%num))

		elif "void inline TFT_iter32" in line:
			code.write(line.replace("TFT_iter32","TFT_iter32_p%i"%num))
			header.write(line.replace("{",";\n").replace("TFT_iter32","TFT_iter32_p%i"%num))
		elif "TFT_iter32" in line:
			code.write(line.replace("TFT_iter32","TFT_iter32_p%i"%num))




		elif "void TFT_Basecase" in line:
			code.write(line.replace("TFT_Basecase","TFT_Basecase_p%i"%num))
			header.write(line.replace("{",";\n").replace("TFT_Basecase","TFT_Basecase_p%i"%num))
		elif "TFT_Basecase" in line:
			code.write(line.replace("TFT_Basecase","TFT_Basecase_p%i"%num))






		elif "void ITFT_DFT" in line:
			code.write(line.replace("ITFT_DFT","ITFT_DFT_p%i"%num))
			header.write(line.replace("{",";\n").replace("ITFT_DFT","ITFT_DFT_p%i"%num))
		elif "ITFT_DFT" in line:
			code.write(line.replace("ITFT_DFT","ITFT_DFT_p%i"%num))


		elif "void ITFT_Core" in line:
			code.write(line.replace("ITFT_Core","ITFT_Core_p%i"%num))
			header.write(line.replace("{",";\n").replace("ITFT_Core","ITFT_Core_p%i"%num))
		elif "ITFT_Core" in line:
			code.write(line.replace("ITFT_Core","ITFT_Core_p%i"%num))
		elif "void ITFT_Wrapper" in line:
			code.write(line.replace("ITFT_Wrapper","ITFT_Wrapper_p%i"%num))
			header.write(line.replace("{",";\n").replace("ITFT_Wrapper","ITFT_Wrapper_p%i"%num))
		elif "ITFT_Wrapper" in line:
			code.write(line.replace("ITFT_Wrapper","ITFT_Wrapper_p%i"%num))
		elif "sfixn normaltomont" in line:
			code.write(line.replace("normaltomont","normaltomont"))
			header.write(line.replace("{",";\n").replace("normaltomont","normaltomont"))
		elif "normaltomont" in line:
			code.write(line.replace("normaltomont","normaltomont"))
		elif "sfixn MontDivMod" in line:
			code.write(line.replace("MontDivMod","MontDivMod"))
			header.write(line.replace("{",";\n").replace("MontDivMod","MontDivMod"))
		elif "MontDivMod" in line:
			code.write(line.replace("MontDivMod","MontDivMod"))

		elif "void InitTFT_Tree" in line:
			code.write(line.replace("InitTFT_Tree","InitTFT_Tree_p%i"%num))
			header.write(line.replace("{",";\n").replace("InitTFT_Tree","InitTFT_Tree_p%i"%num))
		elif "void TFT_Core" in line:
			code.write(line.replace("TFT_Core","TFT_Core_p%i"%num))
			header.write(line.replace("{",";\n").replace("TFT_Core","TFT_Core_p%i"%num))
		elif "namespace TFT_tree" in line:
			code.write(line.replace("TFT_tree","TFT_tree%i"%num))
		elif "void Shuffle1" in line:
			code.write(line.replace("Shuffle1","Shuffle1"))
			header.write(line.replace("{",";\n").replace("Shuffle1","Shuffle1"))
		elif "void TFT_twiddle" in line:
			code.write(line.replace("TFT_twiddle","TFT_twiddle"))
			header.write(line.replace("{",";\n").replace("TFT_twiddle","TFT_twiddle"))
		elif "TFT_Core" in line:
			code.write(line.replace("TFT_Core","TFT_Core_p%i"%num))
		elif "void TFT_New_Core" in line:
			code.write(line.replace("TFT_New_Core","TFT_New_Core_p%i"%num))
			header.write(line.replace("{",";\n").replace("TFT_New_Core","TFT_New_Core_p%i"%num))
		elif "TFT_New_Core" in line:
			code.write(line.replace("TFT_New_Core","TFT_New_Core_p%i"%num))


		elif "DFT_rec" in line:
			code.write(line.replace("DFT_rec","DFT_rec_p%i"%num))
		elif "DFT_eff" in line:
			code.write(line.replace("DFT_eff","DFT_eff_p%i"%num))
		elif "void Shuffle_tft" in line:
			code.write(line.replace("Shuffle_tft","Shuffle_tft"))
			header.write(line.replace("{",";\n").replace("Shuffle_tft","Shuffle_tft"))
		elif "Shuffle_tft" in line:
			code.write(line.replace("Shuffle_tft","Shuffle_tft"))

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
				#code.write("\t\t\t\tPBPAS::ArrayBitReversalSpe%i(A);\n"%tmpH)  #reverse for tft order
				code.write("\t\t\t\tbreak;\n")
				tmpH=tmpH<<1

		elif "void Shuffle" in line:
				header.write(line.replace("{",";\n").replace("Shuffle","Shuffle"))
				code.write(line.replace("Shuffle","Shuffle"))
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
config = open(".config_tft_tree","r")
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
config.close()
