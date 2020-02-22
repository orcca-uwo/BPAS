# This python script builds the appropriate DFT functions for the Finite fields in the filename provided as an argument to this script
import sys # cmd line args
import re # import regular expression module
# mac os --> import psutil
import os # get memory info
from collections import OrderedDict

# print("Number of arguments:")
# print(len(sys.argv))
# print("arguments.")
# print("Argument List:")
# print(str(sys.argv))

if (len(sys.argv)!= 1):
    print("Usage: python dft_config.py")
    sys.exit()

#os.system()
#print os.system(./getprocessorinfo).read()

#file_to_scan = str(sys.argv[1]);
#print("Scanning file:")
#print(file_to_scan)

#avail_memory = None
#mem = psutil.virtual_memory()
#avail_memory = mem.total
#avail_cpu = psutil.cpu_count()
#mem_total_kb = avail_memory / avail_cpu / 1024 ** 2
# print (avail_cpu)
# print (avail_memory)

''' Hardcoded Primes '''
# initial
# N=4096
# big_prime = 4179340454199820289
# omega = 3332048709103214402
# small_prime = 962592769
# omega = 379119492

#N= 4096
#p1 = ((2**63+2**53)**2)+1 -> 85236826359346144956638323529482240001
#o1 = 2631344246756619346648827348049660483
#p2 = ((2**64-2**50)**4)+1 -> 115763822272329310636028559609001827025179711501300126126825041166177555972097
#o2 = 24621278244416973245097080599157415662973930984945449750138591565788235778527

''' Return the information in /proc/meminfo
as a dictionary '''
def meminfo():
    meminfo=OrderedDict()
    with open('/proc/meminfo') as f:
        for line in f:
            meminfo[line.split(':')[0]] = line.split(':')[1].strip()
    return meminfo

''' Return the information in /proc/cpuinfo
as a dictionary in the following format:
    cpu_info['proc0']={...}
    cpu_info['proc1']={...}'''
def cpuinfo():
    cpuinfo=OrderedDict()
    procinfo=OrderedDict()
    nprocs = 0
    with open('/proc/cpuinfo') as f:
        for line in f:
            if not line.strip():
                # end of one processor
                cpuinfo['proc%s' % nprocs] = procinfo
                nprocs=nprocs+1
                # Reset
                procinfo=OrderedDict()
            else:
                if len(line.split(':')) == 2:
                    procinfo[line.split(':')[0].strip()] = line.split(':')[1].strip()
                else:
                    procinfo[line.split(':')[0].strip()] = ''
    return cpuinfo

cpuinfo = cpuinfo()
cachesize = 0
for processor in cpuinfo.keys():
    cz = cpuinfo[processor]['cache size']
    if (cz>cachesize):
        cachesize = cz
## 8kB +12 kuOps trace cache
print('cachesize:',cachesize)

#meminfo = meminfo()
#print(meminfo)
#print('Total memory: {0}'.format(meminfo['MemTotal']))
#print('Free memory: {0}'.format(meminfo['MemFree']))

# try:
#      avail_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
#      meminfo = open('/proc/meminfo').read()
#      matched = re.search(r'^MemTotal:\s+(\d+)', meminfo)
#      if matched:
#          mem_total_kB = int(matched.groups()[0])
# except ValueError:
#      print("Unable to get memory info")
#
# print()

#CREATE BASECASES

LARGESTFIELDELEMENT = 16 # need to set this properly

basecases = []

if (LARGESTFIELDELEMENT * 8 < cachesize):
   	basecases.append(8)
if (LARGESTFIELDELEMENT * 16 < cachesize):
    basecases.append(16)
if (LARGESTFIELDELEMENT * 32 < cachesize):
    basecases.append(32)
if (LARGESTFIELDELEMENT * 64 < cachesize):
    basecases.append(64)


'''
return a list of tokens that are found in group grp
matched by the regular expression in reg_exp
'''
def get_token_list(filename, reg_exp, grp):
    lines = open(filename, 'r').readlines()
    regex = re.compile(reg_exp)
    tokens = []
    for line in lines:
        search_result = regex.search(line)
        if search_result != None:
            tokens.append(search_result.group(grp))
    return tokens

'''
return the line number where the line of text in the file
matches the input text
or -1 if not found
'''
def get_line_numbers(file_name,text):
    #read lines into a list of strings
    lines = open(file_name, 'r').readlines()
    regex = re.compile(text)
    result = []
    line_number=1
    for line in lines:
	       search_result = regex.search(line)
	       if search_result:
  	            result.append(line_number)
	       line_number+=1
    if len(result)>0:
	    return result
    else:
        return -1

'''return n lines from file_name
from start to start+n'''
def get_lines_from_file(file_name,start,n):
    result = []
    lines = open(file_name,'r').readlines()
    for i in range(0,n):
	       result.append(lines[start+i])
    return result

'''overwrite the line of text at line_num with the text in text
in the file named file_name.'''
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

'''overwrite num_lines of contigious lines of text starting at line_num
with the lines of text contained in linesoftext'''
def replace_lines(file_name,line_num,num_lines,linesoftext):
    lines = open(file_name,'r').readlines()
    for i in range(0,num_lines):
	       lines[line_num+i]=linesoftext[i];
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

'''return a list of indices marking the occurances of the substring in fullstring'''
def substring_indices(fullstring, substring):
    results = []
    sub_len = len(substring)
    for i in range(len(fullstring)):  # range returns a list of values from 0 to (len(fullstring) - 1)
            if fullstring[i:i+sub_len] == substring:# this is slice notation; it means take characters i up to (but not including) i + the length of th substring
                results.append(i)
    return results

''' return the number of occurences of substring in string'''
def num_occurences(fullstring, substring):
    results = substring_indices(fullstring,substring)
    return len(results)

''' convert a file into a string '''
def file_to_string(file_name):
    lines = open(file_name,'r').readlines()
    result = ''.join(lines)
    return result


def write_DFT2(fieldtype="none"):
    out="\n"
    if fieldtype=="none":
        out+="inline void DFT2 (long int* a0,long int* a1,long int prime){\n"
        out+="long int sum;\n"
        out+="sum=smallprimefield_add(*a0,*a1,prime);\n"
        out+="*a1=smallprimefield_sub(*a0,*a1,prime);\n"
        out+="*a0=sum;\n"
        out+="}\n"
    else:
        out+="template<class FiniteField>\n"
        out+="inline void DFT2 (FiniteField* a0,FiniteField* a1){\n"
        out+="\tFiniteField sum;\n"
        out+="\tsum=*a0+*a1;\n"
        out+="\t*a1=*a0-*a1;\n"
        out+="\t*a0=sum;\n"
        out+="}\n"
        out+="template void DFT2<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* b);\n"
        out+="template void DFT2<BigPrimeField>(BigPrimeField* a,BigPrimeField* b);\n"
        out+="template void DFT2<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* b);\n"
    out+="\n"
    return out

def write_swap(fieldtype="none"):
    out="\n"
    if fieldtype=="none":
        out+="inline void swap(long int* a,long int* b){\n"
        out+="\tlong int tmp;\n"
        out+="\ttmp = *a;\n"
        out+="\t*a = *b;\n"
        out+="\t*b = tmp;\n"
        out+="}\n"
    else:
        out+="template<class FiniteField>\n"
        out+="inline void swap(FiniteField* a,FiniteField* b){\n"
        out+="\tFiniteField tmp;\n"
        out+="\ttmp = *a;\n"
        out+="\t*a = *b;\n"
        out+="\t*b = tmp;\n"
        out+="\t}\n"
        out+="template void swap<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* b);\n"
        out+="template void swap<BigPrimeField>(BigPrimeField*a,BigPrimeField*b);\n"
        out+="template void swap<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField*a,GeneralizedFermatPrimeField*b);\n"
    out+="\n"
    return out

def write_host_pow():
    out="\n"
    out+="inline long int host_pow(long int x, long int e, long int prime, long long int R, long long int pP){\n"
    out+="\tlong int m=1;\n"
    out+="\tfor(long int i=0; i<e;i++)\n"
    out+="\t\tm=smallprimefield_multi(x, m, prime, R, pP);\n"
    out+="\treturn m;\n"
    out+="}\n"
    return out

def write_precompute_header(fieldtype="none",omega_type="outer"):
    out="\n"
    if omega_type=="outer":
	if fieldtype=="spf" or fieldtype=="bpf":
		out+="template<class FiniteField>\n"
		out+="FiniteField** precomputeOuterPowersOfOmega(int K, int e, FiniteField omega);\n"
	else:
		out+="long int** precomputeOuterPowersOfOmega(int K, int e, long int omega, long int prime, long long int R);\n"
    else:
	if fieldtype=="spf" or fieldtype=="bpf":
		out+="template<class FiniteField>\n"
		out+="FiniteField* precomputeInnerPowersOfOmega(FiniteField omega,int basecase);\n"
	else:
		out+="long int* precomputeInnerPowersOfOmega(long int omega,int basecase,long int prime,long long int R,long long int pP);\n"
    return out

def write_precompute_powers_of_omega(fieldtype="none",omega_type="outer"):
    out="\n"
    if omega_type=="outer":
        if fieldtype=='spf' or fieldtype=='bpf':
            out+="// precompute the powers of omega for the run from N to K\n"
            out+="template<class FiniteField>\n"
            out+="FiniteField** precomputeOuterPowersOfOmega(int K, int e, FiniteField omega){\n"
            out+="// setup and array allocation\n"
            out+="\tFiniteField** precomputed_omegas = new FiniteField*[e];\n"
            out+="\t// set stride\n"
            out+="\tint N=pow(K,e);\n"
            out+="\tint stride=N;\n"
            out+="\t// for each level from e-2 -> 0 allocate space\n"
            out+="\tfor (int l=0;l<=e-2;l++){\n"
            out+="\t\tprecomputed_omegas[l]= new FiniteField[stride];\n"
            out+="\t\tstride/=K; // update stride\n"
            out+="\t\t// Compute D_{K}{K^e-1}=D_{K}{stride}\n"
            out+="\t\tfor(int i=0;i<K;i++){\n"
            out+="\t\t\tfor(int j=0;j<stride;j++){\n"
            out+="\t\t\t\tprecomputed_omegas[l][i*stride+j]=omega^(i*j);\n"
            out+="\t\t\t}\n"
            out+="\t\t}\n"
            out+="\t}\n"
            out+="\treturn precomputed_omegas;\n"
            out+="}\n"
            out+="template SmallPrimeField** precomputeOuterPowersOfOmega<SmallPrimeField>(int K,int e,SmallPrimeField omega);\n"
            out+="template BigPrimeField** precomputeOuterPowersOfOmega<BigPrimeField>(int K,int e,BigPrimeField omega);\n"
        else:
            out+="//precompute the powers of omega for the run from N to K\n"
            out+="long int** precomputeOuterPowersOfOmega(int K, int e, long int omega, long int prime, long long int R){\n"
            out+="\tlong int pP = smallprimefield_getPp(prime,R);\n"
            out+="\tlong int** precomputed_omegas = (long int**)malloc(e*sizeof(long int*));\n"
            out+="\t// set stride\n"
       	    out+="\tint N=pow(K,e);\n"
            out+="\tint stride=N;\n"
            out+="\t// for each level from e-2 -> 0 allocate space\n"
            out+="\tfor (int l=0;l<=e-2;l++){\n"
            out+="\t\tprecomputed_omegas[l]=(long int*)malloc(stride*sizeof(long int));\n"
            out+="\t\tstride/=K; // update stride\n"
            out+="\t\tfor(int i=0;i<K;i++){\n"
    	    out+="\t\t\tfor(int j=0;j<stride;j++){\n"
    	    out+="\t\t\t\tprecomputed_omegas[l][i*stride+j]=host_pow(omega,(long int)(i*j),prime,R,pP);"
    	    out+="\t\t\t}\n"
    	    out+="\t\t}\n"
            out+="\t}\n"
            out+="\treturn precomputed_omegas;\n"
            out+="}\n"
    else:
        if fieldtype=='spf' or fieldtype=='bpf':
            out+="// precomputes the powers of omega for the provided basecase\n"
            out+="template<class FiniteField>\n"
            out+="FiniteField* precomputeInnerPowersOfOmega(FiniteField omega,int basecase){\n"
            out+="\tFiniteField* omegas;\n"
            out+="\tint arr_size;\n"
            out+="\tarr_size = basecase >> 1;\n"
            out+="\tomegas = new FiniteField[arr_size];\n"
            out+="\tomegas[0] = 1;\n"
            out+="\tomegas[1] = omega;\n"
            out+="\tfor (int i=2;i<arr_size;++i){\n"
            out+="\t\tomegas[i]=omegas[i-1]*omega;\n"
            out+="\t}\n"
            out+="\treturn omegas;\n"
            out+="}\n"
            out+="template SmallPrimeField* precomputeInnerPowersOfOmega<SmallPrimeField>(SmallPrimeField omega,int basecase);\n"
            out+="template BigPrimeField* precomputeInnerPowersOfOmega<BigPrimeField>(BigPrimeField omega,int basecase);\n"
        else:
            out+="// precomputes the powers of omega for the provided basecase\n"
            out+="long int* precomputeInnerPowersOfOmega(long int omega,int basecase,long int prime,long long int R,long long int pP){\n"
            out+="\tlong int* omegas;\n"
            out+="\tint arr_size;\n"
            out+="\tarr_size = basecase >> 1;\n"
            out+="\tomegas = (long int*)malloc(sizeof(long int)*arr_size);\n"
            out+="\tomegas[0] = 1;\n"
            out+="\tomegas[1] = omega;\n"
            out+="\tfor (int i=2;i<arr_size;++i){\n"
            out+="\t\tomegas[i]=smallprimefield_multi(omegas[i-1],omega,prime,R,pP);\n"
            out+="\t}\n"
            out+="\treturn omegas;\n"
            out+="}\n"
    out+="\n"
    return out

#def write_precompute_powers_of_omega(omega_type="outer"):
#    out="\n"
#    if omega_type=="outer":
#        out+="// precompute the powers of omega for the run from N to K\n"
#        out+="template<class FiniteField>\n"
#        out+="FiniteField** precomputeOuterPowersOfOmega(int K, int e, FiniteField omega){\n"
#        out+="// setup and array allocation\n"
#        out+="\tFiniteField** precomputed_omegas = new FiniteField*[e];\n"
#        out+="\t// set stride\n"
#        out+="\tint N=pow(K,e);\n"
#        out+="\tint stride=N;\n"
#        out+="\t// for each level from e-2 -> 0 allocate space\n"
#        out+="\tfor (int l=0;l<=e-2;l++){\n"
#        out+="\t\tprecomputed_omegas[l]= new FiniteField[stride];\n"
#        out+="\t\tstride/=K; // update stride\n"
#        out+="\t\t// Compute D_{K}{K^e-1}=D_{K}{stride}\n"
#        out+="\t\tfor(int i=0;i<K;i++){\n"
#        out+="\t\t\tfor(int j=0;j<stride;j++){\n"
#        out+="\t\t\t\tprecomputed_omegas[l][i*stride+j]=omega^(i*j);\n"
#        out+="\t\t\t}\n"
#        out+="\t\t}\n"
#        out+="\t}\n"
#        out+="\treturn precomputed_omegas;\n"
#        out+="}\n"
#        out+="template SmallPrimeField** precomputeOuterPowersOfOmega<SmallPrimeField>(int K,int e,SmallPrimeField omega);\n"
#        out+="template BigPrimeField** precomputeOuterPowersOfOmega<BigPrimeField>(int K,int e,BigPrimeField omega);\n"
#    else:
#        out+="// precomputes the powers of omega for the provided basecase\n"
#        out+="template<class FiniteField>\n"
#        out+="FiniteField* precomputeInnerPowersOfOmega(FiniteField omega,int basecase){\n"
#        out+="\tFiniteField* omegas;\n"
#        out+="\tint arr_size;\n"
#        out+="\tarr_size = basecase >> 1;\n"
#        out+="\tomegas = new FiniteField[arr_size];\n"
#        out+="\tomegas[0] = 1;\n"
#        out+="\tomegas[1] = omega;\n"
#        out+="\tfor (int i=2;i<arr_size;++i){\n"
#        out+="\t\tomegas[i]=omegas[i-1]*omega;\n"
#        out+="\t}\n"
#        out+="\treturn omegas;\n"
#        out+="}\n"
#        out+="template SmallPrimeField* precomputeInnerPowersOfOmega<SmallPrimeField>(SmallPrimeField omega,int basecase);\n"
#        out+="template BigPrimeField* precomputeInnerPowersOfOmega<BigPrimeField>(BigPrimeField omega,int basecase);\n"
#    out+="\n"
#    return out

''' write header code for dft 8,16,32,64
    fieldtype options are spf,bpf,gfpf'''
def write_dft_header(fieldtype="none",basecase="8"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    elif(fieldtype=="bpf"):
        bpf_bool = True
    elif(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool:
        out+="template <class FiniteField>\n"
        out+="extern  FiniteField* DFT_"
        out+=str(basecase)
        out+="(FiniteField* a,FiniteField* omegas);\n"
    elif bpf_bool:
        out+="template <class FiniteField>\n"
        out+="extern  FiniteField* DFT_"
        out+=str(basecase)
        out+="(FiniteField* a,FiniteField* omegas);\n"
    elif gfpf_bool:
        out+="inline GeneralizedFermatPrimeField* DFT_"
        out+=str(basecase)
        out+="(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas);\n"
    else:
	out+="long int* DFT_"
        out+=str(basecase)
        out+="(long int* a,long int* omegas, long int prime,long long int R,long long int pP);\n"

    return out

'''write dft code for basecase size 8
fieldtype options: "spf", "bpf", "gfpf"
default builds c version'''
def write_dft_8(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    elif(fieldtype=="bpf"):
        bpf_bool = True
    elif(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="template<class FiniteField>\n"
        out+="inline FiniteField* DFT_8(FiniteField* a,FiniteField* omegas){\n"
    elif bpf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="template<class FiniteField>\n"
        out+="inline FiniteField* DFT_8(FiniteField* a,FiniteField* omegas){\n"
    elif gfpf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="inline GeneralizedFermatPrimeField* DFT_8(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas){\n"
    else:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="\param prime\n"
        out+="\param R\n"
        out+="\param pP\n"
        out+="*/\n"
        out+="long int* DFT_8(long int* a,long int* omegas, long int prime,long long int R,long long int pP){\n"

    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&a[0],&a[4]); // dft on permutated indexes\n"
        out+="\tDFT2(&a[1],&a[5]);\n"
        out+="\tDFT2(&a[2],&a[6]);\n"
        out+="\tDFT2(&a[3],&a[7]);\n"
    else:
        out+="\tDFT2(&a[0],&a[4],prime); // dft on permutated indexes\n"
        out+="\tDFT2(&a[1],&a[5],prime);\n"
        out+="\tDFT2(&a[2],&a[6],prime);\n"
        out+="\tDFT2(&a[3],&a[7],prime);\n"

    out+="\n"
    if gfpf_bool:
        out+="\ta[6]=a[6].MulPowR(2);   //twiddle\n"
        out+="\ta[7]=a[7].MulPowR(2);\n"
    elif spf_bool or bpf_bool:
        out+="\ta[6]=a[6]*omegas[2];   //twiddle\n"
        out+="\ta[7]=a[7]*omegas[2];\n"
    else:
        out+="\ta[6]=smallprimefield_multi(a[6],omegas[2],prime,R,pP);   // twiddle\n"
        out+="\ta[7]=smallprimefield_multi(a[7],omegas[2],prime,R,pP);\n"

    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&a[0],&a[2]); // dft on permutated indexes\n"
        out+="\tDFT2(&a[4],&a[6]);\n"
        out+="\tDFT2(&a[1],&a[3]);\n"
        out+="\tDFT2(&a[5],&a[7]);\n"
    else:
        out+="\tDFT2(&a[0],&a[2],prime); // dft on permutated indexes\n"
        out+="\tDFT2(&a[1],&a[3],prime);\n"
        out+="\tDFT2(&a[4],&a[6],prime);\n"
        out+="\tDFT2(&a[5],&a[7],prime);\n"

    out+="\n"
    if gfpf_bool:
        out+="\ta[3]=a[3].MulPowR(2); // twiddle\n"
        out+="\ta[5]=a[5].MulPowR(1);\n"
        out+="\ta[7]=a[7].MulPowR(3);\n"
    elif spf_bool or bpf_bool:
        out+="\ta[3]=a[3]*omegas[2];// twiddle\n"
        out+="\ta[5]=a[5]*omegas[1];\n"
        out+="\ta[7]=a[7]*omegas[3];\n"
    else:
        out+="\ta[3]=smallprimefield_multi(a[3],omegas[2],prime,R,pP);   // twiddle\n"
        out+="\ta[5]=smallprimefield_multi(a[5],omegas[1],prime,R,pP);\n"
        out+="\ta[7]=smallprimefield_multi(a[7],omegas[3],prime,R,pP);\n"

    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&a[0],&a[1]);   // dft on permutated indexes\n"
        out+="\tDFT2(&a[2],&a[3]);\n"
        out+="\tDFT2(&a[4],&a[5]);\n"
        out+="\tDFT2(&a[6],&a[7]);\n"
    else:
        out+="\tDFT2(&a[0],&a[1],prime);   // dft on permutated indexes\n"
        out+="\tDFT2(&a[2],&a[3],prime);\n"
        out+="\tDFT2(&a[4],&a[5],prime);\n"
        out+="\tDFT2(&a[6],&a[7],prime);\n"

    out+="\n"
    out+="\tswap(&a[1],&a[4]); // final permutation\n"
    out+="\tswap(&a[3],&a[6]);\n"
    out+="\n"
    out+="\treturn a; // return result\n"
    out+="}\n"
    if spf_bool or bpf_bool:
        out+="template SmallPrimeField* DFT_8<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* omegas);\n"
        out+="template BigPrimeField* DFT_8<BigPrimeField>(BigPrimeField* a,BigPrimeField* omegas);\n"
    return out

'''write dft code for basecase size 16
fieldtype options are spf, bpf, gfpf
default builds c version'''
def write_dft_16(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if (fieldtype=="bpf"):
        bpf_bool = True
    if (fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool or bpf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="template<class FiniteField>\n"
        out+="inline FiniteField* DFT_16(FiniteField* a,FiniteField* omegas){\n"
    elif gfpf_bool:
        #out+="template<class FiniteField>\n"
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="inline GeneralizedFermatPrimeField* DFT_16(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas){\n"
    else:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="\param prime\n"
        out+="\param R\n"
        out+="\param pP\n"
        out+="*/\n"
        out+="long int* DFT_16(long int* a,long int* omegas, long int prime,long long int R,long long int pP){\n"

    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&a[0],&a[8]);\n"
        out+="\tDFT2(&a[1],&a[9]);\n"
        out+="\tDFT2(&a[2],&a[10]);\n"
        out+="\tDFT2(&a[3],&a[11]);\n"
        out+="\tDFT2(&a[4],&a[12]);\n"
        out+="\tDFT2(&a[5],&a[13]);\n"
        out+="\tDFT2(&a[6],&a[14]);\n"
        out+="\tDFT2(&a[7],&a[15]);\n"
    else:
        out+="\tDFT2(&a[0],&a[8],prime);\n"
        out+="\tDFT2(&a[1],&a[9],prime);\n"
        out+="\tDFT2(&a[2],&a[10],prime);\n"
        out+="\tDFT2(&a[3],&a[11],prime);\n"
        out+="\tDFT2(&a[4],&a[12],prime);\n"
        out+="\tDFT2(&a[5],&a[13],prime);\n"
        out+="\tDFT2(&a[6],&a[14],prime);\n"
        out+="\tDFT2(&a[7],&a[15],prime);\n"

    out+="\n"
    if gfpf_bool:
        out+="\ta[12]=a[12].MulPowR(4);\n"
        out+="\ta[13]=a[13].MulPowR(4);\n"
        out+="\ta[14]=a[14].MulPowR(4);\n"
        out+="\ta[15]=a[15].MulPowR(4);\n"
    elif spf_bool or bpf_bool:
        out+="\ta[12]=a[12]*omegas[4];\n"
        out+="\ta[13]=a[13]*omegas[4];\n"
        out+="\ta[14]=a[14]*omegas[4];\n"
        out+="\ta[15]=a[15]*omegas[4];\n"
    else:
        out+="\ta[12]=smallprimefield_multi(a[12],omegas[4],prime,R,pP);\n"
        out+="\ta[13]=smallprimefield_multi(a[13],omegas[4],prime,R,pP);\n"
        out+="\ta[14]=smallprimefield_multi(a[14],omegas[4],prime,R,pP);\n"
        out+="\ta[15]=smallprimefield_multi(a[15],omegas[4],prime,R,pP);\n"

    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&a[0],&a[4]);\n"
        out+="\tDFT2(&a[1],&a[5]);\n"
        out+="\tDFT2(&a[2],&a[6]);\n"
        out+="\tDFT2(&a[3],&a[7]);\n"
        out+="\tDFT2(&a[8],&a[12]);\n"
        out+="\tDFT2(&a[9],&a[13]);\n"
        out+="\tDFT2(&a[10],&a[14]);\n"
        out+="\tDFT2(&a[11],&a[15]);\n"
    else:
        out+="\tDFT2(&a[0],&a[4],prime);\n"
        out+="\tDFT2(&a[1],&a[5],prime);\n"
        out+="\tDFT2(&a[2],&a[6],prime);\n"
        out+="\tDFT2(&a[3],&a[7],prime);\n"
        out+="\tDFT2(&a[8],&a[12],prime);\n"
        out+="\tDFT2(&a[9],&a[13],prime);\n"
        out+="\tDFT2(&a[10],&a[14],prime);\n"
        out+="\tDFT2(&a[11],&a[15],prime);\n"

    out+="\n"
    if gfpf_bool:
        out+="\ta[6]=a[6].MulPowR(4);\n"
        out+="\ta[7]=a[7].MulPowR(4);\n"
        out+="\ta[10]=a[10].MulPowR(2);\n"
        out+="\ta[11]=a[11].MulPowR(2);\n"
        out+="\ta[14]=a[14].MulPowR(6);\n"
        out+="\ta[15]=a[15].MulPowR(6); \n"
    elif spf_bool or bpf_bool:
        out+="\ta[6]=a[6]*omegas[4];\n"
        out+="\ta[7]=a[7]*omegas[4];\n"
        out+="\ta[10]=a[10]*omegas[2];\n"
        out+="\ta[11]=a[11]*omegas[2];\n"
        out+="\ta[14]=a[14]*omegas[6];\n"
        out+="\ta[15]=a[15]*omegas[6];\n"
    else:
        out+="\ta[6]=smallprimefield_multi(a[6],omegas[4],prime,R,pP);\n"
        out+="\ta[7]=smallprimefield_multi(a[7],omegas[4],prime,R,pP);\n"
        out+="\ta[10]=smallprimefield_multi(a[10],omegas[2],prime,R,pP);\n"
        out+="\ta[11]=smallprimefield_multi(a[11],omegas[2],prime,R,pP);\n"
        out+="\ta[14]=smallprimefield_multi(a[14],omegas[6],prime,R,pP);\n"
        out+="\ta[15]=smallprimefield_multi(a[15],omegas[6],prime,R,pP);\n"

    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&a[0],&a[2]);\n"
        out+="\tDFT2(&a[1],&a[3]);\n"
        out+="\tDFT2(&a[4],&a[6]);\n"
        out+="\tDFT2(&a[5],&a[7]);\n"
        out+="\tDFT2(&a[8],&a[10]);\n"
        out+="\tDFT2(&a[9],&a[11]);\n"
        out+="\tDFT2(&a[12],&a[14]);\n"
        out+="\tDFT2(&a[13],&a[15]);\n"
    else:
        out+="\tDFT2(&a[0],&a[2],prime);\n"
        out+="\tDFT2(&a[1],&a[3],prime);\n"
        out+="\tDFT2(&a[4],&a[6],prime);\n"
        out+="\tDFT2(&a[5],&a[7],prime);\n"
        out+="\tDFT2(&a[8],&a[10],prime);\n"
        out+="\tDFT2(&a[9],&a[11],prime);\n"
        out+="\tDFT2(&a[12],&a[14],prime);\n"
        out+="\tDFT2(&a[13],&a[15],prime);\n"

    out+="\n"
    if gfpf_bool:
        out+="\ta[3]=a[3].MulPowR(4);\n"
        out+="\ta[5]=a[5].MulPowR(2);\n"
        out+="\ta[7]=a[7].MulPowR(6);\n"
        out+="\ta[9]=a[9].MulPowR(1);\n"
        out+="\ta[11]=a[11].MulPowR(5);\n"
        out+="\ta[13]=a[13].MulPowR(3);\n"
        out+="\ta[15]=a[15].MulPowR(7);\n"
    elif spf_bool or bpf_bool:
        out+="\ta[3]=a[3]*omegas[4];\n"
        out+="\ta[5]=a[5]*omegas[2];\n"
        out+="\ta[7]=a[7]*omegas[6];\n"
        out+="\ta[9]=a[9]*omegas[1];\n"
        out+="\ta[11]=a[11]*omegas[5];\n"
        out+="\ta[13]=a[13]*omegas[3];\n"
        out+="\ta[15]=a[15]*omegas[7];\n"
    else:
        out+="\ta[3]=smallprimefield_multi(a[3],omegas[4],prime,R,pP);\n"
        out+="\ta[5]=smallprimefield_multi(a[5],omegas[2],prime,R,pP);\n"
        out+="\ta[7]=smallprimefield_multi(a[7],omegas[6],prime,R,pP);\n"
        out+="\ta[9]=smallprimefield_multi(a[9],omegas[1],prime,R,pP);\n"
        out+="\ta[11]=smallprimefield_multi(a[11],omegas[5],prime,R,pP);\n"
        out+="\ta[13]=smallprimefield_multi(a[13],omegas[3],prime,R,pP);\n"
        out+="\ta[15]=smallprimefield_multi(a[15],omegas[7],prime,R,pP);\n"

    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&a[0],&a[1]);\n"
        out+="\tDFT2(&a[2],&a[3]);\n"
        out+="\tDFT2(&a[4],&a[5]);\n"
        out+="\tDFT2(&a[6],&a[7]);\n"
        out+="\tDFT2(&a[8],&a[9]);\n"
        out+="\tDFT2(&a[10],&a[11]);\n"
        out+="\tDFT2(&a[12],&a[13]);\n"
        out+="\t DFT2(&a[14],&a[15]);\n"
    else:
        out+="\tDFT2(&a[0],&a[1],prime);\n"
        out+="\tDFT2(&a[2],&a[3],prime);\n"
        out+="\tDFT2(&a[4],&a[5],prime);\n"
        out+="\tDFT2(&a[6],&a[7],prime);\n"
        out+="\tDFT2(&a[8],&a[9],prime);\n"
        out+="\tDFT2(&a[10],&a[11],prime);\n"
        out+="\tDFT2(&a[12],&a[13],prime);\n"
        out+="\tDFT2(&a[14],&a[15],prime);\n"

    out+="\n"
    #SWAP
    out+="\tswap(&a[1],&a[8]);\n"
    out+="\tswap(&a[2],&a[4]);\n"
    out+="\tswap(&a[3],&a[12]);\n"
    out+="\tswap(&a[5],&a[10]);\n"
    out+="\tswap(&a[7],&a[14]);\n"
    out+="\tswap(&a[11],&a[13]);\n"
    out+="\n"
    out+="\treturn a;\n"
    out+="}\n"
    if spf_bool or bpf_bool:
        out+="template SmallPrimeField* DFT_16<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* omegas);\n"
        out+="template BigPrimeField* DFT_16<BigPrimeField>(BigPrimeField* a,BigPrimeField* omegas);\n"
    return out

'''write dft code for basecase size 32
fieldtype options are spf, bpf, gfpf
default builds c version'''
def write_dft_32(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="template<class FiniteField>\n"
        out+="inline FiniteField* DFT_32(FiniteField* A,FiniteField* omegas){\n"
    elif bpf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="template<class FiniteField>\n"
        out+="inline FiniteField* DFT_32(FiniteField* A,FiniteField* omegas){\n"
    elif gfpf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="inline GeneralizedFermatPrimeField* DFT_32(GeneralizedFermatPrimeField* A,GeneralizedFermatPrimeField* omegas){\n"
    else:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="\param prime\n"
        out+="\param R\n"
        out+="\param pP\n"
        out+="*/\n"
        out+="long int* DFT_32(long int* A,long int* omegas, long int prime,long long int R,long long int pP){\n"

    out+="\n"
    #DFT 1
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[16]);\n"
        out+="\tDFT2(&A[1],&A[17]);\n"
        out+="\tDFT2(&A[2],&A[18]);\n"
        out+="\tDFT2(&A[3],&A[19]);\n"
        out+="\tDFT2(&A[4],&A[20]);\n"
        out+="\tDFT2(&A[5],&A[21]);\n"
        out+="\tDFT2(&A[6],&A[22]);\n"
        out+="\tDFT2(&A[7],&A[23]);\n"
        out+="\tDFT2(&A[8],&A[24]);\n"
        out+="\tDFT2(&A[9],&A[25]);\n"
        out+="\tDFT2(&A[10],&A[26]);\n"
        out+="\tDFT2(&A[11],&A[27]);\n"
        out+="\tDFT2(&A[12],&A[28]);\n"
        out+="\tDFT2(&A[13],&A[29]);\n"
        out+="\tDFT2(&A[14],&A[30]);\n"
        out+="\tDFT2(&A[15],&A[31]);\n"
    else:
        out+="\tDFT2(&A[0],&A[16],prime);\n"
        out+="\tDFT2(&A[1],&A[17],prime);\n"
        out+="\tDFT2(&A[2],&A[18],prime);\n"
        out+="\tDFT2(&A[3],&A[19],prime);\n"
        out+="\tDFT2(&A[4],&A[20],prime);\n"
        out+="\tDFT2(&A[5],&A[21],prime);\n"
        out+="\tDFT2(&A[6],&A[22],prime);\n"
        out+="\tDFT2(&A[7],&A[23],prime);\n"
        out+="\tDFT2(&A[8],&A[24],prime);\n"
        out+="\tDFT2(&A[9],&A[25],prime);\n"
        out+="\tDFT2(&A[10],&A[26],prime);\n"
        out+="\tDFT2(&A[11],&A[27],prime);\n"
        out+="\tDFT2(&A[12],&A[28],prime);\n"
        out+="\tDFT2(&A[13],&A[29],prime);\n"
        out+="\tDFT2(&A[14],&A[30],prime);\n"
        out+="\tDFT2(&A[15],&A[31],prime);\n"

    out+="\n"
    #TWIDDLE 1
    if gfpf_bool:
        out+="\tA[24]=A[24].MulPowR(8);\n"
        out+="\tA[25]=A[25].MulPowR(8);\n"
        out+="\tA[26]=A[26].MulPowR(8);\n"
        out+="\tA[27]=A[27].MulPowR(8);\n"
        out+="\tA[28]=A[28].MulPowR(8);\n"
        out+="\tA[29]=A[29].MulPowR(8);\n"
        out+="\tA[30]=A[30].MulPowR(8);\n"
        out+="\tA[31]=A[31].MulPowR(8);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[24]=A[24]*omegas[8];\n"
        out+="\tA[25]=A[25]*omegas[8];\n"
        out+="\tA[26]=A[26]*omegas[8];\n"
        out+="\tA[27]=A[27]*omegas[8];\n"
        out+="\tA[28]=A[28]*omegas[8];\n"
        out+="\tA[29]=A[29]*omegas[8];\n"
        out+="\tA[30]=A[30]*omegas[8];\n"
        out+="\tA[31]=A[31]*omegas[8];\n"
    else:
        out+="\tA[24]=smallprimefield_multi(A[24],omegas[8],prime,R,pP);\n"
        out+="\tA[25]=smallprimefield_multi(A[25],omegas[8],prime,R,pP);\n"
        out+="\tA[26]=smallprimefield_multi(A[26],omegas[8],prime,R,pP);\n"
        out+="\tA[27]=smallprimefield_multi(A[27],omegas[8],prime,R,pP);\n"
        out+="\tA[28]=smallprimefield_multi(A[28],omegas[8],prime,R,pP);\n"
        out+="\tA[29]=smallprimefield_multi(A[29],omegas[8],prime,R,pP);\n"
        out+="\tA[30]=smallprimefield_multi(A[30],omegas[8],prime,R,pP);\n"
        out+="\tA[31]=smallprimefield_multi(A[31],omegas[8],prime,R,pP);\n"

    out+="\n"
    #DFT 2
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[8]);\n"
        out+="\tDFT2(&A[1],&A[9]);\n"
        out+="\tDFT2(&A[2],&A[10]);\n"
        out+="\tDFT2(&A[3],&A[11]);\n"
        out+="\tDFT2(&A[4],&A[12]);\n"
        out+="\tDFT2(&A[5],&A[13]);\n"
        out+="\tDFT2(&A[6],&A[14]);\n"
        out+="\tDFT2(&A[7],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[24]);\n"
        out+="\tDFT2(&A[17],&A[25]);\n"
        out+="\tDFT2(&A[18],&A[26]);\n"
        out+="\tDFT2(&A[19],&A[27]);\n"
        out+="\tDFT2(&A[20],&A[28]);\n"
        out+="\tDFT2(&A[21],&A[29]);\n"
        out+="\tDFT2(&A[22],&A[30]);\n"
        out+="\tDFT2(&A[23],&A[31]);\n"
    else:
        out+="\tDFT2(&A[0],&A[8],prime);\n"
        out+="\tDFT2(&A[1],&A[9],prime);\n"
        out+="\tDFT2(&A[2],&A[10],prime);\n"
        out+="\tDFT2(&A[3],&A[11],prime);\n"
        out+="\tDFT2(&A[4],&A[12],prime);\n"
        out+="\tDFT2(&A[5],&A[13],prime);\n"
        out+="\tDFT2(&A[6],&A[14],prime);\n"
        out+="\tDFT2(&A[7],&A[15],prime);\n"
        out+="\tDFT2(&A[16],&A[24],prime);\n"
        out+="\tDFT2(&A[17],&A[25],prime);\n"
        out+="\tDFT2(&A[18],&A[26],prime);\n"
        out+="\tDFT2(&A[19],&A[27],prime);\n"
        out+="\tDFT2(&A[20],&A[28],prime);\n"
        out+="\tDFT2(&A[21],&A[29],prime);\n"
        out+="\tDFT2(&A[22],&A[30],prime);\n"
        out+="\tDFT2(&A[23],&A[31],prime);\n"

    out+="\n"
    #TWIDDLE 2
    if gfpf_bool:
        out+="\tA[12]=A[12].MulPowR(8);\n"
        out+="\tA[13]=A[13].MulPowR(8);\n"
        out+="\tA[14]=A[14].MulPowR(8);\n"
        out+="\tA[15]=A[15].MulPowR(8);\n"
        out+="\tA[20]=A[20].MulPowR(4);\n"
        out+="\tA[21]=A[21].MulPowR(4);\n"
        out+="\tA[22]=A[22].MulPowR(4);\n"
        out+="\tA[23]=A[23].MulPowR(4);\n"
        out+="\tA[28]=A[28].MulPowR(12);\n"
        out+="\tA[29]=A[29].MulPowR(12);\n"
        out+="\tA[30]=A[30].MulPowR(12);\n"
        out+="\tA[31]=A[31].MulPowR(12);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[12]=A[12]*omegas[8];\n"
        out+="\tA[13]=A[13]*omegas[8];\n"
        out+="\tA[14]=A[14]*omegas[8];\n"
        out+="\tA[15]=A[15]*omegas[8];\n"
        out+="\tA[20]=A[20]*omegas[4];\n"
        out+="\tA[21]=A[21]*omegas[4];\n"
        out+="\tA[22]=A[22]*omegas[4];\n"
        out+="\tA[23]=A[23]*omegas[4];\n"
        out+="\tA[28]=A[28]*omegas[12];\n"
        out+="\tA[29]=A[29]*omegas[12];\n"
        out+="\tA[30]=A[30]*omegas[12];\n"
        out+="\tA[31]=A[31]*omegas[12];\n"
    else:
        out+="\tA[12]=smallprimefield_multi(A[12],omegas[8],prime,R,pP);\n"
        out+="\tA[13]=smallprimefield_multi(A[13],omegas[8],prime,R,pP);\n"
        out+="\tA[14]=smallprimefield_multi(A[14],omegas[8],prime,R,pP);\n"
        out+="\tA[15]=smallprimefield_multi(A[15],omegas[8],prime,R,pP);\n"
        out+="\tA[20]=smallprimefield_multi(A[20],omegas[4],prime,R,pP);\n"
        out+="\tA[21]=smallprimefield_multi(A[21],omegas[4],prime,R,pP);\n"
        out+="\tA[22]=smallprimefield_multi(A[22],omegas[4],prime,R,pP);\n"
        out+="\tA[23]=smallprimefield_multi(A[23],omegas[4],prime,R,pP);\n"
        out+="\tA[28]=smallprimefield_multi(A[28],omegas[12],prime,R,pP);\n"
        out+="\tA[29]=smallprimefield_multi(A[29],omegas[12],prime,R,pP);\n"
        out+="\tA[30]=smallprimefield_multi(A[30],omegas[12],prime,R,pP);\n"
        out+="\tA[31]=smallprimefield_multi(A[31],omegas[12],prime,R,pP);\n"

    out+="\n"
    #DFT 3
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[4]);\n"
        out+="\tDFT2(&A[1],&A[5]);\n"
        out+="\tDFT2(&A[2],&A[6]);\n"
        out+="\tDFT2(&A[3],&A[7]);\n"
        out+="\tDFT2(&A[8],&A[12]);\n"
        out+="\tDFT2(&A[9],&A[13]);\n"
        out+="\tDFT2(&A[10],&A[14]);\n"
        out+="\tDFT2(&A[11],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[20]);\n"
        out+="\tDFT2(&A[17],&A[21]);\n"
        out+="\tDFT2(&A[18],&A[22]);\n"
        out+="\tDFT2(&A[19],&A[23]);\n"
        out+="\tDFT2(&A[24],&A[28]);\n"
        out+="\tDFT2(&A[25],&A[29]);\n"
        out+="\tDFT2(&A[26],&A[30]);\n"
        out+="\tDFT2(&A[27],&A[31]);\n"
    else:
        out+="\tDFT2(&A[0],&A[4]);\n"
        out+="\tDFT2(&A[1],&A[5]);\n"
        out+="\tDFT2(&A[2],&A[6]);\n"
        out+="\tDFT2(&A[3],&A[7]);\n"
        out+="\tDFT2(&A[8],&A[12]);\n"
        out+="\tDFT2(&A[9],&A[13]);\n"
        out+="\tDFT2(&A[10],&A[14]);\n"
        out+="\tDFT2(&A[11],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[20]);\n"
        out+="\tDFT2(&A[17],&A[21]);\n"
        out+="\tDFT2(&A[18],&A[22]);\n"
        out+="\tDFT2(&A[19],&A[23]);\n"
        out+="\tDFT2(&A[24],&A[28]);\n"
        out+="\tDFT2(&A[25],&A[29]);\n"
        out+="\tDFT2(&A[26],&A[30]);\n"
        out+="\tDFT2(&A[27],&A[31]);\n"

    out+="\n"
    #TWIDDLE 3
    if gfpf_bool:
        out+="\tA[6]=A[6].MulPowR(8);\n"
        out+="\tA[7]=A[7].MulPowR(8);\n"
        out+="\tA[10]=A[10].MulPowR(4);\n"
        out+="\tA[11]=A[11].MulPowR(4);\n"
        out+="\tA[14]=A[14].MulPowR(12);\n"
        out+="\tA[15]=A[15].MulPowR(12);\n"
        out+="\tA[18]=A[18].MulPowR(2);\n"
        out+="\tA[19]=A[19].MulPowR(2);\n"
        out+="\tA[22]=A[22].MulPowR(10);\n"
        out+="\tA[23]=A[23].MulPowR(10);\n"
        out+="\tA[26]=A[26].MulPowR(6);\n"
        out+="\tA[27]=A[27].MulPowR(6);\n"
        out+="\tA[30]=A[30].MulPowR(14);\n"
        out+="\tA[31]=A[31].MulPowR(14);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[6]=A[6]*omegas[8];\n"
        out+="\tA[7]=A[7]*omegas[8];\n"
        out+="\tA[10]=A[10]*omegas[4];\n"
        out+="\tA[11]=A[11]*omegas[4];\n"
        out+="\tA[14]=A[14]*omegas[12];\n"
        out+="\tA[15]=A[15]*omegas[12];\n"
        out+="\tA[18]=A[18]*omegas[2];\n"
        out+="\tA[19]=A[19]*omegas[2];\n"
        out+="\tA[22]=A[22]*omegas[10];\n"
        out+="\tA[23]=A[23]*omegas[10];\n"
        out+="\tA[26]=A[26]*omegas[6];\n"
        out+="\tA[27]=A[27]*omegas[6];\n"
        out+="\tA[30]=A[30]*omegas[14];\n"
        out+="\tA[31]=A[31]*omegas[14];\n"
    else:
        out+="\tA[6]=smallprimefield_multi(A[6],omegas[8],prime,R,pP);\n"
        out+="\tA[7]=smallprimefield_multi(A[7],omegas[8],prime,R,pP);\n"
        out+="\tA[10]=smallprimefield_multi(A[10],omegas[4],prime,R,pP);\n"
        out+="\tA[11]=smallprimefield_multi(A[11],omegas[4],prime,R,pP);\n"
        out+="\tA[14]=smallprimefield_multi(A[14],omegas[12],prime,R,pP);\n"
        out+="\tA[15]=smallprimefield_multi(A[15],omegas[12],prime,R,pP);\n"
        out+="\tA[18]=smallprimefield_multi(A[18],omegas[2],prime,R,pP);\n"
        out+="\tA[19]=smallprimefield_multi(A[19],omegas[2],prime,R,pP);\n"
        out+="\tA[22]=smallprimefield_multi(A[22],omegas[10],prime,R,pP);\n"
        out+="\tA[23]=smallprimefield_multi(A[23],omegas[10],prime,R,pP);\n"
        out+="\tA[26]=smallprimefield_multi(A[26],omegas[6],prime,R,pP);\n"
        out+="\tA[27]=smallprimefield_multi(A[27],omegas[6],prime,R,pP);\n"
        out+="\tA[30]=smallprimefield_multi(A[30],omegas[14],prime,R,pP);\n"
        out+="\tA[31]=smallprimefield_multi(A[31],omegas[14],prime,R,pP);\n"

    out+="\n"
    #DFT 4
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[2]);\n"
        out+="\tDFT2(&A[1],&A[3]);\n"
        out+="\tDFT2(&A[4],&A[6]);\n"
        out+="\tDFT2(&A[5],&A[7]);\n"
        out+="\tDFT2(&A[8],&A[10]);\n"
        out+="\tDFT2(&A[9],&A[11]);\n"
        out+="\tDFT2(&A[12],&A[14]);\n"
        out+="\tDFT2(&A[13],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[18]);\n"
        out+="\tDFT2(&A[17],&A[19]);\n"
        out+="\tDFT2(&A[20],&A[22]);\n"
        out+="\tDFT2(&A[21],&A[23]);\n"
        out+="\tDFT2(&A[24],&A[26]);\n"
        out+="\tDFT2(&A[25],&A[27]);\n"
        out+="\tDFT2(&A[28],&A[30]);\n"
        out+="\tDFT2(&A[29],&A[31]);\n"
    else:
        out+="\tDFT2(&A[0],&A[2],prime);\n"
        out+="\tDFT2(&A[1],&A[3],prime);\n"
        out+="\tDFT2(&A[4],&A[6],prime);\n"
        out+="\tDFT2(&A[5],&A[7],prime);\n"
        out+="\tDFT2(&A[8],&A[10],prime);\n"
        out+="\tDFT2(&A[9],&A[11],prime);\n"
        out+="\tDFT2(&A[12],&A[14],prime);\n"
        out+="\tDFT2(&A[13],&A[15],prime);\n"
        out+="\tDFT2(&A[16],&A[18],prime);\n"
        out+="\tDFT2(&A[17],&A[19],prime);\n"
        out+="\tDFT2(&A[20],&A[22],prime);\n"
        out+="\tDFT2(&A[21],&A[23],prime);\n"
        out+="\tDFT2(&A[24],&A[26],prime);\n"
        out+="\tDFT2(&A[25],&A[27],prime);\n"
        out+="\tDFT2(&A[28],&A[30],prime);\n"
        out+="\tDFT2(&A[29],&A[31],prime);\n"

    out+="\n"
    #TWIDDLE 4
    if gfpf_bool:
        out+="\tA[3]=A[3].MulPowR(8);\n"
        out+="\tA[5]=A[5].MulPowR(4);\n"
        out+="\tA[7]=A[7].MulPowR(12);\n"
        out+="\tA[9]=A[9].MulPowR(2);\n"
        out+="\tA[11]=A[11].MulPowR(10);\n"
        out+="\tA[13]=A[13].MulPowR(6);\n"
        out+="\tA[15]=A[15].MulPowR(14);\n"
        out+="\tA[17]=A[17].MulPowR(1);\n"
        out+="\tA[19]=A[19].MulPowR(9);\n"
        out+="\tA[21]=A[21].MulPowR(5);\n"
        out+="\tA[23]=A[23].MulPowR(13);\n"
        out+="\tA[25]=A[25].MulPowR(3);\n"
        out+="\tA[27]=A[27].MulPowR(11);\n"
        out+="\tA[29]=A[29].MulPowR(7);\n"
        out+="\tA[31]=A[31].MulPowR(15);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[3]=A[3]*omegas[8];\n"
        out+="\tA[5]=A[5]*omegas[4];\n"
        out+="\tA[7]=A[7]*omegas[12];\n"
        out+="\tA[9]=A[9]*omegas[2];\n"
        out+="\tA[11]=A[11]*omegas[10];\n"
        out+="\tA[13]=A[13]*omegas[6];\n"
        out+="\tA[15]=A[15]*omegas[14];\n"
        out+="\tA[17]=A[17]*omegas[1];\n"
        out+="\tA[19]=A[19]*omegas[9];\n"
        out+="\tA[21]=A[21]*omegas[5];\n"
        out+="\tA[23]=A[23]*omegas[13];\n"
        out+="\tA[25]=A[25]*omegas[3];\n"
        out+="\tA[27]=A[27]*omegas[11];\n"
        out+="\tA[29]=A[29]*omegas[7];\n"
        out+="\tA[31]=A[31]*omegas[15];\n"
    else:
        out+="\tA[3]=smallprimefield_multi(A[3],omegas[8],prime,R,pP);\n"
        out+="\tA[5]=smallprimefield_multi(A[5],omegas[4],prime,R,pP);\n"
        out+="\tA[7]=smallprimefield_multi(A[7],omegas[12],prime,R,pP);\n"
        out+="\tA[9]=smallprimefield_multi(A[9],omegas[2],prime,R,pP);\n"
        out+="\tA[11]=smallprimefield_multi(A[11],omegas[10],prime,R,pP);\n"
        out+="\tA[13]=smallprimefield_multi(A[13],omegas[6],prime,R,pP);\n"
        out+="\tA[15]=smallprimefield_multi(A[15],omegas[14],prime,R,pP);\n"
        out+="\tA[17]=smallprimefield_multi(A[17],omegas[1],prime,R,pP);\n"
        out+="\tA[19]=smallprimefield_multi(A[19],omegas[9],prime,R,pP);\n"
        out+="\tA[21]=smallprimefield_multi(A[21],omegas[5],prime,R,pP);\n"
        out+="\tA[23]=smallprimefield_multi(A[23],omegas[13],prime,R,pP);\n"
        out+="\tA[25]=smallprimefield_multi(A[25],omegas[3],prime,R,pP);\n"
        out+="\tA[27]=smallprimefield_multi(A[27],omegas[11],prime,R,pP);\n"
        out+="\tA[29]=smallprimefield_multi(A[29],omegas[7],prime,R,pP);\n"
        out+="\tA[31]=smallprimefield_multi(A[31],omegas[15],prime,R,pP);\n"

    out+="\n"
    #DFT
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[1]);\n"
        out+="\tDFT2(&A[2],&A[3]);\n"
        out+="\tDFT2(&A[4],&A[5]);\n"
        out+="\tDFT2(&A[6],&A[7]);\n"
        out+="\tDFT2(&A[8],&A[9]);\n"
        out+="\tDFT2(&A[10],&A[11]);\n"
        out+="\tDFT2(&A[12],&A[13]);\n"
        out+="\tDFT2(&A[14],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[17]);\n"
        out+="\tDFT2(&A[18],&A[19]);\n"
        out+="\tDFT2(&A[20],&A[21]);\n"
        out+="\tDFT2(&A[22],&A[23]);\n"
        out+="\tDFT2(&A[24],&A[25]);\n"
        out+="\tDFT2(&A[26],&A[27]);\n"
        out+="\tDFT2(&A[28],&A[29]);\n"
        out+="\tDFT2(&A[30],&A[31]);\n"
    else:
        out+="\tDFT2(&A[0],&A[1],prime);\n"
        out+="\tDFT2(&A[2],&A[3],prime);\n"
        out+="\tDFT2(&A[4],&A[5],prime);\n"
        out+="\tDFT2(&A[6],&A[7],prime);\n"
        out+="\tDFT2(&A[8],&A[9],prime);\n"
        out+="\tDFT2(&A[10],&A[11],prime);\n"
        out+="\tDFT2(&A[12],&A[13],prime);\n"
        out+="\tDFT2(&A[14],&A[15],prime);\n"
        out+="\tDFT2(&A[16],&A[17],prime);\n"
        out+="\tDFT2(&A[18],&A[19],prime);\n"
        out+="\tDFT2(&A[20],&A[21],prime);\n"
        out+="\tDFT2(&A[22],&A[23],prime);\n"
        out+="\tDFT2(&A[24],&A[25],prime);\n"
        out+="\tDFT2(&A[26],&A[27],prime);\n"
        out+="\tDFT2(&A[28],&A[29],prime);\n"
        out+="\tDFT2(&A[30],&A[31],prime);\n"

    out+="\n"
    #swap
    out+="\tswap(&A[1],&A[16]);\n"
    out+="\tswap(&A[2],&A[8]);\n"
    out+="\tswap(&A[3],&A[24]);\n"
    out+="\tswap(&A[5],&A[20]);\n"
    out+="\tswap(&A[6],&A[12]);\n"
    out+="\tswap(&A[7],&A[28]);\n"
    out+="\tswap(&A[9],&A[18]);\n"
    out+="\tswap(&A[11],&A[26]);\n"
    out+="\tswap(&A[13],&A[22]);\n"
    out+="\tswap(&A[15],&A[30]);\n"
    out+="\tswap(&A[19],&A[25]);\n"
    out+="\tswap(&A[23],&A[29]);\n"
    out+="\n"
    out+="\treturn A;\n"
    out+="}\n"

    if spf_bool or bpf_bool:
        out+="template SmallPrimeField* DFT_32<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* omegas);\n"
        out+="template BigPrimeField* DFT_32<BigPrimeField>(BigPrimeField* a,BigPrimeField* omegas);\n"
    return out

'''write dft code for basecase size 64
fieldtype options are spf, bpf, gfpf
default builds c version'''
def write_dft_64(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="template<class FiniteField>\n"
        out+="inline FiniteField* DFT_64(FiniteField* A,FiniteField* omegas){\n"
    elif bpf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="template<class FiniteField>\n"
        out+="inline FiniteField* DFT_64(FiniteField* A,FiniteField* omegas){\n"
    elif gfpf_bool:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="*/\n"
        out+="inline GeneralizedFermatPrimeField* DFT_64(GeneralizedFermatPrimeField* A,GeneralizedFermatPrimeField* omegas){\n"
    else:
        out+="/*!\n"
        out+="\param a coefficient vector\n"
        out+="\param omegas (precomputed powers of omega)\n"
        out+="\param prime\n"
        out+="\param R\n"
        out+="\param pP\n"
        out+="*/\n"
        out+="long int* DFT_64(long int* A,long int* omegas, long int prime,long long int R,long long int pP){\n"

    #DFT 1
    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[32]);\n"
        out+="\tDFT2(&A[1],&A[33]);\n"
        out+="\tDFT2(&A[2],&A[34]);\n"
        out+="\tDFT2(&A[3],&A[35]);\n"
        out+="\tDFT2(&A[4],&A[36]);\n"
        out+="\tDFT2(&A[5],&A[37]);\n"
        out+="\tDFT2(&A[6],&A[38]);\n"
        out+="\tDFT2(&A[7],&A[39]);\n"
        out+="\tDFT2(&A[8],&A[40]);\n"
        out+="\tDFT2(&A[9],&A[41]);\n"
        out+="\tDFT2(&A[10],&A[42]);\n"
        out+="\tDFT2(&A[11],&A[43]);\n"
        out+="\tDFT2(&A[12],&A[44]);\n"
        out+="\tDFT2(&A[13],&A[45]);\n"
        out+="\tDFT2(&A[14],&A[46]);\n"
        out+="\tDFT2(&A[15],&A[47]);\n"
        out+="\tDFT2(&A[16],&A[48]);\n"
        out+="\tDFT2(&A[17],&A[49]);\n"
        out+="\tDFT2(&A[18],&A[50]);\n"
        out+="\tDFT2(&A[19],&A[51]);\n"
        out+="\tDFT2(&A[20],&A[52]);\n"
        out+="\tDFT2(&A[21],&A[53]);\n"
        out+="\tDFT2(&A[22],&A[54]);\n"
        out+="\tDFT2(&A[23],&A[55]);\n"
        out+="\tDFT2(&A[24],&A[56]);\n"
        out+="\tDFT2(&A[25],&A[57]);\n"
        out+="\tDFT2(&A[26],&A[58]);\n"
        out+="\tDFT2(&A[27],&A[59]);\n"
        out+="\tDFT2(&A[28],&A[60]);\n"
        out+="\tDFT2(&A[29],&A[61]);\n"
        out+="\tDFT2(&A[30],&A[62]);\n"
        out+="\tDFT2(&A[31],&A[63]);\n"
    else:
        out+="\tDFT2(&A[0],&A[32],prime);\n"
        out+="\tDFT2(&A[1],&A[33],prime);\n"
        out+="\tDFT2(&A[2],&A[34],prime);\n"
        out+="\tDFT2(&A[3],&A[35],prime);\n"
        out+="\tDFT2(&A[4],&A[36],prime);\n"
        out+="\tDFT2(&A[5],&A[37],prime);\n"
        out+="\tDFT2(&A[6],&A[38],prime);\n"
        out+="\tDFT2(&A[7],&A[39],prime);\n"
        out+="\tDFT2(&A[8],&A[40],prime);\n"
        out+="\tDFT2(&A[9],&A[41],prime);\n"
        out+="\tDFT2(&A[10],&A[42],prime);\n"
        out+="\tDFT2(&A[11],&A[43],prime);\n"
        out+="\tDFT2(&A[12],&A[44],prime);\n"
        out+="\tDFT2(&A[13],&A[45],prime);\n"
        out+="\tDFT2(&A[14],&A[46],prime);\n"
        out+="\tDFT2(&A[15],&A[47],prime);\n"
        out+="\tDFT2(&A[16],&A[48],prime);\n"
        out+="\tDFT2(&A[17],&A[49],prime);\n"
        out+="\tDFT2(&A[18],&A[50],prime);\n"
        out+="\tDFT2(&A[19],&A[51],prime);\n"
        out+="\tDFT2(&A[20],&A[52],prime);\n"
        out+="\tDFT2(&A[21],&A[53],prime);\n"
        out+="\tDFT2(&A[22],&A[54],prime);\n"
        out+="\tDFT2(&A[23],&A[55],prime);\n"
        out+="\tDFT2(&A[24],&A[56],prime);\n"
        out+="\tDFT2(&A[25],&A[57],prime);\n"
        out+="\tDFT2(&A[26],&A[58],prime);\n"
        out+="\tDFT2(&A[27],&A[59],prime);\n"
        out+="\tDFT2(&A[28],&A[60],prime);\n"
        out+="\tDFT2(&A[29],&A[61],prime);\n"
        out+="\tDFT2(&A[30],&A[62],prime);\n"
        out+="\tDFT2(&A[31],&A[63],prime);\n"
    #TWIDDLE
    out+="\n"
    if gfpf_bool:
        out+="\tA[48] = A[48].MulPowR(16);\n"
        out+="\tA[49] = A[49].MulPowR(16);\n"
        out+="\tA[50] = A[50].MulPowR(16);\n"
        out+="\tA[51] = A[51].MulPowR(16);\n"
        out+="\tA[52] = A[52].MulPowR(16);\n"
        out+="\tA[53] = A[53].MulPowR(16);\n"
        out+="\tA[54] = A[54].MulPowR(16);\n"
        out+="\tA[55] = A[55].MulPowR(16);\n"
        out+="\tA[56] = A[56].MulPowR(16);\n"
        out+="\tA[57] = A[57].MulPowR(16);\n"
        out+="\tA[58] = A[58].MulPowR(16);\n"
        out+="\tA[59] = A[59].MulPowR(16);\n"
        out+="\tA[60] = A[60].MulPowR(16);\n"
        out+="\tA[61] = A[61].MulPowR(16);\n"
        out+="\tA[62] = A[62].MulPowR(16);\n"
        out+="\tA[63] = A[63].MulPowR(16);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[48] = A[48]*omegas[16];\n"
        out+="\tA[49] = A[49]*omegas[16];\n"
        out+="\tA[50] = A[50]*omegas[16];\n"
        out+="\tA[51] = A[51]*omegas[16];\n"
        out+="\tA[52] = A[52]*omegas[16];\n"
        out+="\tA[53] = A[53]*omegas[16];\n"
        out+="\tA[54] = A[54]*omegas[16];\n"
        out+="\tA[55] = A[55]*omegas[16];\n"
        out+="\tA[56] = A[56]*omegas[16];\n"
        out+="\tA[57] = A[57]*omegas[16];\n"
        out+="\tA[58] = A[58]*omegas[16];\n"
        out+="\tA[59] = A[59]*omegas[16];\n"
        out+="\tA[60] = A[60]*omegas[16];\n"
        out+="\tA[61] = A[61]*omegas[16];\n"
        out+="\tA[62] = A[62]*omegas[16];\n"
        out+="\tA[63] = A[63]*omegas[16];\n"
    else:
        out+="\tA[48] = smallprimefield_multi(A[48],omegas[16],prime,R,pP);\n"
        out+="\tA[49] = smallprimefield_multi(A[49],omegas[16],prime,R,pP);\n"
        out+="\tA[50] = smallprimefield_multi(A[50],omegas[16],prime,R,pP);\n"
        out+="\tA[51] = smallprimefield_multi(A[51],omegas[16],prime,R,pP);\n"
        out+="\tA[52] = smallprimefield_multi(A[52],omegas[16],prime,R,pP);\n"
        out+="\tA[53] = smallprimefield_multi(A[53],omegas[16],prime,R,pP);\n"
        out+="\tA[54] = smallprimefield_multi(A[54],omegas[16],prime,R,pP);\n"
        out+="\tA[55] = smallprimefield_multi(A[55],omegas[16],prime,R,pP);\n"
        out+="\tA[56] = smallprimefield_multi(A[56],omegas[16],prime,R,pP);\n"
        out+="\tA[57] = smallprimefield_multi(A[57],omegas[16],prime,R,pP);\n"
        out+="\tA[58] = smallprimefield_multi(A[58],omegas[16],prime,R,pP);\n"
        out+="\tA[59] = smallprimefield_multi(A[59],omegas[16],prime,R,pP);\n"
        out+="\tA[60] = smallprimefield_multi(A[60],omegas[16],prime,R,pP);\n"
        out+="\tA[61] = smallprimefield_multi(A[61],omegas[16],prime,R,pP);\n"
        out+="\tA[62] = smallprimefield_multi(A[62],omegas[16],prime,R,pP);\n"
        out+="\tA[63] = smallprimefield_multi(A[63],omegas[16],prime,R,pP);\n"
    #DFT 2
    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[16]);\n"
        out+="\tDFT2(&A[1],&A[17]);\n"
        out+="\tDFT2(&A[2],&A[18]);\n"
        out+="\tDFT2(&A[3],&A[19]);\n"
        out+="\tDFT2(&A[4],&A[20]);\n"
        out+="\tDFT2(&A[5],&A[21]);\n"
        out+="\tDFT2(&A[6],&A[22]);\n"
        out+="\tDFT2(&A[7],&A[23]);\n"
        out+="\tDFT2(&A[8],&A[24]);\n"
        out+="\tDFT2(&A[9],&A[25]);\n"
        out+="\tDFT2(&A[10],&A[26]);\n"
        out+="\tDFT2(&A[11],&A[27]);\n"
        out+="\tDFT2(&A[12],&A[28]);\n"
        out+="\tDFT2(&A[13],&A[29]);\n"
        out+="\tDFT2(&A[14],&A[30]);\n"
        out+="\tDFT2(&A[15],&A[31]);\n"
        out+="\tDFT2(&A[32],&A[48]);\n"
        out+="\tDFT2(&A[33],&A[49]);\n"
        out+="\tDFT2(&A[34],&A[50]);\n"
        out+="\tDFT2(&A[35],&A[51]);\n"
        out+="\tDFT2(&A[36],&A[52]);\n"
        out+="\tDFT2(&A[37],&A[53]);\n"
        out+="\tDFT2(&A[38],&A[54]);\n"
        out+="\tDFT2(&A[39],&A[55]);\n"
        out+="\tDFT2(&A[40],&A[56]);\n"
        out+="\tDFT2(&A[41],&A[57]);\n"
        out+="\tDFT2(&A[42],&A[58]);\n"
        out+="\tDFT2(&A[43],&A[59]);\n"
        out+="\tDFT2(&A[44],&A[60]);\n"
        out+="\tDFT2(&A[45],&A[61]);\n"
        out+="\tDFT2(&A[46],&A[62]);\n"
        out+="\tDFT2(&A[47],&A[63]);\n"
    else:
        out+="\tDFT2(&A[0],&A[16],prime);\n"
        out+="\tDFT2(&A[1],&A[17],prime);\n"
        out+="\tDFT2(&A[2],&A[18],prime);\n"
        out+="\tDFT2(&A[3],&A[19],prime);\n"
        out+="\tDFT2(&A[4],&A[20],prime);\n"
        out+="\tDFT2(&A[5],&A[21],prime);\n"
        out+="\tDFT2(&A[6],&A[22],prime);\n"
        out+="\tDFT2(&A[7],&A[23],prime);\n"
        out+="\tDFT2(&A[8],&A[24],prime);\n"
        out+="\tDFT2(&A[9],&A[25],prime);\n"
        out+="\tDFT2(&A[10],&A[26],prime);\n"
        out+="\tDFT2(&A[11],&A[27],prime);\n"
        out+="\tDFT2(&A[12],&A[28],prime);\n"
        out+="\tDFT2(&A[13],&A[29],prime);\n"
        out+="\tDFT2(&A[14],&A[30],prime);\n"
        out+="\tDFT2(&A[15],&A[31],prime);\n"
        out+="\tDFT2(&A[32],&A[48],prime);\n"
        out+="\tDFT2(&A[33],&A[49],prime);\n"
        out+="\tDFT2(&A[34],&A[50],prime);\n"
        out+="\tDFT2(&A[35],&A[51],prime);\n"
        out+="\tDFT2(&A[36],&A[52],prime);\n"
        out+="\tDFT2(&A[37],&A[53],prime);\n"
        out+="\tDFT2(&A[38],&A[54],prime);\n"
        out+="\tDFT2(&A[39],&A[55],prime);\n"
        out+="\tDFT2(&A[40],&A[56],prime);\n"
        out+="\tDFT2(&A[41],&A[57],prime);\n"
        out+="\tDFT2(&A[42],&A[58],prime);\n"
        out+="\tDFT2(&A[43],&A[59],prime);\n"
        out+="\tDFT2(&A[44],&A[60],prime);\n"
        out+="\tDFT2(&A[45],&A[61],prime);\n"
        out+="\tDFT2(&A[46],&A[62],prime);\n"
        out+="\tDFT2(&A[47],&A[63],prime);\n"
    #TWIDDLE
    out+="\n"
    if gfpf_bool:
        out+="\tA[24] = A[24].MulPowR(16);\n"
        out+="\tA[25] = A[25].MulPowR(16);\n"
        out+="\tA[26] = A[26].MulPowR(16);\n"
        out+="\tA[27] = A[27].MulPowR(16);\n"
        out+="\tA[28] = A[28].MulPowR(16);\n"
        out+="\tA[29] = A[29].MulPowR(16);\n"
        out+="\tA[30] = A[30].MulPowR(16);\n"
        out+="\tA[31] = A[31].MulPowR(16);\n"
        out+="\tA[40] = A[40].MulPowR(8);\n"
        out+="\tA[41] = A[41].MulPowR(8);\n"
        out+="\tA[42] = A[42].MulPowR(8);\n"
        out+="\tA[43] = A[43].MulPowR(8);\n"
        out+="\tA[44] = A[44].MulPowR(8);\n"
        out+="\tA[45] = A[45].MulPowR(8);\n"
        out+="\tA[46] = A[46].MulPowR(8);\n"
        out+="\tA[47] = A[47].MulPowR(8);\n"
        out+="\tA[56] = A[56].MulPowR(24);\n"
        out+="\tA[57] = A[57].MulPowR(24);\n"
        out+="\tA[58] = A[58].MulPowR(24);\n"
        out+="\tA[59] = A[59].MulPowR(24);\n"
        out+="\tA[60] = A[60].MulPowR(24);\n"
        out+="\tA[61] = A[61].MulPowR(24);\n"
        out+="\tA[62] = A[62].MulPowR(24);\n"
        out+="\tA[63] = A[63].MulPowR(24);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[24] = A[24]*omegas[16];\n"
        out+="\tA[25] = A[25]*omegas[16];\n"
        out+="\tA[26] = A[26]*omegas[16];\n"
        out+="\tA[27] = A[27]*omegas[16];\n"
        out+="\tA[28] = A[28]*omegas[16];\n"
        out+="\tA[29] = A[29]*omegas[16];\n"
        out+="\tA[30] = A[30]*omegas[16];\n"
        out+="\tA[31] = A[31]*omegas[16];\n"
        out+="\tA[40] = A[40]*omegas[8];\n"
        out+="\tA[41] = A[41]*omegas[8];\n"
        out+="\tA[42] = A[42]*omegas[8];\n"
        out+="\tA[43] = A[43]*omegas[8];\n"
        out+="\tA[44] = A[44]*omegas[8];\n"
        out+="\tA[45] = A[45]*omegas[8];\n"
        out+="\tA[46] = A[46]*omegas[8];\n"
        out+="\tA[47] = A[47]*omegas[8];\n"
        out+="\tA[56] = A[56]*omegas[24];\n"
        out+="\tA[57] = A[57]*omegas[24];\n"
        out+="\tA[58] = A[58]*omegas[24];\n"
        out+="\tA[59] = A[59]*omegas[24];\n"
        out+="\tA[60] = A[60]*omegas[24];\n"
        out+="\tA[61] = A[61]*omegas[24];\n"
        out+="\tA[62] = A[62]*omegas[24];\n"
        out+="\tA[63] = A[63]*omegas[24];\n"
    else:
        out+="\tA[24] = smallprimefield_multi(A[24],omegas[16],prime,R,pP);\n"
        out+="\tA[25] = smallprimefield_multi(A[25],omegas[16],prime,R,pP);\n"
        out+="\tA[26] = smallprimefield_multi(A[26],omegas[16],prime,R,pP);\n"
        out+="\tA[27] = smallprimefield_multi(A[27],omegas[16],prime,R,pP);\n"
        out+="\tA[28] = smallprimefield_multi(A[28],omegas[16],prime,R,pP);\n"
        out+="\tA[29] = smallprimefield_multi(A[29],omegas[16],prime,R,pP);\n"
        out+="\tA[30] = smallprimefield_multi(A[30],omegas[16],prime,R,pP);\n"
        out+="\tA[31] = smallprimefield_multi(A[31],omegas[16],prime,R,pP);\n"
        out+="\tA[40] = smallprimefield_multi(A[40],omegas[8],prime,R,pP);\n"
        out+="\tA[41] = smallprimefield_multi(A[41],omegas[8],prime,R,pP);\n"
        out+="\tA[42] = smallprimefield_multi(A[42],omegas[8],prime,R,pP);\n"
        out+="\tA[43] = smallprimefield_multi(A[43],omegas[8],prime,R,pP);\n"
        out+="\tA[44] = smallprimefield_multi(A[44],omegas[8],prime,R,pP);\n"
        out+="\tA[45] = smallprimefield_multi(A[45],omegas[8],prime,R,pP);\n"
        out+="\tA[46] = smallprimefield_multi(A[46],omegas[8],prime,R,pP);\n"
        out+="\tA[47] = smallprimefield_multi(A[47],omegas[8],prime,R,pP);\n"
        out+="\tA[56] = smallprimefield_multi(A[56],omegas[24],prime,R,pP);\n"
        out+="\tA[57] = smallprimefield_multi(A[57],omegas[24],prime,R,pP);\n"
        out+="\tA[58] = smallprimefield_multi(A[58],omegas[24],prime,R,pP);\n"
        out+="\tA[59] = smallprimefield_multi(A[59],omegas[24],prime,R,pP);\n"
        out+="\tA[60] = smallprimefield_multi(A[60],omegas[24],prime,R,pP);\n"
        out+="\tA[61] = smallprimefield_multi(A[61],omegas[24],prime,R,pP);\n"
        out+="\tA[62] = smallprimefield_multi(A[62],omegas[24],prime,R,pP);\n"
        out+="\tA[63] = smallprimefield_multi(A[63],omegas[24],prime,R,pP);\n"
    #DFT 3
    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[8]);\n"
        out+="\tDFT2(&A[1],&A[9]);\n"
        out+="\tDFT2(&A[2],&A[10]);\n"
        out+="\tDFT2(&A[3],&A[11]);\n"
        out+="\tDFT2(&A[4],&A[12]);\n"
        out+="\tDFT2(&A[5],&A[13]);\n"
        out+="\tDFT2(&A[6],&A[14]);\n"
        out+="\tDFT2(&A[7],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[24]);\n"
        out+="\tDFT2(&A[17],&A[25]);\n"
        out+="\tDFT2(&A[18],&A[26]);\n"
        out+="\tDFT2(&A[19],&A[27]);\n"
        out+="\tDFT2(&A[20],&A[28]);\n"
        out+="\tDFT2(&A[21],&A[29]);\n"
        out+="\tDFT2(&A[22],&A[30]);\n"
        out+="\tDFT2(&A[23],&A[31]);\n"
        out+="\tDFT2(&A[32],&A[40]);\n"
        out+="\tDFT2(&A[33],&A[41]);\n"
        out+="\tDFT2(&A[34],&A[42]);\n"
        out+="\tDFT2(&A[35],&A[43]);\n"
        out+="\tDFT2(&A[36],&A[44]);\n"
        out+="\tDFT2(&A[37],&A[45]);\n"
        out+="\tDFT2(&A[38],&A[46]);\n"
        out+="\tDFT2(&A[39],&A[47]);\n"
        out+="\tDFT2(&A[48],&A[56]);\n"
        out+="\tDFT2(&A[49],&A[57]);\n"
        out+="\tDFT2(&A[50],&A[58]);\n"
        out+="\tDFT2(&A[51],&A[59]);\n"
        out+="\tDFT2(&A[52],&A[60]);\n"
        out+="\tDFT2(&A[53],&A[61]);\n"
        out+="\tDFT2(&A[54],&A[62]);\n"
        out+="\tDFT2(&A[55],&A[63]);\n"
    else:
        out+="\tDFT2(&A[0],&A[8],prime);\n"
        out+="\tDFT2(&A[1],&A[9],prime);\n"
        out+="\tDFT2(&A[2],&A[10],prime);\n"
        out+="\tDFT2(&A[3],&A[11],prime);\n"
        out+="\tDFT2(&A[4],&A[12],prime);\n"
        out+="\tDFT2(&A[5],&A[13],prime);\n"
        out+="\tDFT2(&A[6],&A[14],prime);\n"
        out+="\tDFT2(&A[7],&A[15],prime);\n"
        out+="\tDFT2(&A[16],&A[24],prime);\n"
        out+="\tDFT2(&A[17],&A[25],prime);\n"
        out+="\tDFT2(&A[18],&A[26],prime);\n"
        out+="\tDFT2(&A[19],&A[27],prime);\n"
        out+="\tDFT2(&A[20],&A[28],prime);\n"
        out+="\tDFT2(&A[21],&A[29],prime);\n"
        out+="\tDFT2(&A[22],&A[30],prime);\n"
        out+="\tDFT2(&A[23],&A[31],prime);\n"
        out+="\tDFT2(&A[32],&A[40],prime);\n"
        out+="\tDFT2(&A[33],&A[41],prime);\n"
        out+="\tDFT2(&A[34],&A[42],prime);\n"
        out+="\tDFT2(&A[35],&A[43],prime);\n"
        out+="\tDFT2(&A[36],&A[44],prime);\n"
        out+="\tDFT2(&A[37],&A[45],prime);\n"
        out+="\tDFT2(&A[38],&A[46],prime);\n"
        out+="\tDFT2(&A[39],&A[47],prime);\n"
        out+="\tDFT2(&A[48],&A[56],prime);\n"
        out+="\tDFT2(&A[49],&A[57],prime);\n"
        out+="\tDFT2(&A[50],&A[58],prime);\n"
        out+="\tDFT2(&A[51],&A[59],prime);\n"
        out+="\tDFT2(&A[52],&A[60],prime);\n"
        out+="\tDFT2(&A[53],&A[61],prime);\n"
        out+="\tDFT2(&A[54],&A[62],prime);\n"
        out+="\tDFT2(&A[55],&A[63],prime);\n"
    #TWIDDLE
    out+="\n"
    if gfpf_bool:
        out+="\tA[12] = A[12].MulPowR(16);\n"
        out+="\tA[13] = A[13].MulPowR(16);\n"
        out+="\tA[14] = A[14].MulPowR(16);\n"
        out+="\tA[15] = A[15].MulPowR(16);\n"
        out+="\tA[20] = A[20].MulPowR(8);\n"
        out+="\tA[21] = A[21].MulPowR(8);\n"
        out+="\tA[22] = A[22].MulPowR(8);\n"
        out+="\tA[23] = A[23].MulPowR(8);\n"
        out+="\tA[28] = A[28].MulPowR(24);\n"
        out+="\tA[29] = A[29].MulPowR(24);\n"
        out+="\tA[30] = A[30].MulPowR(24);\n"
        out+="\tA[31] = A[31].MulPowR(24);\n"
        out+="\tA[36] = A[36].MulPowR(4);\n"
        out+="\tA[37] = A[37].MulPowR(4);\n"
        out+="\tA[38] = A[38].MulPowR(4);\n"
        out+="\tA[39] = A[39].MulPowR(4);\n"
        out+="\tA[44] = A[44].MulPowR(20);\n"
        out+="\tA[45] = A[45].MulPowR(20);\n"
        out+="\tA[46] = A[46].MulPowR(20);\n"
        out+="\tA[47] = A[47].MulPowR(20);\n"
        out+="\tA[52] = A[52].MulPowR(12);\n"
        out+="\tA[53] = A[53].MulPowR(12);\n"
        out+="\tA[54] = A[54].MulPowR(12);\n"
        out+="\tA[55] = A[55].MulPowR(12);\n"
        out+="\tA[60] = A[60].MulPowR(28);\n"
        out+="\tA[61] = A[61].MulPowR(28);\n"
        out+="\tA[62] = A[62].MulPowR(28);\n"
        out+="\tA[63] = A[63].MulPowR(28);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[12] = A[12]*omegas[16];\n"
        out+="\tA[13] = A[13]*omegas[16];\n"
        out+="\tA[14] = A[14]*omegas[16];\n"
        out+="\tA[15] = A[15]*omegas[16];\n"
        out+="\tA[20] = A[20]*omegas[8];\n"
        out+="\tA[21] = A[21]*omegas[8];\n"
        out+="\tA[22] = A[22]*omegas[8];\n"
        out+="\tA[23] = A[23]*omegas[8];\n"
        out+="\tA[28] = A[28]*omegas[24];\n"
        out+="\tA[29] = A[29]*omegas[24];\n"
        out+="\tA[30] = A[30]*omegas[24];\n"
        out+="\tA[31] = A[31]*omegas[24];\n"
        out+="\tA[36] = A[36]*omegas[4];\n"
        out+="\tA[37] = A[37]*omegas[4];\n"
        out+="\tA[38] = A[38]*omegas[4];\n"
        out+="\tA[39] = A[39]*omegas[4];\n"
        out+="\tA[44] = A[44]*omegas[20];\n"
        out+="\tA[45] = A[45]*omegas[20];\n"
        out+="\tA[46] = A[46]*omegas[20];\n"
        out+="\tA[47] = A[47]*omegas[20];\n"
        out+="\tA[52] = A[52]*omegas[12];\n"
        out+="\tA[53] = A[53]*omegas[12];\n"
        out+="\tA[54] = A[54]*omegas[12];\n"
        out+="\tA[55] = A[55]*omegas[12];\n"
        out+="\tA[60] = A[60]*omegas[28];\n"
        out+="\tA[61] = A[61]*omegas[28];\n"
        out+="\tA[62] = A[62]*omegas[28];\n"
        out+="\tA[63] = A[63]*omegas[28];\n"
    else:
        out+="\tA[12] = smallprimefield_multi(A[12],omegas[16],prime,R,pP);\n"
        out+="\tA[13] = smallprimefield_multi(A[13],omegas[16],prime,R,pP);\n"
        out+="\tA[14] = smallprimefield_multi(A[14],omegas[16],prime,R,pP);\n"
        out+="\tA[15] = smallprimefield_multi(A[15],omegas[16],prime,R,pP);\n"
        out+="\tA[20] = smallprimefield_multi(A[20],omegas[8],prime,R,pP);\n"
        out+="\tA[21] = smallprimefield_multi(A[21],omegas[8],prime,R,pP);\n"
        out+="\tA[22] = smallprimefield_multi(A[22],omegas[8],prime,R,pP);\n"
        out+="\tA[23] = smallprimefield_multi(A[23],omegas[8],prime,R,pP);\n"
        out+="\tA[28] = smallprimefield_multi(A[28],omegas[24],prime,R,pP);\n"
        out+="\tA[29] = smallprimefield_multi(A[29],omegas[24],prime,R,pP);\n"
        out+="\tA[30] = smallprimefield_multi(A[30],omegas[24],prime,R,pP);\n"
        out+="\tA[31] = smallprimefield_multi(A[31],omegas[24],prime,R,pP);\n"
        out+="\tA[36] = smallprimefield_multi(A[36],omegas[4],prime,R,pP);\n"
        out+="\tA[37] = smallprimefield_multi(A[37],omegas[4],prime,R,pP);\n"
        out+="\tA[38] = smallprimefield_multi(A[38],omegas[4],prime,R,pP);\n"
        out+="\tA[39] = smallprimefield_multi(A[39],omegas[4],prime,R,pP);\n"
        out+="\tA[44] = smallprimefield_multi(A[44],omegas[20],prime,R,pP);\n"
        out+="\tA[45] = smallprimefield_multi(A[45],omegas[20],prime,R,pP);\n"
        out+="\tA[46] = smallprimefield_multi(A[46],omegas[20],prime,R,pP);\n"
        out+="\tA[47] = smallprimefield_multi(A[47],omegas[20],prime,R,pP);\n"
        out+="\tA[52] = smallprimefield_multi(A[52],omegas[12],prime,R,pP);\n"
        out+="\tA[53] = smallprimefield_multi(A[53],omegas[12],prime,R,pP);\n"
        out+="\tA[54] = smallprimefield_multi(A[54],omegas[12],prime,R,pP);\n"
        out+="\tA[55] = smallprimefield_multi(A[55],omegas[12],prime,R,pP);\n"
        out+="\tA[60] = smallprimefield_multi(A[60],omegas[28],prime,R,pP);\n"
        out+="\tA[61] = smallprimefield_multi(A[61],omegas[28],prime,R,pP);\n"
        out+="\tA[62] = smallprimefield_multi(A[62],omegas[28],prime,R,pP);\n"
        out+="\tA[63] = smallprimefield_multi(A[63],omegas[28],prime,R,pP);\n"
        #if idx%4==0 idx+=4
    #DFT 4
    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[4]);\n"
        out+="\tDFT2(&A[1],&A[5]);\n"
        out+="\tDFT2(&A[2],&A[6]);\n"
        out+="\tDFT2(&A[3],&A[7]);\n"
        out+="\tDFT2(&A[8],&A[12]);\n"
        out+="\tDFT2(&A[9],&A[13]);\n"
        out+="\tDFT2(&A[10],&A[14]);\n"
        out+="\tDFT2(&A[11],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[20]);\n"
        out+="\tDFT2(&A[17],&A[21]);\n"
        out+="\tDFT2(&A[18],&A[22]);\n"
        out+="\tDFT2(&A[19],&A[23]);\n"
        out+="\tDFT2(&A[24],&A[28]);\n"
        out+="\tDFT2(&A[25],&A[29]);\n"
        out+="\tDFT2(&A[26],&A[30]);\n"
        out+="\tDFT2(&A[27],&A[31]);\n"
        out+="\tDFT2(&A[32],&A[36]);\n"
        out+="\tDFT2(&A[33],&A[37]);\n"
        out+="\tDFT2(&A[34],&A[38]);\n"
        out+="\tDFT2(&A[35],&A[39]);\n"
        out+="\tDFT2(&A[40],&A[44]);\n"
        out+="\tDFT2(&A[41],&A[45]);\n"
        out+="\tDFT2(&A[42],&A[46]);\n"
        out+="\tDFT2(&A[43],&A[47]);\n"
        out+="\tDFT2(&A[48],&A[52]);\n"
        out+="\tDFT2(&A[49],&A[53]);\n"
        out+="\tDFT2(&A[50],&A[54]);\n"
        out+="\tDFT2(&A[51],&A[55]);\n"
        out+="\tDFT2(&A[56],&A[60]);\n"
        out+="\tDFT2(&A[57],&A[61]);\n"
        out+="\tDFT2(&A[58],&A[62]);\n"
        out+="\tDFT2(&A[59],&A[63]);\n"
    else:
        out+="\tDFT2(&A[0],&A[4],prime);\n"
        out+="\tDFT2(&A[1],&A[5],prime);\n"
        out+="\tDFT2(&A[2],&A[6],prime);\n"
        out+="\tDFT2(&A[3],&A[7],prime);\n"
        out+="\tDFT2(&A[8],&A[12],prime);\n"
        out+="\tDFT2(&A[9],&A[13],prime);\n"
        out+="\tDFT2(&A[10],&A[14],prime);\n"
        out+="\tDFT2(&A[11],&A[15],prime);\n"
        out+="\tDFT2(&A[16],&A[20],prime);\n"
        out+="\tDFT2(&A[17],&A[21],prime);\n"
        out+="\tDFT2(&A[18],&A[22],prime);\n"
        out+="\tDFT2(&A[19],&A[23],prime);\n"
        out+="\tDFT2(&A[24],&A[28],prime);\n"
        out+="\tDFT2(&A[25],&A[29],prime);\n"
        out+="\tDFT2(&A[26],&A[30],prime);\n"
        out+="\tDFT2(&A[27],&A[31],prime);\n"
        out+="\tDFT2(&A[32],&A[36],prime);\n"
        out+="\tDFT2(&A[33],&A[37],prime);\n"
        out+="\tDFT2(&A[34],&A[38],prime);\n"
        out+="\tDFT2(&A[35],&A[39],prime);\n"
        out+="\tDFT2(&A[40],&A[44],prime);\n"
        out+="\tDFT2(&A[41],&A[45],prime);\n"
        out+="\tDFT2(&A[42],&A[46],prime);\n"
        out+="\tDFT2(&A[43],&A[47],prime);\n"
        out+="\tDFT2(&A[48],&A[52],prime);\n"
        out+="\tDFT2(&A[49],&A[53],prime);\n"
        out+="\tDFT2(&A[50],&A[54],prime);\n"
        out+="\tDFT2(&A[51],&A[55],prime);\n"
        out+="\tDFT2(&A[56],&A[60],prime);\n"
        out+="\tDFT2(&A[57],&A[61],prime);\n"
        out+="\tDFT2(&A[58],&A[62],prime);\n"
        out+="\tDFT2(&A[59],&A[63],prime);\n"
    #TWIDDLE
    out+="\n"
    if gfpf_bool:
        out+="\tA[6] = A[6].MulPowR(16);\n"
        out+="\tA[7] = A[7].MulPowR(16);\n"
        out+="\tA[10] = A[10].MulPowR(8);\n"
        out+="\tA[11] = A[11].MulPowR(8);\n"
        out+="\tA[14] = A[14].MulPowR(24);\n"
        out+="\tA[15] = A[15].MulPowR(24);\n"
        out+="\tA[18] = A[18].MulPowR(4);\n"
        out+="\tA[19] = A[19].MulPowR(4);\n"
        out+="\tA[22] = A[22].MulPowR(20);\n"
        out+="\tA[23] = A[23].MulPowR(20);\n"
        out+="\tA[26] = A[26].MulPowR(12);\n"
        out+="\tA[27] = A[27].MulPowR(12);\n"
        out+="\tA[30] = A[30].MulPowR(28);\n"
        out+="\tA[31] = A[31].MulPowR(28);\n"
        out+="\tA[34] = A[34].MulPowR(2);\n"
        out+="\tA[35] = A[35].MulPowR(2);\n"
        out+="\tA[38] = A[38].MulPowR(18);\n"
        out+="\tA[39] = A[39].MulPowR(18);\n"
        out+="\tA[42] = A[42].MulPowR(10);\n"
        out+="\tA[43] = A[43].MulPowR(10);\n"
        out+="\tA[46] = A[46].MulPowR(26);\n"
        out+="\tA[47] = A[47].MulPowR(26);\n"
        out+="\tA[50] = A[50].MulPowR(6);\n"
        out+="\tA[51] = A[51].MulPowR(6);\n"
        out+="\tA[54] = A[54].MulPowR(22);\n"
        out+="\tA[55] = A[55].MulPowR(22);\n"
        out+="\tA[58] = A[58].MulPowR(14);\n"
        out+="\tA[59] = A[59].MulPowR(14);\n"
        out+="\tA[62] = A[62].MulPowR(30);\n"
        out+="\tA[63] = A[63].MulPowR(30);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[6] = A[6]*omegas[16];\n"
        out+="\tA[7] = A[7]*omegas[16];\n"
        out+="\tA[10] = A[10]*omegas[8];\n"
        out+="\tA[11] = A[11]*omegas[8];\n"
        out+="\tA[14] = A[14]*omegas[24];\n"
        out+="\tA[15] = A[15]*omegas[24];\n"
        out+="\tA[18] = A[18]*omegas[4];\n"
        out+="\tA[19] = A[19]*omegas[4];\n"
        out+="\tA[22] = A[22]*omegas[20];\n"
        out+="\tA[23] = A[23]*omegas[20];\n"
        out+="\tA[26] = A[26]*omegas[12];\n"
        out+="\tA[27] = A[27]*omegas[12];\n"
        out+="\tA[30] = A[30]*omegas[28];\n"
        out+="\tA[31] = A[31]*omegas[28];\n"
        out+="\tA[34] = A[34]*omegas[2];\n"
        out+="\tA[35] = A[35]*omegas[2];\n"
        out+="\tA[38] = A[38]*omegas[18];\n"
        out+="\tA[39] = A[39]*omegas[18];\n"
        out+="\tA[42] = A[42]*omegas[10];\n"
        out+="\tA[43] = A[43]*omegas[10];\n"
        out+="\tA[46] = A[46]*omegas[26];\n"
        out+="\tA[47] = A[47]*omegas[26];\n"
        out+="\tA[50] = A[50]*omegas[6];\n"
        out+="\tA[51] = A[51]*omegas[6];\n"
        out+="\tA[54] = A[54]*omegas[22];\n"
        out+="\tA[55] = A[55]*omegas[22];\n"
        out+="\tA[58] = A[58]*omegas[14];\n"
        out+="\tA[59] = A[59]*omegas[14];\n"
        out+="\tA[62] = A[62]*omegas[30];\n"
        out+="\tA[63] = A[63]*omegas[30];\n"
    else:
        out+="\tA[6] = smallprimefield_multi(A[6],omegas[16],prime,R,pP);\n"
        out+="\tA[7] = smallprimefield_multi(A[7],omegas[16],prime,R,pP);\n"
        out+="\tA[10] = smallprimefield_multi(A[10],omegas[8],prime,R,pP);\n"
        out+="\tA[11] = smallprimefield_multi(A[11],omegas[8],prime,R,pP);\n"
        out+="\tA[14] = smallprimefield_multi(A[14],omegas[24],prime,R,pP);\n"
        out+="\tA[15] = smallprimefield_multi(A[15],omegas[24],prime,R,pP);\n"
        out+="\tA[18] = smallprimefield_multi(A[18],omegas[4],prime,R,pP);\n"
        out+="\tA[19] = smallprimefield_multi(A[19],omegas[4],prime,R,pP);\n"
        out+="\tA[22] = smallprimefield_multi(A[22],omegas[20],prime,R,pP);\n"
        out+="\tA[23] = smallprimefield_multi(A[23],omegas[20],prime,R,pP);\n"
        out+="\tA[26] = smallprimefield_multi(A[26],omegas[12],prime,R,pP);\n"
        out+="\tA[27] = smallprimefield_multi(A[27],omegas[12],prime,R,pP);\n"
        out+="\tA[30] = smallprimefield_multi(A[30],omegas[28],prime,R,pP);\n"
        out+="\tA[31] = smallprimefield_multi(A[31],omegas[28],prime,R,pP);\n"
        out+="\tA[34] = smallprimefield_multi(A[34],omegas[2],prime,R,pP);\n"
        out+="\tA[35] = smallprimefield_multi(A[35],omegas[2],prime,R,pP);\n"
        out+="\tA[38] = smallprimefield_multi(A[38],omegas[18],prime,R,pP);\n"
        out+="\tA[39] = smallprimefield_multi(A[39],omegas[18],prime,R,pP);\n"
        out+="\tA[42] = smallprimefield_multi(A[42],omegas[10],prime,R,pP);\n"
        out+="\tA[43] = smallprimefield_multi(A[43],omegas[10],prime,R,pP);\n"
        out+="\tA[46] = smallprimefield_multi(A[46],omegas[26],prime,R,pP);\n"
        out+="\tA[47] = smallprimefield_multi(A[47],omegas[26],prime,R,pP);\n"
        out+="\tA[50] = smallprimefield_multi(A[50],omegas[6],prime,R,pP);\n"
        out+="\tA[51] = smallprimefield_multi(A[51],omegas[6],prime,R,pP);\n"
        out+="\tA[54] = smallprimefield_multi(A[54],omegas[22],prime,R,pP);\n"
        out+="\tA[55] = smallprimefield_multi(A[55],omegas[22],prime,R,pP);\n"
        out+="\tA[58] = smallprimefield_multi(A[58],omegas[14],prime,R,pP);\n"
        out+="\tA[59] = smallprimefield_multi(A[59],omegas[14],prime,R,pP);\n"
        out+="\tA[62] = smallprimefield_multi(A[62],omegas[30],prime,R,pP);\n"
        out+="\tA[63] = smallprimefield_multi(A[63],omegas[30],prime,R,pP);\n"
    #DFT 5
    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[2]);\n"
        out+="\tDFT2(&A[1],&A[3]);\n"
        out+="\tDFT2(&A[4],&A[6]);\n"
        out+="\tDFT2(&A[5],&A[7]);\n"
        out+="\tDFT2(&A[8],&A[10]);\n"
        out+="\tDFT2(&A[9],&A[11]);\n"
        out+="\tDFT2(&A[12],&A[14]);\n"
        out+="\tDFT2(&A[13],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[18]);\n"
        out+="\tDFT2(&A[17],&A[19]);\n"
        out+="\tDFT2(&A[20],&A[22]);\n"
        out+="\tDFT2(&A[21],&A[23]);\n"
        out+="\tDFT2(&A[24],&A[26]);\n"
        out+="\tDFT2(&A[25],&A[27]);\n"
        out+="\tDFT2(&A[28],&A[30]);\n"
        out+="\tDFT2(&A[29],&A[31]);\n"
        out+="\tDFT2(&A[32],&A[34]);\n"
        out+="\tDFT2(&A[33],&A[35]);\n"
        out+="\tDFT2(&A[36],&A[38]);\n"
        out+="\tDFT2(&A[37],&A[39]);\n"
        out+="\tDFT2(&A[40],&A[42]);\n"
        out+="\tDFT2(&A[41],&A[43]);\n"
        out+="\tDFT2(&A[44],&A[46]);\n"
        out+="\tDFT2(&A[45],&A[47]);\n"
        out+="\tDFT2(&A[48],&A[50]);\n"
        out+="\tDFT2(&A[49],&A[51]);\n"
        out+="\tDFT2(&A[52],&A[54]);\n"
        out+="\tDFT2(&A[53],&A[55]);\n"
        out+="\tDFT2(&A[56],&A[58]);\n"
        out+="\tDFT2(&A[57],&A[59]);\n"
        out+="\tDFT2(&A[60],&A[62]);\n"
        out+="\tDFT2(&A[61],&A[63]);\n"
    else:
        out+="\tDFT2(&A[0],&A[2],prime);\n"
        out+="\tDFT2(&A[1],&A[3],prime);\n"
        out+="\tDFT2(&A[4],&A[6],prime);\n"
        out+="\tDFT2(&A[5],&A[7],prime);\n"
        out+="\tDFT2(&A[8],&A[10],prime);\n"
        out+="\tDFT2(&A[9],&A[11],prime);\n"
        out+="\tDFT2(&A[12],&A[14],prime);\n"
        out+="\tDFT2(&A[13],&A[15],prime);\n"
        out+="\tDFT2(&A[16],&A[18],prime);\n"
        out+="\tDFT2(&A[17],&A[19],prime);\n"
        out+="\tDFT2(&A[20],&A[22],prime);\n"
        out+="\tDFT2(&A[21],&A[23],prime);\n"
        out+="\tDFT2(&A[24],&A[26],prime);\n"
        out+="\tDFT2(&A[25],&A[27],prime);\n"
        out+="\tDFT2(&A[28],&A[30],prime);\n"
        out+="\tDFT2(&A[29],&A[31],prime);\n"
        out+="\tDFT2(&A[32],&A[34],prime);\n"
        out+="\tDFT2(&A[33],&A[35],prime);\n"
        out+="\tDFT2(&A[36],&A[38],prime);\n"
        out+="\tDFT2(&A[37],&A[39],prime);\n"
        out+="\tDFT2(&A[40],&A[42],prime);\n"
        out+="\tDFT2(&A[41],&A[43],prime);\n"
        out+="\tDFT2(&A[44],&A[46],prime);\n"
        out+="\tDFT2(&A[45],&A[47],prime);\n"
        out+="\tDFT2(&A[48],&A[50],prime);\n"
        out+="\tDFT2(&A[49],&A[51],prime);\n"
        out+="\tDFT2(&A[52],&A[54],prime);\n"
        out+="\tDFT2(&A[53],&A[55],prime);\n"
        out+="\tDFT2(&A[56],&A[58],prime);\n"
        out+="\tDFT2(&A[57],&A[59],prime);\n"
        out+="\tDFT2(&A[60],&A[62],prime);\n"
        out+="\tDFT2(&A[61],&A[63],prime);\n"
    #TWIDDLE
    out+="\n"
    if gfpf_bool:
        out+="\tA[3] = A[3].MulPowR(16);\n"
        out+="\tA[5] = A[5].MulPowR(8);\n"
        out+="\tA[7] = A[7].MulPowR(24);\n"
        out+="\tA[9] = A[9].MulPowR(4);\n"
        out+="\tA[11] = A[11].MulPowR(20);\n"
        out+="\tA[13] = A[13].MulPowR(12);\n"
        out+="\tA[15] = A[15].MulPowR(28);\n"
        out+="\tA[17] = A[17].MulPowR(2);\n"
        out+="\tA[19] = A[19].MulPowR(18);\n"
        out+="\tA[21] = A[21].MulPowR(10);\n"
        out+="\tA[23] = A[23].MulPowR(26);\n"
        out+="\tA[25] = A[25].MulPowR(6);\n"
        out+="\tA[27] = A[27].MulPowR(22);\n"
        out+="\tA[29] = A[29].MulPowR(14);\n"
        out+="\tA[31] = A[31].MulPowR(30);\n"
        out+="\tA[33] = A[33].MulPowR(1);\n"
        out+="\tA[35] = A[35].MulPowR(17);\n"
        out+="\tA[37] = A[37].MulPowR(9);\n"
        out+="\tA[39] = A[39].MulPowR(25);\n"
        out+="\tA[41] = A[41].MulPowR(5);\n"
        out+="\tA[43] = A[43].MulPowR(21);\n"
        out+="\tA[45] = A[45].MulPowR(13);\n"
        out+="\tA[47] = A[47].MulPowR(29);\n"
        out+="\tA[49] = A[49].MulPowR(3);\n"
        out+="\tA[51] = A[51].MulPowR(19);\n"
        out+="\tA[53] = A[53].MulPowR(11);\n"
        out+="\tA[55] = A[55].MulPowR(27);\n"
        out+="\tA[57] = A[57].MulPowR(7);\n"
        out+="\tA[59] = A[59].MulPowR(23);\n"
        out+="\tA[61] = A[61].MulPowR(15);\n"
        out+="\tA[63] = A[63].MulPowR(31);\n"
    elif spf_bool or bpf_bool:
        out+="\tA[3] = A[3]*omegas[16];\n"
        out+="\tA[5] = A[5]*omegas[8];\n"
        out+="\tA[7] = A[7]*omegas[24];\n"
        out+="\tA[9] = A[9]*omegas[4];\n"
        out+="\tA[11] = A[11]*omegas[20];\n"
        out+="\tA[13] = A[13]*omegas[12];\n"
        out+="\tA[15] = A[15]*omegas[28];\n"
        out+="\tA[17] = A[17]*omegas[2];\n"
        out+="\tA[19] = A[19]*omegas[18];\n"
        out+="\tA[21] = A[21]*omegas[10];\n"
        out+="\tA[23] = A[23]*omegas[26];\n"
        out+="\tA[25] = A[25]*omegas[6];\n"
        out+="\tA[27] = A[27]*omegas[22];\n"
        out+="\tA[29] = A[29]*omegas[14];\n"
        out+="\tA[31] = A[31]*omegas[30];\n"
        out+="\tA[33] = A[33]*omegas[1];\n"
        out+="\tA[35] = A[35]*omegas[17];\n"
        out+="\tA[37] = A[37]*omegas[9];\n"
        out+="\tA[39] = A[39]*omegas[25];\n"
        out+="\tA[41] = A[41]*omegas[5];\n"
        out+="\tA[43] = A[43]*omegas[21];\n"
        out+="\tA[45] = A[45]*omegas[13];\n"
        out+="\tA[47] = A[47]*omegas[29];\n"
        out+="\tA[49] = A[49]*omegas[3];\n"
        out+="\tA[51] = A[51]*omegas[19];\n"
        out+="\tA[53] = A[53]*omegas[11];\n"
        out+="\tA[55] = A[55]*omegas[27];\n"
        out+="\tA[57] = A[57]*omegas[7];\n"
        out+="\tA[59] = A[59]*omegas[23];\n"
        out+="\tA[61] = A[61]*omegas[15];\n"
        out+="\tA[63] = A[63]*omegas[31];\n"
    else:
        out+="\tA[3] = smallprimefield_multi(A[3],omegas[16],prime,R,pP);\n"
        out+="\tA[5] = smallprimefield_multi(A[5],omegas[8],prime,R,pP);\n"
        out+="\tA[7] = smallprimefield_multi(A[7],omegas[24],prime,R,pP);\n"
        out+="\tA[9] = smallprimefield_multi(A[9],omegas[4],prime,R,pP);\n"
        out+="\tA[11] = smallprimefield_multi(A[11],omegas[20],prime,R,pP);\n"
        out+="\tA[13] = smallprimefield_multi(A[13],omegas[12],prime,R,pP);\n"
        out+="\tA[15] = smallprimefield_multi(A[15],omegas[28],prime,R,pP);\n"
        out+="\tA[17] = smallprimefield_multi(A[17],omegas[2],prime,R,pP);\n"
        out+="\tA[19] = smallprimefield_multi(A[19],omegas[18],prime,R,pP);\n"
        out+="\tA[21] = smallprimefield_multi(A[21],omegas[10],prime,R,pP);\n"
        out+="\tA[23] = smallprimefield_multi(A[23],omegas[26],prime,R,pP);\n"
        out+="\tA[25] = smallprimefield_multi(A[25],omegas[6],prime,R,pP);\n"
        out+="\tA[27] = smallprimefield_multi(A[27],omegas[22],prime,R,pP);\n"
        out+="\tA[29] = smallprimefield_multi(A[29],omegas[14],prime,R,pP);\n"
        out+="\tA[31] = smallprimefield_multi(A[31],omegas[30],prime,R,pP);\n"
        out+="\tA[33] = smallprimefield_multi(A[33],omegas[1],prime,R,pP);\n"
        out+="\tA[35] = smallprimefield_multi(A[35],omegas[17],prime,R,pP);\n"
        out+="\tA[37] = smallprimefield_multi(A[37],omegas[9],prime,R,pP);\n"
        out+="\tA[39] = smallprimefield_multi(A[39],omegas[25],prime,R,pP);\n"
        out+="\tA[41] = smallprimefield_multi(A[41],omegas[5],prime,R,pP);\n"
        out+="\tA[43] = smallprimefield_multi(A[43],omegas[21],prime,R,pP);\n"
        out+="\tA[45] = smallprimefield_multi(A[45],omegas[13],prime,R,pP);\n"
        out+="\tA[47] = smallprimefield_multi(A[47],omegas[29],prime,R,pP);\n"
        out+="\tA[49] = smallprimefield_multi(A[49],omegas[3],prime,R,pP);\n"
        out+="\tA[51] = smallprimefield_multi(A[51],omegas[19],prime,R,pP);\n"
        out+="\tA[53] = smallprimefield_multi(A[53],omegas[11],prime,R,pP);\n"
        out+="\tA[55] = smallprimefield_multi(A[55],omegas[27],prime,R,pP);\n"
        out+="\tA[57] = smallprimefield_multi(A[57],omegas[7],prime,R,pP);\n"
        out+="\tA[59] = smallprimefield_multi(A[59],omegas[23],prime,R,pP);\n"
        out+="\tA[61] = smallprimefield_multi(A[61],omegas[15],prime,R,pP);\n"
        out+="\tA[63] = smallprimefield_multi(A[63],omegas[31],prime,R,pP);\n"
    #DFT 6
    out+="\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tDFT2(&A[0],&A[1]);\n"
        out+="\tDFT2(&A[2],&A[3]);\n"
        out+="\tDFT2(&A[4],&A[5]);\n"
        out+="\tDFT2(&A[6],&A[7]);\n"
        out+="\tDFT2(&A[8],&A[9]);\n"
        out+="\tDFT2(&A[10],&A[11]);\n"
        out+="\tDFT2(&A[12],&A[13]);\n"
        out+="\tDFT2(&A[14],&A[15]);\n"
        out+="\tDFT2(&A[16],&A[17]);\n"
        out+="\tDFT2(&A[18],&A[19]);\n"
        out+="\tDFT2(&A[20],&A[21]);\n"
        out+="\tDFT2(&A[22],&A[23]);\n"
        out+="\tDFT2(&A[24],&A[25]);\n"
        out+="\tDFT2(&A[26],&A[27]);\n"
        out+="\tDFT2(&A[28],&A[29]);\n"
        out+="\tDFT2(&A[30],&A[31]);\n"
        out+="\tDFT2(&A[32],&A[33]);\n"
        out+="\tDFT2(&A[34],&A[35]);\n"
        out+="\tDFT2(&A[36],&A[37]);\n"
        out+="\tDFT2(&A[38],&A[39]);\n"
        out+="\tDFT2(&A[40],&A[41]);\n"
        out+="\tDFT2(&A[42],&A[43]);\n"
        out+="\tDFT2(&A[44],&A[45]);\n"
        out+="\tDFT2(&A[46],&A[47]);\n"
        out+="\tDFT2(&A[48],&A[49]);\n"
        out+="\tDFT2(&A[50],&A[51]);\n"
        out+="\tDFT2(&A[52],&A[53]);\n"
        out+="\tDFT2(&A[54],&A[55]);\n"
        out+="\tDFT2(&A[56],&A[57]);\n"
        out+="\tDFT2(&A[58],&A[59]);\n"
        out+="\tDFT2(&A[60],&A[61]);\n"
        out+="\tDFT2(&A[62],&A[63]);\n"
    else:
        out+="\tDFT2(&A[0],&A[1],prime);\n"
        out+="\tDFT2(&A[2],&A[3],prime);\n"
        out+="\tDFT2(&A[4],&A[5],prime);\n"
        out+="\tDFT2(&A[6],&A[7],prime);\n"
        out+="\tDFT2(&A[8],&A[9],prime);\n"
        out+="\tDFT2(&A[10],&A[11],prime);\n"
        out+="\tDFT2(&A[12],&A[13],prime);\n"
        out+="\tDFT2(&A[14],&A[15],prime);\n"
        out+="\tDFT2(&A[16],&A[17],prime);\n"
        out+="\tDFT2(&A[18],&A[19],prime);\n"
        out+="\tDFT2(&A[20],&A[21],prime);\n"
        out+="\tDFT2(&A[22],&A[23],prime);\n"
        out+="\tDFT2(&A[24],&A[25],prime);\n"
        out+="\tDFT2(&A[26],&A[27],prime);\n"
        out+="\tDFT2(&A[28],&A[29],prime);\n"
        out+="\tDFT2(&A[30],&A[31],prime);\n"
        out+="\tDFT2(&A[32],&A[33],prime);\n"
        out+="\tDFT2(&A[34],&A[35],prime);\n"
        out+="\tDFT2(&A[36],&A[37],prime);\n"
        out+="\tDFT2(&A[38],&A[39],prime);\n"
        out+="\tDFT2(&A[40],&A[41],prime);\n"
        out+="\tDFT2(&A[42],&A[43],prime);\n"
        out+="\tDFT2(&A[44],&A[45],prime);\n"
        out+="\tDFT2(&A[46],&A[47],prime);\n"
        out+="\tDFT2(&A[48],&A[49],prime);\n"
        out+="\tDFT2(&A[50],&A[51],prime);\n"
        out+="\tDFT2(&A[52],&A[53],prime);\n"
        out+="\tDFT2(&A[54],&A[55],prime);\n"
        out+="\tDFT2(&A[56],&A[57],prime);\n"
        out+="\tDFT2(&A[58],&A[59],prime);\n"
        out+="\tDFT2(&A[60],&A[61],prime);\n"
        out+="\tDFT2(&A[62],&A[63],prime);\n"
    #swap
    out+="\n"
    out+="\tswap(&A[1],&A[32]);\n"
    out+="\tswap(&A[2],&A[16]);\n"
    out+="\tswap(&A[3],&A[48]);\n"
    out+="\tswap(&A[4],&A[8]);\n"
    out+="\tswap(&A[5],&A[40]);\n"
    out+="\tswap(&A[6],&A[24]);\n"
    out+="\tswap(&A[7],&A[56]);\n"
    out+="\tswap(&A[9],&A[36]);\n"
    out+="\tswap(&A[10],&A[20]);\n"
    out+="\tswap(&A[11],&A[52]);\n"
    out+="\tswap(&A[13],&A[44]);\n"
    out+="\tswap(&A[14],&A[28]);\n"
    out+="\tswap(&A[15],&A[60]);\n"
    out+="\tswap(&A[17],&A[34]);\n"
    out+="\tswap(&A[19],&A[50]);\n"
    out+="\tswap(&A[21],&A[42]);\n"
    out+="\tswap(&A[22],&A[26]);\n"
    out+="\tswap(&A[23],&A[58]);\n"
    out+="\tswap(&A[25],&A[38]);\n"
    out+="\tswap(&A[27],&A[54]);\n"
    out+="\tswap(&A[29],&A[46]);\n"
    out+="\tswap(&A[31],&A[62]);\n"
    out+="\tswap(&A[35],&A[49]);\n"
    out+="\tswap(&A[37],&A[41]);\n"
    out+="\tswap(&A[39],&A[57]);\n"
    out+="\tswap(&A[43],&A[53]);\n"
    out+="\tswap(&A[47],&A[61]);\n"
    out+="\tswap(&A[55],&A[59]);\n"
    out+="\n"
    out+="\treturn A;\n"
    out+="}\n"

    if spf_bool or bpf_bool:
        out+="template SmallPrimeField* DFT_64<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* omegas);\n"
        out+="template BigPrimeField* DFT_64<BigPrimeField>(BigPrimeField* a,BigPrimeField* omegas);\n"

    return out

'''write code for the stride permutation function
    fieldtype options are spf, bpf, gfpf
    default builds c version'''
def write_stride_permutation(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    out+="#ifndef BLOCKSIZE\n"
    out+="#define BLOCKSIZE 16\n"
    out+="#endif\n"

    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void stride_permutation(FiniteField* A,int m, int n){\n"
        out+="\tint blocksize=m^((m^n)&(-(m>n)));\n"
        out+="\tblocksize=BLOCKSIZE^((BLOCKSIZE^blocksize)&(-(BLOCKSIZE>blocksize)));\n"
        out+="\tFiniteField* B = new FiniteField[m*n];\n"
        out+="\tfor (int i = 0; i < n; i += blocksize) {\n"
        out+="\t\tfor (int j = 0; j < m; j += blocksize) {\n"
        out+="\t\t\t// transpose the block beginning at [i,j]\n"
        out+="\t\t\tfor (int k = i; k < i + blocksize; ++k) {\n"
        out+="\t\t\t\tfor (int l = j; l < j + blocksize; ++l) {\n"
        out+="\t\t\t\t\tB[k+l*n] = A[l+k*m];\n"
        out+="\t\t\t\t}\n"
        out+="\t\t\t}\n"
        out+="\t\t}\n"
        out+="\t}\n"
        out+="\tfor (long int i=0;i<m*n;i++)\n"
        out+="\t\tA[i]=B[i];\n"
        out+="\tdelete[] B;\n"
        out+="}\n"
        out+="template void stride_permutation<SmallPrimeField>(SmallPrimeField* A,int m,int n);\n"
        out+="template void stride_permutation<BigPrimeField>(BigPrimeField* A,int m,int n);\n"
        out+="template void stride_permutation<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* A,int m,int n);\n"
    else:
        out+="void stride_permutation(long int* A,int m, int n){\n"
        out+="\tint blocksize=m^((m^n)&(-(m>n)));\n"
        out+="\tblocksize=BLOCKSIZE^((BLOCKSIZE^blocksize)&(-(BLOCKSIZE>blocksize)));\n"
        out+="\tlong int* B = new long int[m*n];\n"
        out+="\tfor (int i = 0; i < n; i += blocksize) {\n"
        out+="\t\tfor (int j = 0; j < m; j += blocksize) {\n"
        out+="\t\t\t// transpose the block beginning at [i,j]\n"
        out+="\t\t\tfor (int k = i; k < i + blocksize; ++k) {\n"
        out+="\t\t\t\tfor (int l = j; l < j + blocksize; ++l) {\n"
        out+="\t\t\t\t\tB[k+l*n] = A[l+k*m];\n"
        out+="\t\t\t\t}\n"
        out+="\t\t\t}\n"
        out+="\t\t}\n"
        out+="\t}\n"
        out+="\tfor (long int i=0;i<m*n;i++)\n"
        out+="\t\tA[i]=B[i];\n"
        out+="\n"
        out+="\tdelete[] B;\n"
        out+="}\n"
    return out

'''write code for twiddle function
    fieldtype options are spf, bpf, gfpf
    default builds c version'''
def write_twiddle(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void precomputed_twiddle(FiniteField* vector,int m, int n, const FiniteField* omegas){\n"
        out+="\tFiniteField t;\n"
        out+="\tfor (int j=0;j<n;j++)\n"
        out+="\t\tfor(int i=0;i<m;i++){\n"
        out+="\t\t\tt = omegas[j*m+i];\n"
        out+="\t\t\tvector[j*m+i]=vector[j*m+i]*t;\n"
        out+="\t\t}\n"
        out+="\t}\n"
        out+="template void precomputed_twiddle<SmallPrimeField>(SmallPrimeField* vector,int m,int n,const SmallPrimeField* omegas);\n"
        out+="template void precomputed_twiddle<BigPrimeField>(BigPrimeField* vector,int m,int n,const BigPrimeField* omegas);\n"
        out+="template void precomputed_twiddle<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* vector,int m,int n,const GeneralizedFermatPrimeField* omegas);\n"
    else:
        out+="void precomputed_twiddle(long int* vector,int m, int n, const long int* omegas,long int prime, long long int R, long long int pP){\n"
        out+="\tlong int t;\n"
        out+="\tfor (int j=0;j<n;j++)\n"
        out+="\t\tfor(int i=0;i<m;i++){\n"
        out+="\t\t\tt = omegas[j*m+i];\n"
        out+="\t\t\tvector[j*m+i]=smallprimefield_multi(vector[j*m+i],t,prime,R,pP);\n"
        out+="\t\t}\n"
        out+="}\n"
    return out

def write_general_hdr(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void precomputed_DFT_general(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner);\n"
    else:
        out+="void precomputed_DFT_general(long int* vector, int K, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner);\n"
    out+="\n"
    return out

'''write code for a general dft function
    fieldtype options are spf, bpf, gfpf
    default builds c version'''
def write_general(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void precomputed_DFT_general(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner){\n"
    else:
        out+="void precomputed_DFT_general(long int* vector, int K, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner){\n"
        out+="\tlong long int pP = smallprimefield_getPp(prime,R);\n"
    #  /*       Step I       */
    out+="\tfor (long int i=0;i<=e-2;i++){\n"
    out+="\t\tfor (long int j=0;j<pow(K,i);j++){\n"
    out+="\t\t\tstride_permutation(&vector[j*(long int)pow(K,e-i)],K,(long int)pow(K,e-i-1));\n"
    out+="\t\t}\n"
    out+="\t}\n"
    #  /*       Step II       */
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tFiniteField omega_at_that_level_b;\n"
        out+="\tFiniteField* omegas;\n"
        out+="\tomega_at_that_level_b = powersofomega_inner[e-1];//POW(omega_w,(long int)(pow(K,e-1)),prime,R,pP);\n"
        out+="\tomegas = precomputeInnerPowersOfOmega(omega_at_that_level_b,K);\n"
    else:
        out+="\tlong int omega_at_that_level_b;\n"
        out+="\tlong int* omegas;\n"
        out+="\tomega_at_that_level_b = powersofomega_inner[e-1];//POW(omega_w,(long int)(pow(K,e-1)),prime,R,pP);\n"
        out+="\tomegas = precomputeInnerPowersOfOmega(omega_at_that_level_b,K,prime,R,pP);\n"

    out+="\tfor (long int j=0;j<pow(K,e-1);j++){\n"
    out+="\t\tint idx = j*K;\n"
    if spf_bool or bpf_bool or gfpf_bool:
        #out+="\t\tif (K==2){\n"
        #out+="\t\t\tDFT_2(&vector[idx],omega_at_that_level_b);\n"
        out+="\t\tif (K==8){\n"
        out+="\t\t\tDFT_8(&vector[idx],omegas);\n"
        out+="\t\t}else if (K==16){\n"
        out+="\t\t\tDFT_16(&vector[idx],omegas);\n"
        out+="\t\t}else if (K==32){\n"
        out+="\t\t\tDFT_32(&vector[idx],omegas);\n"
        out+="\t\t}else if (K==64){\n"
        out+="\t\t\tDFT_64(&vector[idx],omegas);\n"
        out+="\t\t}\n"
    else:
        #out+="\t\tif (K==2){\n"
        #out+="\t\t\tDFT_2(&vector[idx],omega_at_that_level_b,prime,R,pP);\n"
        out+="\t\tif (K==8){\n"
        out+="\t\t\tDFT_8(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t}else if (K==16){\n"
        out+="\t\t\tDFT_16(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t}else if (K==32){\n"
        out+="\t\t\tDFT_32(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t}else if (K==64){\n"
        out+="\t\t\tDFT_64(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t}\n"
    out+="\t}\n"

    out+="\tfor (int i=e-2;i>=0;i--){\n"
        #    /*       Step III       */
    out+="\t\tint stride = pow(K,e-i);\n"
    out+="\t\tfor (int j=0;j<pow(K,i);j++){\n"
    out+="\t\t\tint index = j*stride;\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\t\t\tprecomputed_twiddle(&vector[index],stride/K,K,powersofomega_outer[i]);\n"
    else:
        out+="\t\t\tprecomputed_twiddle(&vector[index],stride/K,K,powersofomega_outer[i],prime,R,pP);\n"

    out+="\t\t\tstride_permutation(&vector[index],stride/K,K);\n"
    out+="\t\t}\n"
        #    /*       Step IV       */
    out+="\t\tfor (int j=0;j<pow(K,e-1);j++){\n"
    out+="\t\t\tint idx = j*K;\n"
    if spf_bool or bpf_bool or gfpf_bool:
        #out+="\t\t\tif (K==2){\n"
        #out+="\t\t\t\tDFT_2(&vector[idx],omega_at_that_level_b);\n"
        out+="\t\t\tif (K==8){\n"
        out+="\t\t\t\tDFT_8(&vector[idx],omegas);\n"
        out+="\t\t\t}else if (K==16){\n"
        out+="\t\t\t\tDFT_16(&vector[idx],omegas);\n"
        out+="\t\t\t}else if (K==32){\n"
        out+="\t\t\t\tDFT_32(&vector[idx],omegas);\n"
        out+="\t\t\t}else if (K==64){\n"
        out+="\t\t\t\tDFT_64(&vector[idx],omegas);\n"
        out+="\t\t\t}\n"
    else:
        #out+="\t\t\tif (K==2){\n"
        #out+="\t\t\t\tDFT_2(&vector[idx],omega_at_that_level_b,prime,R,pP);\n"
        out+="\t\t\tif (K==8){\n"
        out+="\t\t\t\tDFT_8(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t\t}else if (K==16){\n"
        out+="\t\t\t\tDFT_16(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t\t}else if (K==32){\n"
        out+="\t\t\t\tDFT_32(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t\t}else if (K==64){\n"
        out+="\t\t\t\tDFT_64(&vector[idx],omegas,prime,R,pP);\n"
        out+="\t\t\t}\n"
    out+="\t\t}\n"
            #/*       Step V       */
    out+="\t\tfor (int j=0;j<pow(K,i);j++){\n"
    out+="\t\t\tint index = j*stride;\n"
    out+="\t\t\tstride_permutation(&vector[index],K,stride/K);\n"
    out+="\t\t}\n"
    out+="\t}// end of step III IV & V for_loop\n"
    out+="\tdelete[] omegas;\n"
    out+="}// end DFT_general\n"
    if spf_bool or bpf_bool or gfpf_bool:
    	out+="template void precomputed_DFT_general<SmallPrimeField>(SmallPrimeField* vector, int K, int e, SmallPrimeField omega_w, SmallPrimeField** powersofomega_outer, SmallPrimeField* powersofomega_inner);\n"
    	out+="template void precomputed_DFT_general<BigPrimeField>(BigPrimeField* vector, int K, int e, BigPrimeField omega_w, BigPrimeField** powersofomega_outer, BigPrimeField* powersofomega_inner);\n"
    	out+="template void precomputed_DFT_general<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* vector, int K, int e, GeneralizedFermatPrimeField omega_w, GeneralizedFermatPrimeField** powersofomega_outer, GeneralizedFermatPrimeField* powersofomega_inner);\n"
    out+="\n"
    return out

''' writes the code for the header
of the DFT General wrapper function '''
def write_general_wrapper_hdr(fieldtype="none"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False
    if (fieldtype=="spf"):
        spf_bool = True
    if (fieldtype=="bpf"):
        bpf_bool = True
    if (fieldtype=="gfpf"):
        gfpf_bool = True
    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void precomputed_DFT_general_new(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner);\n"
    else:
        out+="void precomputed_DFT_general_new(long int* vector, int K, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner);\n"
    out+="\n"
    return out

''' writes the code for the new dft function '''
def write_dft_wrapper(fieldtype="none"):

    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void precomputed_DFT_general_new(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner){\n"
        out+="\tif (K==8){\n"
        out+="\t\tprecomputed_DFT_general_8(vector,e,omega_w,powersofomega_outer,powersofomega_inner);\n"
        out+="\t}else if (K==16){\n"
        out+="\t\tprecomputed_DFT_general_16(vector,e,omega_w,powersofomega_outer,powersofomega_inner);\n"
        out+="\t}else if (K==36){\n"
        out+="\t\tprecomputed_DFT_general_32(vector,e,omega_w,powersofomega_outer,powersofomega_inner);\n"
        out+="\t}else if (K==64){\n"
        out+="\t\tprecomputed_DFT_general_64(vector,e,omega_w,powersofomega_outer,powersofomega_inner);\n"
        out+="\t}\n"
    else:
        out+="void precomputed_DFT_general_new(long int* vector, int K, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner){\n"
        out+="\tif (K==8){\n"
        out+="\t\tprecomputed_DFT_general_8(vector,e,omega_w, prime, R, powersofomega_outer,powersofomega_inner);\n"
        out+="\t}else if (K==16){\n"
        out+="\t\tprecomputed_DFT_general_16(vector,e,omega_w, prime, R, powersofomega_outer,powersofomega_inner);\n"
        out+="\t}else if (K==36){\n"
        out+="\t\tprecomputed_DFT_general_32(vector,e,omega_w, prime, R, powersofomega_outer,powersofomega_inner);\n"
        out+="\t}else if (K==64){\n"
        out+="\t\tprecomputed_DFT_general_64(vector,e,omega_w, prime, R, powersofomega_outer,powersofomega_inner);\n"
        out+="\t}\n"
    out+="}\n"
    return out

''' writes code for the
'''
def write_general_K_hdr(fieldtype="none", K="8"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool = False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void precomputed_DFT_general_"
        out+=str(K)
        out+="(FiniteField* vector, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner);\n"
    else:
        out+="void precomputed_DFT_general_"
        out+=str(K)
        out+="(long int* vector, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner);\n"

    out+="\n"
    return out

'''write code for a dft
fieldtype options are spf, bpf, gfpf
default builds c version'''
def write_general_K(fieldtype="none",K="8"):
    out="\n"
    spf_bool = False
    bpf_bool = False
    gfpf_bool= False

    if (fieldtype=="spf"):
        spf_bool = True
    if(fieldtype=="bpf"):
        bpf_bool = True
    if(fieldtype=="gfpf"):
        gfpf_bool = True

    if spf_bool or bpf_bool or gfpf_bool:
        out+="template<class FiniteField>\n"
        out+="void precomputed_DFT_general_"
        out+=str(K)
        out+="(FiniteField* vector, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner){\n"
    else:
        out+="void precomputed_DFT_general_"
        out+=str(K)
        out+="(long int* vector, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner){\n"
        out+="\tlong long int pP = smallprimefield_getPp(prime,R);\n"

    #  /*       Step I       */
    out+="\tfor (long int i=0;i<=e-2;i++){\n"
    out+="\t\tfor (long int j=0;j<pow("
    out+=str(K)
    out+=",i);j++){\n"
    out+="\t\t\tstride_permutation(&vector[j*(long int)pow("
    out+=str(K)
    out+=",e-i)],"
    out+=str(K)
    out+=",(long int)pow("
    out+=str(K)
    out+=",e-i-1));\n"
    out+="\t\t}\n"
    out+="\t}\n"

    #  /*       Step II       */
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\tFiniteField omega_at_that_level_b;\n"
        out+="\tFiniteField* omegas;\n"
    else:
        out+="\tlong int omega_at_that_level_b;\n"
        out+="\tlong int* omegas;\n"


    if spf_bool or bpf_bool or gfpf_bool:
	out+="\tomegas = precomputeInnerPowersOfOmega(omega_at_that_level_b,"
        out+=str(K)
        out+=");\n"
    else:
	out+="\tomegas = precomputeInnerPowersOfOmega(omega_at_that_level_b,"
        out+=str(K)
        out+=",prime,R,pP);\n"

    out+="\tfor (long int j=0;j<pow("
    out+=str(K)
    out+=",e-1);j++){\n"
    out+="\t\tint idx = j*"
    out+=str(K)
    out+=";\n"
    if spf_bool or bpf_bool or gfpf_bool:
        if (K=="2"):
            out+="\t\tDFT_2(&vector[idx],omega_at_that_level_b);\n"
        if (K=="8"):
            out+="\t\tDFT_8(&vector[idx],omegas);\n"
        if (K=="16"):
            out+="\t\tDFT_16(&vector[idx],omegas);\n"
        if (K=="32"):
            out+="\t\tDFT_32(&vector[idx],omegas);\n"
        if (K=="64"):
            out+="\t\tDFT_64(&vector[idx],omegas);\n"
    else:
        if (K=="2"):
            out+="\t\tDFT_2(&vector[idx],omega_at_that_level_b,prime,R,pP);\n"
        if (K=="8"):
            out+="\t\tDFT_8(&vector[idx],omegas,prime,R,pP);\n"
        if (K=="16"):
            out+="\t\tDFT_16(&vector[idx],omegas,prime,R,pP);\n"
        if (K=="32"):
            out+="\t\tDFT_32(&vector[idx],omegas,prime,R,pP);\n"
        if (K=="64"):
            out+="\t\tDFT_64(&vector[idx],omegas,prime,R,pP);\n"
    out+="\t}\n"
    out+="\tfor (int i=e-2;i>=0;i--){\n"
    #    /*       Step III       */
    out+="\t\tint stride = pow("
    out+=str(K)
    out+=",e-i);\n"
    out+="\t\tfor (int j=0;j<pow("
    out+=str(K)
    out+=",i);j++){\n"
    out+="\t\t\tint index = j*stride;\n"
    if spf_bool or bpf_bool or gfpf_bool:
        out+="\t\t\tprecomputed_twiddle(&vector[index],stride/"
        out+=str(K)
        out+=","
        out+=str(K)
        out+=",powersofomega_outer[i]);\n"
    else:
        out+="\t\t\tprecomputed_twiddle(&vector[index],stride/"
        out+=str(K)
        out+=","
        out+=str(K)
        out+=",powersofomega_outer[i],prime,R,pP);\n"
    out+="\t\t\tstride_permutation(&vector[index],stride/"
    out+=str(K)
    out+=","
    out+=str(K)
    out+=");\n"
    out+="\t\t}\n"
    #    /*       Step IV       */
    out+="\t\tfor (int j=0;j<pow("
    out+=str(K)
    out+=",e-1);j++){\n"
    out+="\t\t\tint idx = j*"
    out+=str(K)
    out+=";\n"
    if spf_bool or bpf_bool or gfpf_bool:
        if (K=="2"):
            out+="\t\t\t\tDFT_2(&vector[idx],omega_at_that_level_b);\n"
        if (K=="8"):
            out+="\t\t\t\tDFT_8(&vector[idx],omegas);\n"
        if (K=="16"):
            out+="\t\t\t\tDFT_16(&vector[idx],omegas);\n"
        if (K=="32"):
            out+="\t\t\t\tDFT_32(&vector[idx],omegas);\n"
        if (K=="64"):
            out+="\t\t\t\tDFT_64(&vector[idx],omegas);\n"
    else:
        if (K==2):
            out+="\t\t\t\tDFT_2(&vector[idx],omega_at_that_level_b,prime,R,pP);\n"
        if (K==8):
            out+="\t\t\t\tDFT_8(&vector[idx],omegas,prime,R,pP);\n"
        if (K==16):
            out+="\t\t\t\tDFT_16(&vector[idx],omegas,prime,R,pP);\n"
        if (K==32):
            out+="\t\t\t\tDFT_32(&vector[idx],omegas,prime,R,pP);\n"
        if (K==64):
            out+="\t\t\t\tDFT_64(&vector[idx],omegas,prime,R,pP);\n"
    out+="\t\t}\n"
    #/*       Step V       */
    out+="\t\tfor (int j=0;j<pow("
    out+=str(K)
    out+=",i);j++){\n"
    out+="\t\t\tint index = j*stride;\n"
    out+="\t\t\tstride_permutation(&vector[index],"
    out+=str(K)
    out+=",stride/"
    out+=str(K)
    out+=");\n"
    out+="\t\t}\n"
    out+="\t}// end of step III IV & V for_loop\n"
    out+="\tdelete[] omegas;\n"
    out+="}// end DFT_general\n"
    out+="\n"

    return out

''' Appends proc info and cache size to end of file named filename '''
def write_proc_info_to_file(proc_name,proc_cachesize,filename="proc_db.txt"):
    lines = open(filename, 'r').readlines()
    text = "" + proc_name + " " + proc_cachesize + "\n"
    lines.append(text)
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
    return


#### TESTS

#how to use num_occurences
#str = "this is a string that has two occurences of the word this."
#sub = "this"
#print(num_occurences(str,sub))
#sub = "that"
#print(num_occurences(str,sub))

#another example of num_occurences
#sub = "dft("
#str = file_to_string("file_1.txt")
#print(num_occurences(str,sub))

#how to use substring_indices
#inds = substring_indices(str,sub)
#print("indices:")
#print(inds)

#how to use function get_line_numbers
#regex = r"dft\((\w+),(\w+),(\w+),(\w+)\)"
#lns = get_line_numbers("file_1.txt",regex)
#print("line numbers:")
#print(lns)


#regex = r"dft\(([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+)\);"
#NOT USED ANYMORE

# regex = r"dft\(([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+),([0-9a-zA-Z]+)\);"
# print(regex)
# lns = get_line_numbers(file_to_scan,regex)
#
# print("dft line numbers:")
# print(lns)
# if lns!=-1:
#     print(len(lns))
#     dft_bool = True
#
# dft_tokens = get_token_list(file_to_scan,regex,1)
# print("dft tokens:")
# print(dft_tokens)
#
# regex = r"SmallPrimeField[ *][ ]?(\w*)"
# spf_tokens = get_token_list(file_to_scan,regex,1)
# print("spf tokens:")
# print(spf_tokens)
#
# regex = r"BigPrimeField[ *][ ]?(\w*)"
# bpf_tokens = get_token_list(file_to_scan,regex,1)
# print("bpf tokens:")
# print(bpf_tokens)
#
# regex = r"GeneralizedFermatPrimeField[ *][ ]?(\w*)"
# gfpf_tokens = get_token_list(file_to_scan,regex,1)
# print("gfpf tokens:")
# print(gfpf_tokens)
#
#
# spf_bool = False
# bpf_bool = False
# gfpf_bool = False
#
# for tok in dft_tokens:
#     if any(x for x in spf_tokens if x == tok):
#         spf_bool = True
#     if any(x for x in bpf_tokens if x==tok):
#         bpf_bool = True
#     if any(x for x in gfpf_tokens if x==tok):
#         gfpf_bool = True


###CREATE HEADERS
headers = []

#FOR THE SWAP AND DFT2/BUTTERFLY FUNCTIONS
swp=write_swap("spf")
butterfly=write_DFT2("spf")
headers.append(swp)
headers.append(butterfly)
###CREATE KERNALS

kernals = []

###FOR EACH DFT BASECASE FIELDTYPE
for basecase in basecases:

    if (basecase==8):
        hdr = write_dft_header("spf",8)
        headers.append(hdr)
        krnl = write_dft_8("spf")
        kernals.append(krnl)
        #     spf_kernal_8 = get_lines_from_file(DFT_8_TEMPLATE_FILE,DFT_8_SPF_START_LINE,DFT_8_SPF_TOTAL_LINES)
        #     kernals.append(gfpf_kernal_8)
        hdr = write_dft_header("gfpf",8)
        headers.append(hdr)
        krnl = write_dft_8("gfpf")
        kernals.append(krnl)
        #     gfpf_kernal_8 = get_lines_from_file(DFT_8_TEMPLATE_FILE,DFT_8_GFPF_START_LINE,DFT_8_GFPF_TOTAL_LINES)
        #     kernals.append(gfpf_kernal_8)

    if (basecase==16):
        hdr = write_dft_header("spf",16)
        headers.append(hdr)
        krnl = write_dft_16("spf")
        kernals.append(krnl)
        hdr = write_dft_header("gfpf",16)
        headers.append(hdr)
        krnl = write_dft_16("gfpf")
        kernals.append(krnl)

    if (basecase==32):
        hdr = write_dft_header("spf",32)
        headers.append(hdr)
        krnl=write_dft_32("spf")
        kernals.append(krnl)
        krnl=write_dft_32("gfpf")
        hdr = write_dft_header("gfpf",32)
        kernals.append(krnl)
        headers.append(hdr)

    if (basecase==64):
        hdr = write_dft_header("spf",64)
        headers.append(hdr)
        krnl = write_dft_64("spf")
        kernals.append(krnl)
        hdr = write_dft_header("gfpf",64)
        headers.append(hdr)
        krnl = write_dft_64("gfpf")
        kernals.append(krnl)

###CREATE STRIDE, TWIDDLE & GENERAL DFT
stride=write_stride_permutation("spf")
kernals.append(stride)
twiddle=write_twiddle("spf")
kernals.append(twiddle)
pcomi = write_precompute_powers_of_omega("spf","inner")
kernals.append(pcomi)
pcomo = write_precompute_powers_of_omega("spf","outer")
kernals.append(pcomo)
hdr=write_general_hdr("spf")
headers.append(hdr)
gen=write_general("spf")
kernals.append(gen)

### updated general dft
genk=write_general_K("spf","8")
kernals.append(genk)
genk=write_general_K("spf","16")
kernals.append(genk)
genk=write_general_K("spf","32")
kernals.append(genk)
genk=write_general_K("spf","64")
kernals.append(genk)

genkhdr=write_general_wrapper_hdr("spf")
headers.append(genkhdr)
genk=write_dft_wrapper("spf")
kernals.append(genk)

#### DFT C VERSION
c_headers = []
c_kernals = []
swpc=write_swap()
butterflyc=write_DFT2()
powc=write_host_pow()
c_headers.append(swpc)
c_headers.append(butterflyc)
c_headers.append(powc)
# C basecases
krnl = write_dft_8()
hdr = write_dft_header('none',8)
c_headers.append(hdr)
c_kernals.append(krnl)
krnl = write_dft_16()
hdr = write_dft_header('none',16)
c_headers.append(hdr)
c_kernals.append(krnl)
krnl = write_dft_32()
hdr = write_dft_header('none',32)
c_headers.append(hdr)
c_kernals.append(krnl)
krnl = write_dft_64()
hdr = write_dft_header('none',64)
c_headers.append(hdr)
c_kernals.append(krnl)

# C stride_permutation
stride=write_stride_permutation()
c_kernals.append(stride)
# C twiddle
twiddle=write_twiddle()
c_kernals.append(twiddle)
# C precompute_powers_of_omega inner and outer
pci=write_precompute_powers_of_omega("none","inner")
pcihdr=write_precompute_header("none","inner")
c_headers.append(pcihdr)
c_kernals.append(pci)
pco=write_precompute_powers_of_omega("none","outer")
c_kernals.append(pco)
pcohdr=write_precompute_header("none","outer")
c_headers.append(pcohdr)
# C general routine
gen=write_general()
c_kernals.append(gen)
hdr=write_general_hdr()
c_headers.append(hdr)

# general_kernal = get_lines_from_file(DFT_GENERAL_TEMPLATE_FILE,DFT_GENERAL_START_LINE,DFT_GENERAL_TOTAL_LINES)
# kernals.append(general_kernal)

hdr_output=open("../../../include/FFT/src/custom_DFT.h",'w')
#hdr_output=open("custom_DFT.h",'w')
hdr_output.write("#ifndef __CUSTOM_DFT_H \n")
hdr_output.write("#define __CUSTOM_DFT_H \n")
hdr_output.write("#include \"../../../include/FiniteFields/SmallPrimeField_Support.h\"\n")
hdr_output.write("#include \"../../../include/FiniteFields/BigPrimeField_Support.h\"\n")
hdr_output.write("#include \"../../../include/FiniteFields/SmallPrimeField.hpp\"\n")
hdr_output.write("#include \"../../../include/FiniteFields/BigPrimeField.hpp\"\n")
hdr_output.write("#include \"../../../include/FiniteFields/GeneralizedFermatPrimeField.hpp\"\n")


for hdr in headers:
    hdr_output.writelines(hdr)
# add c headers
for c_hdr in c_headers:
    hdr_output.writelines(c_hdr)

hdr_output.write("#endif\n")
hdr_output.close()

output = open("custom_DFT.cpp",'w')
output.write("#include \"../../../include/bpas.h\"\n")
output.write("#include \"../../../include/FFT/src/custom_DFT.h\"\n")

for kernal in kernals:
    output.writelines(kernal)
# add c functions
for c_kernal in c_kernals:
    output.writelines(c_kernal)

output.close()
