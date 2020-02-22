## author: Davood
#######################################
from sympy import *
from random import *
from gmpy import mpz 
import sys
import os

########################################
def getenv(env_var):
	return os.environ[env_var];

#######################################
def get_readable_value(n):
	list=["","k","m","g","t","p"];

	idx=0;
	h=n;
	while h>=1024:
		h/=1024
		idx+=1;

	return str(h)+list[idx]

#######################################
def echo_file(file_name):
	t=open(file_name,"r");
	v=t.read().split("\n")[:];
	t.close();
	return v;

#######################################
class big_prime_field:
	
	def set_prime_name(self,prime_name):

		list=[]

		libraryPath=getenv("CUMODP_HOME");
		projectPath=(libraryPath+"/new_tests/bigfft_u32/");
		prime_path=(projectPath+"/test/primes/"+prime_name);
		radix_path=prime_path+"/radix";
		k_path=prime_path+"/k";
		omega_prefix=(prime_path+"/omega/omega_K");
		
		list.append("prime_name");
		list.append("libraryPath");
		list.append("projectPath");
		list.append("prime_path");
		list.append("omega_prefix");


		k=int(echo_file(k_path)[0]);
		r_bits=echo_file(radix_path);
		r=0;	
		for b in r_bits:
			if b!="":
				r+=(1<<int(b));

		list.append("k");
		list.append("r");	

		K=2*k
		p=r**k+1;	

		list.append("K");
		list.append("p");

		#read data files for x and y
		#stored in form of vectors of k=8 coefficients
		# make sure you have executed "make generateData" in "BigPrimeFieldFFT_3"
		
		for x in list:
			setattr(self, x, eval(x))



P = big_prime_field ()#e=2

#######################################
def extract_bitmap(x):
	m=str(bin(x)).replace("0b","")[::-1]


	idx=0
#	print("m-reverse",m)
	p_bitmap=[]
	for x in m:
		if int(x)==1:
#			print(idx)
			p_bitmap.append(idx)
		idx+=1;
#	print("bitmap",p_bitmap)
	return p_bitmap

#######################################
def log2(n):
	return int(log(n,2));

#######################################
def pow_mod_bigp(x,e,p):
#	print("bitmap-phase");
	bitmap=extract_bitmap(e);
	m=max(bitmap)
	max_pow_2=m;#int(log(m,2))
#	print("max_pow_2",max_pow_2)
	power_vector_size=(1+max_pow_2);
	power_vector=power_vector_size*[mpz(0)]
	power_vector[0]=mpz(x);
#	power_vector[1]=x; # x is the number to be inversed

	
#	print("exp-phase-I");
	for i in range(1,max_pow_2+1):
		power_vector[i]=(power_vector[i-1]*power_vector[i-1])%(p);
	
#	print(power_vector)
#	stored_powers=len(bitmap)*[1]

	cnt=0;
	m=1;
#	print("bitmap-phase-II","len_bitmap="+str(len(bitmap)), (bitmap[0]));
	for b in bitmap:
#		m=(m*power_vector[b])%p
#		print("pow_vec_b",power_vector[b]);
		m=(m*mpz(power_vector[b]))%p##

	return m;
#######################################
def is_prime(p):
	return 0;


########################################
def pow_mod_p(x,n,p):
	m=mpz(x)
	n=mpz(n)
	print(n)
	if n==0:
		m=1
	
	i=0;
#	for i in range(n-1):
	while i<n:
		print(i)
		m=(m*x)%p;
		i+=1
	
	return m;
#######################################
def list_int_factors(n):
	f=factorint(n)
	bases=f.keys();
	exponents=f.values();

	exponent_list=[]
	len_list_n=[]
	for i in range(len(bases)):
		tmp_list=[]
		for j in range((exponents[i]+1)):
			tmp_list+=[mpz(bases[i])**j];
#		print(tmp_list)
		exponent_list.append(tmp_list);
		len_list_n.append(len(tmp_list));
	
#	print(exponent_list);
#	print(len_list_n)

	n_product=mpz(1);
	for x in len_list_n:
		n_product*=x;
		
	value_list=exponent_list

	factor_list=[]
	for e in range(n_product):
		v=mpz(e);
		m=mpz(1);
		for i in range(len(len_list_n)):
			t=v%len_list_n[i]
			m*=value_list[i][t]
			v/=len_list_n[i];
#		print("m",m);
		factor_list.append(m);
	
	return factor_list
#######################################
# find omega s.t. (omega)**n == 1 mod p
def generate_omega(n,p,v=1):
	
	a=mpz(0)
	omega=mpz(0)
	q=(p-1)/n;
	if (p-1)%n!=0:
		print ("n not dividing p-1");
		return -1;
	
	# removing 1 and n from the list of factors
	n_factor_list=list_int_factors(n)[1:-1] 
	if v==1:
		print("len(n_factor_list)=",len(n_factor_list));
	
	while True:	
		if v==1:
			print("...random a ");	
		a=randint(2,p-1);
	
		if v==1:
			print("... omega")
#		omega=(a**q)%p
#		omega=pow_mod_p(a,q,p)##	
		omega=pow_mod_bigp(a,q,p)####
		if omega==1:
			continue;
		omega_okay=1;
		for x in n_factor_list:
#			if v==1:
#				print("in factor_list");
			if (pow_mod_bigp(omega,x,p)==1):
#			if ((mpz(omega)**mpz(x))%p)==1:
				if v==1:
					print ("Error: for factor x|n; omega^x==1!",omega, x)
				omega_okay=0;
				break;
		if omega_okay==1:
			break;
	
	if v==1:
		print("last_check!")
#	if ((omega**n)%p)==1  and omega!=1:
	if (pow_mod_bigp(omega,n,p))==1  and omega!=1:##
		if v==1:
			print ("omega is",omega);
		return omega;


#######################################
def toBigFieldElement(n):
#	print("to bigfield",n)
	
	k=P.k;
	r=P.r;
	p=P.p;
	v=k*[0];
	n=mpz(n%p);
	for i in range(k):
		v[i]=n%r;
		n/=r;
	return v;
#######################################
def write_to_file(lines_list, file_path):
#	print("opening file", file_path);
	f=open(file_path,"w");
	for l in lines_list:
		f.write(str(l)+"\n");
	
	f.close();

#######################################

def test_generate_roots_of_unity(prime_name):
	P.set_prime_name(prime_name);
	for e in range(1,10):
		N=(P.K)**(e);
		print("N,K,e",N,P.K, e)
		omega=generate_omega(N,P.p);	
		current_omega_path=P.omega_prefix+str(e);
		if e==1:
			omega=P.r;
		print(current_omega_path);
		omega_in_bigp=toBigFieldElement(omega);
		write_to_file(omega_in_bigp, current_omega_path);

#######################################
def main():
	argc=len(sys.argv);
	if argc<4:
#		print("nothing");
		return 0;
	
	n=int(sys.argv[1]);
	r=int(sys.argv[2]);
	k=int(sys.argv[3]);	

	print("[n=%d,r=%d,k=%d]"%(n,r,k))

	output_name="current_omega"
	if (argc>4):
		output_name=str(sys.argv[4]);
	p=r**k+1;

	while True:
		w=generate_omega(n,p);
		print("w",w)
		if (pow(w,n/(2*k))%p)==r:
			break;
#		else:
#			print("continue");

#	print("omega",w);
#	print(p);
#	print("w^n=",(w**n)%p);
	w_out=k*[0];
	for i in range(k):
		w_out[i]=w%r;
		w/=r;
	write_to_file(w_out, output_name);

	return w;

#######################################

if __name__ == "__main__":
	main();

