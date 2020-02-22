/**
 * Implementation of 'MulSSA.h'
 *
 * @author Yuzhen Xie
 */

#include "../../../include/IntegerPolynomial/Multiplication/MulSSA.h"

/*----------------------------------------
Our new 2-convolution multiplication
-----------------------------------------*/
void MulSSA::mul2C(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c)
{//---------------------------------------------------------------

  //debug
  //a->writeToFile("a.input");
  //b->writeToFile("b.input");

//  unsigned long long start = __cilkview_getticks();

  cilk_spawn
    a->setCoefficientDigits(a->maxCoefficientSize()); 
    b->setCoefficientDigits(b->maxCoefficientSize());
  cilk_sync;
 
  int digitCount = MAX_COEFF_DIGITS(a, b);
  
  determineGoodN(digitCount, MAX(a->getSize(),b->getSize()));

//  unsigned long long end = __cilkview_getticks();
//  std::cout << "determineGoodN (sec): \t" << (end - start) /1000.f << std::endl;

//  std::cout<<"digitCount, N, K, M, LIMB_BITS, bound: " <<digitCount<<", "<<N<<", "<<K<<", "<<M<<", "<<LIMB_BITS<<", "<<2+log(K*a->getSize())/log(2)+2*M<<std::endl;
 
//  start = __cilkview_getticks();

  BivariatePolynomial *biA = cilk_spawn 
    ToBivarMod(a->getCoefficients(), a->getSize());
  BivariatePolynomial *biB = ToBivarMod(b->getCoefficients(), b->getSize());
  cilk_sync;

//  end = __cilkview_getticks();
//  std::cout << "ToBivarMod (sec): \t" << (end - start) /1000.f << std::endl;

  //debug
  //biA->writeToFile("biA.input");
  //biB->writeToFile("biB.input");

//  start = __cilkview_getticks();

  //cyclic and negacyclic convolution
  int d1 = a->getSize();
  int d2 = b->getSize();
  sfixn *ncc = PBPAS::TwoConvolutionMod(biA->getCoefficients(), biB->getCoefficients(), d1, d2, P, K);

  biA->freeHeap();
  biB->freeHeap();

//  end = __cilkview_getticks();
//  std::cout << "\nTwoConvolutionMod (sec): \t" << (end - start) /1000.f << std::endl;

  /*
  std::cout<<"Negacyclic output \n";
  for(int i=0; i<K*(d1+d2-1); ++i)
    std::cout<<ncc[i]<<", ";
  std::cout<<"\n\n";

  std::cout<<"Cyclic output \n";
  for(int i=K*(d1+d2-1); i<K*(d1+d2-1)*2; ++i)
    std::cout<<ncc[i]<<", ";
  std::cout<<"\n\n";
  */

//  start = __cilkview_getticks();

  RecoverProduct(ncc, d1, d2, c->getCoefficients());
  my_free(ncc);

//  end = __cilkview_getticks();
//  std::cout << "RecoverProduct (sec): \t" << (end - start) /1000.f << std::endl;
}

/*---------------------------------------------------------------------
  for both C^+ (ncc) and C^- (cc), for each coefficient of y^i,
  which is a large integer encoded as a polynomial in x, i.e.
  ncc_i(x), cc_i(x), 
  (1) convert ncc_i and cc_i to GMP, get u_i and v_i
  (2) compute c_i = (u_i+v_i)/2+(-u_i+v_i)/2*2^N

good version
--------------------------------------------------------------------------*/
void MulSSA::RecoverProduct(sfixn *ncc, int d1, int d2, mpz_class *c)
{//------------------------------------------------------------------------
  int N1 = N-1;
  int ys = d1+d2-1;
  sfixn *cc = ncc + K*ys;

#pragma cilk_grainsize = 8192; //good for d=65536,dig=65535
  cilk_for(int i=0; i<ys; ++i){//enough parallelism for 48-core
    int Ks = i*K;
    sfixn *ncc_i = ncc+Ks;
    sfixn *cc_i = cc+Ks;
    
    //evaluate p(x) at 2^M to get big int
    mpz_t u_i_mpz, v_i_mpz;
    //mpz_init(u_i_mpz);
    //mpz_init(v_i_mpz);
    
    ToMPZ(ncc_i, LIMB_BITS, u_i_mpz, K, M);
    ToMPZ(cc_i, LIMB_BITS, v_i_mpz, K, M);
    
    mpz_add(c[i].get_mpz_t(), u_i_mpz, v_i_mpz); //c[i] = (u_i + v_i)/2
    c[i]>>=1;
    
    mpz_sub(v_i_mpz, v_i_mpz, u_i_mpz); //v_i = (-u_i + v_i)/2*2^N
    mpz_mul_2exp(v_i_mpz, v_i_mpz, N1);
    
    mpz_add(c[i].get_mpz_t(), c[i].get_mpz_t(), v_i_mpz);
    
    mpz_clear(u_i_mpz);
    mpz_clear(v_i_mpz);  
  }    
}

void MulSSA::RecoverProduct1(sfixn *ncc, int d1, int d2, mpz_class *c)
{//------------------------------------------------------------------------
  int N1 = N-1;
  int ys = d1+d2-1;
  sfixn *cc = ncc + K*ys;

  for(int i=0; i<ys; ++i){
    int Ks = i*K;
    sfixn *ncc_i = ncc+Ks;
    sfixn *cc_i = cc+Ks;

    //evaluate p(x) at 2^M to get big int
    mpz_t u_i_mpz, v_i_mpz;
    //mpz_init(u_i_mpz);
    //mpz_init(v_i_mpz);

    ToMPZ(ncc_i, LIMB_BITS, u_i_mpz, K, M);
    ToMPZ(cc_i, LIMB_BITS, v_i_mpz, K, M);

    mpz_add(c[i].get_mpz_t(), u_i_mpz, v_i_mpz); //c[i] = (u_i + v_i)/2
    c[i]>>=1;
   
    mpz_sub(v_i_mpz, v_i_mpz, u_i_mpz); //v_i = (-u_i + v_i)/2*2^N
    mpz_mul_2exp(v_i_mpz, v_i_mpz, N1);
    
    mpz_add(c[i].get_mpz_t(), c[i].get_mpz_t(), v_i_mpz);
    
    mpz_clear(u_i_mpz);
    mpz_clear(v_i_mpz);   
  } 
}

void MulSSA::RecoverProduct2(sfixn *ncc, int d1, int d2, mpz_class *c)
{//------------------------------------------------------------------------
  int N1 = N-1;
  int ys = d1+d2-1;
  sfixn *cc = ncc + K*ys;

#pragma cilk_grainsize = 16384;
  cilk_for(int i=0; i<ys; ++i){//enough parallelism for 48-core? 
    int Ks = i*K;
    sfixn *ncc_i = ncc+Ks;
    sfixn *cc_i = cc+Ks;

    //evaluate p(x) at 2^M to get big int
    mpz_t u_i_mpz, v_i_mpz;
    //mpz_init(u_i_mpz);
    //mpz_init(v_i_mpz);

    ToMPZ(ncc_i, LIMB_BITS, u_i_mpz, K, M);
    ToMPZ(cc_i, LIMB_BITS, v_i_mpz, K, M);

    mpz_add(c[i].get_mpz_t(), u_i_mpz, v_i_mpz); //c[i] = (u_i + v_i)/2
    c[i]>>=1;
   
    mpz_sub(v_i_mpz, v_i_mpz, u_i_mpz); //v_i = (-u_i + v_i)/2*2^N
    mpz_mul_2exp(v_i_mpz, v_i_mpz, N1);
    
    mpz_add(c[i].get_mpz_t(), c[i].get_mpz_t(), v_i_mpz);
    
    mpz_clear(u_i_mpz);
    mpz_clear(v_i_mpz);  
  }    
}


/*---------------------------------------------------------------------
  for both C^+ (ncc) and C^- (cc), for each coefficient of y^i,
  which is a large integer encoded as a polynomial in x, i.e.
  ncc_i(x), cc_i(x), 
  (1) convert ncc_i and cc_i to GMP, get u_i and v_i
  (2) compute c_i = (u_i+v_i)/2+(-u_i+v_i)/2*2^N
--------------------------------------------------------------------------*/
void MulSSA::ToGMP_RecoverProduct(sfixn *ncc, int d1, int d2, mpz_class *c)
{//------------------------------------------------------------------------
  int N1 = N-1;
  int Ke = K-1;
  int K2 = K-2;
  int ys = d1+d2-1;

  sfixn *cc = ncc + K*ys;
  cilk_for(int i=0; i<ys; ++i){//enough parallelism for 48-core?
 
    int Ks = i*K;
    sfixn *ncc_i = ncc+Ks;
    sfixn *cc_i = cc+Ks;
    //evaluate p(x) at 2^M by Horner's rule

    mpz_class ncc_i_gmp = ncc_i[Ke] ;
    if (ncc_i[Ke] > HALF_P) //bring back to center
      ncc_i_gmp = ncc_i[Ke] - P;
    //if (ncc_i[Ke] < N_HALF_P) //safe check
    //  ncc_i_gmp = ncc_i[Ke] + P;
  
    mpz_class cc_i_gmp = cc_i[Ke];
    if (cc_i[Ke] > HALF_P)
      cc_i_gmp = cc_i[Ke] - P;
    //if (cc_i[Ke] < N_HALF_P)
    //  cc_i_gmp = cc_i[Ke] + P;

    for(int j=K2; j>=0; --j){
      if(ncc_i_gmp!=0)
	ncc_i_gmp <<= M; 
      if (ncc_i[j] > HALF_P)
	ncc_i_gmp = ncc_i_gmp + ncc_i[j] - P;
      //else if (ncc_i[j] < N_HALF_P)
      //  ncc_i_gmp = ncc_i_gmp + ncc_i[j] + P;
      else
	ncc_i_gmp += ncc_i[j];

      if(cc_i_gmp!=0)
	cc_i_gmp <<= M; 
      if (cc_i[j] > HALF_P)
	cc_i_gmp = cc_i_gmp + cc_i[j] - P;
      //else if (cc_i[j] < N_HALF_P)
      //  cc_i_gmp = cc_i_gmp + cc_i[j] + P;
      else
	cc_i_gmp += cc_i[j];
    }
    //c=a*b
    mpz_class c1 = ncc_i_gmp+cc_i_gmp;   
    c1 >>= 1;   
    mpz_class c2 = cc_i_gmp-ncc_i_gmp;    
    c2 <<=N1;    
    c[i] = c1 + c2; 
  }
}


/**
 * Convert a univariate polynomial with large coefficients to 
 * a bivariate polynomial representation with relatively 
 * small coefficients by a base of 2^M and mod a prime number p
 * 
 * univariate --> a1    * y ^ 0 + a2    * y ^ 1 + ... + ad    * y ^ (d-1)
 * bivariate  --> A1(x) * y ^ 0 + A2(x) * y ^ 1 + ... + Ad(x) * y ^ (d-1)
 *
 * Ai(x) = b1 * x ^ 0 + b2 * x ^ 1 + ... + bK * x ^ (K-1)
 *
 * @param coeff the coefficients of large integers
 * @param d partial degree of the univariate polynomial plus one
 *
 * good version
 */
BivariatePolynomial * MulSSA::ToBivarMod(mpz_class *coeff, int d)
{//------------------------------------------------------------------
  BivariatePolynomial *biA = new BivariatePolynomial(d, K, M);
  
  //convert a large int coeff[i] to a univar poly of degree K-1
  //by a base of 2^M
  //store the coefficients of the univar poly to the coeff of y^i
  //setup grain size, need to tune
#pragma cilk_grainsize = 8192; //good for d=65536,dig=65535
  cilk_for(int i=0; i<d; ++i){
    if (coeff[i]!=0) {
      sfixn *Bi = biA->getCoefficients() + i*K;
      mpzToPolyMod((coeff[i]).get_mpz_t(), M, P, Bi);
    }
  }
  
  return biA;
}

//---------------------------------------------------
//just for binary representation mpz
//precondition: z!=0, number of bits of z / M <= K
//K = size of X 
void MulSSA::mpzToPolyMod(mpz_t z, int M, sfixn prime, sfixn *X) 
{//----------------------------------------------------------  
  int s = z->_mp_size; //abs(s)>=1
  int isNeg = 0;
  if (s<0) {
    isNeg = 1;
    s = -s;
  }

  //std::cout<<"s: "<<s<<std::endl;

  mp_ptr limbs = z->_mp_d;
  mp_limb_t carry = 0;
  mp_size_t carry_size = 0;

  int e = s-1; 
  int index = 0;

  //all the limbs except for the last one are full-size (64 bits)
  int b = GMP_LIMB_BITS;
  //size_t b = 8*sizeof(long unsigned int); //8-bytes*8 = #bits
  
  int i = 0;
  while (i<e){//1st to the 2nd last limb
    //std::cout<<"i: "<<i<<std::endl;
    //std::cout<<"index: "<<index<<std::endl;
      
    long unsigned int w = limbs[i];
    //std::cout<<"in w: "<<w<<std::endl;
    
    i++;  
    size_t j = 0; //current position in w
    if (carry_size == 0){
      while (b-j>=M){
	long unsigned int q = w>>M;
	
	X[index] = ((w) ^ (q<<M));
	if (isNeg==1)
	  X[index] = -X[index]+prime;
	
	index++;
	w = q;
	j += M; 
	//std::cout<<"w: "<<w<<std::endl;
      }
    }else{
      //std::cout<<"carry_size: "<<carry_size<<std::endl;
      //std::cout<<"carry: "<<carry<<std::endl;
      size_t a = M - carry_size;    
      while (b-j>=a){ 	
	long unsigned int q = (w>>a); //quotient = w / 2^a
	//std::cout<<"q: "<<q<<std::endl;
	//long unsigned int r = (w) ^ (q<<a);
	//std::cout<<"r: "<<r<<std::endl;
	X[index] = (((w) ^ (q<<a)) <<carry_size) + carry; //w % 2^a + carry, 
	if (isNeg==1)
	  X[index] = -X[index]+prime;
	
	w = q;
	j += a;
	index++;
	//std::cout<<"w: "<<w<<std::endl;
	
	a = M;
	carry_size = 0;
	carry = 0;
      }  
    }
    carry_size = b - j;
    carry = w;  
  }
  
  //last limb, its bit size may be smaller than long
  long unsigned int w = limbs[e];
  //std::cout<<"last limb w: "<<w<<std::endl;
  //b = PBPAS::ceil_log2_long(w);
  b = log2(w)+1; //num of bits
  //std::cout<<"log2+1, b: "<<b<<std::endl;
  
  size_t j = 0; //current position in w    
  if (carry_size == 0){
    while (b-j>=M){
      long unsigned int q = w>>M;     
      X[index] = ((w) ^ (q<<M));
      if (isNeg==1)
	X[index] = -X[index]+prime;
      
      index++;
      w = q;
      j += M; 
      //std::cout<<"w: "<<w<<std::endl;
    }
    if (w!=0){ //carry
      X[index] = w;
      if (isNeg==1)
	X[index] = -X[index]+prime;
    }
  }else{
    //std::cout<<"carry_size: "<<carry_size<<std::endl;
    //std::cout<<"carry: "<<carry<<std::endl;

    size_t a = M - carry_size;  
    while (b-j>=a){ 
      long unsigned int q = (w>>a); //quotient = w / 2^a
      //std::cout<<"q: "<<q<<std::endl;
      X[index] = (((w) ^ (q<<a)) <<carry_size) + carry; //w % 2^a + carry, 
      if (isNeg==1)
	X[index] = -X[index]+prime;
      
      w = q;
      j += a;
      index++;
      //std::cout<<"w: "<<w<<std::endl;
      
      a = M;
      carry_size = 0;
      carry = 0;
    } 
    if (w!=0){
      X[index] = (w <<carry_size) + carry; 
      if (isNeg==1)
	X[index] = -X[index]+prime;
    }
  }  
}       
 

BivariatePolynomial * MulSSA::ToBivarMod0(mpz_class *coeff, int d)
{//------------------------------------------------------------------
  BivariatePolynomial *biA = new BivariatePolynomial(d, K, M);
  
  //convert a large int coeff[i] to a univar poly of degree K-1
  //by a base of 2^M
  //store the coefficients of the univar poly to the coeff of y^i
  cilk_for(int i=0; i<d; ++i){

    if (coeff[i]!=0){     
      mpz_class absValue = abs(coeff[i]);
      bool isNegative = coeff[i] < 0;
      
      std::string str = absValue.get_str(2);
      int start = str.size() % M;
      int count =  ceil((float)str.size() / M);
      if(start == 0) start = M;
      
      sfixn *Bi = biA->getCoefficients() + i*K;
      
      int c1 = count-1;

      sfixn tmp = mpz_class(str.substr(0, start), 2).get_ui();
      Bi[c1] = isNegative ? (-tmp+P) : tmp;
      
      for(int i=1; i<=c1; ++i){
	tmp = mpz_class(str.substr((i-1)*M+start, M), 2).get_ui();
	int ci = c1-i;
	Bi[ci] = isNegative ? (-tmp+P) : tmp;		
      }
    }
  }
  
  return biA;
}

//K should be changed to sfixn
//mpz_init2(mpz_t x, mp_bitcnt_t n)
//typedef unsigned long int	mp_bitcnt_t in gmp.h
//mpz_clear(mpz_t x)
//precondition: GMP_LIMB_BITS>=M
//Z is the mpz rep of q with base 2^M
//K = number of q
//mpz_t Z; mpz_init(Z);
//limb_bits: bit bound on ui and vi + 1
//assume q[i]<2^M
void MulSSA::ToMPZ0(sfixn *q, unsigned long limb_bits, mpz_t Zp,
		    int K1, int M1)
{//----------------------------------------------------------
  mpz_init2(Zp, limb_bits);
  Zp->_mp_size = 0;  
  mpz_t Zn;
  mpz_init2(Zn, limb_bits);
  Zn->_mp_size = 0;

  sfixn i     = 0;
  sfixn index = 0;
  sfixn sh    = 0;
  //int M2      = GMP_LIMB_BITS; //= 64

  mp_limb_t Tn = 0; //typedef unsigned long int	mp_limb_t;
  mp_limb_t Tp = 0;

  //precondition: GMP_LIMB_BITS>=M1
  //treat q[0]
  if (q[0]>HALF_P) //negative
    Tn = P - q[0];
  else
    Tp = q[0];

  std::cout<<"Tp: "<<Tp<<std::endl;
  i++;
  std::cout<<"GMP_LIMB_BITS: "<<GMP_LIMB_BITS<<"\n";
  int delta = GMP_LIMB_BITS - M1;
  if (delta == 0) {
    Zp->_mp_d[index] = Tp;    
    Tp = 0;
    Zn->_mp_d[index] = Tn;     
    Tn = 0;

    index++;
    Zp->_mp_size = index;
    Zn->_mp_size = index;

    sh += GMP_LIMB_BITS;
    delta = GMP_LIMB_BITS;
  }
  
  sfixn D = 0;
  while (i<K1){
    D += M1;
    std::cout<<"D: "<<D<<"\n";
    if (delta>=M1) {
      if (q[i]!=0){
	if (q[i]>HALF_P) {//negative, HALF_P = (P+1)/2 or (P-1)/2 ????
	  std::cout<<"P - q[i]: "<<P - q[i]<<std::endl;
	  Tn += ((P - q[i])<<(D-sh));
	  std::cout<<"Tn: "<<Tn<<std::endl; std::cout<<"Tp: "<<Tp<<std::endl;
	}else//positive
	  Tp += (q[i]<<(D-sh));
	//if (i==1) 
	{std::cout<<"i, Tp: "<<i<<", "<<Tp<<std::endl;} //30
      }//else put M number of 0s into both
      delta -= M1;
    }else{//delta<M
      if (q[i]!=0){
	if (q[i]>HALF_P) {//negative,
	  mp_limb_t v = P - q[i];
	  mp_limb_t q_delta_h = v >> delta;
	  mp_limb_t q_delta   = v << (GMP_LIMB_BITS-delta); //(v ^ (q_delta_h << delta));
	  
	  std::cout<<"Tn, q_delta: "<<Tn <<", " <<q_delta <<"\n";
	  Zn->_mp_d[index] = Tn + q_delta;
	  std::cout<<"index, Zn->_mp_d[index]: "<<index<<", "<<Zn->_mp_d[index]<<"\n";
	  Tn = q_delta_h;

	  Zp->_mp_d[index] = Tp;
	  //std::cout<<"index, Zp->_mp_d[index]: "<<index<<", "<<Zp->_mp_d[index]<<"\n";
	  Tp = 0;

	  index++;
	  Zp->_mp_size = index;
	  Zn->_mp_size = index;
	  sh += GMP_LIMB_BITS;
	  delta = GMP_LIMB_BITS - (M1-delta);
	}else{//positive
	  mp_limb_t q_delta_h = (q[i]) >> delta;
	  mp_limb_t q_delta   = (q[i])<<(GMP_LIMB_BITS-delta); //((q[i]) ^ (q_delta_h << delta));

	  //std::cout<<"q_delta_h, q_delta: "<< q_delta_h<<", "<<q_delta<<"\n";
	  Zp->_mp_d[index] = Tp + q_delta;
	  //std::cout<<"Zp->_mp_d[index]: "<<Zp->_mp_d[index]<<"\n";
	  Tp = q_delta_h;
	  Zn->_mp_d[index] = Tn; 
	  Tn = 0;
	  index++;
	  Zp->_mp_size = index;
	  Zn->_mp_size = index;
	  sh += GMP_LIMB_BITS;
	  delta = GMP_LIMB_BITS - (M1-delta);
	}      
      }else{ //q[i]=0
	Zp->_mp_d[index] = Tp;
	Tp = 0;
	Zn->_mp_d[index] = Tn; 
	Tn = 0;
	index++;
	Zp->_mp_size = index;
	Zn->_mp_size = index;
	sh += GMP_LIMB_BITS;
	delta = GMP_LIMB_BITS - (M1-delta);
      }
    }
    if (delta == 0) {
      Zp->_mp_d[index] = Tp;
      Tp = 0;
      Zn->_mp_d[index] = Tn; 
      Tn = 0;
      index++;
      Zp->_mp_size = index;
      Zn->_mp_size = index;
      sh += GMP_LIMB_BITS;
      delta = GMP_LIMB_BITS;
    }
    i++;
  }
  if (Tp!=0) {
    Zp->_mp_d[index] = Tp;
    index++;
    Zp->_mp_size = index;
    Zn->_mp_size = index;
  }
  if (Tn!=0){
    Zn->_mp_d[index] = Tn;
    index++;
    Zp->_mp_size = index;
    Zn->_mp_size = index;
  }

  if (Zn!=0){
    mpz_sub(Zp, Zp, Zn);
    mpz_clear(Zn);
  }else
    mpz_clear(Zn);
  
} //end of ToMPZ0


void MulSSA::ToMPZ(sfixn *q, unsigned long limb_bits, mpz_t Zp,
		   int K1, int M1)
{//----------------------------------------------------------

  mpz_init2(Zp, limb_bits);
  Zp->_mp_size = 0;  
  mpz_t Zn;
  mpz_init2(Zn, limb_bits);
  Zn->_mp_size = 0;

  unsigned __int128 Tn = 0;
  unsigned __int128 Tp = 0;

  //precondition: GMP_LIMB_BITS>=M
  //treat q[0]
  if (q[0]>HALF_P) //negative
    Tn = P - q[0];
  else
    Tp = q[0];

  sfixn i     = 1;
  sfixn index = 0;
  sfixn sh    = 0;
  sfixn D     = 0;
  //GMP_LIMB_BITS is 64
  //mp_limb_t is unsigned long int 
  while (i<K1){
    D += M1; 

    if (D-sh >= GMP_LIMB_BITS){//write out one limb      
      Zp->_mp_d[index] = (mp_limb_t)Tp;
      Tp >>= GMP_LIMB_BITS;

      Zn->_mp_d[index] = (mp_limb_t)Tn;
      Tn >>= GMP_LIMB_BITS;

      index++;
      Zp->_mp_size = index;
      Zn->_mp_size = index;
      sh += GMP_LIMB_BITS;

      if (q[i]!=0){
	if (q[i]>HALF_P) {//negative, HALF_P = (P-1)/2 
	  Tn += (((unsigned __int128)(P-q[i])) << (D-sh));
	}else//positive
	  Tp = Tp + (((unsigned __int128)q[i])<<(D-sh));
      }//else add D-sh number of 0s into both           
    }else{
      if (q[i]!=0){
	if (q[i]>HALF_P) {//negative
	  Tn += (((unsigned __int128)(P-q[i])) << (D-sh));
	}else//positive
	  Tp += (((unsigned __int128)q[i])<<(D-sh));
      }//else put D-sh number of 0s into both      
    }
    i++;
  }
  if (Tp!=0) {
    Zp->_mp_d[index] = (mp_limb_t)Tp;
    Tp >>= GMP_LIMB_BITS;
    Zp->_mp_size += 1;
 
    if (Tp!=0) {
      Zp->_mp_d[(Zp->_mp_size)] = (mp_limb_t)Tp;
      Zp->_mp_size += 1;
    }
  }
  if (Tn!=0){
    Zn->_mp_d[index] = (mp_limb_t)Tn;
    Tn >>= GMP_LIMB_BITS;
    Zn->_mp_size += 1;

    if (Tn!=0){
      Zn->_mp_d[(Zn->_mp_size)] = (mp_limb_t)Tn;
      Zn->_mp_size += 1;
    }
  }
  if (Zn!=0){
    mpz_sub(Zp, Zp, Zn);
    mpz_clear(Zn);
  }else{
    mpz_clear(Zn);
  }

} //end of ToMPZ


int MulSSA::binarySearch(int n, int start, int end){
  int mid = (start + end) / 2;

  if (n == goodN[mid])
    return mid;
  else if(n < goodN[mid] && mid > 0 && n > goodN[mid - 1]) 
    return mid;
  else if(n > goodN[mid] && n <= goodN[end])
    return binarySearch(n, mid, end);
  else if(n < goodN[mid] && n >= goodN[start]) 
    return binarySearch(n, start, mid);
  else{
    std::cout<<"Cannot find good N, K, M.\n";
    exit(1);
  }
}

//d: max(d1,d2)
void MulSSA::determineGoodN(int n, sfixn d){
  int index = binarySearch(n, 0, size-1); 
  N = goodN[index];
  K = 1 << goodK[index];
  M = goodM[index];
  //u_i, v_i has bit bound 2dK2^(M+N) 
  //+1 for the extra limb recommended by mpz_t
  LIMB_BITS = (1 + goodK[index] + logceiling(d) + M + N) + 64; //tight bound
}
