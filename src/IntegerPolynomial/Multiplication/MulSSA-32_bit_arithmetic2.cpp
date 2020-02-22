/**
 * Implementation of 'MulSSA.h'
 *
 * @author Ning Xie and Yuzhen Xie
 */

#include "../../../include/IntegerPolynomial/Multiplication/MulSSA.h"

/*----------------------------------------
Our new 2-convolution multiplication
-----------------------------------------*/
void MulSSA::mul2C(Polynomial *a, Polynomial *b, Polynomial *c) {
	mul2C2(a, b, c);
}

void MulSSA::mul2C2(Polynomial *a, Polynomial *b, Polynomial *c)
{//---------------------------------------------------------------

  //debug
  //a->writeToFile("a.input");
  //b->writeToFile("b.input");
  unsigned long long start = __cilkview_getticks();
  a->coefficientDigits = cilk_spawn a->getBitCount(); 
  b->coefficientDigits = b->getBitCount();
  cilk_sync;
 
  int digitCount = MAX_COEFF_DIGITS(a, b);
  //cout << "mul2C, digitCount: " << digitCount << endl;
  
  sfixn d1 = a->degree;
  sfixn d2 = b->degree;
  sfixn d = (d1 > d2)? d1 : d2;
  determineGoodN(digitCount, d);
  unsigned long long end = __cilkview_getticks();

  cout << "determineGoodN (sec): \t" << (end - start) /1000.f << endl;

  cout<<"d, digitCount, N, K, M, size: " <<d<<", "<<digitCount<<", "<<N<<", "<<K<<", "<<M<<", "<<2+log(K*a->degree)/log(2)+2*M<<endl;

  start = __cilkview_getticks();
  BivariatePolynomial *biA = cilk_spawn 
    ToBivarTwoMod(a->coefficients, d1, prime1, prime2);
  BivariatePolynomial *biB = ToBivarTwoMod(b->coefficients, d2, prime1, prime2);
  cilk_sync;
  end = __cilkview_getticks();
  cout << "ToBivarTwoMod (sec): \t" << (end - start) /1000.f << endl;

  //denug
  //biA->writeToFile("biA.input");
  //biB->writeToFile("biB.input");

  start = __cilkview_getticks();
  // cyclic and negacyclic convolution
  sfixn *ncc1 = PBPAS::TwoConvolutionMod(biA->coefficients, biB->coefficients, biA->coefficients+biA->size, biB->coefficients+biB->size, d1, d2, prime1, prime2, biA->K);
  end = __cilkview_getticks();
  cout << "TwoConvolutionMod (sec): \t" << (end - start) /1000.f << endl;

  biA->freeHeap();
  biB->freeHeap();

  start = __cilkview_getticks();
  CRT_ToGMP_Recovering2(ncc1, d1, d2, c->coefficients);
  end = __cilkview_getticks();
  cout << "CRT_ToGMP_Recovering (sec): \t" << (end - start) /1000.f << endl;

  my_free(ncc1);
}

/**
 * Convert from (a mod p1, b mod p2) to GMP class
 *
 * @size: the size of array a & b
 **/
void MulSSA::CRTtoMPZ(mpz_t zp, sfixn* a, sfixn* b, int size) {
  mpz_t zn;
  mpz_init2(zp, LIMB_BITS);
  mpz_init2(zn, LIMB_BITS);
  
  int idx = 0;
  __int128 Tp = 0, Tn = 0;
  sfixn shifter = 0;
  for (int i = 0; i < size; ++i) {
    long int elem = (long int) ((b[i] - a[i]) % prime2) * U1 % prime2;
    elem = a[i] + elem * prime1;
    if (elem > HALF_P1_P2) { elem -= P1_P2; }
    else if (elem < N_HALF_P1_P2) { elem += P1_P2; }
    if (elem < 0)
      Tn += (__int128) (-elem) << shifter;
    else if (elem)
      Tp += (__int128) elem << shifter;
    shifter += M;
    
    if (shifter >= GMP_LIMB_BITS) {
      if (Tp > 0) {
	zp->_mp_d[idx] = (mp_limb_t) Tp;
	zp->_mp_size = idx + 1;
      }
      else { zp->_mp_d[idx] = 0; }
      Tp >>= GMP_LIMB_BITS;
      if (Tn > 0) {
	zn->_mp_d[idx] = (mp_limb_t) Tn;
	zn->_mp_size = idx + 1;
      }
      else { zn->_mp_d[idx] = 0; }
      Tn >>= GMP_LIMB_BITS;
      shifter -= GMP_LIMB_BITS;
      idx++;
    }
  }
  if (Tp > 0) {
    zp->_mp_d[idx] = (mp_limb_t) Tp;
    zp->_mp_size = idx + 1;
  }
  if (Tn > 0) {
    zn->_mp_d[idx] = (mp_limb_t) Tn;
    zn->_mp_size = idx + 1;
  }
  if (zn > 0) { mpz_sub (zp, zp, zn); }
  mpz_clear(zn);
}

void MulSSA::recoverMPZ(mpz_t res, mpz_t u, mpz_t v) {
	if (!u->_mp_size && !v->_mp_size)
		return;

	mpz_init2(res, 2*N);

	bool carry1 = 0, carry2 = 1, pl;
	int usize = abs(u->_mp_size), vsize = abs(v->_mp_size);
	int maxs = (usize > vsize)? usize : vsize;

	int q = (N - 1) / GMP_LIMB_BITS;
	int r = (N - 1) % GMP_LIMB_BITS;
	res->_mp_d[q] = 0;

	int isxor = mpz_cmpabs(v, u);	// 0: v = u; >0: v > u; <0: v < u;
	mp_limb_t MAXNUM = 0xFFFFFFFFFFFFFFFF;

	// u and v always have the same sign?
	//if (u->_mp_size * v->_mp_size < 0)
	//	cout << "Signs are different.." << endl;

	mpz_t avg;
	mpz_init2(avg, LIMB_BITS);
	for (int i = 0; i < maxs; ++i) {
		// Initialize (u+v) >> 1 and (v-u) << (N-1)
		mp_limb_t a = 0, b = 0;
		if (i < usize) { a = u->_mp_d[i]; }
		if (i < vsize) { b = v->_mp_d[i]; }

		// Compute (u+v) >> 1
		mp_limb_t uv = a + b + carry1;
		carry1 = uv < a;
		uv >>= 1;

		if (i < usize-1) { pl = u->_mp_d[i+1] & 1; }
		else { pl = 0; }
		if (i < vsize-1) { pl = pl ^ (v->_mp_d[i+1] & 1); }
		else { pl = pl ^ 0; }
		pl = pl ^ carry1;
		if (pl)
			uv = uv | ((mp_limb_t) pl << (GMP_LIMB_BITS - 1));

		avg->_mp_d[i] = uv;
		if (avg->_mp_d[i]) { avg->_mp_size = i + 1; }

		// Compute (v-u) << (N-1)
		if (isxor) {
			if (isxor > 0) { a = a ^ MAXNUM; }
			else { b = b ^ MAXNUM; }
			uv = a + b + carry2;
			carry2 = uv < a;

			if (!r) {
				res->_mp_d[q+i] = uv;
				if (res->_mp_d[q+i]) { res->_mp_size = q + i + 1; }
			}
			else {
				res->_mp_d[q+i] = res->_mp_d[q+i] | (uv << r);
				res->_mp_d[q+i+1] = uv >> (GMP_LIMB_BITS - r);
				if (res->_mp_d[q+i+1]) { res->_mp_size = q + i + 2; }
				else if (res->_mp_d[q+i]) { res->_mp_size = q + i + 1; }
			}
		}
	}
	if (avg->_mp_size && u->_mp_size < 0)
		avg->_mp_size = -avg->_mp_size;
	if (res->_mp_size) {
		for (int i = 0; i < q; ++i)
			res->_mp_d[i] = 0;
		if (u->_mp_size * isxor < 0) { res->_mp_size = -res->_mp_size; }
	}

	mpz_add (res, res, avg);
	mpz_clear(avg);
}

/*---------------------------------------------------------------------
  for both C^+ (ncc) and C^- (cc), for each coefficient of y^i,
  which is a large integer encoded as a polynomial in x, i.e.
  ncc1_i(x), ncc2_i(x), cc1_i(x), cc2_i(x)
  (1) for each coefficient of x, apply CRT for the two prime numbers, 
      get ncc_i and cc_i
  (2) convert ncc_i and cc_i to GMP, get u_i and v_i
  (3) compute c_i = (u_i+v_i)/2+(-u_i+v_i)/2*2^N
------------------------------------------------------------------------*/
void MulSSA::CRT_ToGMP_Recovering2(sfixn *ncc1, int d1, int d2, mpz_class *c) {
  int csize = d1 + d2 - 1;
  int fullsize = K * csize;
  sfixn *ncc2 = ncc1 + fullsize;
  sfixn *cc1 = ncc2 + fullsize;
  sfixn *cc2 = cc1 + fullsize;
  
#pragma cilk_grainsize = 1024;
  cilk_for (int i = 0; i < csize; ++i) {
    int Ks = i * K;
    mpz_t u, v;
    CRTtoMPZ(u, ncc1+Ks, ncc2+Ks, K);
    CRTtoMPZ(v, cc1+Ks, cc2+Ks, K);

    mpz_add (c[i].get_mpz_t(), u, v);
    c[i] >>= 1;
    mpz_sub (v, v, u);
    c[i] += mpz_class (v) << (N - 1);

    //recoverMPZ(c[i].get_mpz_t(), u, v);

    mpz_clear(u);
    mpz_clear(v);
  }
}

void MulSSA::CRT_ToGMP_Recovering(sfixn *ncc1, int d1, int d2, mpz_class *c)
{//------------------------------------------------------------------------
  int N1 = N-1;
  int Ke = K-1;
  int K2 = K-2;

  int ys = d1+d2-1;
  int xys = K * ys;
  int *ncc2 = ncc1 + xys;
  int *cc1 = ncc2 + xys;
  int *cc2 = cc1 + xys;

  /*
  ofstream ofs("ncc1.out", ofstream::out);	
  for(int j = 0; j < xys; ++j) 
    ofs << ncc1[j] << " ";
  ofs << "\n\n";
  for(int j = 0; j < xys; ++j) 
    ofs << ncc2[j] << " ";
  ofs << "\n\n";
  for(int j = 0; j < xys; ++j) 
    ofs << cc1[j] << " ";
  ofs << "\n\n";
  for(int j = 0; j < xys; ++j) 
    ofs << cc2[j] << " ";
  ofs << "\n\n";
  ofs.close();
  */

  //ofstream ofs1("crt_gmp.out", ofstream::out);
  //ofstream ofs2("crt.out", ofstream::out);
 
  /*
  std::cout<<"prime1, prime1, P1_P2: "<<prime1<<", "<<prime2<<", "<<P1_P2<<std::endl;
  std::cout<<"HALF_P1_P2, N_HALF_P1_P2: "<<HALF_P1_P2<<", "<<N_HALF_P1_P2<<std::endl;
  std::cout<<"P1_U1, P2_U2: "<<P1_U1<<", "<<P2_U2<<std::endl;
  */

  cilk_for(int i=0; i<ys; ++i){//enough parallelism for 48-core?
    int Ks = i*K;

    //CRT and then to GMP by evaluation via Horner's rule    
    //ncc_i
    int *ncc1_i = ncc1+Ks;
    int *ncc2_i = ncc2+Ks;
    long int ncc_i_gmp = 0;
    long int crt_i = 0;
    if (ncc1_i[Ke]!=0 || ncc2_i[Ke]!=0){
      //crt_i = (ncc1_i[Ke]+(ncc2_i[Ke]-ncc1_i[Ke])*U1_P1) % P1_P2; 
      mpz_class tmp = (ncc1_i[Ke]*P2_U2 + ncc2_i[Ke]*P1_U1) % P1_P2;
      crt_i = tmp.get_si(); 
      if(crt_i > HALF_P1_P2) 
      	crt_i -= P1_P2;
      if(crt_i < N_HALF_P1_P2)
	crt_i += P1_P2;

      //ofs2 << crt_i <<" ";

      ncc_i_gmp = crt_i;
    }
    for(int j=K2; j>=0; --j){
      //if(ncc_i_gmp!=0)
      ncc_i_gmp <<= M; //
      if (ncc1_i[j]!=0 || ncc2_i[j]!=0){
	//crt_i = (ncc1_i[j]+(ncc2_i[j]-ncc1_i[j])*U1_P1)%P1_P2;
	mpz_class tmp = (ncc1_i[j]*P2_U2 + ncc2_i[j]*P1_U1) % P1_P2;
	crt_i = tmp.get_si(); 

	if(crt_i > HALF_P1_P2) 
	  crt_i -= P1_P2;
	if(crt_i < N_HALF_P1_P2) 
	  crt_i += P1_P2;

	//ofs2 << crt_i <<" ";

	ncc_i_gmp += crt_i;
      }
    }

    //cc_i
    int *cc1_i = cc1+Ks;
    int *cc2_i = cc2+Ks;
    mpz_class cc_i_gmp = 0;
    if (cc1_i[Ke]!=0 || cc2_i[Ke]!=0){
      //crt_i = (cc1_i[Ke]+(cc2_i[Ke]-cc1_i[Ke])*U1_P1)%P1_P2; 
      mpz_class tmp = (cc1_i[Ke]*P2_U2 + cc2_i[Ke]*P1_U1) % P1_P2;
      crt_i = tmp.get_si(); 

      if(crt_i > HALF_P1_P2) 
      	crt_i -= P1_P2;    
      if(crt_i < N_HALF_P1_P2) 
      	crt_i += P1_P2;  

      //ofs2 << crt_i <<" ";

      cc_i_gmp = crt_i;
    }
    for(int j=K2; j>=0; --j){
      //if(cc_i_gmp!=0)
      cc_i_gmp <<= M; 
      if (cc1_i[j]!=0 || cc2_i[j]!=0){
	//crt_i = (cc1_i[j]+(cc2_i[j]-cc1_i[j])*U1_P1)%P1_P2;
	mpz_class tmp = (cc1_i[j]*P2_U2 + cc2_i[j]*P1_U1) % P1_P2;
	crt_i = tmp.get_si(); 

	if(crt_i > HALF_P1_P2) 
	  crt_i -= P1_P2;
	if(crt_i < N_HALF_P1_P2) 
	  crt_i += P1_P2; 

	//ofs2 << crt_i <<" ";

	cc_i_gmp += crt_i;
      }
    }
    //ofs2 << "\n";

   //ToDo: combine the two for-loops? compare the timing!

    //c=a*b
    mpz_class c1 = ncc_i_gmp+cc_i_gmp;   
    c1 >>= 1;   
    mpz_class c2 = cc_i_gmp-ncc_i_gmp;    
    c2 <<=N1;    
    c[i] = c1 + c2; 

    //ofs1 << ncc_i_gmp << " " <<ncc_i_gmp<<" " <<c1<<" " <<c2<<" "<<c[i] <<"\n";

    //gmp error
    //c[i] = (ncc_i_gmp+cc_i_gmp)<<1 + (cc_i_gmp-ncc_i_gmp)>>N1; 
  }
  //ofs2.close();
  //ofs1.close();
  
}

/**
 * Convert a univariate polynomial with large coefficients to 
 * a bivariate polynomial representation with relatively 
 * small coefficients.
 *
 * univariate --> a1    * y ^ 0 + a2    * y ^ 1 + ... + ad    * y ^ (d-1)
 * bivariate  --> A1(x) * y ^ 0 + A2(x) * y ^ 1 + ... + Ad(x) * y ^ (d-1)
 *
 * Ai(x) = b1 * x ^ 0 + b2 * x ^ 1 + ... + bK * x ^ (K-1)
 *
 * @param coeff the coefficients of large integers
 * @param d partial degree of the univariate polynomial plus one
 */
BivariatePolynomial * MulSSA::ToBivarTwoMod(mpz_class *coeff, int d, 
					    int p1, int p2)
{//-------------------------------------------------------------------
  BivariatePolynomial *biA = new BivariatePolynomial(d, K, M, 2);//double size
  int *A2 = biA->coefficients + biA->size; 
  // biA->coefficients has length of K*d*2

  //convert a large int coeff[i] to a univar poly of degree K-1
  //by a base of 2^M
  //store the coefficients of the univar poly to the coeff of y^i
  #pragma cilk_grainsize = 8192;
  cilk_for(int i=0; i<d; ++i){
    if (coeff[i]!=0){  
      sfixn ci = i*K;
      sfixn *Ai = biA->coefficients + ci;
      sfixn *A2i = A2 + ci;
      mpzToPolyTwoMod((coeff[i]).get_mpz_t(), M, p1, p2, Ai, A2i);
    }
  }
  
  return biA;
}

//---------------------------------------------------
//just for binary representation mpz
//precondition: z!=0, number of bits of z / M <= K
//K = size of X1 
void MulSSA::mpzToPolyTwoMod(mpz_t z, int M, 
			     sfixn prime1, sfixn prime2,
			     sfixn *X1, sfixn *X2) 
{//----------------------------------------------------------  
  int s = z->_mp_size; //abs(s)>=1
  int isNeg = 0;
  if (s<0) {
    isNeg = 1;
    s = -s;
  }

  //cout<<"s: "<<s<<endl;

  mp_ptr limbs = z->_mp_d;
  mp_limb_t carry = 0;
  mp_size_t carry_size = 0;

  int e = s-1; 
  int index = 0;

  //all the limbs except for the last one are full-size (64 bits)
  size_t b = GMP_LIMB_BITS;
  //size_t b = 8*sizeof(long unsigned int); //8-bytes*8 = #bits
  
  int i = 0;
  while (i<e){//1st to the 2nd last limb
    //cout<<"i: "<<i<<endl;
    //cout<<"index: "<<index<<endl;
    // mp_limb_t
    long unsigned int w = limbs[i];
    //cout<<"in w: "<<w<<endl;
    
    i++;  
    size_t j = 0; //current position in w
    if (carry_size == 0){
      while (b-j>=M){
	long unsigned int q = w>>M;
	X1[index] = (sfixn) ((w) ^ (q<<M));	
	if (isNeg==1){
	  X1[index] = -X1[index];
	  X2[index] = X1[index]+prime2;
	  X1[index] = X1[index]+prime1;	  
	}else
	  X2[index] = X1[index];

	index++;
	w = q;
	j += M; 
	//cout<<"w: "<<w<<endl;
      }
    }else{
      //cout<<"carry_size: "<<carry_size<<endl;
      //cout<<"carry: "<<carry<<endl;
      size_t a = M - carry_size;    
      while (b-j>=a){ 	
	long unsigned int q = (w>>a); //quotient = w / 2^a
	//cout<<"q: "<<q<<endl;
	//long unsigned int r = (w) ^ (q<<a);
	//cout<<"r: "<<r<<endl;
	X1[index] = (sfixn) ((((w) ^ (q<<a)) <<carry_size) + carry); //w % 2^a + carry, 
	if (isNeg==1){
	  X1[index] = -X1[index];
	  X2[index] = X1[index]+prime2;
	  X1[index] = X1[index]+prime1;	  
	}else
	  X2[index] = X1[index];
	
	w = q;
	j += a;
	index++;
	//cout<<"w: "<<w<<endl;
	
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
  //cout<<"last limb w: "<<w<<endl;
  //b = PBPAS::ceil_log2_long(w);
  b = log2(w)+1; //num of bits
  //cout<<"log2+1, b: "<<b<<endl;
  
  size_t j = 0; //current position in w    
  if (carry_size == 0){
    while (b-j>=M){
      long unsigned int q = w>>M;  
 
      X1[index] = (sfixn) ((w) ^ (q<<M));
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];

      index++;
      w = q;
      j += M; 
      //cout<<"w: "<<w<<endl;
    }
    if (w!=0) {//carry
      X1[index] = (sfixn) w;
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];
    }
 
  }else{
    //cout<<"carry_size: "<<carry_size<<endl;
    //cout<<"carry: "<<carry<<endl;

    size_t a = M - carry_size;  
    while (b-j>=a){ 
      long unsigned int q = (w>>a); //quotient = w / 2^a
      //cout<<"q: "<<q<<endl;
      X1[index] = (sfixn) ((((w) ^ (q<<a)) <<carry_size) + carry); //w % 2^a + carry, 
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];

      w = q;
      j += a;
      index++;
      //cout<<"w: "<<w<<endl;
      
      a = M;
      carry_size = 0;
      carry = 0;
    } 

    if (w!=0){
      X1[index] = (sfixn) ((w <<carry_size) + carry); 
      if (isNeg==1){
	X1[index] = -X1[index];
	X2[index] = X1[index]+prime2;
	X1[index] = X1[index]+prime1;	  
      }else
	X2[index] = X1[index];
    }
  } 
}       
 

BivariatePolynomial * MulSSA::ToBivarTwoMod0(mpz_class *coeff, int d, 
					     int p1, int p2)
{//-------------------------------------------------------------------
  BivariatePolynomial *biA = new BivariatePolynomial(d, K, M, 2);
  int *B2 = biA->coefficients + biA->size; 
  // biA->coefficients has length of K*d*2

  //convert a large int coeff[i] to a univar poly of degree K-1
  //by a base of 2^M
  //store the coefficients of the univar poly to the coeff of y^i
  cilk_for(int i=0; i<d; ++i){
    if (coeff[i]!=0){  
      mpz_class absValue = abs(coeff[i]);
      bool isNegative = coeff[i] < 0;
      
      string str = absValue.get_str(2);
      int start = str.size() % M;
      int count =  ceil((float)str.size() / M);
      if(start == 0) start = M;

      int *B1_i = biA->coefficients + i*K;
      int *B2_i = B2 + i*K;
      int c1 = count-1;

      int tmp = mpz_class(str.substr(0, start), 2).get_ui();
      B1_i[c1] = isNegative ? (-tmp+p1) : tmp;
      B2_i[c1] = isNegative ? (-tmp+p2) : tmp;
      
      for(int i=1; i<=c1; ++i){
	tmp = mpz_class(str.substr((i-1)*M+start, M), 2).get_ui();
	int ci = c1-i;
	B1_i[ci] = isNegative ? (-tmp+p1) : tmp;
	B2_i[ci] = isNegative ? (-tmp+p2) : tmp;	
      }
    }
  }
  
  return biA;
}

int MulSSA::binarySearch(int n, int start, int end) {
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

void MulSSA::determineGoodN(int n, sfixn d){
  int index = binarySearch(n, 0, size-1); 
  N = goodN[index];
  K = 1 << goodK[index];
  M = goodM[index];	
  LIMB_BITS = (1 + goodK[index] + logceiling(d) + M + N) + 64;
}
