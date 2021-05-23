#include <gmpxx.h>
#include <bpas.h>
#include <iostream>
#include <time.h>

#include "../../../include/ModularPolynomial/DUSP_Support.h"
#include "../../../include/Utils/Unix_Timer.h"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

using namespace NTL;


void convertToZZ_pX (duspoly_t* a, ZZ_pX& f, Prime_ptr* Pptr)
{

    ZZ p = NTL::ZZ(Pptr->prime);
    ZZ_p::init(p);

    if (isZero_spX (a)) {
      SetCoeff(f, 0, ZZ_p(0));
      return;
    }

    for (long i = 0; i <= POLY_LT(a); i++) {
	ZZ_p elm = ZZ_p (a->elems[i]);
	SetCoeff (f, i, elm);
    }
}


void test_mul1 (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr)
{

  int ntrials = 10;
  unsigned long long start;
  float elaspedSerial, sum;

  duspoly_t* a2 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b2 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* c2 = NULL;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    plainMulPolynomialsInForm_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_plainMul] " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    KaratsubaMulPolynomialsInForm_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_karMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  // sum = 0;
  // for (int i = 0; i < ntrials; ++i) {
  //   elaspedSerial = 0;
  //   _startTimer(&start);
  //   fftMulPolynomialsInForm_spX (a2, b2, &c2, Pptr);
  //   _stopTimer(&start, &elaspedSerial);
  //   sum += elaspedSerial;

  //   freePolynomial_spX (&c2); c2 = NULL;
  // }
  // cout << "[dusp_fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  // << "\t\ttime= " << sum/ntrials << endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    // void _fftMulPolynomialsInForm_4step_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr);
    _fftMulPolynomialsInForm_4step_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_4fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    _fftMulPolynomialsInForm_6step_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_6fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;



  ZZ_pX A;
  ZZ_pX B;
  ZZ_pX C;

  duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
  duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);
  convertToZZ_pX (aa2, A, Pptr);
  convertToZZ_pX (bb2, B, Pptr);
  freePolynomial_spX (&aa2);
  freePolynomial_spX (&bb2);

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    mul (C, A, B);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;   
  }
  cout << "[ntl_Mul]       " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;


  freePolynomial_spX (&a2);
  freePolynomial_spX (&b2);
}

void test_mul2 (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr)
{

  int ntrials = 2;
  unsigned long long start;
  float elaspedSerial, sum;

  duspoly_t* a2 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b2 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* c2 = NULL;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    plainMulPolynomialsInForm_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_plainMul] " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    KaratsubaMulPolynomialsInForm_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_karMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  // sum = 0;
  // for (int i = 0; i < ntrials; ++i) {
  //   elaspedSerial = 0;
  //   _startTimer(&start);
  //   fftMulPolynomialsInForm_spX (a2, b2, &c2, Pptr);
  //   _stopTimer(&start, &elaspedSerial);
  //   sum += elaspedSerial;

  //   freePolynomial_spX (&c2); c2 = NULL;
  // }
  // cout << "[dusp_fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  // << "\t\ttime= " << sum/ntrials << endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    // void _fftMulPolynomialsInForm_4step_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr);
    _fftMulPolynomialsInForm_4step_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_4fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    _fftMulPolynomialsInForm_6step_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_6fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;


  ZZ_pX A;
  ZZ_pX B;
  ZZ_pX C;

  duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
  duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);
  convertToZZ_pX (aa2, A, Pptr);
  convertToZZ_pX (bb2, B, Pptr);
  freePolynomial_spX (&aa2);
  freePolynomial_spX (&bb2);

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    mul (C, A, B);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;   
  }
  cout << "[ntl_Mul]       " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  freePolynomial_spX (&a2);
  freePolynomial_spX (&b2);
}



void test_mul_fft (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr)
{

  int ntrials = 10;
  unsigned long long start;
  float elaspedSerial, sum;

  duspoly_t* a2 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b2 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* c2 = NULL;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    // void _fftMulPolynomialsInForm_4step_spX (const duspoly_t* a, const duspoly_t* b, duspoly_t** c, const Prime_ptr* Pptr);
    _fftMulPolynomialsInForm_4step_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_4fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    _fftMulPolynomialsInForm_6step_spX (a2, b2, &c2, Pptr);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;

    freePolynomial_spX (&c2); c2 = NULL;
  }
  cout << "[dusp_6fftMul]   " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;


  ZZ_pX A;
  ZZ_pX B;
  ZZ_pX C;

  duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
  duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);
  convertToZZ_pX (aa2, A, Pptr);
  convertToZZ_pX (bb2, B, Pptr);
  freePolynomial_spX (&aa2);
  freePolynomial_spX (&bb2);

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    elaspedSerial = 0;
    _startTimer(&start);
    mul (C, A, B);
    _stopTimer(&start, &elaspedSerial);
    sum += elaspedSerial;   
  }
  cout << "[ntl_Mul]        " << "[type 1] " << "sz1= " << sz1 << "\tsz2= " << sz2 << "\tpr= " << Pptr->prime 
  << "\t\ttime= " << sum/ntrials << endl;


  freePolynomial_spX (&a2);
  freePolynomial_spX (&b2);
}


int main(int argc, char const *argv[])
{
    polysize_t sz1 = 10;
    polysize_t sz2 = 10;
    prime_t pr =  4179340454199820289; // 4179340454199820289; // 2485986994308513793; // 180143985094819841;
    // 9223372036854602819 => bad prime for bpas-fft :(
    Prime_ptr* Pptr = smallprimefield_get_prime_constants (pr);

    if (argc == 3) {
      sz1 = atoi(argv[1]);
      sz2 = atoi(argv[2]);
    }

  test_mul1 (sz1, sz2, Pptr);
  // test_mul_fft (sz1, sz2, Pptr);

  return 0;
}
