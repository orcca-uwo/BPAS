#include <gmpxx.h>
#include <bpas.h>
#include <iostream>
#include <time.h>
#include <fstream>

#include "../../../include/ModularPolynomial/DUSP_Support.h"
#include "../../../include/ModularPolynomial/DUSP_FFT_Support.h"
#include "../../../include/ModularPolynomial/DUSP_Support_Test.h"
#include "../../../include/ModularPolynomial/DUSP_Support_Factoring.h"
#include "../../../include/ModularPolynomial/DUSP_NTL_Support.h"
#include "../../../include/IntegerPolynomial/DUZP_Support.h"
#include "../../../include/RegularChain/regularchain.hpp"
#include "../../../include/ModularPolynomial/DUSP_Parallel.hpp"

#include "../../../include/Utils/Parallel/ExecutorThreadPool.hpp"
#include "../../../include/Utils/Unix_Timer.h"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ_p.h> 
#include <NTL/pair_ZZ_pX_long.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace NTL;
using namespace std;

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


void convertToZZX (duspoly_t* a, ZZX& f)
{

    if (isZero_spX (a)) {
      SetCoeff(f, 0, 0);
      return;
    }

    for (long i = 0; i <= POLY_LT(a); i++) {
	    SetCoeff (f, i, a->elems[i]);
    }
}

///////////////////////////////////////////////////////////////////////////////////////

void test_freePoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);
    duspoly_t* c = makePolynomial_spX (sz1);
    duspoly_t* d = makePolynomial_spX (sz2);

    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);
    freePolynomial_spX (&d);

    if (isZero_spX (a) + isZero_spX (b) + isZero_spX (c) + isZero_spX (d) == -4) {
        fprintf (stderr, "[DUSP] freePolynomial... PASS\n");
    } else {

        fprintf(stderr, "[DUSP] isZero_spX (a) = %d\n", isZero_spX (a));
        fprintf(stderr, "[DUSP] isZero_spX (b) = %d\n", isZero_spX (b));
        fprintf(stderr, "[DUSP] isZero_spX (c) = %d\n", isZero_spX (c));
        fprintf(stderr, "[DUSP] isZero_spX (d) = %d\n", isZero_spX (d));

        fprintf (stderr, "[DUSP] a->alloc = %lu\n", a->alloc);
        fprintf (stderr, "[DUSP] b->alloc = %lu\n", b->alloc);
        fprintf (stderr, "[DUSP] c->alloc = %lu\n", c->alloc);
        fprintf (stderr, "[DUSP] d->alloc = %lu\n", d->alloc);

        fprintf (stderr, "[DUSP] a->lt = %lu\n", a->lt);
        fprintf (stderr, "[DUSP] b->lt = %lu\n", b->lt);
        fprintf (stderr, "[DUSP] c->lt = %lu\n", c->lt);
        fprintf (stderr, "[DUSP] d->lt = %lu\n", d->lt);

        fprintf (stderr, "[DUSP] a = ");
        printPolynomial_spX (a, Pptr);
        fprintf (stderr, "[DUSP] b = ");
        printPolynomial_spX (b, Pptr);
        fprintf (stderr, "[DUSP] c = ");
        printPolynomial_spX (c, Pptr);
        fprintf (stderr, "[DUSP] d = ");
        printPolynomial_spX (d, Pptr);

        fprintf (stderr, "[DUSP] freePolynomial... FAIL\n");
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////////////

void test_deepCopyPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag) {

    duspoly_t* a = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);

    duspoly_t* acp = deepCopyPolynomial_spX (a);
    duspoly_t* bcp = deepCopyPolynomial_spX (b);

    if (!isEqual_spX (a, acp)) {
        fprintf (stderr, "[DUSP] a = ");
        printPolynomial_spX (a, Pptr);
        fprintf (stderr, "[DUSP] dcp(a) = ");
        printPolynomial_spX (acp, Pptr);

        fprintf (stderr, "[DUSP] deepCopyPolynomial... FAIL\n");
        exit(1);
    }

    if (!isEqual_spX (b, bcp)) {
        fprintf (stderr, "[DUSP] b = ");
        printPolynomial_spX (b, Pptr);
        fprintf (stderr, "[DUSP] dcp(b) = ");
        printPolynomial_spX (bcp, Pptr);

        fprintf (stderr, "[DUSP] deepCopyPolynomial... FAIL\n");
        exit(1);
    }

    freePolynomial_spX (&a);
    freePolynomial_spX (&b);

    if (isZero_spX (acp) || isZero_spX (bcp)) {

        fprintf (stderr, "[DUSP] a = ");
        printPolynomial_spX (a, Pptr);
        fprintf (stderr, "[DUSP] dcp(a) = ");
        printPolynomial_spX (acp, Pptr);

        fprintf (stderr, "[DUSP] b = ");
        printPolynomial_spX (b, Pptr);
        fprintf (stderr, "[DUSP] dcp(b) = ");
        printPolynomial_spX (bcp, Pptr);

        fprintf (stderr, "[DUSP] deepCopyPolynomial... FAIL\n");
        exit(1);
    }

    freePolynomial_spX (&acp);
    freePolynomial_spX (&bcp);
    fprintf (stderr, "[DUSP] deepCopyPolynomial... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

void test_makePoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag) {

    duspoly_t* zero = zeroPolynomial_spX ();

    if (isZero_spX (zero) != 1) {
        fprintf (stderr, "[DUSP] zero = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] zeroPolynomial... FAIL\n");
        fprintf (stderr, "[DUSP] isZero... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&zero);

    if (!isZero_spX (zero)) {
        fprintf (stderr, "[DUSP] zero = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] free(zero)... FAIL\n");

        exit(1);
    }

    do {
        zero = constPolynomialInForm_spX (sz1*sz2+ (int) rand()%7, Pptr);
        if (isZero_spX (zero)) {
            freePolynomial_spX (&zero);
        } else { break; }
    } while (1);

    if (!isConstant_spX (zero)) {

        fprintf (stderr, "[DUSP] const = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] constPolynomial... FAIL\n");
        fprintf (stderr, "[DUSP] isConstant... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&zero);

    do {
        zero = constPolynomial_spX (sz1*sz2 + (int) rand()%7, Pptr);
        if (isZero_spX (zero)) {
            freePolynomial_spX (&zero);
        } else { break; }
    } while (1);

    if (!isConstant_spX (zero)) {

        fprintf (stderr, "[DUSP] const = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] constPolynomialInForm... FAIL\n");
        fprintf (stderr, "[DUSP] isConstant... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&zero);

    zero = constPolynomialInForm_spX (1, Pptr);

    if (!isOneInForm_spX (zero, Pptr)) {

        fprintf (stderr, "[DUSP] one = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] isOneInForm... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&zero);

    zero = constPolynomial_spX (1, Pptr);

    if (!isOne_spX (zero)) {

        fprintf (stderr, "[DUSP] one = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] isOne... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&zero);

    zero = XpolynomialInForm_spX (Pptr);

    if (!isXInForm_spX (zero, Pptr)) {

        fprintf (stderr, "[DUSP] XPoly = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] Xpolynomial... FAIL\n");
        fprintf (stderr, "[DUSP] isXInForm... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&zero);

    zero = Xpolynomial_spX ();

    if (!isX_spX (zero)) {

        fprintf (stderr, "[DUSP] XPoly = ");
        printPolynomial_spX (zero, Pptr);

        fprintf (stderr, "[DUSP] XpolynomialInForm... FAIL\n");
        fprintf (stderr, "[DUSP] isX... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&zero);

    fprintf (stderr, "[DUSP] zeroPolynomial... PASS\n");
    fprintf (stderr, "[DUSP] isZero... PASS\n");
    fprintf (stderr, "[DUSP] constPolynomial... PASS\n");
    fprintf (stderr, "[DUSP] constPolynomialInForm... PASS\n");
    fprintf (stderr, "[DUSP] isConstant... PASS\n");
    fprintf (stderr, "[DUSP] isOne... PASS\n");
    fprintf (stderr, "[DUSP] isOneInForm... PASS\n");
    fprintf (stderr, "[DUSP] Xpolynomial... PASS\n");
    fprintf (stderr, "[DUSP] XpolynomialInForm... PASS\n");
    fprintf (stderr, "[DUSP] isX... PASS\n");
    fprintf (stderr, "[DUSP] isXInForm... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

void test_normalize (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = makePolynomial_spX (sz1);
    a->elems[0] = 1;
    a->lt = 0;
    duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);


    normalize_spX (&a);

    if (a->alloc > 16) {

      fprintf(stderr, "a->alloc = %lu\n", a->alloc);
      fprintf (stderr, "[DUSP] normalize_spX... FAIL\n");
      exit (1);
    }

    freePolynomial_spX (&a);

    a = makePolynomial_spX (10000);
    for (polysize_t i = 0; i < 1000; i++) {
      a->elems[i] = i+1;
    }

    setLT_spX (&a, 0);

    if (POLY_LT(a) != 999) {

      fprintf(stderr, "a = \n");
      printPolynomial_spX (a, Pptr);

      fprintf(stderr, "a->alloc = %lu\n", a->alloc);
      fprintf(stderr, "a->lt = %lu\n", a->lt);

      fprintf (stderr, "[DUSP] setLT_spX... FAIL\n");
      exit (1);
    }

    normalize_spX (&a);

    if (POLY_ALLOC(a) != 1000) {

      fprintf(stderr, "a = \n");
      printPolynomial_spX (a, Pptr);

      fprintf (stderr, "[DUSP] normalize_spX... FAIL\n");
      exit (1);
    }

    for (polysize_t i = 0; i < b->alloc; i++) {
      b->elems[i] = 0;
    }

    setLT_spX (&b, 1);

    if (!isZero_spX (b)) {

      fprintf(stderr, "b = \n");
      printPolynomial_spX (b, Pptr);

      fprintf (stderr, "[DUSP] normalize_spX... FAIL\n");
      fprintf (stderr, "[DUSP] setLT_spX... FAIL\n");
      exit (1);
    }

    freePolynomial_spX (&a);
    freePolynomial_spX (&b);

    fprintf(stderr, "[DUSP] reallocatePolynomial... PASS\n");
    fprintf (stderr, "[DUSP] normalize... PASS\n");
    fprintf (stderr, "[DUSP] setLT... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

void test_convert (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  duspoly_t* a = NULL;
  duspoly_t* a_in = convertPolynomialToMontgomery_spX (a, Pptr);
  duspoly_t* a_out = convertPolynomialFromMontgomery_spX (a_in, Pptr);

  if (a_in != NULL) {
    fprintf(stderr, "[DUSP] (1) convertPolynomialToMontgomery... FAIL\n");
    exit (1);
  }

  if (a_out != NULL) {
    fprintf(stderr, "[DUSP] (1) convertPolynomialFromMontgomery... FAIL\n");
    exit (1);
  }

  a = randomPolynomialInForm_spX (sz1, Pptr);
  freePolynomial_spX (&a);
  a_in = convertPolynomialToMontgomery_spX (a, Pptr);
  a_out = convertPolynomialFromMontgomery_spX (a, Pptr);

  if (a_in != NULL) {
    fprintf(stderr, "[DUSP] (2) convertPolynomialToMontgomery... FAIL\n");
    exit (1);
  }

  if (a_out != NULL) {
    fprintf(stderr, "[DUSP] (2) convertPolynomialFromMontgomery... FAIL\n");
    exit (1);
  }


  a = randomPolynomialInForm_spX (sz1, Pptr);
  a_in = convertPolynomialToMontgomery_spX (a, Pptr);
  a_out = convertPolynomialFromMontgomery_spX (a_in, Pptr);

  if (!isEqual_spX (a, a_out)) {

    fprintf(stderr, "a = ");
    printPolynomial_spX (a, Pptr);

    fprintf(stderr, "a_in = ");
    printPolynomial_spX (a_in, Pptr);

    fprintf(stderr, "a_out = ");
    printPolynomial_spX (a_out, Pptr);

    fprintf(stderr, "[DUSP] (3) convertPolynomialToMontgomery... FAIL\n");
    fprintf(stderr, "[DUSP] (3) convertPolynomialFromMontgomery... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a);
  freePolynomial_spX (&a_in);
  freePolynomial_spX (&a_out);


  a = randomFullyDensePolynomial_spX (sz2, Pptr);
  a_out = deepCopyPolynomial_spX (a);
  convertPolynomialToMontgomery_spX_inp (&a_out, Pptr);
  convertPolynomialFromMontgomery_spX_inp (&a_out, Pptr);

  if (!isEqual_spX (a, a_out)) {
    fprintf(stderr, "[DUSP] (4) convertPolynomialToMontgomery_inp... FAIL\n");
    fprintf(stderr, "[DUSP] (4) convertPolynomialFromMontgomery_inp... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a);
  freePolynomial_spX (&a_out);

  convertPolynomialToMontgomery_spX_inp (&a, Pptr);
  if (a != NULL) {
    fprintf(stderr, "[DUSP] (5) convertPolynomialToMontgomery_inp... FAIL\n");
    exit (1);
  }

  convertPolynomialFromMontgomery_spX_inp (&a, Pptr);
  if (a != NULL) {
    fprintf(stderr, "[DUSP] (5) convertPolynomialFromMontgomery_inp... FAIL\n");
    exit (1);
  }


  fprintf (stderr, "[DUSP] convertPolynomialToMontgomery... PASS\n");
  fprintf (stderr, "[DUSP] convertPolynomialFromMontgomery... PASS\n");


}

///////////////////////////////////////////////////////////////////////////////////////

void test_shiftPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = NULL;
    duspoly_t* a_in = NULL;
    elem_t a_elems[5] = {0,0,0,0,1}; // x^4
    elem_t one = smallprimefield_convert_in (1, Pptr);

    elem_t a_in_elems[5] = {0,0,0,0,one}; // x^4

    setCoefs_spX (&a, a_elems, 5, Pptr);
    setCoefsInForm_spX (&a_in, a_in_elems, 5);

    if (!isEqual_spX (a, a_in)) {

  fprintf (stderr, "a = ");
  printPolynomial_spX (a, Pptr);

  fprintf (stderr, "a_in = ");
  printPolynomial_spX (a_in, Pptr);

  fprintf (stderr, "[DUSP] setCoefs... FAIL\n");
  fprintf (stderr, "[DUSP] setCoefsInForm... FAIL\n");

  exit (1);
    }

    freePolynomial_spX (&a); a = NULL;

    duspoly_t* sa = NULL;
    rightShiftPolynomial_spX (a_in, &sa, 4);

    if (!isOneInForm_spX (sa, Pptr)) {
  fprintf (stderr, "a_in = ");
  printPolynomial_spX (a_in, Pptr);

  fprintf (stderr, "(a >> 4) = ");
  printPolynomial_spX (sa, Pptr);

  fprintf (stderr, "[DUSP] rightShiftPolynomial... FAIL\n");

  exit(1);

    }

    leftShiftPolynomial_spX (sa, &a, 4);


    if (!isEqual_spX (a, a_in)) {
  fprintf (stderr, "a_in = ");
  printPolynomial_spX (a_in, Pptr);

  fprintf (stderr, "(a_in << 4) = ");
  printPolynomial_spX (a, Pptr);

  fprintf (stderr, "[DUSP] leftShiftPolynomial... FAIL\n");

  exit(1);

    }

    freePolynomial_spX (&a); a = NULL;
    freePolynomial_spX (&sa); sa = NULL;

    rightShiftPolynomial_spX (a_in, &sa, 10);

    if (!isZero_spX (sa)) {
  fprintf (stderr, "a_in = ");
  printPolynomial_spX (a_in, Pptr);

  fprintf (stderr, "(a >> 10) = ");
  printPolynomial_spX (sa, Pptr);

  fprintf (stderr, "[DUSP] rightShiftPolynomial... FAIL\n");

  exit(1);

    }

    freePolynomial_spX (&sa); sa = NULL;


    fprintf (stderr, "[DUSP] setCoefs... PASS\n");
    fprintf (stderr, "[DUSP] setCoefsInForm... PASS\n");
    fprintf (stderr, "[DUSP] rightShiftPolynomial... PASS\n");
    fprintf (stderr, "[DUSP] leftShiftPolynomial... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

void test_splitPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b = NULL;

    duspoly_t* a_l = NULL;
    duspoly_t* a_r = NULL;
    duspoly_t* a_rb = NULL;

    leftSplitPolynomial_spX  (a, &a_l, sz1/2);
    rightShiftPolynomial_spX (a, &a_r, sz1/2);
    leftShiftPolynomial_spX  (a_r, &a_rb, sz1/2);

    addPolynomialsInForm_spX (a_l, a_rb, &b, Pptr);

    if (!isEqual_spX (a, b)) {
  fprintf (stderr, "a = ");
  printPolynomial_spX (a, Pptr);

  fprintf (stderr, "left(a) = ");
  printPolynomial_spX (a_l, Pptr);

  fprintf (stderr, "right(a) = ");
  printPolynomial_spX (a_rb, Pptr);

  fprintf (stderr, "a_l + a_r = ");
  printPolynomial_spX (b, Pptr);

  fprintf (stderr, "[DUSP] leftSplitPolynomial... FAIL\n");

  exit(1);
    }

    freePolynomial_spX (&a);
    freePolynomial_spX (&a_l);
    freePolynomial_spX (&a_r);
    freePolynomial_spX (&a_rb);
    freePolynomial_spX (&b);

    fprintf (stderr, "[DUSP] leftSplitPolynomial... PASS\n");
}


///////////////////////////////////////////////////////////////////////////////////////

void test_addPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = NULL;
    duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);
    duspoly_t* c = NULL;

    duspoly_t* aa;
    duspoly_t* bb;
    duspoly_t* cc;


    addPolynomialsInForm_spX (a, b, &c, Pptr);

    ZZ_pX A;
    ZZ_pX B;
    ZZ_pX C;
    ZZ_pX D;

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    add (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (1) addPolynomialsInForm... FAIL\n");
    }


    //////////////////////////////
    a = randomPolynomialInForm_spX (sz1, Pptr);
    freePolynomial_spX (&c);

    addPolynomialsInForm_spX (a, b, &c, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    add (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (2) addPolynomialsInForm... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = randomPolynomialInForm_spX (2, Pptr);
    b = randomPolynomialInForm_spX (sz2, Pptr);

    addPolynomialsInForm_spX (a, b, &c, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    add (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (3) addPolynomialsInForm... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = NULL;
    b = randomPolynomialInForm_spX (sz2, Pptr);
    c = NULL;

    addPolynomialsInForm_spX_inp (&c, b, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    add (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (4) addPolynomialsInForm_inp... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = randomPolynomialInForm_spX (sz1, Pptr);
    b = randomPolynomialInForm_spX (sz2, Pptr);
    c = deepCopyPolynomial_spX (a);

    addPolynomialsInForm_spX_inp (&c, b, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    add (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (5) addPolynomialsInForm_inp... FAIL\n");
    }


    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = randomPolynomialInForm_spX (sz1, Pptr);
    b = randomPolynomialInForm_spX (2, Pptr);
    c = deepCopyPolynomial_spX (a);

    addPolynomialsInForm_spX_inp (&c, b, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    add (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (6) addPolynomialsInForm_inp... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    duspoly_t* a7 = randomPolynomialInForm_spX (sz1, Pptr);
    duspoly_t* b7 = NULL;
    duspoly_t* c7 = deepCopyPolynomial_spX (a7);

    addPolynomialsInForm_spX_inp (&c7, b7, Pptr);

    duspoly_t* cc7 = convertPolynomialFromMontgomery_spX (c7, Pptr);
    duspoly_t* aa7 = convertPolynomialFromMontgomery_spX (a7, Pptr);
    duspoly_t* bb7 = convertPolynomialFromMontgomery_spX (b7, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc7, C, Pptr);
    convertToZZ_pX (aa7, A, Pptr);
    convertToZZ_pX (bb7, B, Pptr);

    freePolynomial_spX (&cc7);
    freePolynomial_spX (&aa7);
    freePolynomial_spX (&bb7);

    add (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (7) addPolynomialsInForm_inp... FAIL\n");
    }
    freePolynomial_spX (&a7);
    freePolynomial_spX (&b7);
    freePolynomial_spX (&c7);

  fprintf(stderr, "[DUSP] addPolynomialsInForm... PASS\n");
}

///////////////////////////////////////////////////////////////////////////////////////

void test_subPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = NULL;
    duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);
    duspoly_t* c = NULL;

    duspoly_t* aa;
    duspoly_t* bb;
    duspoly_t* cc;


    subPolynomialsInForm_spX (a, b, &c, Pptr);

    ZZ_pX A;
    ZZ_pX B;
    ZZ_pX C;
    ZZ_pX D;

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    sub (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (1) subPolynomialsInForm_spX... FAIL\n");
    }


    //////////////////////////////
    a = randomPolynomialInForm_spX (sz1, Pptr);
    freePolynomial_spX (&c);

    subPolynomialsInForm_spX (a, b, &c, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    sub (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (2) subPolynomialsInForm_spX... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = randomPolynomialInForm_spX (2, Pptr);
    b = randomPolynomialInForm_spX (sz2, Pptr);

    subPolynomialsInForm_spX (a, b, &c, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    sub (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (3) subPolynomialsInForm... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = NULL;
    b = randomPolynomialInForm_spX (sz2, Pptr);
    c = NULL;

    subPolynomialsInForm_spX_inp (&c, b, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    sub (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (4) subPolynomialsInForm_inp... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = randomPolynomialInForm_spX (sz1, Pptr);
    b = randomPolynomialInForm_spX (sz2, Pptr);
    c = deepCopyPolynomial_spX (a);

    subPolynomialsInForm_spX_inp (&c, b, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    sub (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (5) subPolynomialsInForm_inp... FAIL\n");
    }


    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    a = randomPolynomialInForm_spX (sz1, Pptr);
    b = randomPolynomialInForm_spX (2, Pptr);
    c = deepCopyPolynomial_spX (a);

    subPolynomialsInForm_spX_inp (&c, b, Pptr);

    cc = convertPolynomialFromMontgomery_spX (c, Pptr);
    aa = convertPolynomialFromMontgomery_spX (a, Pptr);
    bb = convertPolynomialFromMontgomery_spX (b, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc, C, Pptr);
    convertToZZ_pX (aa, A, Pptr);
    convertToZZ_pX (bb, B, Pptr);

    freePolynomial_spX (&cc);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&bb);

    sub (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (6) subPolynomialsInForm_inp... FAIL\n");
    }

    //////////////////////////////
    freePolynomial_spX (&a);
    freePolynomial_spX (&b);
    freePolynomial_spX (&c);

    duspoly_t* a7 = randomPolynomialInForm_spX (sz1, Pptr);
    duspoly_t* b7 = NULL;
    duspoly_t* c7 = deepCopyPolynomial_spX (a7);

    subPolynomialsInForm_spX_inp (&c7, b7, Pptr);

    duspoly_t* cc7 = convertPolynomialFromMontgomery_spX (c7, Pptr);
    duspoly_t* aa7 = convertPolynomialFromMontgomery_spX (a7, Pptr);
    duspoly_t* bb7 = convertPolynomialFromMontgomery_spX (b7, Pptr);

    C = 0;
    A = 0;
    B = 0;
    D = 0;

    convertToZZ_pX (cc7, C, Pptr);
    convertToZZ_pX (aa7, A, Pptr);
    convertToZZ_pX (bb7, B, Pptr);

    freePolynomial_spX (&cc7);
    freePolynomial_spX (&aa7);
    freePolynomial_spX (&bb7);

    sub (D, A, B);

    if (D != C) {
      std::cerr << "A = " << A << std::endl;
      std::cerr << "B = " << B << std::endl;

      std::cerr << "DUSP: A+B = " << C << std::endl;
      std::cerr << "NTL : A+B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (7) subPolynomialsInForm_inp... FAIL\n");
    }
    freePolynomial_spX (&a7);
    freePolynomial_spX (&b7);
    freePolynomial_spX (&c7);

  fprintf(stderr, "[DUSP] subPolynomialsInForm... PASS\n");
}

///////////////////////////////////////////////////////////////////////////////////////

void test_negPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = randomPolynomialInForm_spX (sz1, Pptr);
    duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);

    duspoly_t* aa = deepCopyPolynomial_spX (a);
    duspoly_t* bb = deepCopyPolynomial_spX (b);

    negPolynomialInForm_spX (aa, Pptr);
    negPolynomialInForm_spX (bb, Pptr);

    duspoly_t* aaa = NULL;
    duspoly_t* bbb = NULL;

    addPolynomialsInForm_spX (aa, a, &aaa, Pptr);

    if (!isZero_spX (aaa)) {
      fprintf(stderr, "a = ");
      printPolynomialOutForm_spX (a, Pptr);

      fprintf(stderr, "-a = ");
      printPolynomialOutForm_spX (aa, Pptr);

      fprintf(stderr, "a-a = ");
      printPolynomialOutForm_spX (aaa, Pptr);

      fprintf(stderr, "[DUSP] negPolynomialInForm... FAIL\n");
      exit(1);
    }

    addPolynomialsInForm_spX (b, bb, &bbb, Pptr);
    if (!isZero_spX (bbb)) {
      fprintf(stderr, "b = ");
      printPolynomialOutForm_spX (b, Pptr);

      fprintf(stderr, "-b = ");
      printPolynomialOutForm_spX (bb, Pptr);

      fprintf(stderr, "b-b = ");
      printPolynomialOutForm_spX (bbb, Pptr);

      fprintf(stderr, "[DUSP] negPolynomialInForm... FAIL\n");
      exit(1);
    }

    freePolynomial_spX (&a);
    freePolynomial_spX (&aa);
    freePolynomial_spX (&aaa);

    fprintf(stderr, "[DUSP] negPolynomialInForm... PASS\n");
}

///////////////////////////////////////////////////////////////////////////////////////

void test_monicPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a = randomPolynomialInForm_spX (sz1, Pptr);

    duspoly_t* ma = NULL;

    elem_t ca;

    monicPolynomialInForm_spX (a, &ma, &ca, Pptr);

    if (leadingCoeffInForm_spX  (ma) != smallprimefield_convert_in (1, Pptr)) {

        fprintf(stderr, "a = \n");
        printPolynomialOutForm_spX (a, Pptr);

        fprintf(stderr, "monic(a) = \n");
        printPolynomialOutForm_spX (ma, Pptr);

        fprintf(stderr, "[DUSP] monicPolynomialInForm... FAIL\n");

        exit (1);
    }

    duspoly_t* aa = NULL;

    scalarMulPolynomialInForm_spX (ma, ca, &aa, Pptr);

    if (!isEqual_spX (a, aa)) {
        fprintf(stderr, "a = \n");
        printPolynomialOutForm_spX (a, Pptr);

        fprintf(stderr, "monic(a) = \n");
        printPolynomialOutForm_spX (ma, Pptr);

        fprintf(stderr, "cont(a).monic(a) = \n");
        printPolynomialOutForm_spX (aa, Pptr);

        fprintf(stderr, "[DUSP] scalarMulPolynomialInForm... FAIL\n");
        exit (1);
    }

    fprintf(stderr, "[DUSP] scalarMulPolynomialInForm... PASS\n");
    fprintf(stderr, "[DUSP] monicPolynomialInForm... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

void test_mulPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag) {

  duspoly_t* a1 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b1 = randomPolynomialInForm_spX (sz2, Pptr);
  duspoly_t* c1 = NULL;

  plainMulPolynomialsInForm_spX (a1, b1, &c1, Pptr);

  duspoly_t* cc1 = convertPolynomialFromMontgomery_spX (c1, Pptr);
  duspoly_t* aa1 = convertPolynomialFromMontgomery_spX (a1, Pptr);
  duspoly_t* bb1 = convertPolynomialFromMontgomery_spX (b1, Pptr);

  ZZ_pX C;
  ZZ_pX A;
  ZZ_pX B;
  ZZ_pX D;

  convertToZZ_pX (cc1, C, Pptr);
  convertToZZ_pX (aa1, A, Pptr);
  convertToZZ_pX (bb1, B, Pptr);

  freePolynomial_spX (&cc1);
  freePolynomial_spX (&aa1);
  freePolynomial_spX (&bb1);

  mul (D, A, B);

  if (D != C) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;

    std::cerr << "DUSP: A*B = " << C << std::endl;
    std::cerr << "NTL : A*B = " << D << std::endl;

    fprintf(stderr, "[DUSP] (1) plainMulPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a1);
  freePolynomial_spX (&b1);
  freePolynomial_spX (&c1);

  ////////////////////////////////////////
  duspoly_t* a2 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b2 = NULL;
  duspoly_t* c2 = NULL;

  plainMulPolynomialsInForm_spX (a2, b2, &c2, Pptr);

  duspoly_t* cc2 = convertPolynomialFromMontgomery_spX (c2, Pptr);
  duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
  duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);

  C = 0;
  A = 0;
  B = 0;
  D = 0;

  convertToZZ_pX (cc2, C, Pptr);
  convertToZZ_pX (aa2, A, Pptr);
  convertToZZ_pX (bb2, B, Pptr);

  freePolynomial_spX (&cc2);
  freePolynomial_spX (&aa2);
  freePolynomial_spX (&bb2);

  mul (D, A, B);

  if (D != C) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;

    std::cerr << "DUSP: A*B = " << C << std::endl;
    std::cerr << "NTL : A*B = " << D << std::endl;

    fprintf(stderr, "[DUSP] (2) plainMulPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a2);
  freePolynomial_spX (&b2);
  freePolynomial_spX (&c2);

  ////////////////////////////////////////
  duspoly_t* a3 = NULL;
  duspoly_t* b3 = randomPolynomialInForm_spX (sz2, Pptr);
  duspoly_t* c3 = NULL;

  plainMulPolynomialsInForm_spX (a3, b3, &c3, Pptr);

  duspoly_t* cc3 = convertPolynomialFromMontgomery_spX (c3, Pptr);
  duspoly_t* aa3 = convertPolynomialFromMontgomery_spX (a3, Pptr);
  duspoly_t* bb3 = convertPolynomialFromMontgomery_spX (b3, Pptr);

  C = 0;
  A = 0;
  B = 0;
  D = 0;

  convertToZZ_pX (cc3, C, Pptr);
  convertToZZ_pX (aa3, A, Pptr);
  convertToZZ_pX (bb3, B, Pptr);

  freePolynomial_spX (&cc3);
  freePolynomial_spX (&aa3);
  freePolynomial_spX (&bb3);

  mul (D, A, B);

  if (D != C) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;

    std::cerr << "DUSP: A*B = " << C << std::endl;
    std::cerr << "NTL : A*B = " << D << std::endl;

    fprintf(stderr, "[DUSP] (3) plainMulPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a3);
  freePolynomial_spX (&b3);
  freePolynomial_spX (&c3);

  ////////////////////////////////////////
  duspoly_t* a4 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b4 = randomPolynomialInForm_spX (sz2, Pptr);
  duspoly_t* c4 = NULL;
  duspoly_t* d4 = NULL;

  plainSqrPolynomialInForm_spX  (a4, &c4, Pptr);
  plainSqrPolynomialInForm_spX  (b4, &d4, Pptr);

  duspoly_t* cc4 = convertPolynomialFromMontgomery_spX (c4, Pptr);
  duspoly_t* dd4 = convertPolynomialFromMontgomery_spX (d4, Pptr);
  duspoly_t* aa4 = convertPolynomialFromMontgomery_spX (a4, Pptr);
  duspoly_t* bb4 = convertPolynomialFromMontgomery_spX (b4, Pptr);

  ZZ_pX C0;
  ZZ_pX D0;
  ZZ_pX C1;
  ZZ_pX D1;
  A = 0;
  B = 0;

  convertToZZ_pX (cc4, C0, Pptr);
  convertToZZ_pX (dd4, D0, Pptr);
  convertToZZ_pX (aa4, A, Pptr);
  convertToZZ_pX (bb4, B, Pptr);

  freePolynomial_spX (&cc4);
  freePolynomial_spX (&dd4);
  freePolynomial_spX (&aa4);
  freePolynomial_spX (&bb4);

  sqr (C1, A);
  sqr (D1, B);

  if (C0 != C1) {
    std::cerr << "A = " << A << std::endl;

    std::cerr << "DUSP: A^2 = " << C0 << std::endl;
    std::cerr << "NTL : A^2 = " << C1 << std::endl;

    fprintf(stderr, "[DUSP] plainSqrPolynomialInForm... FAIL\n");
    exit (1);
  }

  if (D0 != D1) {
    std::cerr << "B = " << B << std::endl;

    std::cerr << "DUSP: B^2 = " << D0 << std::endl;
    std::cerr << "NTL : B^2 = " << D1 << std::endl;

    fprintf(stderr, "[DUSP] plainSqrPolynomialInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a4);
  freePolynomial_spX (&b4);
  freePolynomial_spX (&c4);
  freePolynomial_spX (&d4);

  fprintf(stderr, "[DUSP] plainSqrPolynomialInForm... PASS\n");

  ////////////////////////////////////////5
  duspoly_t* a5 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b5 = randomPolynomialInForm_spX (sz2, Pptr);
  duspoly_t* c5 = NULL;

  KaratsubaMulPolynomialsInForm_spX (a5, b5, &c5, Pptr);

  duspoly_t* cc5 = convertPolynomialFromMontgomery_spX (c5, Pptr);
  duspoly_t* aa5 = convertPolynomialFromMontgomery_spX (a5, Pptr);
  duspoly_t* bb5 = convertPolynomialFromMontgomery_spX (b5, Pptr);

  C = 0;
  A = 0;
  B = 0;
  D = 0;

  convertToZZ_pX (cc5, C, Pptr);
  convertToZZ_pX (aa5, A, Pptr);
  convertToZZ_pX (bb5, B, Pptr);

  freePolynomial_spX (&cc5);
  freePolynomial_spX (&aa5);
  freePolynomial_spX (&bb5);

  mul (D, A, B);

  if (D != C) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;

    std::cerr << "DUSP: A*B = " << C << std::endl;
    std::cerr << "NTL : A*B = " << D << std::endl;

    fprintf(stderr, "[DUSP] (5) KaratsubaMulPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a5);
  freePolynomial_spX (&b5);
  freePolynomial_spX (&c5);

  fprintf(stderr, "[DUSP] KaratsubaMulPolynomialsInForm... PASS\n");

   ////////////////////////////////////////6
  if (Pptr->prime == 4179340454199820289) {
      duspoly_t* a6 = randomPolynomialInForm_spX (sz1, Pptr);
      duspoly_t* b6 = randomPolynomialInForm_spX (sz2, Pptr);
      duspoly_t* c6 = NULL;

      fastMulPolynomialInForm_FFT_spX (a6, b6, &c6, Pptr);

      // std::cerr << "c6:" << std::endl;
      // printPolynomial_spX (c6, Pptr);

      duspoly_t* cc6 = convertPolynomialFromMontgomery_spX (c6, Pptr);
      duspoly_t* aa6 = convertPolynomialFromMontgomery_spX (a6, Pptr);
      duspoly_t* bb6 = convertPolynomialFromMontgomery_spX (b6, Pptr);

      C = 0;
      A = 0;
      B = 0;
      D = 0;

      convertToZZ_pX (cc6, C, Pptr);
      convertToZZ_pX (aa6, A, Pptr);
      convertToZZ_pX (bb6, B, Pptr);

      freePolynomial_spX (&cc6);
      freePolynomial_spX (&aa6);
      freePolynomial_spX (&bb6);

      mul (D, A, B);

      if (D != C) {
        std::cerr << "A = " << A << std::endl;
        std::cerr << "B = " << B << std::endl;

        std::cerr << "DUSP: A*B = " << C << std::endl;
        std::cerr << "NTL : A*B = " << D << std::endl;

        fprintf(stderr, "[DUSP] (6) fastMulPolynomialInForm_FFT... FAIL\n");
        exit (1);
      }

      freePolynomial_spX (&a6);
      freePolynomial_spX (&b6);
      freePolynomial_spX (&c6);

      fprintf(stderr, "[DUSP] fastMulPolynomialInForm_FFT... PASS\n");
  }

  ////////////////////////////////////////7 (Bug Report in KarMul! (24 Sep, 2019))

  elem_t elem_a7[19] = {1, 0, 0, 0, 2, 0, 2, 0, 1, 0, 2, 1, 2, 2, 0, 2, 2, 0, 2};
  elem_t elem_b7[1] = {2};

  duspoly_t* a7 = NULL;
  duspoly_t* b7 = NULL;

  setCoefsInForm_spX (&a7, elem_a7, 19);
  setCoefsInForm_spX (&b7, elem_b7, 1);

  // printPolynomial_spX (a7, Pptr);
  // printPolynomial_spX (b7, Pptr);

  duspoly_t* c7 = NULL;

  mulPolynomialsInForm_spX (a7, b7, &c7, Pptr);

  duspoly_t* cc7 = convertPolynomialFromMontgomery_spX (c7, Pptr);
  duspoly_t* aa7 = convertPolynomialFromMontgomery_spX (a7, Pptr);
  duspoly_t* bb7 = convertPolynomialFromMontgomery_spX (b7, Pptr);

  C = 0;
  A = 0;
  B = 0;
  D = 0;

  convertToZZ_pX (cc7, C, Pptr);
  convertToZZ_pX (aa7, A, Pptr);
  convertToZZ_pX (bb7, B, Pptr);

  freePolynomial_spX (&cc7);
  freePolynomial_spX (&aa7);
  freePolynomial_spX (&bb7);

  mul (D, A, B);

  if (D != C) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;

    std::cerr << "DUSP: A*B = " << C << std::endl;
    std::cerr << "NTL : A*B = " << D << std::endl;

    fprintf(stderr, "[DUSP] (7) mulPolynomialsInForm_spX... FAIL\n");
    exit (1);
  }


  freePolynomial_spX (&a7);
  freePolynomial_spX (&b7);
  freePolynomial_spX (&c7);

////////////////////////////////////////8 (Bug Report in KarMul! (24 Sep, 2019))

  elem_t elem_a8[3] = {0, 1, 2};
  elem_t elem_b8[2] = {1,2};

  duspoly_t* a8 = NULL;
  duspoly_t* b8 = NULL;

  setCoefsInForm_spX (&a8, elem_a8, 3);
  setCoefsInForm_spX (&b8, elem_b8, 2);

  // printPolynomial_spX (a8, Pptr);
  // printPolynomial_spX (b8, Pptr);

  duspoly_t* c8 = NULL;

  mulPolynomialsInForm_spX (a8, b8, &c8, Pptr);

  duspoly_t* cc8 = convertPolynomialFromMontgomery_spX (c8, Pptr);
  duspoly_t* aa8 = convertPolynomialFromMontgomery_spX (a8, Pptr);
  duspoly_t* bb8 = convertPolynomialFromMontgomery_spX (b8, Pptr);

  C = 0;
  A = 0;
  B = 0;
  D = 0;

  convertToZZ_pX (cc8, C, Pptr);
  convertToZZ_pX (aa8, A, Pptr);
  convertToZZ_pX (bb8, B, Pptr);

  freePolynomial_spX (&cc8);
  freePolynomial_spX (&aa8);
  freePolynomial_spX (&bb8);

  mul (D, A, B);

  if (D != C) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;

    std::cerr << "DUSP: A*B = " << C << std::endl;
    std::cerr << "NTL : A*B = " << D << std::endl;

    fprintf(stderr, "[DUSP] (8) mulPolynomialsInForm_spX... FAIL\n");
    exit (1);
  }


  freePolynomial_spX (&a8);
  freePolynomial_spX (&b8);
  freePolynomial_spX (&c8);

  fprintf(stderr, "[DUSP] mulPolynomialsInForm_spX... PASS\n");


}

///////////////////////////////////////////////////////////////////////////////////////

void test_divPoly (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  duspoly_t* a1 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b1 = deepCopyPolynomial_spX (a1);

  duspoly_t* q1 = NULL;
  duspoly_t* r1 = NULL;

  plainDivPolynomialsInForm_spX (a1, b1, &r1, &q1, Pptr);

  if (r1 != NULL || !isOneInForm_spX (q1, Pptr)) {
    fprintf (stderr, "a1 = ");
    printPolynomial_spX (a1, Pptr);

    fprintf(stderr, "b1 = ");
    printPolynomial_spX (b1, Pptr);

    fprintf(stderr, "r1 = ");
    printPolynomial_spX (r1, Pptr);

    fprintf(stderr, "q1 = ");
    printPolynomial_spX (q1, Pptr);

    fprintf(stderr, "[DUSP] (1) plainDivPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a1);
  freePolynomial_spX (&b1);
  freePolynomial_spX (&q1);
  freePolynomial_spX (&r1);

  ////////////////////////2
  duspoly_t* a2 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b2 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* q2 = NULL;
  duspoly_t* r2 = NULL;

  plainDivPolynomialsInForm_spX (a2, b2, &r2, &q2, Pptr);

  duspoly_t* qq2 = convertPolynomialFromMontgomery_spX (q2, Pptr);
  duspoly_t* rr2 = convertPolynomialFromMontgomery_spX (r2, Pptr);
  duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
  duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);

  ZZ_pX A;
  ZZ_pX B;
  ZZ_pX R;
  ZZ_pX Q;

  convertToZZ_pX (qq2, Q, Pptr);
  convertToZZ_pX (rr2, R, Pptr);
  convertToZZ_pX (aa2, A, Pptr);
  convertToZZ_pX (bb2, B, Pptr);

  freePolynomial_spX (&qq2);
  freePolynomial_spX (&rr2);
  freePolynomial_spX (&aa2);
  freePolynomial_spX (&bb2);

  ZZ_pX RR;
  ZZ_pX QQ;

  DivRem (QQ, RR, A, B);

  if (QQ != Q || RR != R) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: R = " << R << std::endl;
    std::cerr << "NTL:  R = " << RR << std::endl;
    std::cerr << "DUSP: Q = " << Q << std::endl;
    std::cerr << "NTL : Q = " << QQ << std::endl;


    fprintf(stderr, "[DUSP] (2) plainDivPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a2);
  freePolynomial_spX (&b2);
  freePolynomial_spX (&q2);
  freePolynomial_spX (&r2);

  fprintf(stderr, "[DUSP] plainDivPolynomialsInForm... PASS\n");

  ////////////////////////3
  duspoly_t* a3 = randomPolynomialInForm_spX (sz2, Pptr);
  duspoly_t* b3 = randomPolynomialInForm_spX (sz1, Pptr);

  duspoly_t* r3 = NULL;

  plainRemPolynomialsInForm_spX (a3, b3, &r3, Pptr);

  duspoly_t* rr3 = convertPolynomialFromMontgomery_spX (r3, Pptr);
  duspoly_t* aa3 = convertPolynomialFromMontgomery_spX (a3, Pptr);
  duspoly_t* bb3 = convertPolynomialFromMontgomery_spX (b3, Pptr);

  A = 0;
  B = 0;
  R = 0;

  convertToZZ_pX (rr3, R, Pptr);
  convertToZZ_pX (aa3, A, Pptr);
  convertToZZ_pX (bb3, B, Pptr);

  freePolynomial_spX (&rr3);
  freePolynomial_spX (&aa3);
  freePolynomial_spX (&bb3);

  RR = 0;

  rem (RR, A, B);

  if (RR != R) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: R = " << R << std::endl;
    std::cerr << "NTL:  R = " << RR << std::endl;

    fprintf(stderr, "[DUSP] (3) plainRemPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a3);
  freePolynomial_spX (&b3);
  freePolynomial_spX (&r3);

 ////////////////////////4
  duspoly_t* a4 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b4 = NULL;
  exponentiatePolynomialInForm_spX (a4, 4, &b4, Pptr);

  duspoly_t* r4 = NULL;

  plainRemPolynomialsInForm_spX (a4, b4, &r4, Pptr);

  duspoly_t* rr4 = convertPolynomialFromMontgomery_spX (r4, Pptr);
  duspoly_t* aa4 = convertPolynomialFromMontgomery_spX (a4, Pptr);
  duspoly_t* bb4 = convertPolynomialFromMontgomery_spX (b4, Pptr);

  A = 0;
  B = 0;
  R = 0;

  convertToZZ_pX (rr4, R, Pptr);
  convertToZZ_pX (aa4, A, Pptr);
  convertToZZ_pX (bb4, B, Pptr);

  freePolynomial_spX (&rr4);
  freePolynomial_spX (&aa4);
  freePolynomial_spX (&bb4);

  RR = 0;

  rem (RR, A, B);

  if (RR != R) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: R = " << R << std::endl;
    std::cerr << "NTL:  R = " << RR << std::endl;

    fprintf(stderr, "[DUSP] (4) plainRemPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a4);
  freePolynomial_spX (&b4);
  freePolynomial_spX (&r4);

 ////////////////////////5
  duspoly_t* a5 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b5 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* r5 = deepCopyPolynomial_spX (a5);

  plainRemPolynomialsInForm_spX_inp (&r5, b5, Pptr);

  duspoly_t* rr5 = convertPolynomialFromMontgomery_spX (r5, Pptr);
  duspoly_t* aa5 = convertPolynomialFromMontgomery_spX (a5, Pptr);
  duspoly_t* bb5 = convertPolynomialFromMontgomery_spX (b5, Pptr);

  A = 0;
  B = 0;
  R = 0;

  convertToZZ_pX (rr5, R, Pptr);
  convertToZZ_pX (aa5, A, Pptr);
  convertToZZ_pX (bb5, B, Pptr);

  freePolynomial_spX (&rr5);
  freePolynomial_spX (&aa5);
  freePolynomial_spX (&bb5);

  RR = 0;

  rem (RR, A, B);

  if (RR != R) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: R = " << R << std::endl;
    std::cerr << "NTL:  R = " << RR << std::endl;

    fprintf(stderr, "[DUSP] (5) plainRemPolynomialsInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a5);
  freePolynomial_spX (&b5);
  freePolynomial_spX (&r5);

  fprintf(stderr, "[DUSP] plainRemPolynomialsInForm... PASS\n");

  /////////////////////////7

  duspoly_t* a7 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b7 = deepCopyPolynomial_spX (a7);

  duspoly_t* q7 = NULL;
  duspoly_t* r7 = NULL;

  fastDivPolynomialInForm_wPSInv_spX (a7, b7, &r7, &q7, Pptr);

  if (r7 != NULL || !isOneInForm_spX (q7, Pptr)) {
    fprintf (stderr, "a7 = ");
    printPolynomial_spX (a7, Pptr);

    fprintf(stderr, "b7 = ");
    printPolynomial_spX (b7, Pptr);

    fprintf(stderr, "r7 = ");
    printPolynomial_spX (r7, Pptr);

    fprintf(stderr, "q7 = ");
    printPolynomial_spX (q7, Pptr);

    fprintf(stderr, "[DUSP] (7) fastDivPolynomialInForm_wPSInv... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a7);
  freePolynomial_spX (&b7);
  freePolynomial_spX (&q7);
  freePolynomial_spX (&r7);

  ////////////////////////8
  duspoly_t* a8 = randomPolynomialInForm_spX (sz1*2, Pptr);
  duspoly_t* b8 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* q8 = NULL;
  duspoly_t* r8 = NULL;

  fastDivPolynomialInForm_wPSInv_spX (a8, b8, &r8, &q8, Pptr);

  duspoly_t* qq8 = convertPolynomialFromMontgomery_spX (q8, Pptr);
  duspoly_t* rr8 = convertPolynomialFromMontgomery_spX (r8, Pptr);
  duspoly_t* aa8 = convertPolynomialFromMontgomery_spX (a8, Pptr);
  duspoly_t* bb8 = convertPolynomialFromMontgomery_spX (b8, Pptr);

   A = 0;
   B = 0;
   R = 0;
   Q = 0;

  convertToZZ_pX (qq8, Q, Pptr);
  convertToZZ_pX (rr8, R, Pptr);
  convertToZZ_pX (aa8, A, Pptr);
  convertToZZ_pX (bb8, B, Pptr);

  freePolynomial_spX (&qq8);
  freePolynomial_spX (&rr8);
  freePolynomial_spX (&aa8);
  freePolynomial_spX (&bb8);

   RR = 0;
   QQ = 0;

  DivRem (QQ, RR, A, B);

  if (QQ != Q || RR != R) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: R = " << R << std::endl;
    std::cerr << "NTL:  R = " << RR << std::endl;
    std::cerr << "DUSP: Q = " << Q << std::endl;
    std::cerr << "NTL : Q = " << QQ << std::endl;


    fprintf(stderr, "[DUSP] (8) fastDivPolynomialInForm_wPSInv... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a8);
  freePolynomial_spX (&b8);
  freePolynomial_spX (&q8);
  freePolynomial_spX (&r8);

  fprintf(stderr, "[DUSP] fastDivPolynomialInForm... PASS\n");

   ////////////////////////6

  duspoly_t* a6 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b6 = NULL;
  exponentiatePolynomialInForm_spX (a6, 4, &b6, Pptr);

  if (!isDividablePolys_spX (b6, a6, Pptr)) {
    fprintf(stderr, "a6 = ");
    printPolynomial_spX (a6, Pptr);

    fprintf(stderr, "b6 = ");
    printPolynomial_spX (b6, Pptr);

    fprintf(stderr, "[DUSP] (6) isDivisiblePolys... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a6);
  freePolynomial_spX (&b6);

  fprintf(stderr, "[DUSP] isDivisiblePolys... PASS\n");

  ///////////////////// 9

  if (Pptr->prime == 4179340454199820289) {

  duspoly_t* a9 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b9 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* r9 = NULL;
  duspoly_t* q9 = NULL;

  fastDivPolynomialInForm_wPSInv_FFT_spX (a9, b9, &r9, &q9, Pptr);

  duspoly_t* qq9 = convertPolynomialFromMontgomery_spX (q9, Pptr);
  duspoly_t* rr9 = convertPolynomialFromMontgomery_spX (r9, Pptr);
  duspoly_t* aa9 = convertPolynomialFromMontgomery_spX (a9, Pptr);
  duspoly_t* bb9 = convertPolynomialFromMontgomery_spX (b9, Pptr);

   A = 0;
   B = 0;
   R = 0;
   Q = 0;

  convertToZZ_pX (qq9, Q, Pptr);
  convertToZZ_pX (rr9, R, Pptr);
  convertToZZ_pX (aa9, A, Pptr);
  convertToZZ_pX (bb9, B, Pptr);

  freePolynomial_spX (&qq9);
  freePolynomial_spX (&rr9);
  freePolynomial_spX (&aa9);
  freePolynomial_spX (&bb9);

   RR = 0;
   QQ = 0;

  DivRem (QQ, RR, A, B);

  if (QQ != Q || RR != R) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: R = " << R << std::endl;
    std::cerr << "NTL:  R = " << RR << std::endl;
    std::cerr << "DUSP: Q = " << Q << std::endl;
    std::cerr << "NTL : Q = " << QQ << std::endl;


    fprintf(stderr, "[DUSP] fastDivPolynomialInForm_wPSInv_FFT... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a9);
  freePolynomial_spX (&b9);
  freePolynomial_spX (&q9);
  freePolynomial_spX (&r9);

  fprintf(stderr, "[DUSP] fastDivPolynomialInForm_FFT... PASS\n");

}

}

///////////////////////////////////////////////////////////////////////////////////////

void test_GCD (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  duspoly_t* a1 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b1 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* g1 = NULL;

  plainGCDInForm_spX (a1, b1, &g1, Pptr);

  duspoly_t* gg1 = convertPolynomialFromMontgomery_spX (g1, Pptr);
  duspoly_t* aa1 = convertPolynomialFromMontgomery_spX (a1, Pptr);
  duspoly_t* bb1 = convertPolynomialFromMontgomery_spX (b1, Pptr);

  ZZ_pX A;
  ZZ_pX B;
  ZZ_pX G;

  convertToZZ_pX (gg1, G, Pptr);
  convertToZZ_pX (aa1, A, Pptr);
  convertToZZ_pX (bb1, B, Pptr);

  freePolynomial_spX (&gg1);
  freePolynomial_spX (&aa1);
  freePolynomial_spX (&bb1);

  ZZ_pX GG;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;


    fprintf(stderr, "[DUSP] (1) plainGCDInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a1);
  freePolynomial_spX (&b1);
  freePolynomial_spX (&g1);

  /////////////////////////////////2
  duspoly_t* a2 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b2 = NULL;
  exponentiatePolynomialInForm_spX (a2, 5, &b2, Pptr);

  duspoly_t* g2 = NULL;

  plainGCDInForm_spX (a2, b2, &g2, Pptr);

  duspoly_t* gg2 = convertPolynomialFromMontgomery_spX (g2, Pptr);
  duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
  duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);

  A = 0;
  B = 0;
  G = 0;

  convertToZZ_pX (gg2, G, Pptr);
  convertToZZ_pX (aa2, A, Pptr);
  convertToZZ_pX (bb2, B, Pptr);

  freePolynomial_spX (&gg2);
  freePolynomial_spX (&aa2);
  freePolynomial_spX (&bb2);

  GG = 0;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;


    fprintf(stderr, "[DUSP] (2) plainGCDInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a2);
  freePolynomial_spX (&b2);
  freePolynomial_spX (&g2);


  /////////////////////////////////3
  duspoly_t* a3 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b3 = NULL;

  duspoly_t* g3 = NULL;

  plainGCDInForm_spX (a3, b3, &g3, Pptr);

  duspoly_t* gg3 = convertPolynomialFromMontgomery_spX (g3, Pptr);
  duspoly_t* aa3 = convertPolynomialFromMontgomery_spX (a3, Pptr);
  duspoly_t* bb3 = convertPolynomialFromMontgomery_spX (b3, Pptr);

  A = 0;
  B = 0;
  G = 0;

  convertToZZ_pX (gg3, G, Pptr);
  convertToZZ_pX (aa3, A, Pptr);
  convertToZZ_pX (bb3, B, Pptr);

  freePolynomial_spX (&gg3);
  freePolynomial_spX (&aa3);
  freePolynomial_spX (&bb3);

  GG = 0;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;

    fprintf(stderr, "[DUSP] (3) plainGCDInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a3);
  freePolynomial_spX (&b3);
  freePolynomial_spX (&g3);


  /////////////////////////////////4
  duspoly_t* a4 = randomPolynomialInForm_spX (2, Pptr);
  duspoly_t* b4 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* g4 = NULL;

  plainGCDInForm_spX (a4, b4, &g4, Pptr);

  duspoly_t* gg4 = convertPolynomialFromMontgomery_spX (g4, Pptr);
  duspoly_t* aa4 = convertPolynomialFromMontgomery_spX (a4, Pptr);
  duspoly_t* bb4 = convertPolynomialFromMontgomery_spX (b4, Pptr);

  A = 0;
  B = 0;
  G = 0;

  convertToZZ_pX (gg4, G, Pptr);
  convertToZZ_pX (aa4, A, Pptr);
  convertToZZ_pX (bb4, B, Pptr);

  freePolynomial_spX (&gg4);
  freePolynomial_spX (&aa4);
  freePolynomial_spX (&bb4);

  GG = 0;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;

    fprintf(stderr, "[DUSP] (4) plainGCDInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a4);
  freePolynomial_spX (&b4);
  freePolynomial_spX (&g4);

  fprintf(stderr, "[DUSP] plainGCDInForm... PASS\n");


  ///////////////// 5
  duspoly_t* a5 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b5 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* g5 = NULL;

  GCDInForm_spX (a5, b5, &g5, Pptr);

  duspoly_t* gg5 = convertPolynomialFromMontgomery_spX (g5, Pptr);
  duspoly_t* aa5 = convertPolynomialFromMontgomery_spX (a5, Pptr);
  duspoly_t* bb5 = convertPolynomialFromMontgomery_spX (b5, Pptr);

   A = 0;
   B = 0;
   G = 0;

  convertToZZ_pX (gg5, G, Pptr);
  convertToZZ_pX (aa5, A, Pptr);
  convertToZZ_pX (bb5, B, Pptr);

  freePolynomial_spX (&gg5);
  freePolynomial_spX (&aa5);
  freePolynomial_spX (&bb5);

  GG = 0;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;

    fprintf(stderr, "[DUSP] (5) halfGCDMatrixInForm... FAIL\n");
    fprintf(stderr, "[DUSP] (5) GCDInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a5);
  freePolynomial_spX (&b5);
  freePolynomial_spX (&g5);

  /////////////////////////////////6
  duspoly_t* a6 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b6 = NULL;
  exponentiatePolynomialInForm_spX (a6, 5, &b6, Pptr);

  duspoly_t* g6 = NULL;

  GCDInForm_spX (a6, b6, &g6, Pptr);

  duspoly_t* gg6 = convertPolynomialFromMontgomery_spX (g6, Pptr);
  duspoly_t* aa6 = convertPolynomialFromMontgomery_spX (a6, Pptr);
  duspoly_t* bb6 = convertPolynomialFromMontgomery_spX (b6, Pptr);

  A = 0;
  B = 0;
  G = 0;

  convertToZZ_pX (gg6, G, Pptr);
  convertToZZ_pX (aa6, A, Pptr);
  convertToZZ_pX (bb6, B, Pptr);

  freePolynomial_spX (&gg6);
  freePolynomial_spX (&aa6);
  freePolynomial_spX (&bb6);

  GG = 0;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;


    fprintf(stderr, "[DUSP] (6) halfGCDMatrixInForm... FAIL\n");
    fprintf(stderr, "[DUSP] (6) GCDInForm... FAIL\n");
     exit (1);
  }

  freePolynomial_spX (&a6);
  freePolynomial_spX (&b6);
  freePolynomial_spX (&g6);


  /////////////////////////////////7
  duspoly_t* a7 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b7 = NULL;

  duspoly_t* g7 = NULL;

  GCDInForm_spX (a7, b7, &g7, Pptr);

  duspoly_t* gg7 = convertPolynomialFromMontgomery_spX (g7, Pptr);
  duspoly_t* aa7 = convertPolynomialFromMontgomery_spX (a7, Pptr);
  duspoly_t* bb7 = convertPolynomialFromMontgomery_spX (b7, Pptr);

  A = 0;
  B = 0;
  G = 0;

  convertToZZ_pX (gg7, G, Pptr);
  convertToZZ_pX (aa7, A, Pptr);
  convertToZZ_pX (bb7, B, Pptr);

  freePolynomial_spX (&gg7);
  freePolynomial_spX (&aa7);
  freePolynomial_spX (&bb7);

  GG = 0;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;

    fprintf(stderr, "[DUSP] (7) halfGCDMatrixInForm... FAIL\n");
    fprintf(stderr, "[DUSP] (7) GCDInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a7);
  freePolynomial_spX (&b7);
  freePolynomial_spX (&g7);


  /////////////////////////////////8
  duspoly_t* a8 = randomPolynomialInForm_spX (2, Pptr);
  duspoly_t* b8 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* g8 = NULL;

  GCDInForm_spX (a8, b8, &g8, Pptr);

  duspoly_t* gg8 = convertPolynomialFromMontgomery_spX (g8, Pptr);
  duspoly_t* aa8 = convertPolynomialFromMontgomery_spX (a8, Pptr);
  duspoly_t* bb8 = convertPolynomialFromMontgomery_spX (b8, Pptr);

  A = 0;
  B = 0;
  G = 0;

  convertToZZ_pX (gg8, G, Pptr);
  convertToZZ_pX (aa8, A, Pptr);
  convertToZZ_pX (bb8, B, Pptr);

  freePolynomial_spX (&gg8);
  freePolynomial_spX (&aa8);
  freePolynomial_spX (&bb8);

  GG = 0;

  GCD (GG, A, B);

  if (G != GG) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
    std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;

    fprintf(stderr, "[DUSP] (8) halfGCDMatrixInForm... FAIL\n");
    fprintf(stderr, "[DUSP] (8) GCDInForm... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a8);
  freePolynomial_spX (&b8);
  freePolynomial_spX (&g8);

  fprintf(stderr, "[DUSP] halfGCDMatrixInForm... PASS\n");
  fprintf(stderr, "[DUSP] GCDInForm... PASS\n");

  /////////////////////// 9
  if (Pptr->prime == 4179340454199820289) {

      duspoly_t* a9 = randomPolynomialInForm_spX (sz1, Pptr);
      duspoly_t* b9 = randomPolynomialInForm_spX (sz2, Pptr);

      duspoly_t* g9 = NULL;

      GCDInForm_FFT_spX (a9, b9, &g9, Pptr);

      duspoly_t* gg9 = convertPolynomialFromMontgomery_spX (g9, Pptr);
      duspoly_t* aa9 = convertPolynomialFromMontgomery_spX (a9, Pptr);
      duspoly_t* bb9 = convertPolynomialFromMontgomery_spX (b9, Pptr);

      ZZ_pX A;
      ZZ_pX B;
      ZZ_pX G;

      convertToZZ_pX (gg9, G, Pptr);
      convertToZZ_pX (aa9, A, Pptr);
      convertToZZ_pX (bb9, B, Pptr);

      freePolynomial_spX (&gg9);
      freePolynomial_spX (&aa9);
      freePolynomial_spX (&bb9);

      ZZ_pX GG;

      GCD (GG, A, B);

      if (G != GG) {
        std::cerr << "A = " << A << std::endl;
        std::cerr << "B = " << B << std::endl;
        std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
        std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;


        fprintf(stderr, "[DUSP] (9) GCDInForm_FFT... FAIL\n");
        exit (9);
      }

      freePolynomial_spX (&a9);
      freePolynomial_spX (&b9);
      freePolynomial_spX (&g9);

      fprintf(stderr, "[DUSP] halfGCDMatrixInForm_FFT... PASS\n");
      fprintf(stderr, "[DUSP] GCDInForm_FFT... PASS\n");

  }

}

///////////////////////////////////////////////////////////////////////////////////////

void test_extGCD (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    elem_t a1_elems[4] = {-4,8,11,3};
    elem_t b1_elems[5] = {-54, 243,-315,25,125};

    duspoly_t* a1 = NULL;
    duspoly_t* b1 = NULL;
    duspoly_t* u1 = NULL;
    duspoly_t* v1 = NULL;
    duspoly_t* g1 = NULL;

    setCoefs_spX (&a1, a1_elems, 4, Pptr);
    setCoefs_spX (&b1, b1_elems, 5, Pptr);

    // plainDivPolynomialsInForm_spX (b1, a1, &r1, &q1, Pptr);
    plainExtGCDInForm_spX (a1, b1, &u1, &v1, &g1, Pptr);

    ZZ_pX A1;
    ZZ_pX B1;
    ZZ_pX U1, V1, G1;
    ZZ_pX UU1, VV1, GG1;

    duspoly_t* aa1 = convertPolynomialFromMontgomery_spX (a1, Pptr);
    duspoly_t* bb1 = convertPolynomialFromMontgomery_spX (b1, Pptr);
    duspoly_t* uu1 = convertPolynomialFromMontgomery_spX (u1, Pptr);
    duspoly_t* vv1 = convertPolynomialFromMontgomery_spX (v1, Pptr);
    duspoly_t* gg1 = convertPolynomialFromMontgomery_spX (g1, Pptr);

    convertToZZ_pX (aa1, A1, Pptr);
    convertToZZ_pX (bb1, B1, Pptr);
    convertToZZ_pX (uu1, UU1, Pptr);
    convertToZZ_pX (vv1, VV1, Pptr);
    convertToZZ_pX (gg1, GG1, Pptr);

    freePolynomial_spX (&aa1);
    freePolynomial_spX (&bb1);
    freePolynomial_spX (&uu1);
    freePolynomial_spX (&vv1);
    freePolynomial_spX (&gg1);

    PlainXGCD (G1, U1, V1, A1, B1);

    if (G1 != GG1 || U1 != UU1 || V1 != VV1) {

        std::cerr << "A1 = " << A1 << std::endl;
        std::cerr << "B1 = " << B1 << std::endl;
        std::cerr << "U1 = " << U1 << std::endl;
        std::cerr << "UU1 = " << UU1 << std::endl;
        std::cerr << "V1 = " << V1 << std::endl;
        std::cerr << "VV1 = " << VV1 << std::endl;
        std::cerr << "G1 = " << G1 << std::endl;
        std::cerr << "GG1 = " << GG1 << std::endl;

        fprintf (stderr, "[DUSP] plainExtGCD... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&a1);
    freePolynomial_spX (&b1);
    freePolynomial_spX (&g1);
    freePolynomial_spX (&u1);
    freePolynomial_spX (&v1);

    ////////////////////////////// 2
    elem_t a2_elems[10] = {1,-3,3,2,-6,3,3,-3,0,1};
    elem_t b2_elems[4] = {1,-1,0,1};

    duspoly_t* a2 = NULL;
    duspoly_t* b2 = NULL;
    duspoly_t* u2 = NULL;
    duspoly_t* v2 = NULL;
    duspoly_t* g2 = NULL;

    setCoefs_spX (&a2, a2_elems, 10, Pptr);
    setCoefs_spX (&b2, b2_elems, 4, Pptr);

    // plainDivPolynomialsInForm_spX (b1, a1, &r1, &q1, Pptr);
    plainExtGCDInForm_spX (a2, b2, &u2, &v2, &g2, Pptr);

    ZZ_pX A2;
    ZZ_pX B2;
    ZZ_pX U2, V2, G2;
    ZZ_pX UU2, VV2, GG2;

    duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
    duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);
    duspoly_t* uu2 = convertPolynomialFromMontgomery_spX (u2, Pptr);
    duspoly_t* vv2 = convertPolynomialFromMontgomery_spX (v2, Pptr);
    duspoly_t* gg2 = convertPolynomialFromMontgomery_spX (g2, Pptr);

    convertToZZ_pX (aa2, A2, Pptr);
    convertToZZ_pX (bb2, B2, Pptr);
    convertToZZ_pX (uu2, UU2, Pptr);
    convertToZZ_pX (vv2, VV2, Pptr);
    convertToZZ_pX (gg2, GG2, Pptr);

    freePolynomial_spX (&aa2);
    freePolynomial_spX (&bb2);
    freePolynomial_spX (&uu2);
    freePolynomial_spX (&vv2);
    freePolynomial_spX (&gg2);

    PlainXGCD (G2, U2, V2, A2, B2);

    if (G2 != GG2 || U2 != UU2 || V2 != VV2) {

        std::cerr << "A2 = " << A2 << std::endl;
        std::cerr << "B2 = " << B2 << std::endl;
        std::cerr << "U2 = " << U2 << std::endl;
        std::cerr << "UU2 = " << UU2 << std::endl;
        std::cerr << "V2 = " << V2 << std::endl;
        std::cerr << "VV2 = " << VV2 << std::endl;
        std::cerr << "G2 = " << G2 << std::endl;
        std::cerr << "GG2 = " << GG2 << std::endl;

        fprintf (stderr, "[DUSP] (2) plainExtGCD... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&a2);
    freePolynomial_spX (&b2);
    freePolynomial_spX (&g2);
    freePolynomial_spX (&u2);
    freePolynomial_spX (&v2);

    /////////////////////// 3
    duspoly_t* a3 = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b3 = randomFullyDensePolynomial_spX (sz2, Pptr);
    duspoly_t* u3 = NULL;
    duspoly_t* v3 = NULL;
    duspoly_t* g3 = NULL;

    plainExtGCDInForm_spX (a3, b3, &u3, &v3, &g3, Pptr);

    ZZ_pX A3;
    ZZ_pX B3;
    ZZ_pX U3, V3, G3;
    ZZ_pX UU3, VV3, GG3;


    duspoly_t* aa3 = convertPolynomialFromMontgomery_spX (a3, Pptr);
    duspoly_t* bb3 = convertPolynomialFromMontgomery_spX (b3, Pptr);
    duspoly_t* uu3 = convertPolynomialFromMontgomery_spX (u3, Pptr);
    duspoly_t* vv3 = convertPolynomialFromMontgomery_spX (v3, Pptr);
    duspoly_t* gg3 = convertPolynomialFromMontgomery_spX (g3, Pptr);

    convertToZZ_pX (aa3, A3, Pptr);
    convertToZZ_pX (bb3, B3, Pptr);
    convertToZZ_pX (uu3, UU3, Pptr);
    convertToZZ_pX (vv3, VV3, Pptr);
    convertToZZ_pX (gg3, GG3, Pptr);

    freePolynomial_spX (&aa3);
    freePolynomial_spX (&bb3);
    freePolynomial_spX (&uu3);
    freePolynomial_spX (&vv3);
    freePolynomial_spX (&gg3);

    PlainXGCD (G3, U3, V3, A3, B3);

    if (G3 != GG3 || U3 != UU3 || V3 != VV3) {

      std::cerr << "A3 = " << A3 << std::endl;
      std::cerr << "B3 = " << B3 << std::endl;
      std::cerr << "U3 = " << U3 << std::endl;
      std::cerr << "UU3 = " << UU3 << std::endl;
      std::cerr << "V3 = " << V3 << std::endl;
      std::cerr << "VV3 = " << VV3 << std::endl;
      std::cerr << "G3 = " << G3 << std::endl;
      std::cerr << "GG3 = " << GG3 << std::endl;

      fprintf (stderr, "[DUSP] (3) plainExtGCD... FAIL\n");

      exit(1);
    }

    freePolynomial_spX (&a3);
    freePolynomial_spX (&b3);
    freePolynomial_spX (&g3);
    freePolynomial_spX (&u3);
    freePolynomial_spX (&v3);


    /////////////////////// 4
    duspoly_t* a4 = randomFullyDensePolynomial_spX (2, Pptr);
    duspoly_t* b4 = randomFullyDensePolynomial_spX (sz2, Pptr);
    duspoly_t* u4 = NULL;
    duspoly_t* v4 = NULL;
    duspoly_t* g4 = NULL;

    plainExtGCDInForm_spX (a4, b4, &u4, &v4, &g4, Pptr);

    ZZ_pX A4;
    ZZ_pX B4;
    ZZ_pX U4, V4, G4;
    ZZ_pX UU4, VV4, GG4;


    duspoly_t* aa4 = convertPolynomialFromMontgomery_spX (a4, Pptr);
    duspoly_t* bb4 = convertPolynomialFromMontgomery_spX (b4, Pptr);
    duspoly_t* uu4 = convertPolynomialFromMontgomery_spX (u4, Pptr);
    duspoly_t* vv4 = convertPolynomialFromMontgomery_spX (v4, Pptr);
    duspoly_t* gg4 = convertPolynomialFromMontgomery_spX (g4, Pptr);

    convertToZZ_pX (aa4, A4, Pptr);
    convertToZZ_pX (bb4, B4, Pptr);
    convertToZZ_pX (uu4, UU4, Pptr);
    convertToZZ_pX (vv4, VV4, Pptr);
    convertToZZ_pX (gg4, GG4, Pptr);

    freePolynomial_spX (&aa4);
    freePolynomial_spX (&bb4);
    freePolynomial_spX (&uu4);
    freePolynomial_spX (&vv4);
    freePolynomial_spX (&gg4);

    PlainXGCD (G4, U4, V4, A4, B4);

    if (G4 != GG4 || U4 != UU4 || V4 != VV4) {

      std::cerr << "A4 = " << A4 << std::endl;
      std::cerr << "B4 = " << B4 << std::endl;
      std::cerr << "U4 = " << U4 << std::endl;
      std::cerr << "UU4 = " << UU4 << std::endl;
      std::cerr << "V4 = " << V4 << std::endl;
      std::cerr << "VV4 = " << VV4 << std::endl;
      std::cerr << "G4 = " << G4 << std::endl;
      std::cerr << "GG4 = " << GG4 << std::endl;

      fprintf (stderr, "[DUSP] (4) plainExtGCD... FAIL\n");

      exit(1);
    }

    freePolynomial_spX (&a4);
    freePolynomial_spX (&b4);
    freePolynomial_spX (&g4);
    freePolynomial_spX (&u4);
    freePolynomial_spX (&v4);


    /////////////////////// 5
    duspoly_t* a5 = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b5 = randomFullyDensePolynomial_spX (1, Pptr);
    duspoly_t* u5 = NULL;
    duspoly_t* v5 = NULL;
    duspoly_t* g5 = NULL;

    plainExtGCDInForm_spX (a5, b5, &u5, &v5, &g5, Pptr);

    ZZ_pX A5;
    ZZ_pX B5;
    ZZ_pX U5, V5, G5;
    ZZ_pX UU5, VV5, GG5;


    duspoly_t* aa5 = convertPolynomialFromMontgomery_spX (a5, Pptr);
    duspoly_t* bb5 = convertPolynomialFromMontgomery_spX (b5, Pptr);
    duspoly_t* uu5 = convertPolynomialFromMontgomery_spX (u5, Pptr);
    duspoly_t* vv5 = convertPolynomialFromMontgomery_spX (v5, Pptr);
    duspoly_t* gg5 = convertPolynomialFromMontgomery_spX (g5, Pptr);

    convertToZZ_pX (aa5, A5, Pptr);
    convertToZZ_pX (bb5, B5, Pptr);
    convertToZZ_pX (uu5, UU5, Pptr);
    convertToZZ_pX (vv5, VV5, Pptr);
    convertToZZ_pX (gg5, GG5, Pptr);

    freePolynomial_spX (&aa5);
    freePolynomial_spX (&bb5);
    freePolynomial_spX (&uu5);
    freePolynomial_spX (&vv5);
    freePolynomial_spX (&gg5);

    PlainXGCD (G5, U5, V5, A5, B5);

    if (!ExtEuclideanTestInForm (aa5, bb5, uu5, vv5, gg5, Pptr) || G5 != GG5) {
        if (G5 != GG5 || U5 != UU5 || V5 != VV5) {

          std::cerr << "A5 = " << A5 << std::endl;
          std::cerr << "B5 = " << B5 << std::endl;
          std::cerr << "U5 = " << U5 << std::endl;
          std::cerr << "UU5 = " << UU5 << std::endl;
          std::cerr << "V5 = " << V5 << std::endl;
          std::cerr << "VV5 = " << VV5 << std::endl;
          std::cerr << "G5 = " << G5 << std::endl;
          std::cerr << "GG5 = " << GG5 << std::endl;

          fprintf (stderr, "[DUSP] (5) plainExtGCD... FAIL\n");

          exit(1);
        }
    }

    freePolynomial_spX (&a5);
    freePolynomial_spX (&b5);
    freePolynomial_spX (&g5);
    freePolynomial_spX (&u5);
    freePolynomial_spX (&v5);


    fprintf (stderr, "[DUSP] plainExtGCDInForm... PASS\n");

    ////////////////////////// 6

        elem_t a6_elems[4] = {-4,8,11,3};
    elem_t b6_elems[5] = {-54, 243,-315,25,125};

    duspoly_t* a6 = NULL;
    duspoly_t* b6 = NULL;
    duspoly_t* u6 = NULL;
    duspoly_t* v6 = NULL;
    duspoly_t* g6 = NULL;

    setCoefs_spX (&a6, a6_elems, 4, Pptr);
    setCoefs_spX (&b6, b6_elems, 5, Pptr);

    extGCDInForm_spX (a6, b6, &u6, &v6, &g6, Pptr);

    ZZ_pX A6;
    ZZ_pX B6;
    ZZ_pX U6, V6, G6;
    ZZ_pX UU6, VV6, GG6;

    duspoly_t* aa6 = convertPolynomialFromMontgomery_spX (a6, Pptr);
    duspoly_t* bb6 = convertPolynomialFromMontgomery_spX (b6, Pptr);
    duspoly_t* uu6 = convertPolynomialFromMontgomery_spX (u6, Pptr);
    duspoly_t* vv6 = convertPolynomialFromMontgomery_spX (v6, Pptr);
    duspoly_t* gg6 = convertPolynomialFromMontgomery_spX (g6, Pptr);

    convertToZZ_pX (aa6, A6, Pptr);
    convertToZZ_pX (bb6, B6, Pptr);
    convertToZZ_pX (uu6, UU6, Pptr);
    convertToZZ_pX (vv6, VV6, Pptr);
    convertToZZ_pX (gg6, GG6, Pptr);

    freePolynomial_spX (&aa6);
    freePolynomial_spX (&bb6);
    freePolynomial_spX (&uu6);
    freePolynomial_spX (&vv6);
    freePolynomial_spX (&gg6);

    PlainXGCD (G6, U6, V6, A6, B6);

    if (G6 != GG6 || U6 != UU6 || V6 != VV6) {

        std::cerr << "A6 = " << A6 << std::endl;
        std::cerr << "B6 = " << B6 << std::endl;
        std::cerr << "U6 = " << U6 << std::endl;
        std::cerr << "UU6 = " << UU6 << std::endl;
        std::cerr << "V6 = " << V6 << std::endl;
        std::cerr << "VV6 = " << VV6 << std::endl;
        std::cerr << "G6 = " << G6 << std::endl;
        std::cerr << "GG6 = " << GG6 << std::endl;

        fprintf (stderr, "[DUSP] extGCDInForm... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&a6);
    freePolynomial_spX (&b6);
    freePolynomial_spX (&g6);
    freePolynomial_spX (&u6);
    freePolynomial_spX (&v6);

    ////////////////////////////// 7
    elem_t a7_elems[10] = {1,-3,3,2,-6,3,3,-3,0,1};
    elem_t b7_elems[4] = {1,-1,0,1};

    duspoly_t* a7 = NULL;
    duspoly_t* b7 = NULL;
    duspoly_t* u7 = NULL;
    duspoly_t* v7 = NULL;
    duspoly_t* g7 = NULL;

    setCoefs_spX (&a7, a7_elems, 10, Pptr);
    setCoefs_spX (&b7, b7_elems, 4, Pptr);

    extGCDInForm_spX (a7, b7, &u7, &v7, &g7, Pptr);

    ZZ_pX A7;
    ZZ_pX B7;
    ZZ_pX U7, V7, G7;
    ZZ_pX UU7, VV7, GG7;

    duspoly_t* aa7 = convertPolynomialFromMontgomery_spX (a7, Pptr);
    duspoly_t* bb7 = convertPolynomialFromMontgomery_spX (b7, Pptr);
    duspoly_t* uu7 = convertPolynomialFromMontgomery_spX (u7, Pptr);
    duspoly_t* vv7 = convertPolynomialFromMontgomery_spX (v7, Pptr);
    duspoly_t* gg7 = convertPolynomialFromMontgomery_spX (g7, Pptr);

    convertToZZ_pX (aa7, A7, Pptr);
    convertToZZ_pX (bb7, B7, Pptr);
    convertToZZ_pX (uu7, UU7, Pptr);
    convertToZZ_pX (vv7, VV7, Pptr);
    convertToZZ_pX (gg7, GG7, Pptr);

    freePolynomial_spX (&aa7);
    freePolynomial_spX (&bb7);
    freePolynomial_spX (&uu7);
    freePolynomial_spX (&vv7);
    freePolynomial_spX (&gg7);

    PlainXGCD (G7, U7, V7, A7, B7);

    if (G7 != GG7 || U7 != UU7 || V7 != VV7) {

        std::cerr << "A7 = " << A7 << std::endl;
        std::cerr << "B7 = " << B7 << std::endl;
        std::cerr << "U7 = " << U7 << std::endl;
        std::cerr << "UU7 = " << UU7 << std::endl;
        std::cerr << "V7 = " << V7 << std::endl;
        std::cerr << "VV7 = " << VV7 << std::endl;
        std::cerr << "G7 = " << G7 << std::endl;
        std::cerr << "GG7 = " << GG7 << std::endl;

        fprintf (stderr, "[DUSP] (7) plainExtGCD... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&a7);
    freePolynomial_spX (&b7);
    freePolynomial_spX (&g7);
    freePolynomial_spX (&u7);
    freePolynomial_spX (&v7);

    ///////////////////////8
    duspoly_t* a8 = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b8 = randomFullyDensePolynomial_spX (sz2, Pptr);
    duspoly_t* u8 = NULL;
    duspoly_t* v8 = NULL;
    duspoly_t* g8 = NULL;

    extGCDInForm_spX (a8, b8, &u8, &v8, &g8, Pptr);

    ZZ_pX A8;
    ZZ_pX B8;
    ZZ_pX U8, V8, G8;
    ZZ_pX UU8, VV8, GG8;


    duspoly_t* aa8 = convertPolynomialFromMontgomery_spX (a8, Pptr);
    duspoly_t* bb8 = convertPolynomialFromMontgomery_spX (b8, Pptr);
    duspoly_t* uu8 = convertPolynomialFromMontgomery_spX (u8, Pptr);
    duspoly_t* vv8 = convertPolynomialFromMontgomery_spX (v8, Pptr);
    duspoly_t* gg8 = convertPolynomialFromMontgomery_spX (g8, Pptr);

    convertToZZ_pX (aa8, A8, Pptr);
    convertToZZ_pX (bb8, B8, Pptr);
    convertToZZ_pX (uu8, UU8, Pptr);
    convertToZZ_pX (vv8, VV8, Pptr);
    convertToZZ_pX (gg8, GG8, Pptr);

    freePolynomial_spX (&aa8);
    freePolynomial_spX (&bb8);
    freePolynomial_spX (&uu8);
    freePolynomial_spX (&vv8);
    freePolynomial_spX (&gg8);

    PlainXGCD (G8, U8, V8, A8, B8);

    if (G8 != GG8 || U8 != UU8 || V8 != VV8) {

      std::cerr << "A8 = " << A8 << std::endl;
      std::cerr << "B8 = " << B8 << std::endl;
      std::cerr << "U8 = " << U8 << std::endl;
      std::cerr << "UU8 = " << UU8 << std::endl;
      std::cerr << "V8 = " << V8 << std::endl;
      std::cerr << "VV8 = " << VV8 << std::endl;
      std::cerr << "G8 = " << G8 << std::endl;
      std::cerr << "GG8 = " << GG8 << std::endl;

      fprintf (stderr, "[DUSP] (8) extGCDInForm... FAIL\n");

      exit(1);
    }

    freePolynomial_spX (&a8);
    freePolynomial_spX (&b8);
    freePolynomial_spX (&g8);
    freePolynomial_spX (&u8);
    freePolynomial_spX (&v8);


    ///////////////////////9
    duspoly_t* a9 = randomFullyDensePolynomial_spX (2, Pptr);
    duspoly_t* b9 = randomFullyDensePolynomial_spX (sz2, Pptr);
    duspoly_t* u9 = NULL;
    duspoly_t* v9 = NULL;
    duspoly_t* g9 = NULL;

    extGCDInForm_spX (a9, b9, &u9, &v9, &g9, Pptr);

    ZZ_pX A9;
    ZZ_pX B9;
    ZZ_pX U9, V9, G9;
    ZZ_pX UU9, VV9, GG9;


    duspoly_t* aa9 = convertPolynomialFromMontgomery_spX (a9, Pptr);
    duspoly_t* bb9 = convertPolynomialFromMontgomery_spX (b9, Pptr);
    duspoly_t* uu9 = convertPolynomialFromMontgomery_spX (u9, Pptr);
    duspoly_t* vv9 = convertPolynomialFromMontgomery_spX (v9, Pptr);
    duspoly_t* gg9 = convertPolynomialFromMontgomery_spX (g9, Pptr);

    convertToZZ_pX (aa9, A9, Pptr);
    convertToZZ_pX (bb9, B9, Pptr);
    convertToZZ_pX (uu9, UU9, Pptr);
    convertToZZ_pX (vv9, VV9, Pptr);
    convertToZZ_pX (gg9, GG9, Pptr);

    freePolynomial_spX (&aa9);
    freePolynomial_spX (&bb9);
    freePolynomial_spX (&uu9);
    freePolynomial_spX (&vv9);
    freePolynomial_spX (&gg9);

    PlainXGCD (G9, U9, V9, A9, B9);

    if (!ExtEuclideanTestInForm (a9, b9, u9, v9, g9, Pptr) || G9 != GG9) {

        if (G9 != GG9 || U9 != UU9 || V9 != VV9) {

          std::cerr << "A9 = " << A9 << std::endl;
          std::cerr << "B9 = " << B9 << std::endl;
          std::cerr << "U9 = " << U9 << std::endl;
          std::cerr << "UU9 = " << UU9 << std::endl;
          std::cerr << "V9 = " << V9 << std::endl;
          std::cerr << "VV9 = " << VV9 << std::endl;
          std::cerr << "G9 = " << G9 << std::endl;
          std::cerr << "GG9 = " << GG9 << std::endl;

          fprintf (stderr, "[DUSP] (9) extGCDInForm... FAIL\n");

          exit(1);
        }
    }

    freePolynomial_spX (&a9);
    freePolynomial_spX (&b9);
    freePolynomial_spX (&g9);
    freePolynomial_spX (&u9);
    freePolynomial_spX (&v9);


    /////////////////////// 10
    duspoly_t* a10 = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b10 = randomFullyDensePolynomial_spX (1, Pptr);
    duspoly_t* u10 = NULL;
    duspoly_t* v10 = NULL;
    duspoly_t* g10 = NULL;

    extGCDInForm_spX (a10, b10, &u10, &v10, &g10, Pptr);

    ZZ_pX A10;
    ZZ_pX B10;
    ZZ_pX U10, V10, G10;
    ZZ_pX UU10, VV10, GG10;


    duspoly_t* aa10 = convertPolynomialFromMontgomery_spX (a10, Pptr);
    duspoly_t* bb10 = convertPolynomialFromMontgomery_spX (b10, Pptr);
    duspoly_t* uu10 = convertPolynomialFromMontgomery_spX (u10, Pptr);
    duspoly_t* vv10 = convertPolynomialFromMontgomery_spX (v10, Pptr);
    duspoly_t* gg10 = convertPolynomialFromMontgomery_spX (g10, Pptr);

    convertToZZ_pX (aa10, A10, Pptr);
    convertToZZ_pX (bb10, B10, Pptr);
    convertToZZ_pX (uu10, UU10, Pptr);
    convertToZZ_pX (vv10, VV10, Pptr);
    convertToZZ_pX (gg10, GG10, Pptr);

    freePolynomial_spX (&aa10);
    freePolynomial_spX (&bb10);
    freePolynomial_spX (&uu10);
    freePolynomial_spX (&vv10);
    freePolynomial_spX (&gg10);

    PlainXGCD (G10, U10, V10, A10, B10);

    if (!ExtEuclideanTestInForm (a10, b10, u10, v10, g10, Pptr) || G10 != GG10) {

        if (G10 != GG10 || U10 != UU10 || V10 != VV10) {

          std::cerr << "A10 = " << A10 << std::endl;
          std::cerr << "B10 = " << B10 << std::endl;
          std::cerr << "U10 = " << U10 << std::endl;
          std::cerr << "UU10 = " << UU10 << std::endl;
          std::cerr << "V10 = " << V10 << std::endl;
          std::cerr << "VV10 = " << VV10 << std::endl;
          std::cerr << "G10 = " << G10 << std::endl;
          std::cerr << "GG10 = " << GG10 << std::endl;

          fprintf (stderr, "[DUSP] (10) extGCDInForm... FAIL\n");

          exit(1);
        }
    }

    freePolynomial_spX (&a10);
    freePolynomial_spX (&b10);
    freePolynomial_spX (&g10);
    freePolynomial_spX (&u10);
    freePolynomial_spX (&v10);


    fprintf (stderr, "[DUSP] extHalfGCDMatrixInForm... PASS\n");
    fprintf (stderr, "[DUSP] extGCDInForm... PASS\n");

   /////////////////////// 11
    if (Pptr->prime == 4179340454199820289) {

        duspoly_t* a11 = randomFullyDensePolynomial_spX (sz1, Pptr);
        duspoly_t* b11 = randomFullyDensePolynomial_spX (sz2, Pptr);
        duspoly_t* u11 = NULL;
        duspoly_t* v11 = NULL;
        duspoly_t* g11 = NULL;

        extGCDInForm_FFT_spX (a11, b11, &u11, &v11, &g11, Pptr);

        ZZ_pX A11;
        ZZ_pX B11;
        ZZ_pX U11, V11, G11;
        ZZ_pX UU11, VV11, GG11;


        duspoly_t* aa11 = convertPolynomialFromMontgomery_spX (a11, Pptr);
        duspoly_t* bb11 = convertPolynomialFromMontgomery_spX (b11, Pptr);
        duspoly_t* uu11 = convertPolynomialFromMontgomery_spX (u11, Pptr);
        duspoly_t* vv11 = convertPolynomialFromMontgomery_spX (v11, Pptr);
        duspoly_t* gg11 = convertPolynomialFromMontgomery_spX (g11, Pptr);

        convertToZZ_pX (aa11, A11, Pptr);
        convertToZZ_pX (bb11, B11, Pptr);
        convertToZZ_pX (uu11, UU11, Pptr);
        convertToZZ_pX (vv11, VV11, Pptr);
        convertToZZ_pX (gg11, GG11, Pptr);

        freePolynomial_spX (&aa11);
        freePolynomial_spX (&bb11);
        freePolynomial_spX (&uu11);
        freePolynomial_spX (&vv11);
        freePolynomial_spX (&gg11);

        PlainXGCD (G11, U11, V11, A11, B11);

        if (G11 != GG11 || U11 != UU11 || V11 != VV11) {

        std::cerr << "A11 = " << A11 << std::endl;
        std::cerr << "B11 = " << B11 << std::endl;
        std::cerr << "U11 = " << U11 << std::endl;
        std::cerr << "UU11 = " << UU11 << std::endl;
        std::cerr << "V11 = " << V11 << std::endl;
        std::cerr << "VV11 = " << VV11 << std::endl;
        std::cerr << "G11 = " << G11 << std::endl;
        std::cerr << "GG11 = " << GG11 << std::endl;

        fprintf (stderr, "[DUSP] (11) extGCDInForm_FFT... FAIL\n");

        exit(1);
        }

        freePolynomial_spX (&a11);
        freePolynomial_spX (&b11);
        freePolynomial_spX (&g11);
        freePolynomial_spX (&u11);
        freePolynomial_spX (&v11);


        fprintf (stderr, "[DUSP] extHalfGCDMatrixInForm_FFT... PASS\n");
        fprintf (stderr, "[DUSP] extGCDInForm_FFT... PASS\n");

    }
}

///////////////////////////////////////////////////////////////////////////////////////

void test_diff (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  duspoly_t* a1 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* a2 = randomPolynomialInForm_spX (sz2, Pptr);
  duspoly_t* a3 = NULL;
  duspoly_t* a4 = randomPolynomialInForm_spX (2, Pptr);

  duspoly_t* da1 = NULL;
  duspoly_t* da2 = deepCopyPolynomial_spX (a2);
  duspoly_t* da3 = NULL;
  duspoly_t* da4 = NULL;

  derivativePolyInForm_spX_inp (&da2, Pptr);
  derivativePolyInForm_spX (a1, &da1, Pptr);
  derivativePolyInForm_spX (a3, &da3, Pptr);
  derivativePolyInForm_spX (a4, &da4, Pptr);

  ZZ_pX A1, DA1, DD1;
  ZZ_pX A2, DA2, DD2;
  ZZ_pX A3, DA3, DD3;
  ZZ_pX A4, DA4, DD4;

  duspoly_t* tmp1 = convertPolynomialFromMontgomery_spX (a1, Pptr);
  duspoly_t* tmp2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
  duspoly_t* tmp3 = convertPolynomialFromMontgomery_spX (a3, Pptr);
  duspoly_t* tmp4 = convertPolynomialFromMontgomery_spX (a4, Pptr);

  duspoly_t* tmp11 = convertPolynomialFromMontgomery_spX (da1, Pptr);
  duspoly_t* tmp22 = convertPolynomialFromMontgomery_spX (da2, Pptr);
  duspoly_t* tmp33 = convertPolynomialFromMontgomery_spX (da3, Pptr);
  duspoly_t* tmp44 = convertPolynomialFromMontgomery_spX (da4, Pptr);

  convertToZZ_pX (tmp1, A1, Pptr);
  convertToZZ_pX (tmp2, A2, Pptr);
  convertToZZ_pX (tmp3, A3, Pptr);
  convertToZZ_pX (tmp4, A4, Pptr);

  convertToZZ_pX (tmp11, DD1, Pptr);
  convertToZZ_pX (tmp22, DD2, Pptr);
  convertToZZ_pX (tmp33, DD3, Pptr);
  convertToZZ_pX (tmp44, DD4, Pptr);

  freePolynomial_spX (&a1);
  freePolynomial_spX (&a2);
  freePolynomial_spX (&a3);
  freePolynomial_spX (&a4);

  freePolynomial_spX (&tmp1);
  freePolynomial_spX (&tmp2);
  freePolynomial_spX (&tmp3);
  freePolynomial_spX (&tmp4);

  freePolynomial_spX (&tmp11);
  freePolynomial_spX (&tmp22);
  freePolynomial_spX (&tmp33);
  freePolynomial_spX (&tmp44);

  diff (DA1, A1);
  diff (DA2, A2);
  diff (DA3, A3);
  diff (DA4, A4);

  if (DA1 != DD1) {
    std::cerr << "A1 = " << A1 << std::endl;
    std::cerr << "DUSP: diff(A1) = " << DD1 << std::endl;
    std::cerr << "NTL:  diff(A1) = " << DA1 << std::endl;

    fprintf (stderr, "[DUSP] (1) derivativePolyInForm... FAIL\n");
    exit (1);
  }

  if (DA2 != DD2) {
    std::cerr << "A2 = " << A2 << std::endl;
    std::cerr << "DUSP: diff(A2) = " << DD2 << std::endl;
    std::cerr << "NTL:  diff(A2) = " << DA2 << std::endl;

    fprintf (stderr, "[DUSP] (2) derivativePolyInForm... FAIL\n");
    exit (1);
  }

  if (DA3 != DD3) {
    std::cerr << "A3 = " << A3 << std::endl;
    std::cerr << "DUSP: diff(A3) = " << DD3 << std::endl;
    std::cerr << "NTL:  diff(A3) = " << DA3 << std::endl;

    fprintf (stderr, "[DUSP] (3) derivativePolyInForm... FAIL\n");
    exit (1);
  }

  if (DA4 != DD4) {
    std::cerr << "A4 = " << A4 << std::endl;
    std::cerr << "DUSP: diff(A4) = " << DD4 << std::endl;
    std::cerr << "NTL:  diff(A4) = " << DA4 << std::endl;

    fprintf (stderr, "[DUSP] (4) derivativePolyInForm... FAIL\n");
    exit (1);
  }

    fprintf (stderr, "[DUSP] derivativePolyInForm... PASS\n");
}

void test_Mat4 (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  dusmat4_t* A = NULL;
  freeMat4_spX (&A);

  A = (dusmat4_t*) malloc (sizeof (dusmat4_t));
  A->polys[0] = randomPolynomialInForm_spX (sz1, Pptr);
  A->polys[1] = randomPolynomialInForm_spX (sz2, Pptr);
  A->polys[2] = randomPolynomialInForm_spX (sz1, Pptr);
  A->polys[3] = randomPolynomialInForm_spX (sz2, Pptr);

  freeMat4_spX (&A);

  if (A != NULL) {
    fprintf(stderr, "[DUSP] freeMat4... FAIL\n");
    exit (1);
  }

    fprintf(stderr, "[DUSP] freeMat4... PASS\n");

  dusmat4_t* B = NULL;

  if (!isEqualMat4_spX (A, B)) {
    fprintf(stderr, "[DUSP] (1) isEqualMat4... FAIL\n");
    exit (1);
  }

  A = (dusmat4_t*) malloc (sizeof (dusmat4_t));
  A->polys[0] = randomPolynomialInForm_spX (sz1, Pptr);
  A->polys[1] = randomPolynomialInForm_spX (sz2, Pptr);
  A->polys[2] = randomPolynomialInForm_spX (sz1, Pptr);
  A->polys[3] = randomPolynomialInForm_spX (sz2, Pptr);

  B = (dusmat4_t*) malloc (sizeof (dusmat4_t));
  B->polys[0] = deepCopyPolynomial_spX (A->polys[0]);
  B->polys[1] = NULL;
  B->polys[2] = deepCopyPolynomial_spX (A->polys[2]);
  B->polys[3] = deepCopyPolynomial_spX (A->polys[3]);

  if (isEqualMat4_spX (A, B)) {
    fprintf(stderr, "[DUSP] (2) isEqualMat4... FAIL\n");
    exit (1);
  }

  B->polys[1] = deepCopyPolynomial_spX (A->polys[1]);

  if (!isEqualMat4_spX (A, B)) {
    fprintf(stderr, "[DUSP] (3) isEqualMat4... FAIL\n");
    exit (1);
  }

  fprintf(stderr, "[DUSP] isEqualMat4... PASS\n");

  freeMat4_spX (&A);
  freeMat4_spX (&B);

  identityMat4InForm_spX (&A, Pptr);
  identityMat4InForm_spX (&B, Pptr);

  if (!isEqualMat4_spX (A, B) || !isIdentityMat4InForm_spX (A, Pptr) || !isIdentityMat4InForm_spX (B, Pptr)) {

      fprintf(stderr, "A = \n");
      printMat4OutForm_spX (A, Pptr);

      fprintf(stderr, "B = \n");
      printMat4OutForm_spX (B, Pptr);

      fprintf(stderr, "[DUSP] identityMat4InForm... FAIL\n");
      fprintf(stderr, "[DUSP] isIdentityMat4InForm... FAIL\n");
      exit (1);
  }

  freeMat4_spX (&A);
  freeMat4_spX (&B);


  fprintf(stderr, "[DUSP] identityMat4InForm... PASS\n");
  fprintf(stderr, "[DUSP] isIdentityMat4InForm... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

void test_res (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    elem_t a1_elems[4] = {-4,8,11,3};
    elem_t b1_elems[5] = {-54, 243,-315,25,125};

    duspoly_t* a1 = NULL;
    duspoly_t* b1 = NULL;
    duspoly_t* r1 = NULL;

    setCoefs_spX (&a1, a1_elems, 4, Pptr);
    setCoefs_spX (&b1, b1_elems, 5, Pptr);

    sylvResultantInForm_spX (a1, b1, &r1, Pptr);

    duspoly_t* hr1 = NULL;
    hgcdResultantInForm_spX (a1, b1, &hr1, Pptr);
    if (!isEqual_spX (r1, hr1)) {
        cerr << "r1 := ";
        printPolynomial_spX (r1, Pptr);
        cerr << "hgcd-r1 := ";
        printPolynomial_spX (hr1, Pptr);
        fprintf (stderr, "[DUSP] hgcdResultantInForm... FAIL\n");
        exit(1);
    }
    // else {
    //     fprintf (stderr, "[DUSP] hgcdResultantInForm... PASS\n");
    // }
    freePolynomial_spX (&hr1);

    ZZ_pX A1;
    ZZ_pX B1;
    ZZ_pX RR1;
    ZZ_p R1;

    duspoly_t* aa1 = convertPolynomialFromMontgomery_spX (a1, Pptr);
    duspoly_t* bb1 = convertPolynomialFromMontgomery_spX (b1, Pptr);
    duspoly_t* rr1 = convertPolynomialFromMontgomery_spX (r1, Pptr);

    convertToZZ_pX (aa1, A1, Pptr);
    convertToZZ_pX (bb1, B1, Pptr);
    convertToZZ_pX (rr1, RR1, Pptr);

    freePolynomial_spX (&aa1);
    freePolynomial_spX (&bb1);
    freePolynomial_spX (&rr1);

    R1 = resultant (A1, B1);

    if (R1 != LeadCoeff (RR1)) {

        std::cerr << "A1 = " << A1 << std::endl;
        std::cerr << "B1 = " << B1 << std::endl;
        std::cerr << "DUSP: res(A1, B1) = " << RR1 << std::endl;
        std::cerr << "NTL:  res(A1, B1) = " << R1 << std::endl;

        fprintf (stderr, "[DUSP] sylvResultantInForm... FAIL\n");

        exit(1);
    }

    freePolynomial_spX (&a1);
    freePolynomial_spX (&b1);
    freePolynomial_spX (&r1);

    ////////////////////////////// 2
    elem_t a2_elems[10] = {1,-3,3,2,-6,3,3,-3,0,1};
    elem_t b2_elems[4] = {1,-1,0,1};

    duspoly_t* a2 = NULL;
    duspoly_t* b2 = NULL;
    duspoly_t* r2 = NULL;

    setCoefs_spX (&a2, a2_elems, 10, Pptr);
    setCoefs_spX (&b2, b2_elems, 4, Pptr);

    sylvResultantInForm_spX (a2, b2, &r2, Pptr);

    duspoly_t* hr2 = NULL;
    hgcdResultantInForm_spX (a2, b2, &hr2, Pptr);
    if (!isEqual_spX (r2, hr2)) {
        cerr << "r2 := ";
        printPolynomial_spX (r2, Pptr);
        cerr << "hgcd-r2 := ";
        printPolynomial_spX (hr2, Pptr);
        fprintf (stderr, "[DUSP] hgcdResultantInForm (2)... FAIL\n");
        exit(1);
    }
    // else {
    //     fprintf (stderr, "[DUSP] hgcdResultantInForm (2)... PASS\n");
    // }
    freePolynomial_spX (&hr2);

    ZZ_pX A2;
    ZZ_pX B2;
    ZZ_pX RR2;
    ZZ_p R2;

    duspoly_t* aa2 = convertPolynomialFromMontgomery_spX (a2, Pptr);
    duspoly_t* bb2 = convertPolynomialFromMontgomery_spX (b2, Pptr);
    duspoly_t* rr2 = convertPolynomialFromMontgomery_spX (r2, Pptr);

    convertToZZ_pX (aa2, A2, Pptr);
    convertToZZ_pX (bb2, B2, Pptr);
    convertToZZ_pX (rr2, RR2, Pptr);

    freePolynomial_spX (&aa2);
    freePolynomial_spX (&bb2);
    freePolynomial_spX (&rr2);

    R2 = resultant (A2, B2);

    if (R2 != LeadCoeff (RR2)) {

        std::cerr << "A2 = " << A2 << std::endl;
        std::cerr << "B2 = " << B2 << std::endl;
        std::cerr << "DUSP: res(A2, B2) = " << RR2 << std::endl;
        std::cerr << "NTL:  res(A2, B2) = " << R2 << std::endl;

        fprintf (stderr, "[DUSP] sylvResultantInForm... FAIL\n");

        exit(2);
    }

    freePolynomial_spX (&a2);
    freePolynomial_spX (&b2);
    freePolynomial_spX (&r2);

    /////////////////////// 3
    duspoly_t* a3 = randomFullyDensePolynomial_spX (sz1, Pptr);
    duspoly_t* b3 = randomFullyDensePolynomial_spX (sz2, Pptr);
    duspoly_t* r3 = NULL;

    sylvResultantInForm_spX (a3, b3, &r3, Pptr);

    duspoly_t* hr3 = NULL;
    hgcdResultantInForm_spX (a3, b3, &hr3, Pptr);
    if (!isEqual_spX (r3, hr3)) {
        cerr << "r3 := ";
        printPolynomial_spX (r3, Pptr);
        cerr << "hgcd-r3 := ";
        printPolynomial_spX (hr3, Pptr);
        fprintf (stderr, "[DUSP] hgcdResultantInForm (3)... FAIL\n");
        exit(1);
    }
    // else {
    //     fprintf (stderr, "[DUSP] hgcdResultantInForm (3)... PASS\n");
    // }
    freePolynomial_spX (&hr3);

    ZZ_pX A3;
    ZZ_pX B3;
    ZZ_pX RR3;
    ZZ_p R3;

    duspoly_t* aa3 = convertPolynomialFromMontgomery_spX (a3, Pptr);
    duspoly_t* bb3 = convertPolynomialFromMontgomery_spX (b3, Pptr);
    duspoly_t* rr3 = convertPolynomialFromMontgomery_spX (r3, Pptr);

    convertToZZ_pX (aa3, A3, Pptr);
    convertToZZ_pX (bb3, B3, Pptr);
    convertToZZ_pX (rr3, RR3, Pptr);

    freePolynomial_spX (&aa3);
    freePolynomial_spX (&bb3);
    freePolynomial_spX (&rr3);

    R3 = resultant (A3, B3);

    if (R3 != LeadCoeff (RR3)) {

        std::cerr << "A3 = " << A3 << std::endl;
        std::cerr << "B3 = " << B3 << std::endl;
        std::cerr << "DUSP: res(A3, B3) = " << RR3 << std::endl;
        std::cerr << "NTL:  res(A3, B3) = " << R3 << std::endl;

        fprintf (stderr, "[DUSP] sylvResultantInForm... FAIL\n");

        exit(3);
    }

    freePolynomial_spX (&a3);
    freePolynomial_spX (&b3);
    freePolynomial_spX (&r3);

    /////////////////////// 4
    duspoly_t* a4 = randomFullyDensePolynomial_spX (2, Pptr);
    duspoly_t* b4 = randomFullyDensePolynomial_spX (sz2, Pptr);
    duspoly_t* r4 = NULL;

    sylvResultantInForm_spX (a4, b4, &r4, Pptr);

    duspoly_t* hr4 = NULL;
    hgcdResultantInForm_spX (a4, b4, &hr4, Pptr);
    if (!isEqual_spX (r4, hr4)) {
        cerr << "r4 := ";
        printPolynomial_spX (r4, Pptr);
        cerr << "hgcd-r4 := ";
        printPolynomial_spX (hr4, Pptr);
        fprintf (stderr, "[DUSP] hgcdResultantInForm (4)... FAIL\n");
        exit(1);
    }
    // else {
    //     fprintf (stderr, "[DUSP] hgcdResultantInForm (2)... PASS\n");
    // }
    freePolynomial_spX (&hr2);

    ZZ_pX A4;
    ZZ_pX B4;
    ZZ_pX RR4;
    ZZ_p R4;

    duspoly_t* aa4 = convertPolynomialFromMontgomery_spX (a4, Pptr);
    duspoly_t* bb4 = convertPolynomialFromMontgomery_spX (b4, Pptr);
    duspoly_t* rr4 = convertPolynomialFromMontgomery_spX (r4, Pptr);

    convertToZZ_pX (aa4, A4, Pptr);
    convertToZZ_pX (bb4, B4, Pptr);
    convertToZZ_pX (rr4, RR4, Pptr);

    freePolynomial_spX (&aa4);
    freePolynomial_spX (&bb4);
    freePolynomial_spX (&rr4);

    R4 = resultant (A4, B4);

    if (R4 != LeadCoeff (RR4)) {

        std::cerr << "A4 = " << A4 << std::endl;
        std::cerr << "B4 = " << B4 << std::endl;
        std::cerr << "DUSP: res(A4, B4) = " << RR4 << std::endl;
        std::cerr << "NTL:  res(A4, B4) = " << R4 << std::endl;

        fprintf (stderr, "[DUSP] sylvResultantInForm... FAIL\n");

        exit(4);
    }

    freePolynomial_spX (&a4);
    freePolynomial_spX (&b4);
    freePolynomial_spX (&r4);


    fprintf (stderr, "[DUSP] sylvResultantInForm... PASS\n");
    fprintf (stderr, "[DUSP] hgcdResultantInForm... PASS\n");
    fprintf (stderr, "[DUSP] sylvSubResultantInForm... PASS\n");
}

///////////////////////////////////////////////////////////////////////////////////////

void test_FactBasics (polysize_t sz1, polysize_t  sz2,Prime_ptr* Pptr, polysize_t flag)
{
  factors_t* f1 = NULL;

  freeFactors_spX (&f1);

  if (f1 != NULL) {
    fprintf(stderr, "[DUSP] freeFactors for f1... FAIL\n");
    exit(1);
  }

  f1 = initFactors_spX (sz1);

  if (f1->alloc != sz1) {
    fprintf(stderr, "[DUSP] initFactors for f1... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] initFactors... PASS\n");

  freeFactors_spX (&f1);
  freeFactors_spX (&f1);
  freeFactors_spX (&f1);

  if (!isZeroFactors_spX (f1)){
    fprintf(stderr, "[DUSP] freeFactors for f11... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] freeFactors... PASS\n");

  f1 = initFactors_spX (sz1);

  for (polysize_t i = 0; i < f1->alloc; i++){
    f1->polys[i] = randomFullyDensePolynomial_spX (sz2, Pptr);
    f1->exps[i] = i;
  }

  factors_t* f2 = deepCopyFactors_spX (f1);

  if (f1->alloc != f2->alloc) {
    fprintf(stderr, "[DUSP] deepCopyFactors for f1/f2->alloc ... FAIL\n");
    exit(1);
  }

  for (polysize_t i = 0; i < f1->alloc; i++) {
    if (!isEqual_spX (f1->polys[i], f2->polys[i])) {

      fprintf(stderr, "[DUSP] f1->polys[%lu] = \n", i);
      printPolynomial_spX (f1->polys[i], Pptr);

      fprintf(stderr, "[DUSP] f2->polys[%lu] = \n", i);
      printPolynomial_spX (f2->polys[i], Pptr);

      fprintf(stderr, "[DUSP] deepCopyFactors for f1/f2->polys[%lu] ... FAIL\n", i);
      exit(1);
    }
  }

  fprintf(stderr, "[DUSP] deepCopyFactors... PASS\n");

  if (!isEqualFactors_spX (f1, f2)) {
    fprintf(stderr, "[DUSP] isEqualFactors for f1/f2... FAIL\n");
    exit(1);
  }

  freeFactors_spX (&f1);

  if (!isZeroFactors_spX (f1)) {
    fprintf(stderr, "[DUSP] isZeroFactors for f1... FAIL\n");
    exit(1);
  }

  freeFactors_spX (&f2);

  if (!isZeroFactors_spX (f2)) {
    fprintf(stderr, "[DUSP] isZeroFactors for f2... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] isZeroFactors... PASS\n");

  factors_t* f3 = initFactors_spX (sz1);
  factors_t* f4 = initFactors_spX (sz1);

  if (!isEqualFactors_spX (f3, f4)) {
    fprintf(stderr, "[DUSP] isEqualFactors for f3/f4... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] isEqualFactors... PASS\n");


  // ///////////////////////////

  factsll_t* ff1 = NULL;

  freeFactsll_spX (&ff1);

  if (ff1 != NULL) {
    fprintf(stderr, "[DUSP] freeFactsll for ff1... FAIL\n");
    exit(1);
  }

  ff1 = initFactsll_spX (sz1);

  if (ff1->exp != sz1) {
    fprintf(stderr, "[DUSP] initFactsll for ff1... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] initFactsll... PASS\n");

  freeFactsll_spX (&ff1);
  freeFactsll_spX (&ff1);
  freeFactsll_spX (&ff1);

  if (!isZeroFactsll_spX (ff1)){
    fprintf(stderr, "[DUSP] freeFactsll for ff1... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] freeFactsll... PASS\n");

  ff1 = initFactsll_spX (sz1);
  ff1->poly = randomPolynomialInForm_spX (sz1, Pptr);

  factsll_t* ff2 = deepCopyFactsll_spX (ff1);

  if (ff1->exp != ff2->exp || !isEqual_spX (ff1->poly, ff2->poly)) {
    fprintf(stderr, "[DUSP] deepCopyFactsll for ff1/ff2->exp ... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] deepCopyFactsll... PASS\n");

  if (!isEqualFactsll_spX (ff1, ff2)) {
    fprintf(stderr, "[DUSP] isEqualFactsll for ff1/ff2... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] isEqualFactsll... PASS\n");

  freeFactsll_spX (&ff1);

  if (!isZeroFactsll_spX (ff1)) {
    fprintf(stderr, "[DUSP] isZeroFactsll for ff1... FAIL\n");
    exit(1);
  }

  freeFactsll_spX (&ff2);

  if (!isZeroFactsll_spX (ff2)) {
    fprintf(stderr, "[DUSP] isZeroFactsll for ff2... FAIL\n");
    exit(1);
  }

  fprintf(stderr, "[DUSP] isZeroFactsll... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

void test_SFD (polysize_t sz1, polysize_t  sz2,Prime_ptr* Pptr, polysize_t flag)
{
  elem_t lc1 = 0;
  duspoly_t* a1 = randomPolynomialInForm_spX (sz1, Pptr);
  monicPolynomialInForm_spX_inp (&a1, &lc1, Pptr);

  factsll_t* f1 = NULL;

  squareFreeFactorizationInFormll_spX (a1, &f1, Pptr);

  duspoly_t* aa1 = convertPolynomialFromMontgomery_spX (a1, Pptr);

  ZZ_pX A1;

  convertToZZ_pX (aa1, A1, Pptr);

  freePolynomial_spX (&aa1);

  vec_pair_ZZ_pX_long U1;
  SquareFreeDecomp (U1, A1);

  ZZ_pX F1_i;

  factsll_t* cur1 = f1;
  duspoly_t* ff1_i = NULL;
  int isThere1;

  while (cur1 != NULL) {
    ff1_i = convertPolynomialFromMontgomery_spX (cur1->poly, Pptr);
    convertToZZ_pX (ff1_i, F1_i, Pptr);


    if (isOne_spX (ff1_i)) {

      freePolynomial_spX (&ff1_i);
      F1_i = 0;

      cur1 = cur1->next;
      continue;
    }

    isThere1 = 0;
    for (polysize_t i = 0; i < U1.length(); i++) {
      if (U1[i].b == cur1->exp) {
          isThere1 = 1;
          if (U1[i].a != F1_i) {
            std::cerr << "A1 = ";
            printPolynomialOutForm_spX (a1, Pptr);

            std::cerr << "Cursor Exp = " << cur1->exp << std::endl;
            std::cerr << "NTL:  U1["<<i<<"] = " << U1[i].a << std::endl;
            std::cerr << "DUSP: F1_i = " << F1_i << std::endl;

            std::cerr << "NTL: U1 = ";
            for (long j = 0; j < U1.length(); j++) {
              std::cerr << "(" << U1[j].a << ")^" << U1[j].b << " * ";
            }
            std::cerr << std::endl;

            std::cerr << "DUSP F1 = ";
            purePrintFactsllOutForm_spX (f1, Pptr);

            fprintf(stderr, "[DUSP] squareFreeFactorizationInFormll... FAIL\n");
            exit (1);
          }
        break;
      }
    }

    if (!isThere1) {
      std::cerr << "A1 = ";
      printPolynomialOutForm_spX (a1, Pptr);

      std::cerr << "Cursor Exp = " << cur1->exp << std::endl;
      std::cerr << "DUSP: F1_i = " << F1_i << std::endl;
      fprintf(stderr, " ******* F1_i couldn't find in NTL output! *******\n");
      std::cerr << std::endl;

      std::cerr << "NTL: U1 = ";
      for (long j = 0; j < U1.length(); j++) {
        std::cerr << "(" << U1[j].a << ")^" << U1[j].b << " * ";
      }
      std::cerr << std::endl;

      std::cerr << "DUSP F1 = ";
      purePrintFactsllOutForm_spX (f1, Pptr);

      fprintf(stderr, "[DUSP] squareFreeFactorizationInFormll... FAIL\n");
      exit (1);
    }

    freePolynomial_spX (&ff1_i);
    F1_i = 0;
    cur1 = cur1->next;
  }

      fprintf(stderr, "[DUSP] squareFreeFactorizationInFormll... PASS\n");

}

///////////////////////////////////////////////////////////////////////////////////////

// void test_Fact (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
// {

//   prime_t pr =  Pptr->prime;
//   // Pptr = smallprimefield_get_prime_constants (pr);
//   // sz1 = 20;

//   elem_t lc1 = 0;
//   duspoly_t* a1 = randomPolynomialInForm_spX (sz1, Pptr);
//   monicPolynomialInForm_spX_inp (&a1, &lc1, Pptr);

//   factsll_t* f1 = NULL;

//   // squareFreeFactorizationInFormll_spX (a1, &f1, Pptr);
//   std::cerr << "modFactoring starts..." << std::endl;
//   modFactorizationInFormll1_spX (a1, &f1, Pptr);
//   std::cerr << "modFactoring ends..." << std::endl;

//   duspoly_t* aa1 = convertPolynomialFromMontgomery_spX (a1, Pptr);

//   ZZ_pX A1;

//   convertToZZ_pX (aa1, A1, Pptr);

//   freePolynomial_spX (&aa1);

//   vec_pair_ZZ_pX_long U1;
//   berlekamp (U1, A1);

//   ZZ_pX F1_i;

//   factsll_t* cur1 = f1;
//   duspoly_t* ff1_i = NULL;
//   int isThere1;

//   while (0) { // TODO: cur1 != NULL
//     ff1_i = convertPolynomialFromMontgomery_spX (cur1->poly, Pptr);
//     convertToZZ_pX (ff1_i, F1_i, Pptr);


//     if (isOne_spX (ff1_i)) {

//       freePolynomial_spX (&ff1_i);
//       F1_i = 0;

//       cur1 = cur1->next;
//       continue;
//     }

//     isThere1 = 0;
//     for (polysize_t i = 0; i < U1.length(); i++) {
//       if (U1[i].b == cur1->exp) {
//           isThere1 = 1;
//           if (U1[i].a != F1_i) {
//             std::cerr << "A1 = ";
//             printPolynomialOutForm_spX (a1, Pptr);

//             std::cerr << "Cursor Exp = " << cur1->exp << std::endl;
//             std::cerr << "NTL:  U1["<<i<<"] = " << U1[i].a << std::endl;
//             std::cerr << "DUSP: F1_i = " << F1_i << std::endl;

//             std::cerr << "NTL: U1 = ";
//             for (long j = 0; j < U1.length(); j++) {
//               std::cerr << "(" << U1[j].a << ")^" << U1[j].b << " * ";
//             }
//             std::cerr << std::endl;

//             std::cerr << "DUSP F1 = ";
//             purePrintFactsllOutForm_spX (f1, Pptr);

//             fprintf(stderr, "[DUSP] (p= %lld) modFactorizationInFormll... FAIL\n", pr);
//             exit (1);
//           }
//         break;
//       }
//     }

//     if (!isThere1) {
//       std::cerr << "A1 = ";
//       printPolynomialOutForm_spX (a1, Pptr);

//       std::cerr << "Cursor Exp = " << cur1->exp << std::endl;
//       std::cerr << "DUSP: F1_i = " << F1_i << std::endl;
//       fprintf(stderr, " ******* F1_i couldn't find in NTL output! *******\n");
//       std::cerr << std::endl;

//       std::cerr << "NTL: U1 = ";
//       for (long j = 0; j < U1.length(); j++) {
//         std::cerr << "(" << U1[j].a << ")^" << U1[j].b << " * ";
//       }
//       std::cerr << std::endl;

//       std::cerr << "DUSP F1 = ";
//       purePrintFactsllOutForm_spX (f1, Pptr);

//       fprintf(stderr, "[DUSP] (p= %lld) modFactorizationInFormll... FAIL\n", pr);
//       exit (1);
//     }

//     freePolynomial_spX (&ff1_i);
//     F1_i = 0;
//     cur1 = cur1->next;
//   }

//   if (modFactorVerificationInFormll1_spX (a1, f1, Pptr)) {
//         fprintf(stderr, "[DUSP] distinctDegFactorizationInFormll... PASS\n");
//         fprintf(stderr, "[DUSP] equalDegFactorizationInFromll... PASS\n");
//       fprintf(stderr, "[DUSP] modFactorizationInFormll... PASS\n");
//   } else {
//       std::cerr << "NTL: U1 = ";
//       for (long j = 0; j < U1.length(); j++) {
//         std::cerr << "(" << U1[j].a << ")^" << U1[j].b << " * ";
//       }
//       std::cerr << std::endl;

//       std::cerr << "DUSP F1 = ";
//       purePrintFactsllOutForm_spX (f1, Pptr);

//       fprintf(stderr, "[DUSP] (p= %lld) modFactorizationInFormll... FAIL\n", pr);
//       exit (1);
//   }
// }

void test_EDFs (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{


  prime_t pr =  3;
  Pptr = smallprimefield_get_prime_constants (pr);
  sz1 = 20;

  elem_t lc1 = 0;
  duspoly_t* a1 = randomPolynomialInForm_spX (3, Pptr);
  monicPolynomialInForm_spX_inp (&a1, &lc1, Pptr);

  elem_t lc2 = 0;
  duspoly_t* a2 = randomPolynomialInForm_spX (5, Pptr);
  monicPolynomialInForm_spX_inp (&a2, &lc2, Pptr);

  elem_t lc3 = 0;
  duspoly_t* a3 = randomPolynomialInForm_spX (10, Pptr);
  monicPolynomialInForm_spX_inp (&a3, &lc3, Pptr);


  factsll_t* G1 = NULL;
  factsll_t* G2 = NULL;
  factsll_t* G3 = NULL;

  distinctDegFactorizationInFormll_spX (a1, &G1, Pptr);
  distinctDegFactorizationInFormll_spX (a2, &G2, Pptr);
  distinctDegFactorizationInFormll_spX (a3, &G3, Pptr);

  freePolynomial_spX (&a1);
  freePolynomial_spX (&a2);
  freePolynomial_spX (&a3);

  std::cerr << "G1 = ";
  purePrintFactsllOutForm_spX (G1, Pptr);
  std::cerr << "G2 = ";
  purePrintFactsllOutForm_spX (G2, Pptr);
  std::cerr << "G3 = ";
  purePrintFactsllOutForm_spX (G3, Pptr);

  factsll_t* G1G2 = concatFactsll_spX (G1, G2);
  factsll_t* G1G2G3 = concatFactsll_spX (G1G2, G3);

  std::cerr << "G1G2 = ";
  purePrintFactsllOutForm_spX (G1G2, Pptr);

  std::cerr << "G1G2G3 = ";
  purePrintFactsllOutForm_spX (G1G2G3, Pptr);

  freeFactsll_spX (&G1);
  freeFactsll_spX (&G2);
  freeFactsll_spX (&G3);

  freeFactsll_spX (&G1G2);

  factsll_t* D = NULL;

  equalDegFactorsInFormll_spX (G1G2G3, &D, Pptr);

  std::cerr << "D = ";
  purePrintFactsllOutForm_spX (D, Pptr);

  freeFactsll_spX (&G1G2G3);

  factors_t* Facts = convertToFactorsList_spX (D, Pptr);

  freeFactsll_spX (&D);

  purePrintFactorsOutForm_spX (Facts, Pptr);

  freeFactors_spX (&Facts);

}



///////////////////////////////////////////////////////////////////////////////////////

// void test_NTLBas (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
// {

//   duspoly_t* a = randomPolynomialInForm_spX (sz1, Pptr);
//   duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);

//   duspoly_t *c1, *c2;

//   mulPolynomialsInForm_spX (a, b, &c1, Pptr);
//   ntlMulPolynomialsInForm_spX (a, b, &c2, Pptr);

//   if (!isEqual_spX (c1, c2)) {
//       fprintf(stderr, "a = \n");
//       printPolynomial_spX  (a, Pptr);

//       fprintf(stderr, "b = \n");
//       printPolynomial_spX (b, Pptr);

//       fprintf(stderr, "[DUSP] a*b = \n");
//       printPolynomial_spX (c1, Pptr);

//       fprintf(stderr, "[DNTL] a*b = \n");
//       printPolynomial_spX (c2, Pptr);

//       fprintf(stderr, "[DUSP] ntlMulPolynomialsInForm... FAIL\n");
//       exit (1);
//   }

//   fprintf(stderr, "[DUSP] ntlMulPolynomialsInForm... PASS\n");


//   freePolynomial_spX (&c1);
//   freePolynomial_spX (&c2);

//   sqrPolynomialInForm_spX (a, &c1, Pptr);
//   ntlSqrPolynomialInForm_spX (a, &c2, Pptr);

//   if (!isEqual_spX (c1, c2)) {
//       fprintf(stderr, "a = \n");
//       printPolynomial_spX  (a, Pptr);

//       fprintf(stderr, "[DUSP] a^2 = \n");
//       printPolynomial_spX (c1, Pptr);

//       fprintf(stderr, "[DNTL] a^2 = \n");
//       printPolynomial_spX (c2, Pptr);

//       fprintf(stderr, "[DUSP] ntlSqrPolynomialInForm... FAIL\n");
//       exit (1);
//   }

//   fprintf(stderr, "[DUSP] ntlSqrPolynomialInForm... PASS\n");

//   freePolynomial_spX (&c1);
//   freePolynomial_spX (&c2);

//   duspoly_t *r1, *r2;
//   duspoly_t *q1, *q2;

//   plainDivPolynomialsInForm_spX (a, b, &r1, &q1, Pptr);
//   ntlDivPolynomialsInForm_spX (a, b, &r2, &q2, Pptr);

//   if (!isEqual_spX (r1, r2)) {
//       fprintf(stderr, "a = \n");
//       printPolynomial_spX  (a, Pptr);

//       fprintf(stderr, "b = \n");
//       printPolynomial_spX  (b, Pptr);

//       fprintf(stderr, "[DUSP] r1 = \n");
//       printPolynomial_spX (r1, Pptr);

//       fprintf(stderr, "[DNTL] r2 = \n");
//       printPolynomial_spX (r2, Pptr);

//       fprintf(stderr, "[DUSP] ntlDivPolynomialInForm... FAIL\n");
//       exit (1);
//   }

//   if (!isEqual_spX (q1, q2)) {
//       fprintf(stderr, "a = \n");
//       printPolynomial_spX  (a, Pptr);

//       fprintf(stderr, "b = \n");
//       printPolynomial_spX  (b, Pptr);

//       fprintf(stderr, "[DUSP] q1 = \n");
//       printPolynomial_spX (q1, Pptr);

//       fprintf(stderr, "[DNTL] q2 = \n");
//       printPolynomial_spX (q2, Pptr);

//       fprintf(stderr, "[DUSP] ntlDivPolynomialInForm... FAIL\n");
//       exit (1);
//   }

//   fprintf(stderr, "[DUSP] ntlDivPolynomialInForm... PASS\n");

//   freePolynomial_spX (&r1);
//   freePolynomial_spX (&r2);

//   freePolynomial_spX (&q1);
//   freePolynomial_spX (&q2);

//   plainGCDInForm_spX (a, b, &q1, Pptr);
//   ntlGCDInForm_spX (a, b, &q2, Pptr);


//   if (!isEqual_spX (q1, q2)) {
//     fprintf(stderr, "a = ");
//     printPolynomial_spX (a,Pptr);

//     fprintf(stderr, "b = ");
//     printPolynomial_spX (b, Pptr);

//     fprintf(stderr, "[DUSP] gcd =  ");
//     printPolynomial_spX (q1, Pptr);

//     fprintf(stderr, "[DNTL] gcd = ");
//     printPolynomial_spX (q2, Pptr);

//       fprintf(stderr, "[DUSP] ntlGCDInForm... FAIL\n");
//       exit (1);
//   }

//   fprintf(stderr, "[DUSP] ntlGCDInForm... PASS\n");

//   freePolynomial_spX (&q1);
//   freePolynomial_spX (&q2);

//   duspoly_t *u1, *v1, *g1;
//   duspoly_t *u2, *v2, *g2;

//   extGCDInForm_spX (a, b, &u1, &v1, &g1, Pptr);
//   ntlExtGCDInForm_spX (a, b, &u2, &v2, &g2, Pptr);

//   if (!isEqual_spX (u1, u2) || !isEqual_spX (v1, v2) || !isEqual_spX (g1, g2)) {
//     fprintf(stderr, "a = ");
//     printPolynomial_spX (a, Pptr);

//     fprintf(stderr, "b = ");
//     printPolynomial_spX (b, Pptr);

//     fprintf(stderr, "g1 = ");
//     printPolynomial_spX (g1, Pptr);

//     fprintf(stderr, "g2 = ");
//     printPolynomial_spX (g2, Pptr);

//     fprintf(stderr, "u1 = ");
//     printPolynomial_spX (u1, Pptr);

//     fprintf(stderr, "u2 = ");
//     printPolynomial_spX (u2, Pptr);

//     fprintf(stderr, "v1 = ");
//     printPolynomial_spX (v1, Pptr);

//     fprintf(stderr, "v2 = ");
//     printPolynomial_spX (v2, Pptr);

//     fprintf(stderr, "[DUSP] ntlExtGCDInForm... FAIL\n");
//     exit (1);
//   }

//   freePolynomial_spX (&a);
//   freePolynomial_spX (&b);
//   freePolynomial_spX (&g1);
//   freePolynomial_spX (&g2);
//   freePolynomial_spX (&u1);
//   freePolynomial_spX (&u2);
//   freePolynomial_spX (&v1);
//   freePolynomial_spX (&v2);

//   fprintf(stderr, "[DUSP] ntlExtGCDInForm... PASS\n");
// }

// void test_NTLFact (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
// {
//     duspoly_t* a = randomPolynomialInForm_spX (sz1, Pptr);
//     elem_t lc;
//     monicPolynomialInForm_spX_inp (&a, &lc, Pptr);

//     factors_t* facts = NULL;

//     ntlFactoringInForm_spX (a, &facts, Pptr);

//     // purePrintFactorsOutForm_spX (facts, Pptr);
//     fprintf(stderr, "[DUSP] ntlFactoringInForm... PASS\n");

//     freePolynomial_spX (&a);
//     freeFactors_spX (&facts);
// }

// void test_interpolation(polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
// {
//   duspoly_t* a = randomPolynomialInForm_spX(sz1, Pptr);

//   int n = a->lt;
//   elem_t* t = (elem_t*) malloc(sizeof(elem_t)*n+1);
//   elem_t* v = (elem_t*) malloc(sizeof(elem_t)*n+1);

//   duspoly_t div;
//   elem_t divElems[2];
//   divElems[1] = smallprimefield_convert_in(1, Pptr);
//   div.elems = divElems;
//   div.lt = 1;
//   div.alloc = 2;

//   elem_t neg1 = smallprimefield_convert_in(-1, Pptr);
//   int curPoint = 5;
//   duspoly_t* r = NULL;
//   for (int i = 0; i <= n; ++i) {
//     curPoint += (rand() % 100) + 1;
//     t[i] = smallprimefield_convert_in(curPoint, Pptr);
//     divElems[0] = smallprimefield_mul(t[i], neg1, Pptr);
//     plainRemPolynomials_spX(a, &div, &r, Pptr);
//     v[i] = r->elems[0];
//     freePolynomial_spX(&r);
//   }

//   duspoly_t* p = interpolatePolynomialInForm_spX(t, v, n+1, Pptr);

//   for (int i = 0; i <= n; ++i) {
//     if (a->elems[i] != p->elems[i]) {
//       fprintf(stderr, "[DUSP] interpolatePolynomialInForm... FAIL\n");
//       fprintf(stderr, "a: ");
//       printPolynomialOutForm_spX(a, Pptr);
//       fprintf(stderr, "p: ");
//       printPolynomialOutForm_spX(p, Pptr);
//       exit(1);
//     }
//   }

//   fprintf(stderr, "[DUSP] interpolatePolynomialInForm... PASS\n");

// }

void test_memmove (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{

  elem_t a[5] = {1,2,3,4,5};
  for (int i = 0; i < 5; i++){
    std::cout << "a["<<i<<"] = " << a[i] << std::endl;
  }

  elem_t* b = (elem_t*) calloc (3, sizeof (elem_t));

  memmove (b, a+2, 3*sizeof(elem_t));


}


void test_duspolyA (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  duspolysA_t* aa = NULL;
  freePolysA_spX (aa);

  aa = makePolysA_spX (sz1);
  freePolysA_spX (aa);

  aa = makePolysA_spX (sz2);
  for (polysize_t i = 0; i < sz2; i++) {
    // aa->polys[i] = randomFullyDensePolynomial_spX (sz2+1, Pptr);
    aa->polys[i] = randomPolynomialInForm_spX (i+2, Pptr);
  }

  for (polysize_t i = 0; i < sz2; i++) {
    cerr << "aa["<<i<<"] = ";
    printPolynomial_spX (aa->polys[i], Pptr);
  }

  duspolysA_t* aa_T = transposePolysA_spX (aa, sz2);

  for (polysize_t i = 0; i < aa_T->size; i++) {
    cerr << "aa_T["<<i<<"] = ";
    printPolynomial_spX (aa_T->polys[i], Pptr);
  }

  duspolysA_t* aa_TT = transposePolysA_spX (aa_T, sz2-1);

  for (polysize_t i = 0; i < aa_TT->size; i++) {
    cerr << "aa_TT["<<i<<"] = ";
    printPolynomial_spX (aa_TT->polys[i], Pptr);
  }

  fprintf (stderr, "isEqual(aa, aa_TT) == %d\n", isPolysAEqual_spX (aa, aa_TT));


  freePolysA_spX (aa_TT);
  freePolysA_spX (aa_T);
  freePolysA_spX (aa);


  // aa = makePolysA_spX (sz2);
  // DUZP_t* a = modularCRTs_DUZP (aa, sz1);
}

void test_MPZ (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  int n = sz1;
  mpz_t m; // product (m_i, i=1,...n)
	mpz_init_set_ui (m, 1l);
	// mpz_t mpz_pr;
	for (int i = 0; i < n; i++) {
		mpz_mul_si (m, m,  (&prime64_ptr[i])->prime);

    fprintf (stderr, "m[%d] = %lld\n", i, (&prime64_ptr[i])->prime);
	}
  gmp_printf ("\n\nMul(m) = %Zd\n\n", m);

	mpz_t divm[n];
	for (int i = 0; i < n; i++) {
		mpz_init (divm[i]);
		mpz_divexact_ui (divm[i], m, (unsigned long) (&prime64_ptr[i])->prime);
	}

  for (int i = 0; i < n; i++) {
    gmp_printf ("divm[%d] = %Zd\n", i, divm[i]);
  }
}

void test_modularRes (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{

    int degA = 9;
    DUZP_t* a = makePolynomial_DUZP (degA+1);
    a->lt = degA;
    mpz_init_set_si (a->coefs[0], 20);
    mpz_init_set_si (a->coefs[1], 3);
    mpz_init_set_si (a->coefs[2], 3);
    mpz_init_set_si (a->coefs[3], 3);
    mpz_init_set_si (a->coefs[4], 3);
    mpz_init_set_si (a->coefs[5], 3);
    mpz_init_set_si (a->coefs[6], 3);
    mpz_init_set_si (a->coefs[7], 3);
    mpz_init_set_si (a->coefs[8], 3);
    mpz_init_set_si (a->coefs[9], 3);

    int degB = 8;
    DUZP_t* b = makePolynomial_DUZP (degB+1);
    b->lt = degB;
    mpz_init_set_si (b->coefs[0], 20);
    mpz_init_set_si (b->coefs[1], 5);
    mpz_init_set_si (b->coefs[2], 7);
    mpz_init_set_si (b->coefs[3], 2);
    mpz_init_set_si (b->coefs[4], 17);
    mpz_init_set_si (b->coefs[5], 234);
    mpz_init_set_si (b->coefs[6], 111);
    mpz_init_set_si (b->coefs[7], 2);
    mpz_init_set_si (b->coefs[8], 12);


  // DUZP_t* a = makePolynomial_DUZP (3);
  // mpz_init_set_si (a->coefs[0], 1000000000);
  // // mpz_init_set_si (a->coefs[1], 9223372036854602819);
  // mpz_init_set_si (a->coefs[1], 20000000);
  // mpz_init_set_si (a->coefs[2], 300000000);
  // a->lt = 2;

  // DUZP_t* b = makePolynomial_DUZP (2);
  // mpz_init_set_si (b->coefs[0], 10000000000);
  // mpz_init_set_si (b->coefs[1], 5000000000000);
  // b->lt = 1;

  char sym[1] = {'x'};
  cerr << "a = ";
  printPoly_DUZP (a, sym);
  cerr << "b = ";
  printPoly_DUZP (b, sym);

  mpz_t uBound;
  mpz_init (uBound);

  DUZP_t* cc = modularResultant_DUZP (a, b, uBound, 1, 0);
  DUZP_t* c = modularResultant_DUZP (a, b, NULL, 0, 0);

  cerr << "res(a,b) = ";
  printPoly_DUZP (cc, sym);

  cerr << "res1(a,b) = ";
  printPoly_DUZP (c, sym);

  gmp_printf ("uBound = %Zd\n", uBound);

  freePolynomial_DUZP (a);
  freePolynomial_DUZP (b);
  freePolynomial_DUZP (c);
  freePolynomial_DUZP (cc);

}

void test_NewFFTMul (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    duspoly_t* a6 = randomPolynomialInForm_spX (sz1, Pptr);
    duspoly_t* b6 = randomPolynomialInForm_spX (sz2, Pptr);
    duspoly_t* c6 = NULL;

    // fastMulPolynomialInForm_FFT_spX (a6, b6, &c6, Pptr);
    fftMulPolynomialsInForm_spX (a6, b6, &c6, Pptr);

    duspoly_t* cc6 = convertPolynomialFromMontgomery_spX (c6, Pptr);
    duspoly_t* aa6 = convertPolynomialFromMontgomery_spX (a6, Pptr);
    duspoly_t* bb6 = convertPolynomialFromMontgomery_spX (b6, Pptr);

    ZZ_pX C(0), A(0), B(0), D(0);

    convertToZZ_pX (cc6, C, Pptr);
    convertToZZ_pX (aa6, A, Pptr);
    convertToZZ_pX (bb6, B, Pptr);

    mul (D, A, B);

    if (D != C) {
      // std::cerr << "A = " << A << std::endl;
      // std::cerr << "B = " << B << std::endl;

      std::cerr << "A = ";
      printPolynomial_spX (aa6, Pptr);
      std::cerr << "B = ";
      printPolynomial_spX (bb6, Pptr);
      std::cerr << "DUSP: A*B = ";
      printPolynomial_spX (cc6, Pptr);

      // std::cerr << "DUSP: A*B = " << C << std::endl;
      std::cerr << "NTL : A*B = " << D << std::endl;

      fprintf(stderr, "[DUSP] (6) genFastMulPolynomialInForm_FFT... FAIL\n");
      exit (1);
    }

    freePolynomial_spX (&cc6);
    freePolynomial_spX (&aa6);
    freePolynomial_spX (&bb6);
    freePolynomial_spX (&a6);
    freePolynomial_spX (&b6);
    freePolynomial_spX (&c6);

    fprintf(stderr, "[DUSP] genFastMulPolynomialInForm_FFT... PASS\n");
}
void test_memleak (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  ////////////////////////////// SUMMARY
  ////////////////////////////// => initFactors works correctly with NO LEAKS in combination with freeFactors_spX
  ////////////////////////////// => initFacts works correctly with NO LEAKS in combniation with freeFactsll_spX
  ////////////////////////////// => computePolyToPowerPInv_spX(_inp) NO LEAKS
  ////////////////////////////// => Missing freePoly... r0, r1 at the end of algorithm ... FIXED ... Now there is NO LEAKS.
  ////////////////////////////// => Missing freePoly... h2, h21 at the end .... FIXED ... NO LEAKS
  ////////////////////////////// => Updatng plainGCD to GCD in exponentiatePoly, distinctDegreeFact, equalDegreeSplitting...
  ////////////////////////////// => Missing freePoly.... hmx, t, f1, h in distinctDegreeFactor.... FIXED ... NO LEAKS

  // factors_t* facts = initFactors_spX (10);
  // freeFactors_spX (&facts);


  // factsll_t* fact = initFactsll_spX (10);
  // freeFactsll_spX (&fact);

  // duspoly_t* a = randomPolynomialInForm_spX (sz1, Pptr);
  // // duspoly_t* b = NULL;
  // duspoly_t* aa = NULL;

  // computePolyToPowerPInv_spX (a, &aa, Pptr);
  // computePolyToPowerPInv_spX_inp (&a, Pptr);

  // freePolynomial_spX(&a);
  // freePolynomial_spX(&aa);

  // a = randomPolynomialInForm_spX (sz1, Pptr);
  // monicPolynomialInForm_spX_inp (&a, NULL, Pptr);
  // b = randomPolynomialInForm_spX (sz2, Pptr);

  // mulPolynomialsInForm_spX (a, b, &aa, Pptr);
  // squareFreeFactorizationInForm_spX (a, &facts, Pptr);
  // squareFreeFactorizationInFormll_spX (a, &fact, Pptr);
  // distinctDegFactorizationInFormll_spX(a, &fact, Pptr);
  // modFactorizationInFormll_spX (a, &fact, Pptr);
  // // isIrreducibleInForm_spX (a, Pptr);

  // freePolynomial_spX (&a);
  // // freePolynomial_spX (&b);
  // // freePolynomial_spX (&aa);

  // // freeFactors_spX (&facts);

  // freeFactsll_spX (&fact);
  //
}


void test_modularSubres (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    // duspoly_t* a = randomPolynomialInForm_spX (sz1, Pptr);
    // duspoly_t* b = randomPolynomialInForm_spX (sz2, Pptr);

    // cerr << "a=";
    // printPolynomial_spX (a, Pptr);
    // cerr << "b=";
    // printPolynomial_spX (b, Pptr);

    // duspolysA_t* subr = makePolysA_spX(MIN_spX(a->lt, b->lt)+0);
    // polysize_t sz = 0;
    // // subresultantChainInForm_spX (a, b, &subr, &sz, Pptr);
    // _subresultantChainInForm_spX (a, b, &subr, 0, &sz, Pptr);

    // cerr << "sz=" << sz << endl;
    // cerr << "subr->sz=" << subr->size << endl;

    // for (int i = 0; i < subr->size; i++) {
    //   cerr << "subr["<<i<<"]=";
    //   printPolynomial_spX (subr->polys[i], Pptr);
    // }


    // int degA = 9;
    // DUZP_t* a = makePolynomial_DUZP (degA+1);
    // a->lt = degA;
    // mpz_init_set_si (a->coefs[0], 20);
    // mpz_init_set_si (a->coefs[1], 3);
    // mpz_init_set_si (a->coefs[2], 3);
    // mpz_init_set_si (a->coefs[3], 3);
    // mpz_init_set_si (a->coefs[4], 3);
    // mpz_init_set_si (a->coefs[5], 3);
    // mpz_init_set_si (a->coefs[6], 3);
    // mpz_init_set_si (a->coefs[7], 3);
    // mpz_init_set_si (a->coefs[8], 3);
    // mpz_init_set_si (a->coefs[9], 3);

    // int degB = 8;
    // DUZP_t* b = makePolynomial_DUZP (degB+1);
    // b->lt = degB;
    // mpz_init_set_si (b->coefs[0], 20);
    // mpz_init_set_si (b->coefs[1], 1);
    // mpz_init_set_si (b->coefs[2], 1);
    // mpz_init_set_si (b->coefs[3], 1);
    // mpz_init_set_si (b->coefs[4], 1);
    // mpz_init_set_si (b->coefs[5], 1);
    // mpz_init_set_si (b->coefs[6], 1);
    // mpz_init_set_si (b->coefs[7], 1);
    // mpz_init_set_si (b->coefs[8], 1);

  // DUZP_t* a = makePolynomial_DUZP (3);
  // mpz_init_set_si (a->coefs[0], 1000000000);
  // // mpz_init_set_si (a->coefs[1], 9223372036854602819);
  // mpz_init_set_si (a->coefs[1], 20000000);
  // mpz_init_set_si (a->coefs[2], 300000000);
  // a->lt = 2;

  // DUZP_t* b = makePolynomial_DUZP (2);
  // mpz_init_set_si (b->coefs[0], 10000000000);
  // mpz_init_set_si (b->coefs[1], 5000000000000);
  // b->lt = 1;

  int coefBound = 50;
  float sparsity = 0.7;
  int includeNeg = 1;

  DUZP_t* a = buildRandomPoly_DUZP (sz1, coefBound, sparsity, includeNeg);
  DUZP_t* b = buildRandomPoly_DUZP (sz2, coefBound, sparsity, includeNeg);

  AltArrZ_t* aZ = convertToAltArrZ_DUZP (a);
  AltArrZ_t* bZ = convertToAltArrZ_DUZP (b);

  char sym[1] = {'x'};
  cerr << "a = ";
  printPoly_DUZP (a, sym);
  cerr << "b = ";
  printPoly_DUZP (b, sym);

  mpz_t uBound;
  mpz_init (uBound);
  DUZP_t** subr = NULL;
  int chain_size = 0;

  subr = modularSubresultantChain_DUZP (a, b, &chain_size, uBound, 0);

  // for (int i = 0; i<chain_size; i++) {
  //   cerr <<"  subr["<<i<<"]=";
  //   printPoly_DUZP (subr[i], sym);
  // }

  AltArrsZ_t* SC = NULL;
  int len = 0;
  DucosSubresultantChainZ (aZ, bZ, &SC, &len);


  AltArrsZ_t* cur = SC;
  const char* syms[1] = {"x1"};
  for (int i = 0; i < MIN_spX (chain_size, len); i++) {
    AltArrZ_t* tmp = convertToAltArrZ_DUZP (subr[i]);
    if (!isExactlyEqual_AAZ_unpk (tmp, cur->poly)) {
      fprintf (stderr, "NOT EQUAL at idx=%d\nsubr[%d] := ", i, i);
      printPoly_AAZ_unpk (stderr, tmp, syms, 1);
      fprintf (stderr, "\ncur->poly[%d] := ", i);
      printPoly_AAZ_unpk (stderr, cur->poly, syms, 1);
      fprintf (stderr, "\n\n");
    }
    cur = cur->next;
    freePolynomial_AAZ (tmp);
  }

  freeAltArrsZ (SC);

  for (int i = 0; i < chain_size; i++) {
    freePolynomial_DUZP (subr[i]);
  }
  if (subr != NULL) {
    free (subr);
  }


  subr = modularSubresultantAtDeg_DUZP (a, b, 0, &chain_size, uBound, 0, 1);
  for (int i = 0; i<chain_size; i++) {
    cerr <<"  HGCDsubr["<<i<<"]=";
    printPoly_DUZP (subr[i], sym);
  }


  // cerr << "\n\n\n";

  // const Prime_ptr* Ppptr = smallprimefield_get_prime_constants (9223372036854602819);
  // duspoly_t* aP = convertToDUSP_DUZP (a, Ppptr);
  // duspoly_t* bP = convertToDUSP_DUZP (b, Ppptr);

  //   duspolysA_t* subresA1;
  //   polysize_t size1;
  //   subresultantChainInForm_spX (aP, bP, &subresA1, &size1, Ppptr);

  //   for (int i = 0; i < size1; i++) {
  //     cerr << "subresA1["<<i<<"] := ";
  //     printPolynomialOutForm_spX (subresA1->polys[i], Ppptr);
  //   }

  //   cerr << "\n";

  //   duspolysA_t* subresA2;
  //   polysize_t size2;

  //   for (int k = 0; k <= 5; k++) {
  //       hgcdSubResultantInFormA_spX (aP, bP, k, &subresA2, &size2, Ppptr);

  //       for (int i = 0; i < size2; i++) {
  //         cerr << "subresA2["<<k<<"] := ";
  //         printPolynomialOutForm_spX (subresA2->polys[i], Ppptr);
  //       }

  //       freePolysA_spX (subresA2);
  //       size2 = 0;
  //   }




}

// void test_REMREM (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
// {
//   duspoly_t* a4 = randomPolynomialInForm_spX (sz1, Pptr);
//   duspoly_t* b4 = randomPolynomialInForm_spX (sz2, Pptr);

//   duspoly_t* g4 = NULL;
//   elem_t* lcrems = NULL;

//   _plainGCDInForm_spX (a4, b4, &g4, &lcrems, Pptr);

//   duspoly_t* g5 = NULL;
//   elem_t* lcrems5 = NULL;
//   _GCDInForm_spX (a4, b4, &g5, &lcrems5, Pptr);

//   duspoly_t* gg4 = convertPolynomialFromMontgomery_spX (g4, Pptr);
//   duspoly_t* aa4 = convertPolynomialFromMontgomery_spX (a4, Pptr);
//   duspoly_t* bb4 = convertPolynomialFromMontgomery_spX (b4, Pptr);

//   ZZ_pX A;
//   ZZ_pX B;
//   ZZ_pX G;

//   convertToZZ_pX (gg4, G, Pptr);
//   convertToZZ_pX (aa4, A, Pptr);
//   convertToZZ_pX (bb4, B, Pptr);

//   freePolynomial_spX (&gg4);
//   freePolynomial_spX (&aa4);
//   freePolynomial_spX (&bb4);

//   ZZ_pX GG;

//   GCD (GG, A, B);

//   if (1) { // G != GG
//     std::cerr << "A = " << A << std::endl;
//     std::cerr << "B = " << B << std::endl;
//     std::cerr << "DUSP: gcd(A,B) = " << G << std::endl;
//     std::cerr << "NTL:  gcd(A,B) = " << GG << std::endl;

//     int count_neq = 0;
//     if (lcrems != NULL) {
//       for (int i = 0; i < MIN_spX (a4->lt, b4->lt); i++) {
//         cerr << "lc["<<i<<"] = " <<  lcrems[i] << "\t , \t" << lcrems5[i] << endl;
//         if (lcrems[i] != lcrems5[i]) {
//           count_neq++;
//         }
//       }
//       cerr << "\nnumber of non-equals = " <<  count_neq << endl;
//     } else {
//       cerr << "lcrems = NULL = lcrems5" << endl;
//     }

//     // cerr << "\n\n\n";

//     // if (lcrems5 != NULL) {
//     //   for (int i = 0; i < MIN_spX (a4->lt, b4->lt); i++) {
//     //     cerr << "lc5["<<i<<"] = " <<  lcrems5[i] << endl;
//     //   }
//     // } else {
//     //   cerr << "lcrems5 = NULL" << endl;
//     // }
//     // fprintf(stderr, "[DUSP] (4) plainGCDInForm... FAIL\n");
//   }

//   freePolynomial_spX (&a4);
//   freePolynomial_spX (&b4);
//   freePolynomial_spX (&g4);
// // fprintf(stderr, "[DUSP] plainGCDInForm... PASS\n");
// }


void _test_halfsubresZp (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{

    duspoly_t* a3 = NULL;
    duspoly_t* b3 = NULL;

    // 1:
    // a3 :=  1 + 5*x^1 + 9*x^2 + 1*x^3 + 6*x^4 + 9*x^5 + 4*x^6 + 4*x^7 + 2*x^9 in Z_{11}[x];
    // b3 :=  1*x^1 + 3*x^2 + 8*x^3 + 5*x^4 + 5*x^5 + 7*x^6 + 2*x^7 in Z_{11}[x];
    // FIXED :)
    // polysize_t a3_bad_size = 10;
    // polysize_t b3_bad_size = 8;
    // elem_t a3_bad[a3_bad_size] = {1, 5, 9, 1, 6, 9, 4, 4, 0, 2};
    // elem_t b3_bad[b3_bad_size]  = {0, 1, 3, 8, 5, 5, 7, 2};
    // setCoefs_spX (&a3, a3_bad, a3_bad_size, Pptr);
    // setCoefs_spX (&b3, b3_bad, b3_bad_size, Pptr);

    // 2:
    // a3 := 2 + 7*x^1 + 5*x^2 + 6*x^3 + 3*x^4  in Z_{11}[x];
    // b3 :=  9 + 8*x^1 + 10*x^2 + 9*x^3 + 10*x^4 in Z_{11}[x];
    // FIXED :)
    // polysize_t a3_bad_size = 5;
    // polysize_t b3_bad_size = 5;
    // elem_t a3_bad[a3_bad_size] = {2,7,5,6,3};
    // elem_t b3_bad[b3_bad_size]  = {9,8,10,9,10};
    // setCoefs_spX (&a3, a3_bad, a3_bad_size, Pptr);
    // setCoefs_spX (&b3, b3_bad, b3_bad_size, Pptr);



    a3 = randomFullyDensePolynomial_spX (sz1, Pptr);
    b3 = randomFullyDensePolynomial_spX (sz2, Pptr);

    cerr << "a3 := ";
    printPolynomialOutForm_spX (a3, Pptr);
    cerr << "b3 := ";
    printPolynomialOutForm_spX (b3, Pptr);


    duspolysA_t* subresA1;
    polysize_t size1;
    subresultantChainInForm_spX (a3, b3, &subresA1, &size1, Pptr);

    for (int i = 0; i < size1; i++) {
      cerr << "subresA1["<<i<<"] := ";
      printPolynomialOutForm_spX (subresA1->polys[i], Pptr);
    }

    cerr << "\n\n";

    duspolysA_t* subresA2;
    polysize_t size2;

    for (int k = 0; k <= 5; k++) {
        hgcdSubResultantInFormA_spX (a3, b3, k, &subresA2, &size2, Pptr);

        for (int i = 0; i < size2; i++) {
          cerr << "subresA2["<<k<<"] := ";
          printPolynomialOutForm_spX (subresA2->polys[i], Pptr);
        }

        freePolysA_spX (subresA2);
        size2 = 0;
    }

  freePolynomial_spX (&a3);
  freePolynomial_spX (&b3);
  freePolysA_spX (subresA1);
  // freePolysA_spX (subresA2);
}

void test_fftMul (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  int prime_idx;
  unsigned long long min_prime=1945555039024054273UL, max_prime=0;
  for (prime_idx=0; prime_idx < 10000; prime_idx++) {
    if (fourier_primes_u64_table[prime_idx] > max_prime) {
      max_prime = fourier_primes_u64_table[prime_idx];
    }
    if (fourier_primes_u64_table[prime_idx] < min_prime) {
      min_prime = fourier_primes_u64_table[prime_idx];
    }
  }

  fprintf (stderr, "min: %llu  ||| max: %llu\n", min_prime, max_prime);


}

void test_biSubresPr (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  // polysize_t* ad = (polysize_t*) malloc (2*sizeof(polysize_t));
  // polysize_t* bd = (polysize_t*) malloc (2*sizeof(polysize_t));

  // elem_t* a = (elem_t*) calloc (6, sizeof(elem_t));
  // elem_t* b = (elem_t*) calloc (4, sizeof(elem_t));
  // // y*x^2+1
  //   a[0] = smallprimefield_convert_in (1, Pptr);
  //   a[5] = smallprimefield_convert_in (1, Pptr);
  // // b=(y+2)x+1
  //   b[0] = smallprimefield_convert_in (1, Pptr);
  //   b[1] = smallprimefield_convert_in (2, Pptr);
  //   b[2] = smallprimefield_convert_in (0, Pptr);
  //   b[3] = smallprimefield_convert_in (1, Pptr);

  // elem_t* a = (elem_t*) calloc (16, sizeof(elem_t));
  // elem_t* b = (elem_t*) calloc (6, sizeof(elem_t));
  // // // y^3*x^3+1
  //   a[12] = smallprimefield_convert_in (1, Pptr);
  //   a[15] = smallprimefield_convert_in (1, Pptr);
  // // // b=(y+2)x^2+1
  //   b[0] = smallprimefield_convert_in (2, Pptr);
  //   b[2] = smallprimefield_convert_in (2, Pptr);
  //   b[3] = smallprimefield_convert_in (1, Pptr);
  //   b[5] = smallprimefield_convert_in (1, Pptr);
  // ad[0] = 3;   ad[1] = 3;
  // bd[0] = 1;   bd[1] = 2;

  // elem_t* a = (elem_t*) calloc (6*6, sizeof(elem_t));
  // elem_t* b = (elem_t*) calloc (2*5, sizeof(elem_t));
  // // y^5 * x^5
  //   a[6*5+5] = smallprimefield_convert_in (1, Pptr);
  // // b=(y+2)*x^4
  //   b[4] = smallprimefield_convert_in (2, Pptr);
  //   b[5+4] = smallprimefield_convert_in (1, Pptr);
  // ad[0] = 5;   ad[1] = 5;
  // bd[0] = 1;   bd[1] = 4;


  // elem_t* a = (elem_t*) calloc (6*6, sizeof(elem_t));
  // elem_t* b = (elem_t*) calloc (2*5, sizeof(elem_t));
  // // y^5 * x^5 + x^2
  //   a[2] = smallprimefield_convert_in (1, Pptr);
  //   a[6*5+5] = smallprimefield_convert_in (1, Pptr);
  // // b=(y+2)*x^4 -1
  //   b[0] = smallprimefield_convert_in (-1, Pptr);
  //   b[4] = smallprimefield_convert_in (2, Pptr);
  //   b[5+4] = smallprimefield_convert_in (1, Pptr);
  // ad[0] = 5;   ad[1] = 5;
  // bd[0] = 1;   bd[1] = 4;

  int maxDeg[2];
  maxDeg[0] = sz1;  maxDeg[1] = sz2;
  unsigned long ulBound = 5;
  float fsparsity = 0.5;
  int inNeg = 0;
  AltArrZ_t* a = buildRandomZPolyFromMax (2, maxDeg, ulBound, fsparsity, inNeg);
  maxDeg[0] = sz1/2 + 1; maxDeg[1] = sz2/2 + 1;
  AltArrZ_t* b = buildRandomZPolyFromMax (2, maxDeg, ulBound, fsparsity, inNeg);

  const char* sym[2] = {"x", "y"}; // TEST
	fprintf (stderr, "a := ");
	printPoly_AAZ (stderr, a, sym, a->nvar);
	fprintf (stderr, "\nb := ");
	printPoly_AAZ (stderr, b, sym, b->nvar);
	fprintf (stderr, "\n");

  polysize_t ad[2];
  elem_t* aZ = convertFromAltArrZ_DBSP (a, ad, Pptr);
  // std::cerr << "apd = " << apd[0] << " , " << apd[1] << std::endl;
  std::cerr << "converted(a) := " << std::endl;
	printDBSP_AAZ (aZ, ad, Pptr);
  polysize_t bd[2];
  elem_t* bZ = convertFromAltArrZ_DBSP (b, bd, Pptr);
  // std::cerr << "bpd = " << bd[0] << " , " << bpd[1] << std::endl;
  std::cerr << "converted(b) := " << std::endl;
  printDBSP_AAZ (bZ, bd, Pptr);

  // for (polysize_t i = 0; i < (ad[0]+1)*(ad[1]+1); i++)
  //   std::cerr << "a["<<i<<"]="<< smallprimefield_convert_out (a[i], Pptr) << std::endl;
  // std::cerr << std::endl;
  // for (polysize_t i = 0; i < (bd[0]+1)*(bd[1]+1); i++)
  //   std::cerr << "b["<<i<<"]="<< smallprimefield_convert_out (b[i], Pptr) << std::endl;

  // exit (1);

  biSubresPr_t* subres = NULL;
  biSylvSubResultantInForm_spX (aZ, ad, bZ, bd, &subres, Pptr);

  std::cerr << "\n\n";
  for (polysize_t i = 0; i < subres->n; i++) {
    std::cerr << "S["<<i<<"]->deg\t"<< subres->deg[i][0] << "  and  " << subres->deg[i][1] << std::endl;
    for (polysize_t j = 0; j < subres->size[i]; j++) {
      std::cerr << "S["<<i<<"]["<<j<<"]="<<smallprimefield_convert_out (subres->coefs[i][j], Pptr) <<std::endl;
    }
  }

  for (polysize_t i = 0; i < subres->n; i++) {
    AltArrZ_t* t = convertToAltArrZ_DBSP (subres->coefs[i], subres->deg[i], Pptr);
    std::cerr << "S["<<i<<"] := " << std::endl;
    printPoly_AAZ (stderr, t, sym, t->nvar);
    fprintf (stderr, "\n");
    freePolynomial_AAZ (t);
  }

  // free (a);
  // free (b);
  // free (ad);
  // free (bd);
  // freeBiSubresultantInForm_spX (subres);
}

void test_conertDBZP (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{

  int maxDeg[2];
  maxDeg[0] = sz1;  maxDeg[1] = sz2;
  unsigned long ulBound = 4;
  float fsparsity = 0.5;
  int inNeg = 0;
  AltArrZ_t* a = buildRandomZPolyFromMax (2, maxDeg, ulBound, fsparsity, inNeg);
  maxDeg[0] = sz1/2 + 1; maxDeg[1] = sz2/2 + 1;
  AltArrZ_t* b = buildRandomZPolyFromMax (2, maxDeg, ulBound, fsparsity, inNeg);

  const char* sym[2] = {"x", "y"}; // TEST
	fprintf (stderr, "a := ");
	printPoly_AAZ (stderr, a, sym, a->nvar);
	fprintf (stderr, "\nb := ");
	printPoly_AAZ (stderr, b, sym, b->nvar);
	fprintf (stderr, "\n");

  polysize_t apd[2];
  mpz_t* mpz_a = convertFromAltArrZ_DBZP (a, apd);
  std::cerr << "apd = " << apd[0] << " , " << apd[1] << std::endl;
  std::cerr << "converted(a) := " << std::endl;
	// for (int i = 0; i < (apd[0]+1)*(apd[1]+1); i++) {
	// 	gmp_fprintf (stderr, "ret_a[%d] = %Zd\n", i, mpz_a[i]);
	// }
  printDBZP_AAZ (mpz_a, apd);
  AltArrZ_t* aa = convertToAltArrZ_DBZP (mpz_a, apd);
	fprintf (stderr, "aa := ");
	printPoly_AAZ (stderr, aa, sym, aa->nvar);
	fprintf (stderr, "\n");

  polysize_t bpd[2];
  mpz_t* mpz_b = convertFromAltArrZ_DBZP (b, bpd);
  std::cerr << "bpd = " << bpd[0] << " , " << bpd[1] << std::endl;
  std::cerr << "converted(b) := " << std::endl;
  printDBZP_AAZ (mpz_b, bpd);
	// for (int i = 0; i < (bpd[0]+1)*(bpd[1]+1); i++) {
	// 	gmp_fprintf (stderr, "ret_b[%d] = %Zd\n", i, mpz_b[i]);
	// }
  AltArrZ_t* bb = convertToAltArrZ_DBZP (mpz_b, bpd);
	fprintf (stderr, "bb := ");
	printPoly_AAZ (stderr, bb, sym, bb->nvar);
	fprintf (stderr, "\n");


  freePolynomial_AAZ (a);
  freePolynomial_AAZ (b);
  freePolynomial_AAZ (aa);
  freePolynomial_AAZ (bb);
}



void test_biSubresZ (polysize_t sz1, polysize_t sz2, polysize_t flag, unsigned long ulBound = 10, int nthreads = -1)
{
  int maxDeg[2];
  maxDeg[0] = sz1;  maxDeg[1] = sz2;
  float fsparsity = 0.0;
  int inNeg = 0;
  AltArrZ_t* a = buildRandomZPolyFromMax (2, maxDeg, ulBound, fsparsity, inNeg);
  maxDeg[0] = sz1/2 + 1; maxDeg[1] = sz2/2 + 1;
  AltArrZ_t* b = buildRandomZPolyFromMax (2, maxDeg, ulBound, fsparsity, inNeg);

  const char* sym[2] = {"x", "y"}; // TEST
	fprintf (stderr, "a := ");
	printPoly_AAZ (stderr, a, sym, a->nvar);
	fprintf (stderr, "\nb := ");
	printPoly_AAZ (stderr, b, sym, b->nvar);
	fprintf (stderr, "\n");

  // polysize_t ad[2];
  // elem_t* aZ = convertFromAltArrZ_DBSP (a, ad, Pptr);
  // // std::cerr << "apd = " << apd[0] << " , " << apd[1] << std::endl;
  // std::cerr << "converted(a) := " << std::endl;
	// printDBSP_AAZ (aZ, ad, Pptr);
  // polysize_t bd[2];
  // elem_t* bZ = convertFromAltArrZ_DBSP (b, bd, Pptr);
  // // std::cerr << "bpd = " << bd[0] << " , " << bpd[1] << std::endl;
  // std::cerr << "converted(b) := " << std::endl;
  // printDBSP_AAZ (bZ, bd, Pptr);
  // biSubresPr_t* subres = NULL;
  // biSylvSubResultantInForm_spX (aZ, ad, bZ, bd, &subres, Pptr);
  // for (polysize_t i = 0; i < subres->n; i++) {
  //   AltArrZ_t* t = convertToAltArrZ_DBSP (subres->coefs[i], subres->deg[i], Pptr);
  //   std::cerr << "S["<<i<<"] := " << std::endl;
  //   printPoly_AAZ (stderr, t, sym, t->nvar);
  //   fprintf (stderr, "\n");
  //   freePolynomial_AAZ (t);
  // }


  unsigned long long start;
  float elaspedSerial, elaspedParallel;

  AltArrsZ_t* Subres = NULL;
  int chain_size;
  int num_primes = 0;
  _startTimer(&start);
  biModularSubresultantChainZ (a, b, &Subres, &chain_size, &num_primes);
  // DUSP::biModularSubresultantChainZ_serial (a, b, &Subres, &chain_size, &num_primes);
  _stopTimer(&start, &elaspedSerial);
  fprintf(stderr, "Serial   Time: %.5f,        NumPrimes: %d\n", elaspedSerial, num_primes);

#if !(defined(SERIAL) && SERIAL)
  ExecutorThreadPool::getThreadPool(); //wake up the pool before timings
  AltArrsZ_t* Subres2 = NULL;
  int chain_size2;
  int num_primes2 = 0;
  _startTimer(&start);
  DUSP::biModularSubresultantChainZ (a, b, &Subres2, &chain_size2, &num_primes2);
  _stopTimer(&start, &elaspedParallel);
  AltArrsZ_t* cur = Subres;
  AltArrsZ_t* cur2 = Subres2;
  int idx = 0;
  while (cur != NULL && cur2 != NULL) {
    if (!isExactlyEqual_AAZ(cur->poly, cur2->poly)) {
      cerr << "Bivar Modular Subresultant Test FAILED!" << endl;
      exit(1);
    }
    cur = cur->next;
    cur2 = cur2->next;
    idx++;
  }
  cerr << "chain_size  := " << chain_size << std::endl;
  cerr << "chain_size2 := " << chain_size2 << std::endl;
  fprintf(stderr, "Serial   Time: %.5f,        NumPrimes: %d\n", elaspedSerial, num_primes);
  fprintf(stderr, "Parallel Time: %.5f,        NumPrimes: %d\n", elaspedParallel, num_primes2);
#endif
  cerr << "Bivar Modular Subresultant Test PASSED!" << endl;

  // AltArrsZ_t* cur = Subres;
  // int idx = 0;
  // while (cur != NULL) {
  //   cerr << "Subres["<< idx << "] := ";
  // 	printPoly_AAZ (stderr, cur->poly, sym, 2);
  //   cerr << "; \n";
  //   cur = cur->next;
  //   idx++;
  // }
  // cerr << "chain_size := " << chain_size << std::endl;

  freePolynomial_AAZ (a);
  freePolynomial_AAZ (b);
  // free (aZ);
  // free (bZ);
  // freeBiSubresultantInForm_spX (subres);
}


void test_SRC (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
 
  // elem_t a_el[2] = {0, 1291970143101910932}; 
  // elem_t b_el[1] = {645985071550955466};
  // duspoly_t* a = NULL;
  // setCoefs_spX (&a, a_el, 2, Pptr);
  // duspoly_t* b = NULL;
  // setCoefs_spX (&b, b_el, 1, Pptr);

  duspoly_t *a = randomPolynomialInForm_spX (sz1+1, Pptr);
  duspoly_t *b = randomPolynomialInForm_spX (sz2+1, Pptr);

    duspolysA_t* SC = NULL;
    polysize_t sz;
    subresultantChainInForm_spX (a, b, &SC, &sz, Pptr);
          // std::cerr << "a := "; printPolynomialOutForm_spX (a, Pptr);
          // std::cerr << "b := "; printPolynomialOutForm_spX (b, Pptr);
          // for (int i = 0; i < sz; i++) {
          //   std::cerr << "S["<<i<<"] := "; printPolynomialOutForm_spX (SC->polys[i], Pptr);
          // }
    
    for (int i = 0; i < MIN_spX(sz1, sz2)-1; ++i) {
        duspolysA_t* ssrc = NULL;
        // std::cerr << "i = " << i << std::endl;
        plainSpeculativesubresultantChainInForm_spX (a, b, i, &ssrc, Pptr);
        if (!isEqual_spX(SC->polys[i+1], ssrc->polys[0]) || !isEqual_spX(SC->polys[i], ssrc->polys[1])) {
          std::cerr << "a := "; printPolynomialOutForm_spX (a, Pptr);
          std::cerr << "b := "; printPolynomialOutForm_spX (b, Pptr);
          for (int i = 0; i < sz; i++) {
            std::cerr << "S["<<i<<"] := "; printPolynomialOutForm_spX (SC->polys[i], Pptr);
          }
          if (ssrc->size == 2) {
            std::cerr << "SS["<<i+1<<"] := "; printPolynomialOutForm_spX (ssrc->polys[0], Pptr);
            std::cerr << "SS["<<i<<"] := ";   printPolynomialOutForm_spX (ssrc->polys[1], Pptr);
          } else {
            std::cerr << "Spec. Subres of a, b for k = " << i << " is a list of len < 2!" << std::endl;
          }
          std::cerr << "[DUSP] subresultantChainInForm_spX... FAIL" << std::endl;
          std::cerr << "[DUSP] plainSpeculativesubresultantChainInForm_spX... FAIL" << std::endl;
          exit (1);
        }
        if (ssrc != NULL)
          freePolysA_spX (ssrc);
    }

  std::cerr << "[DUSP] subresultantChainInForm_spX... PASS" << std::endl;
  std::cerr << "[DUSP] plainSpeculativesubresultantChainInForm_spX... PASS" << std::endl;

  freePolynomial_spX(&a);
  freePolynomial_spX(&b);
  freePolysA_spX (SC);
}

void test_RegGCDSRC (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  duspoly_t *a = randomPolynomialInForm_spX (sz1+1, Pptr);
  duspoly_t *b = randomPolynomialInForm_spX (sz2+1, Pptr);

    duspolysA_t* SC = NULL;
    polysize_t sz;
    subresultantChainInForm_spX (a, b, &SC, &sz, Pptr);
    std::cerr << "a := "; printPolynomialOutForm_spX (a, Pptr);
    std::cerr << "b := "; printPolynomialOutForm_spX (b, Pptr);
    for (int i = 0; i < sz; i++) {
      std::cerr << "S["<<i<<"] := "; printPolynomialOutForm_spX (SC->polys[i], Pptr);
    }

    std::cerr << std::endl;

    duspolysA_t *results = NULL;
    for (int i = 0; i < MIN_spX(sz1, sz2)-1; ++i) {
        // std::cerr << "i = " << i << std::endl;
        hgcdSubResultantInFormA_spX (a, b, i, &results, &sz, Pptr);
        std::cerr << "[old] S["<<i+1<< "] := "; printPolynomialOutForm_spX (results->polys[0], Pptr);
        std::cerr << "[old] S["<<i<< "] := "; printPolynomialOutForm_spX (results->polys[1], Pptr);
        freePolysA_spX (results); results = NULL; sz = 0;
    }

    std::cerr << std::endl;

    // duspolysA_t *results = NULL;
    specA_spX_t *A = createSpecA_spX (sz1);
    // cerr << "A: size := " << A->size << "  alloc := " << A->alloc << endl;
    specQ_spX_t *Q = createSpecQ_spX (sz1);
    // cerr << "Q: size := " << Q->size << "  alloc := " << Q->alloc << endl;
    specR_spX_t *R = createSpecR_spX (sz1);
    // cerr << "R: size := " << R->size << "  alloc := " << R->alloc << endl;

    for (int i = 0; i < MIN_spX(sz1, sz2)-1; ++i) {
        polysize_t degs[2] = {0, 0};
        // std::cerr << "i = " << i << std::endl;
        regularGCDSpecSRC_spX (a, b, i, &results, degs, 1, A, Q, R, Pptr);
        std::cerr << "degs[0] = " << degs[0] << " degs[1] = " << degs[1] << std::endl;
        std::cerr << "[new] SS["<<i+1<< "] := "; printPolynomialOutForm_spX (results->polys[0], Pptr);
        std::cerr << "[new] SS["<<i<< "] := "; printPolynomialOutForm_spX (results->polys[1], Pptr);
        freePolysA_spX (results); results = NULL;
    }

  freePolynomial_spX(&a);
  freePolynomial_spX(&b);
  freePolysA_spX (SC);
}



void _test_subresultant (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
   // 4*x^10 + 9223372036854602799*x^9 + 131262*x^8 + 9223372036854077891*x^7 + 1722761168*x^6 + 9223372031688156479*x^5 + 11303897220108*x^4 + 9223349437670294127*x^3 + 37080740769906054*x^2 + 9186302594815864203*x + 2532334113279498244
    elem_t* ae = (elem_t*) malloc (11*sizeof(elem_t));
    ae[10] = 4;
    ae[9] = 9223372036854602799;
    ae[8] = 131262;
    ae[7] = 9223372036854077891;
    ae[6] = 1722761168;
    ae[5] = 9223372031688156479;
    ae[4] = 11303897220108;
    ae[3] = 9223349437670294127;
    ae[2] = 37080740769906054;
    ae[1] = 9186302594815864203;
    ae[0] = 2532334113279498244;

    // 2*x^8 + 9223372036854602811*x^7 + 52500*x^6 + 9223372036854445347*x^5 + 516744362*x^4 + 9223372035821376539*x^3 + 2260297252514*x^2 + 9223369777073937201*x + 3707170181956171
    elem_t* be = (elem_t*) malloc (9*sizeof(elem_t));
    be[8] = 2;
    be[7] = 9223372036854602811;
    be[6] = 52500;
    be[5] = 9223372036854445347;
    be[4] = 516744362;
    be[3] = 9223372035821376539;
    be[2] = 2260297252514;
    be[1] = 9223369777073937201;
    be[0] = 3707170181956171;

    duspoly_t* a = NULL;
    setCoefs_spX (&a, ae, 11, Pptr);
    duspoly_t* b = NULL;
    setCoefs_spX (&b, be, 9, Pptr);

    cerr << "a := ";
    printPolynomialOutForm_spX (a, Pptr);
    cerr << "b := ";
    printPolynomialOutForm_spX (b, Pptr);


    duspolysA_t* SC = NULL;
    polysize_t sz;
    subresultantChainInForm_spX (a, b, &SC, &sz, Pptr);

    for (polysize_t i = 0; i < sz; i++) {
      cerr << "SC["<<i<<"] := ";
      printPolynomialOutForm_spX (SC->polys[i], Pptr);
    }
    
}

#include "../../../include/FFT/src/fft_iter_main.h"
void _temp_FFT_Mul (duspoly_t* a, duspoly_t* b, duspoly_t** c, Prime_ptr* Pptr)
{
    polysize_t max_deg = MAX_spX (degPolynomial_spX (a), degPolynomial_spX (b));

    polysize_t r = Log2_spX (max_deg)+2;
    polysize_t n = 1<<r;
    polysize_t n2 = n<<1;
    // std::cout << "r = " << r << "\tn = " << n << "\tn2 = " << n2 << std::endl;
    elem_t* f = (elem_t*) calloc (n, sizeof (elem_t ));
    elem_t* g = (elem_t*) calloc (n, sizeof (elem_t));
    for (long i = 0; i <= degPolynomial_spX (a); ++i) {
      f[i] = (long) a->elems[i];
    }
    for (long i = 0; i <= degPolynomial_spX (b); ++i) {
      g[i] = (long) b->elems[i];
    }

    elem_t *SB1, *Omega, *OmegaInv, *pwm_fg;
    Omega = (elem_t*) calloc (n2<<1, sizeof (elem_t));
    OmegaInv = (elem_t*) calloc (n2<<1, sizeof (elem_t));
    RootsTable2_spf(n2, r+1, Omega, OmegaInv, Pptr);

    int H = 1024, j = 0;

    elem_t* _omega = Omega;
    elem_t* _omegaInv = OmegaInv;
    Omega += n2;
    OmegaInv += n2;
    SB1 = (elem_t*) calloc (n, sizeof (elem_t));

    DFT_eff(n, r, f, Omega, Pptr, H, SB1);
    DFT_eff(n, r, g, Omega, Pptr, H, SB1);
    // PBPAS::DFT_eff(n, r, f, Omega, pPtr, H, NULL, SB1, 1);
    // PBPAS::DFT_eff(n, r, g, Omega, pPtr, H, NULL, SB1, 1);

    pwm_fg = (elem_t*) calloc (n, sizeof (elem_t));
    for (j = 0; j < n; j++) {
      pwm_fg[j] = smallprimefield_mul (f[j], g[j], Pptr);
    }

    free (f);
    free (g);
    elem_t invn = smallprimefield_inv(smallprimefield_convert_in(n, Pptr), Pptr);
    InvDFT_eff(n, r, pwm_fg, OmegaInv, Pptr, H, SB1, invn);
    // PBPAS::InvDFT_eff(n, r, pwm_fg, OmegaInv, pPtr, H, NULL, SB1, invn, 1);

    setCoefsInForm_spX (c, pwm_fg, n);
    normalize_spX (c);

    free(_omega);
    free(_omegaInv);
    free (SB1);
}

void test_fft_iter_main(int sz1, int sz2, Prime_ptr* Pptr, polysize_t flag) {
  duspoly_t* a6 = randomPolynomialInForm_spX (sz1, Pptr);
  duspoly_t* b6 = randomPolynomialInForm_spX (sz2, Pptr);
  duspoly_t* c6 = NULL;

// fastMulPolynomialInForm_FFT_spX (a6, b6, &c6, Pptr);
// fftMulPolynomialsInForm_spX (a6, b6, &c6, Pptr);
  _temp_FFT_Mul(a6, b6, &c6, Pptr);

  duspoly_t* cc6 = convertPolynomialFromMontgomery_spX (c6, Pptr);
  duspoly_t* aa6 = convertPolynomialFromMontgomery_spX (a6, Pptr);
  duspoly_t* bb6 = convertPolynomialFromMontgomery_spX (b6, Pptr);

  NTL::ZZ_pX C(0), A(0), B(0), D(0);

  convertToZZ_pX (cc6, C, Pptr);
  convertToZZ_pX (aa6, A, Pptr);
  convertToZZ_pX (bb6, B, Pptr);

  NTL::mul (D, A, B);

  if (D != C) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;

    std::cerr << "A = ";
    printPolynomial_spX (aa6, Pptr);
    std::cerr << "B = ";
    printPolynomial_spX (bb6, Pptr);
    std::cerr << "DUSP: A*B = ";
    printPolynomial_spX (cc6, Pptr);

    std::cerr << "DUSP: A*B = " << C << std::endl;
        std::cerr << "NTL : A*B = " << D << std::endl;

    fprintf(stderr, "[DUSP] test_fft_iter_main... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&cc6);
  freePolynomial_spX (&aa6);
  freePolynomial_spX (&bb6);
  freePolynomial_spX (&a6);
  freePolynomial_spX (&b6);
  freePolynomial_spX (&c6);

  fprintf(stderr, "[DUSP] test_fft_iter_main... PASS\n");
}


void test_newFastDiv (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
  duspoly_t* a8 = randomPolynomialInForm_spX (sz1*2, Pptr);
  duspoly_t* b8 = randomPolynomialInForm_spX (sz2, Pptr);

  duspoly_t* q8 = NULL;
  duspoly_t* r8 = NULL;

  fastDivPolynomialInForm_wPSInv_4stepFFT_spX (a8, b8, &r8, &q8, Pptr);

  // r8 = NULL; q8 = NULL;
  // fastDivPolynomialInForm_wPSInv_FFT_spX (a8, b8, &r8, &q8, Pptr);
  // exit(1);

  duspoly_t* qq8 = convertPolynomialFromMontgomery_spX (q8, Pptr);
  duspoly_t* rr8 = convertPolynomialFromMontgomery_spX (r8, Pptr);
  duspoly_t* aa8 = convertPolynomialFromMontgomery_spX (a8, Pptr);
  duspoly_t* bb8 = convertPolynomialFromMontgomery_spX (b8, Pptr);

  ZZ_pX A;
  ZZ_pX B;
  ZZ_pX R;
  ZZ_pX Q;

  convertToZZ_pX (qq8, Q, Pptr);
  convertToZZ_pX (rr8, R, Pptr);
  convertToZZ_pX (aa8, A, Pptr);
  convertToZZ_pX (bb8, B, Pptr);

  freePolynomial_spX (&qq8);
  freePolynomial_spX (&rr8);
  freePolynomial_spX (&aa8);
  freePolynomial_spX (&bb8);

  ZZ_pX RR;
  ZZ_pX QQ;

  DivRem (QQ, RR, A, B);

  if (QQ != Q || RR != R) {
    std::cerr << "A = " << A << std::endl;
    std::cerr << "B = " << B << std::endl;
    std::cerr << "DUSP: R = " << R << std::endl;
    std::cerr << "NTL:  R = " << RR << std::endl;
    std::cerr << "DUSP: Q = " << Q << std::endl;
    std::cerr << "NTL : Q = " << QQ << std::endl;


    fprintf(stderr, "[DUSP] fastDivPolynomialInForm_wPSInv_4stepFFT... FAIL\n");
    exit (1);
  }

  freePolynomial_spX (&a8);
  freePolynomial_spX (&b8);
  freePolynomial_spX (&q8);
  freePolynomial_spX (&r8);

  fprintf(stderr, "[DUSP] fastDivPolynomialInForm_4stepFFT... PASS\n");
}

void test_newDUZPSRC (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    int coefBoundBits = 32;
    float sparsity = 0;
    int includeNeg = 0;
    DUZP_t *a = buildRandomPoly_DUZP(sz1, coefBoundBits, sparsity, includeNeg);
    DUZP_t *b = buildRandomPoly_DUZP(sz2, coefBoundBits, sparsity, includeNeg);

    char sym[1] = {'x'};
    fprintf(stderr, "a := \n");
    printPoly_DUZP(a, sym);
    fprintf(stderr, "b := \n");
    printPoly_DUZP(b, sym);

    int csize; 
    mpz_t uBound;
    mpz_init(uBound);
    DUZP_t ** ss = modularSubresultantChain_DUZP (a, b, &csize, uBound, 0);

    for (int i = 0; i < csize; ++i) {
      fprintf(stderr, "ss[%d] := ", i);
      printPoly_DUZP(ss[i], sym);
    }

    polysize_t degs[2] = {0, 0};
    int k = 0;
    int csizek, csizek2, info_size = 0; 
    specAQRArray_spX_t *uspecInfo = NULL;
    DUZP_t ** sss1;
    DUZP_t ** sss2; 

    unsigned long long int start;
    float elapsed; 

    for (int k = 0; k < b->lt; ++k) {
      _startTimer(&start);
      sss1 = modularSubresultantAtDeg_DUZP (a, b, k, &csizek, uBound, 0, 1);
      _stopTimer(&start, &elapsed);
      std::cerr << "[old] [k="<<k<<"]" << "Time: " << elapsed << endl;
      // for (int i = 0; i < csizek; ++i) {
      //   fprintf(stderr, "[k=%d] sss1[%d] := ", k, i);
      //   printPoly_DUZP(sss1[i], sym);
      //   freePolynomial_DUZP (sss1[i]);
      // }
      // cerr << endl;
      _startTimer(&start);
      sss2 = regularGCDUnivariateSpeculativeSRC_DUZP (a, b, k, &csizek2, degs, uBound, &uspecInfo, &info_size, 1);
      _stopTimer(&start, &elapsed);
      std::cerr << "[new] [k="<<k<<"]" << "Time: " << elapsed << endl;
      // for (int i = 0; i < csizek2; ++i) {
      //   fprintf(stderr, "[k=%d] [deg=%ld] sss2[%d] := ", k, degs[i], i);
      //   printPoly_DUZP(sss2[i], sym);
      //   freePolynomial_DUZP (sss2[i]);
      // }
      // cerr << endl;
      csizek = 0; csizek2 = 0;
    }
    


    freePolynomial_DUZP (a);
    freePolynomial_DUZP (b);
    for (int i = 0; i < csize; ++i) {
      freePolynomial_DUZP (ss[i]);
    }

 }

void test_newDBZPSRC (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    int coefBoundBits = 5;
    float sparsity = 0;
    int includeNeg = 0;
    int maxDegs[2] = {(int) sz1, (int) sz2};
    int maxDegs_b[2] = {(int) sz1-1, (int) sz2-1};
    AltArrZ_t *a = buildRandomZPolyFromMax (2, maxDegs, coefBoundBits, sparsity, includeNeg);
    AltArrZ_t *b = buildRandomZPolyFromMax (2, maxDegs_b, coefBoundBits, sparsity, includeNeg);

    const char* sym[10] = {"x", "y", "z", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST
    fprintf (stderr, "a := ");
    printPoly_AAZ (stderr, a, sym, a->nvar);
    fprintf (stderr, "\nb := ");
    printPoly_AAZ (stderr, b, sym, b->nvar);

    int csize = 0, no_prime = 0; 
    AltArrsZ_t *Subres = NULL;
    int isGoodOmega;
    // isGoodOmega = biModularFFTSubresultantChainZ (a, b, &Subres, &csize, &no_prime);

    AltArrsZ_t *cur = Subres;
    // for (int i = 0; i < csize; ++i) {
    //   fprintf(stderr, "\nS[%d] := ", i);
    //   printPoly_AAZ (stderr, cur->poly, sym, a->nvar);
    //   cur = cur->next;
    // }
    // std::cerr << std::endl;

    AltArrsZ_t *sss1=NULL, *sss2=NULL;
    int csize1=0, csize2=0;
    unsigned long long int start; 
    float elapsed;
    specAQRArray_spXY_t *bspecInfo = NULL;
    polysize_t mdegs[2];
    int info_size = 0;
    int results_mode = 0;
    for (int k = 0; k < mainDegree_AAZ (a)-1; ++k) {
        _startTimer(&start);
        isGoodOmega = hgcdBiModularFFTSubresultantChainZ (a, b, k, &sss1, &csize1);
        _stopTimer(&start, &elapsed);
        std::cerr << "[old] Time[k="<< k <<"]: " << elapsed << std::endl; 
        // cur = sss1;
        // for (int i = 0; i < csize1; ++i) {
        //   fprintf(stderr, "\n[old] S[%d][%d] := ", k, i);
        //   printPoly_AAZ (stderr, cur->poly, sym, a->nvar);
        //   cur = cur->next;
        // }
        std::cerr << std::endl;
        _startTimer(&start);
        isGoodOmega = regularGCDBiModularFFTSubresultantChainZ (a, b, k, &sss2, &csize2, 
                        mdegs, &bspecInfo, &info_size, results_mode);
        _stopTimer(&start, &elapsed);
        std::cerr << "[new] Time[k="<< k <<"]: " << elapsed << std::endl; 

        // print bspecInfo: 
        std::cerr << "\n\nPrint bspecInfo for k="<< k << " :"<< std::endl;
        std::cerr << "size: " << bspecInfo->size << std::endl;
        std::cerr << "W array: " << std::endl;
        for (int i = 0; i < bspecInfo->size; ++i) {
          std::cerr << bspecInfo->bspecArray[i]->W << std::endl;
        }

        // cur = sss2; 
        // for (int i = 0; i < csize1; ++i) {
          // fprintf(stderr, "\n[new] S[%d][%d] := ", k, i);
          // printPoly_AAZ (stderr, cur->poly, sym, a->nvar);
          // cur = cur->next;
        // }
        std::cerr << std::endl;
    }

    freePolynomial_AAZ (a);
    freePolynomial_AAZ (b);
    cur = Subres;
    for (int i = 0; i < csize; ++i) {
        if (cur->poly != NULL) 
          freePolynomial_AAZ (cur->poly);
        cur = cur->next;
    }
 }

void test_smqpSRC1 (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag) {
  unsigned long long start;
  float elasped;

  SMQP a = SMQP ("x^10 + x^4 + x^2 + 1");
  SMQP b = SMQP ("x^8 + x^5 + 1");

  std::cerr << "a := " << a << std::endl;
  std::cerr << "b := " << b << std::endl;
  int idx = 0;

    // std::vector<SMQP> S1 = a.subresultantChainAtIdx(b, Symbol("x"), idx);
    // for (int i=0; i < S1.size(); ++i) {
    //   std::cerr << "S1[" << i <<"] := " << S1.at(i) << std::endl;
    // }

    specSRC_AAZ *lazyInfo = NULL;
    _startTimer(&start);
    std::vector<SMQP> S2 = a.subresultantAtIdx (b, Symbol("x"), 0, &lazyInfo);
    _stopTimer(&start, &elasped);
    std::cerr << "[0-call] Time: " << elasped << std::endl;
    for (int i=0; i < S2.size(); ++i) {
      std::cerr << "S2[ idx=0 ][ i= " << i <<" ] := " << S2.at(i) << std::endl;
    }

    for (idx = 1; idx < b.degree()-1; idx++) {
        _startTimer(&start);
        S2 = a.subresultantAtIdx (b, Symbol("x"), idx, &lazyInfo);
        _stopTimer(&start, &elasped);
        std::cerr << "[i-call] Time: " << elasped << std::endl;
        for (int i=0; i < S2.size(); ++i) {
          std::cerr << "S2[ idx=" << idx << " ][ i= "<< i <<" ] := " << S2.at(i) << std::endl;
        }
    }

    cerr << endl;
    Integer degs1, degs2;
    degree_t degs[2];
    std::vector<SMQP> S3 = a.subresultantInitialAtIdx (b, Symbol("x"), idx, 
                                                      degs1, degs2, &lazyInfo);
    degs[0] = degs1.get_d(); degs[1] = degs2.get_d();
    for (int i=0; i < S3.size(); ++i) {
      std::cerr << "[deg= "<< degs[i] << "] S3[" << i <<"] := " << S3.at(i) << std::endl;
    }

}

void test_smqpSRC2 (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag) {
  // SMQP a = SMQP ("(y^2 + 6)*(x_1 - 1) - y*(x_1^2 + 1)");
  // SMQP b = SMQP ("(x_1^2 + 6)*(y - 1) - x_1*(y^2 + 1)");
  SMQP a = SMQP ("4*x_1^10 - 20*x_1^9 + 20*x_1^8*x_2^2 + 42*x_1^8 - 80*x_1^7*x_2^2 - 48*x_1^7 + 40*x_1^6*x_2^4 + 136*x_1^6*x_2^2 + 32*x_1^6 - 120*x_1^5*x_2^4 - 128*x_1^5*x_2^2 - 12*x_1^5 + 40*x_1^4*x_2^6 + 156*x_1^4*x_2^4 + 72*x_1^4*x_2^2 - 80*x_1^3*x_2^6 - 112*x_1^3*x_2^4 - 24*x_1^3*x_2^2 + 4*x_1^3 + 20*x_1^2*x_2^8 + 72*x_1^2*x_2^6 + 48*x_1^2*x_2^4 - 6*x_1^2 - 20*x_1*x_2^8 - 32*x_1*x_2^6 - 12*x_1*x_2^4 + 4*x_1*x_2^2 + 4*x_1 + 4*x_2^10 + 10*x_2^8 + 8*x_2^6 - 2*x_2^2 - 1");
  SMQP b = SMQP ("2*x_1^8 - 8*x_1^7 + 8*x_1^6*x_2^2 + 12*x_1^6 - 24*x_1^5*x_2^2 - 8*x_1^5 + 12*x_1^4*x_2^4 + 28*x_1^4*x_2^2 + 2*x_1^4 - 24*x_1^3*x_2^4 - 16*x_1^3*x_2^2 + 8*x_1^2*x_2^6 + 20*x_1^2*x_2^4 + 4*x_1^2*x_2^2 + 2*x_1^2 - 8*x_1*x_2^6 - 8*x_1*x_2^4 - 2*x_1 + 2*x_2^8 + 4*x_2^6 + 2*x_2^4 + 2*x_2^2 + 1");

  std::cerr << "a := " << a << std::endl;
  std::cerr << "b := " << b << std::endl;
  int idx = 0;

  // for (int idx = 0;  idx < b.degree(); idx++) {
    // std::vector<SMQP> S1 = a.subresultantChainAtIdx(b, Symbol("x_1"), idx);
    // for (int i=0; i < S1.size(); ++i) {
    //   std::cerr << "S1[" << i <<"] := " << S1.at(i) << std::endl;
    // }

    specSRC_AAZ *lazyInfo = NULL;
    std::vector<SMQP> S2 = a.subresultantAtIdx (b, Symbol("x_1"), 0, &lazyInfo);
    for (int i=0; i < S2.size(); ++i) {
      std::cerr << "S2[0][" << i <<"] := " << S2.at(i) << std::endl;
    }

    S2 = a.subresultantAtIdx (b, Symbol("x_1"), 1, &lazyInfo);
    for (int i=0; i < S2.size(); ++i) {
      std::cerr << "S2[1][" << i <<"] := " << S2.at(i) << std::endl;
    }

    S2 = a.subresultantAtIdx (b, Symbol("x_1"), 2, &lazyInfo);
    for (int i=0; i < S2.size(); ++i) {
      std::cerr << "S2[2][" << i <<"] := " << S2.at(i) << std::endl;
    }

    S2 = a.subresultantAtIdx (b, Symbol("x_1"), 4, &lazyInfo);
    for (int i=0; i < S2.size(); ++i) {
      std::cerr << "S2[4][" << i <<"] := " << S2.at(i) << std::endl;
    }
  // }

    cerr << endl;
    degree_t degs[2];
    Integer degs1, degs2;
    std::vector<SMQP> S3 = a.subresultantInitialAtIdx (b, Symbol("x_1"), 0, degs1, degs2, &lazyInfo);
    degs[0] = degs1.get_d(); degs[1] = degs2.get_d();
    for (int i=0; i < S3.size(); ++i) {
      std::cerr << "[deg= "<< degs[i] << "] S3[0][" << i <<"] := " << S3.at(i) << std::endl;
    }

    S3 = a.subresultantInitialAtIdx (b, Symbol("x_1"), 1, degs1, degs2, &lazyInfo);
    degs[0] = degs1.get_d(); degs[1] = degs2.get_d();
    for (int i=0; i < S3.size(); ++i) {
      std::cerr << "[deg= "<< degs[i] << "] S3[1][" << i <<"] := " << S3.at(i) << std::endl;
    }

    S3 = a.subresultantInitialAtIdx (b, Symbol("x_1"), 2, degs1, degs2, &lazyInfo);
    degs[0] = degs1.get_d(); degs[1] = degs2.get_d();
    for (int i=0; i < S3.size(); ++i) {
      std::cerr << "[deg= "<< degs[i] << "] S3[2][" << i <<"] := " << S3.at(i) << std::endl;
    }

    S3 = a.subresultantInitialAtIdx (b, Symbol("x_1"), 4, degs1, degs2, &lazyInfo);
    degs[0] = degs1.get_d(); degs[1] = degs2.get_d();
    for (int i=0; i < S3.size(); ++i) {
      std::cerr << "[deg= "<< degs[i] << "] S3[4][" << i <<"] := " << S3.at(i) << std::endl;
    }

}


void test_smqpSRC3 (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag) {
  unsigned long long start;
  float elasped;

  SMQP a = SMQP ("- 13/6*x^10 + 1/15*x^8 - 13/4*x^7 - 9/7*x^6 - 11/9*x^3");
  SMQP b = SMQP ("- 65/3*x^9 + 8/15*x^7 - 91/4*x^6 - 54/7*x^5 - 11/3*x^2");

  std::cerr << "a := " << a << std::endl;
  std::cerr << "b := " << b << std::endl;
  int idx = 0;

    std::vector<SMQP> S1 = a.subresultantChain(b, Symbol("x"));
    for (int i=0; i < S1.size(); ++i) {
      std::cerr << "S1[" << i <<"] := " << S1.at(i) << std::endl;
    }

    // specSRC_AAZ *lazyInfo = NULL;
    // _startTimer(&start);
    // std::vector<SMQP> S2 = a.subresultantAtIdx (b, Symbol("x"), 0, &lazyInfo);
    // _stopTimer(&start, &elasped);
    // std::cerr << "[0-call] Time: " << elasped << std::endl;
    // for (int i=0; i < S2.size(); ++i) {
    //   std::cerr << "S2[ idx=0 ][ i= " << i <<" ] := " << S2.at(i) << std::endl;
    // }

    // for (idx = 1; idx < b.degree()-1; idx++) {
    //     _startTimer(&start);
    //     S2 = a.subresultantAtIdx (b, Symbol("x"), idx, &lazyInfo);
    //     _stopTimer(&start, &elasped);
    //     std::cerr << "[i-call] Time: " << elasped << std::endl;
    //     for (int i=0; i < S2.size(); ++i) {
    //       std::cerr << "S2[ idx=" << idx << " ][ i= "<< i <<" ] := " << S2.at(i) << std::endl;
    //     }
    // }

    // cerr << endl;
    // Integer degs1, degs2;
    // degree_t degs[2];
    // std::vector<SMQP> S3 = a.subresultantInitialAtIdx (b, Symbol("x"), idx, 
    //                                                   degs1, degs2, &lazyInfo);
    // degs[0] = degs1.get_d(); degs[1] = degs2.get_d();
    // for (int i=0; i < S3.size(); ++i) {
    //   std::cerr << "[deg= "<< degs[i] << "] S3[" << i <<"] := " << S3.at(i) << std::endl;
    // }

}


void time_mul_DUSP_vs_DUZP (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag) 
{

  int coefBoundBits = flag;
  int includeNeg = 1;
  float sparsity = 0;
  
  DUZP_t *A = buildRandomPoly_DUZP (sz1, coefBoundBits, sparsity, includeNeg);
  duspoly_t *a = convertToDUSP_DUZP (A, Pptr);
  duspoly_t *a_out = convertPolynomialFromMontgomery_spX (a, Pptr);
  
  DUZP_t *B = buildRandomPoly_DUZP (sz2, coefBoundBits, sparsity, includeNeg);
  duspoly_t *b = convertToDUSP_DUZP (B, Pptr);
  duspoly_t *b_out = convertPolynomialFromMontgomery_spX (b, Pptr);

  duspoly_t *c;

  int ntrials = 10;
  unsigned long long int start;
  float elapsed, sum;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    _startTimer (&start);
    multiplyPolynomials_DUZP (A, B);
    _stopTimer (&start, &elapsed);
    sum += elapsed;
  }
  std::cout << "[DUZP] coefBound= " << coefBoundBits << " sz1= " << sz1 << " sz2= " << sz2 << " p= " << Pptr->prime << " time= " << sum/ntrials << std::endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    _startTimer (&start);
    plainMulPolynomials_spX (a, b, &c, Pptr);
    _stopTimer (&start, &elapsed);
    sum += elapsed;
    freePolynomial_spX (&c);
  }
  std::cout << "[DUSP-inForm] coefBound= " << coefBoundBits << " sz1= " << sz1 << " sz2= " << sz2 << " p= " << Pptr->prime << " time= " << sum/ntrials << std::endl;

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    _startTimer (&start);
    plainMulPolynomials_outform_spX (a_out, b_out, &c, Pptr);
    _stopTimer (&start, &elapsed);
    sum += elapsed;
    freePolynomial_spX (&c);
  }
  std::cout << "[DUSP-CInt] coefBound= " << coefBoundBits << " sz1= " << sz1 << " sz2= " << sz2 << " p= " << Pptr->prime << " time= " << sum/ntrials << std::endl;

  ZZ_pX An;
  ZZ_pX Bn;
  ZZ_pX Cn;
  convertToZZ_pX (a_out, An, Pptr);
  convertToZZ_pX (b_out, Bn, Pptr);

  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    _startTimer (&start);
    PlainMul (Cn, An, Bn);
    _stopTimer (&start, &elapsed);
    sum += elapsed;
  }
  std::cout << "[NTL-ZZpX] coefBound= " << coefBoundBits << " sz1= " << sz1 << " sz2= " << sz2 << " p= " << Pptr->prime << " time= " << sum/ntrials << std::endl;

  ZZX AZ;
  ZZX BZ;
  ZZX CZ;
  convertToZZX (a_out, AZ);
  convertToZZX (a_out, BZ);
  
  sum = 0;
  for (int i = 0; i < ntrials; ++i) {
    _startTimer (&start);
    PlainMul (CZ, AZ, BZ);
    _stopTimer (&start, &elapsed);
    sum += elapsed;
  }
  std::cout << "[NTL-ZZX] coefBound= " << coefBoundBits << " sz1= " << sz1 << " sz2= " << sz2 << " p= " << Pptr->prime << " time= " << sum/ntrials << std::endl;


  freePolynomial_DUZP (A);
  freePolynomial_DUZP (B);
  freePolynomial_spX (&a);
  freePolynomial_spX (&b);
  freePolynomial_spX (&a_out);
  freePolynomial_spX (&b_out);
}

void test_bivarSysSolver (polysize_t sz1, polysize_t sz2, Prime_ptr* Pptr, polysize_t flag)
{
    int coefBoundBits = 1;
    float sparsity = 0.1;
    int includeNeg = 1;

    int d = sz1 + sz2;
    // std::vector<int> maxDegs = {sz1/2, sz2};
    // std::vector<int> maxDegs_b = {(sz1-1)/2, sz2-1};
    std::vector<int> maxDegs = {d, (int) sz1};
    std::vector<int> maxDegs_b = {d-1, (int) sz2};
    // std::vector<int> maxDegs = {(int) sz1, d};
    // std::vector<int> maxDegs_b = {(int) sz2, d-1};
    
    SMQP A;
    SMQP B;
    A.randomPolynomial(maxDegs, coefBoundBits,sparsity, includeNeg);
    B.randomPolynomial(maxDegs_b, coefBoundBits,sparsity, includeNeg);

    // std::ofstream out("polys.txt");
    // std::streambuf *coutbuf = std::cout.rdbuf();
    // std::cout.rdbuf(out.rdbuf());
    cout << "A := " << A << endl;
    cout << "B := " << B << endl;

    RegularChain<RationalNumber, SMQP> R(A.ringVariables()); 
    std::vector<SMQP> sys = {A, B};

    unsigned long long start;
    float elapsed;

    _startTimer(&start);
    R.triangularize(sys);
    _stopTimer(&start, &elapsed);

    cerr << "sz1: " << sz1 << " sz2: " << sz2 << " Triangularize Time: " << elapsed << endl;

    // for (auto res: R.triangularize(sys)) {
    //   cout << res << endl;
    // }

    // const char* sym[10] = {"x", "y", "z", "x4", "x5", "x6", "x7", "x8", "x9", "x10"}; // TEST
    // fprintf (stderr, "a := ");
    // printPoly_AAZ (stderr, a, sym, a->nvar);
    // fprintf (stderr, "\nb := ");
    // printPoly_AAZ (stderr, b, sym, b->nvar);

    // freePolynomial_AAZ (a);
    // freePolynomial_AAZ (b);
 }

int main (int argc, char** argv) {

    polysize_t sz1 = 10;
    polysize_t sz2 = 9;
    polysize_t flag = 1;
    prime_t pr =  4179340454199820289; // 9211409350344572929 // 2485986994308513793; // 180143985094819841;
    Prime_ptr* Pptr = smallprimefield_get_prime_constants (pr);

    // if (argc == 6 && strcmp(argv[1], "src")  == 0) {
    //   time_src(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), Pptr, flag);
    //   return 0;
    // }

    if (argc == 3) {
        sz1 = (polysize_t) atoi (argv[1]);
        sz2 = (polysize_t) atoi (argv[2]);
    }

    if (argc == 4) {
        sz1 = (polysize_t) atoi (argv[1]);
        sz2 = (polysize_t) atoi (argv[2]);
        // pr = fourier_primes_u64_table[atoi(argv[3])]; //
        pr = (prime_t) atoi (argv[3]);
        free (Pptr);
        Pptr = smallprimefield_get_prime_constants (pr);
    }

    unsigned long bound = 5;
    if (argc > 3) {
      bound = atoi(argv[3]);
    }
    int nthreads = -1;
    if (argc > 4) {
      nthreads = atoi(argv[4]);
    }

    ////////////////// Time
    // time_mul_DUSP_vs_DUZP (sz1, sz2, Pptr, 16);
    // std::cout << std::endl;
    // time_mul_DUSP_vs_DUZP (sz1, sz2, Pptr, 32);
    // std::cout << std::endl;
    // time_mul_DUSP_vs_DUZP (sz1, sz2, Pptr, 64);

   
    ////////////////// Tests:
    // fprintf(stderr, "\n----------------------------\npr = %lld \n----------------------------\n", Pptr->prime);
    // fprintf(stderr, "rsquare = %lld\n", Pptr->rsquare);
    // fprintf(stderr, "rsquare = %lld\n", Pptr->prime_inv);

    test_newFastDiv (sz1, sz2, Pptr, flag);
    test_freePoly     (sz1, sz2, Pptr, flag);
    test_deepCopyPoly (sz1, sz2, Pptr, flag);
    test_makePoly     (sz1, sz2, Pptr, flag);
    test_normalize    (sz1, sz2, Pptr, flag);
    test_convert      (sz1, sz2, Pptr, flag);
    test_shiftPoly    (sz1, sz2, Pptr, flag);
    test_splitPoly    (sz1, sz2, Pptr, flag);
    test_addPoly      (sz1, sz2, Pptr, flag);
    test_subPoly      (sz1, sz2, Pptr, flag);
    test_negPoly      (sz1, sz2, Pptr, flag);
    test_monicPoly    (sz1, sz2, Pptr, flag);
    test_mulPoly      (sz1, sz2, Pptr, flag);
    // // test_NewFFTMul    (sz1, sz2, Pptr, flag);
    test_divPoly      (sz1, sz2, Pptr, flag);
    test_GCD          (sz1, sz2, Pptr, flag);
    test_extGCD       (sz1, sz2, Pptr, flag);
    test_diff         (sz1, sz2, Pptr, flag);
    test_Mat4         (sz1, sz2, Pptr, flag);
    test_res          (sz1, sz2, Pptr, flag);
    test_FactBasics   (sz1, sz2, Pptr, flag);
    test_SFD          (sz1, sz2, Pptr, flag);
    // test_NTLBas       (sz1, sz2, Pptr, flag);
    // test_NTLFact      (sz1, sz2, Pptr, flag);
    // test_interpolation   (sz1, sz2, Pptr, flag);
    test_SRC          (sz1, sz2, Pptr, flag);
    // test_smqpSRC1      (sz1, sz2, Pptr, flag);
    // test_smqpSRC2      (sz1, sz2, Pptr, flag);
    // test_newDBZPSRC    (sz1, sz2, Pptr, flag);
    // test_newDUZPSRC    (sz1, sz2, Pptr, flag); 
    // test_bivarSysSolver (sz1, sz2, Pptr, flag);
    // test_smqpSRC3(sz1, sz2, Pptr, flag);

    // Other Test functions
    // test_RegGCDSRC    (sz1, sz2, Pptr, flag);
    // test_newDUZPSRC     (sz1, sz2, Pptr, flag);
    // test_newDBZPSRC     (sz1, sz2, Pptr, flag);
    // test_Fact         (sz1, sz2, Pptr, flag);
    // test_memmove          (sz1, sz2, Pptr, flag);
    // test_memleak          (sz1, sz2, Pptr, flag);
    // test_duspolyA         (sz1, sz2, Pptr, flag);
    // test_MPZ              (sz1, sz2, Pptr, flag);
    // test_modularRes       (sz1, sz2, Pptr, flag);
    // test_NewFFTMul       (sz1, sz2, Pptr, flag);
    // test_modularSubres   (sz1, sz2, Pptr, flag);
    // test_REMREM    (sz1, sz2, Pptr, flag);
    // _test_halfsubresZp    (sz1, sz2, Pptr, flag);
    // test_fftMul           (sz1, sz2, Pptr, flag);
    // test_biSubresPr            (sz1, sz2, Pptr, flag);
    // test_conertDBZP             (sz1, sz2, Pptr, flag);
    // test_biSubresZ            (sz1, sz2, flag, bound, nthreads);
    // test_fft_iter_main(sz1, sz2, Pptr, flag);
    // _test_subresultant

    fprintf (stderr, " \n");
    free(Pptr);
    return 0;
}
// 