
#ifndef _MAPLE_INTERFACE_H_
#define _MAPLE_INTERFACE_H_

#include "../IntegerPolynomial/SMZP_Support.h"

#include <maplec.h>
#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif


/********************************
 * The MapleInterface is an interface to a running maple kernel from C.
 * It can be used by the program to communicate with maple, execute things 
 * on the kernel, and generally pass data back and forth.
 * 
 * WARNING: Since the maple kernel is a singleton, it not thread-safe. Calling 
 * the kernel from many threads will lead to trouble.
 *
 ********************************/




/**
 * Get a reference to the maple kernel vector. 
 * This pointer is a singleton that exists for the lifetime of the program.
 * 
 * Note: Maple kernels cannot be stopped and subsequently restarted via the C interface.
 *       Thus, this is the one and only kernel for maple during the entire program. 
 *
 * returns a pointer to the maple kernel vector.
 */
MKernelVector* getMapleKVSingleton();

/**
 * Start up the maple kernel using some arguments.
 * argc: the number of arguments.
 * argv: the arguments to pass to kernel initialization. 
 */
void startMapleKernel(int argc, char* argv[]);

/**
 * Stop the maple kernel.
 * Note: once stopped, the kernel can never be started again during the lifetime
 *       of a program.
 */ 
void stopMapleKernel();

/**
 * Restart the maple kernel.
 * Note: As of maple 2017, this has no effect on the maple kernel via the C interface.
 *       The kernel is *not* restarted.
 */
// void restartMapleKernel();

/**
 * Convert an multivariate integer polynomial (AltArrZ_t) to a maple ALGEB object.
 * aa: the polynomial to convert.
 * syms: the symbols representing each variable in the polynomial.
 *
 * returns the ALGEB constructed from aa. 
 */
ALGEB convertToMaple_AAZ(const AltArrZ_t* aa, const char** syms);

/**
 * Conver a maple ALGEB object to a string. The returned string is dynamically
 * allocated and must be freed by the caller.
 *
 * in: The ALGEB to convert.
 *
 * returns the string representation of the ALGEB in.
 */
char* algebToString_MplInt(ALGEB in);

/**
 * Get the GCD of two multivariate integer polynomials using Maple. Both polynomials
 * should belong to the same ring (i.e. have the same variable in the same order).
 *
 * aa: the first polynomial.
 * ba: the second polynomial.
 * syms: the symbols representing each variable in the polynomials.
 *
 * returns the GCD of the two input polynomials.
 */
AltArrZ_t* gcd_MplInt_AAZ(const AltArrZ_t* aa, const AltArrZ_t* ba, const char** syms);

/**
 * Get the GCD of two algebraic types by calling Maple's gcd. 
 *
 * f: the first algebraic object as a string.
 * g: the second algebraic object as a string.
 *
 * returns the GCD of the two inputs.
 */
char* gcd_MplInt_string(const char* f, const char* g);

/**
 * Factor a polynomial over the integers.
 * A common integer factor is returned in numericFact. Irreducible polynomial factors
 * are returned in retFactsPtr with a parallel array retExps encoding the
 * exponents of each factor. 
 * 
 * aa: the polynomial to factor
 * symsL the symbols representing each variable in the polynomial.
 * numfacts: the number of factors returned in retFactsPtr and retExps
 * retFactsPtr: a pointer to the dynamically-allocated return array of factors
 * retExps: a pointer to the dynamically-allocated return array of exponents
 * numericFact: an integer factor (or 1 if no integer factor)
 * 
 */
void factorPolynomial_MplInt_AAZ(const AltArrZ_t* aa, const char** syms, int* numFacts, AltArrZ_t*** retFactsPtr, int** retExps, mpz_t numericFact);

/**
 * Factors a polynomial by calling maple's factors command.
 * A common integer factor is returned in numericFact. Irreducible polynomial factors
 * are returned in retFactsPtr with a parallel array retExps encoding the
 * exponents of each factor. 
 * 
 * poly: the polynomial to factor, the last character should be a colon (:).
 * numfacts: the number of factors returned in retFactsPtr and retExps
 * retFactsPtr: a pointer to the dynamically-allocated return array of factors
 * retExps: a pointer to the dynamically-allocated return array of exponents
 * numericFact: an integer factor (or 1 if no integer factor)
 * 
 */
void factorPolynomial_MplInt_string(const char* poly, int* numFacts, char*** retFactsPtr, char*** retExps, char** numericFact);

/**
 * Determine if the the result of triangularize command is correct.
 * inputs: {<F>, <RC>, <vars>, <results>}
 * nInputs: number of inputs (i.e. 4)
 * isLazard: Check the results in the Lazard sense?
 *
 * returns true iff the result is correct according to RegularChains package.
 */
int triangularizeValidate_MplInt(const char** inputs, int nInputs, int isLazard);


#ifdef __cplusplus
}
#endif


#endif