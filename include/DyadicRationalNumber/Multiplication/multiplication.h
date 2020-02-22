#ifndef _UNIPOLYMULTIPLICATION_H_
#define _UNIPOLYMULTIPLICATION_H_

/**
 * Function Interface of Univariate Rational Polynomial Multiplication
 * The denom of rational coefficients must be the power of 2
 **/

#include "../globals.h"
#include "../../IntegerPolynomial/Multiplication/Multiplication.h"

/**
 * Get the max number of bits of denom of rational coefficients
 * 
 * Input:
 * @a: Rational coefficients
 * @n: Size of a
 *
 * Return:
 * 0: Coefficients are integer
 * >0: The max number of bits of denom
 **/
unsigned long long getMaxBits(DyadicRationalNumber* a, int n);
unsigned long long getMaxBits(mpq_class* a, int n);

/**
 * Univariate polynomial multiplication
 * 
 * Output:
 * @mul: =a*b
 * Input:
 * @a: An input univariate polynomial, with rational coefficients
 * @n: Size of a
 * @b: An input univariate polynomial, with rational coefficients
 * @m: Size of b
 **/
void univariateMultiplication(lfixq* mul, DyadicRationalNumber* a, int n, DyadicRationalNumber* b, int m);
void univariateMultiplication(lfixq* mul, mpq_class* a, int n, mpq_class* b, int m);

/**
 * Univariate polynomial multiplication
 * 
 * Output:
 * @mul: =a*b
 * Input:
 * @a: An input univariate polynomial, with rational coefficients
 * @n: Size of a
 * @b: An input univariate polynomial, with integer coefficients
 * @m: Size of b
 **/
void univariateMultiplication(lfixq* mul, DyadicRationalNumber* a, int n, lfixz* b, int m);
void univariateMultiplication(lfixq* mul, mpq_class* a, int n, lfixz* b, int m);

/**
 * Univariate polynomial multiplication
 * Calling from integer univariate polynomial multiplication
 *
 * Output:
 * @mul: =a*b
 * Input:
 * @a: An input univariate polynomial, with integer coefficients
 * @n: Size of a
 * @b: An input univariate polynomial, with integer coefficients
 * @m: Size of b
 **/
void univariateMultiplication(lfixz* mul, lfixz* a, int n, lfixz* b, int m);

/**
 * Naive univariate multiplication
 **/
void naiveUnivariateMultiplication(lfixz* mul, lfixz* a, int n, lfixz* b, int m);

#endif
