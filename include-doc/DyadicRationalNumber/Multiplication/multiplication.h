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
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


