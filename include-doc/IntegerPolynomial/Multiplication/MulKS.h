#pragma once
/**
	Implementation of polynomial multiplication using Kronecker Substitution trick and GMP Integer Multiplication.

	@author Farnam Mansouri
*/

#include "Mul.h"

class MulKS: public Mul{
	public:
		virtual UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
			UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
		        setBitCounts(a, b, &c);
		        mulBigIntSer(a, b, &c, 0);
			return c;
		};
		virtual void multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){
                        setBitCounts(a, b, c);
                        mulBigIntSer(a, b, c, 0);
                };

		/**
		 * Polynomial Multiplication using Kronecker substitution (also known as Bruno Salvy's trick).
		 * This trick has 3 parts:
		 *      1. Convert representation of polynomials to big integers. (With Zero Padding)
		 *      2. Do the integer multiplication.
		 *      3. Convert back to the polynomial representation.
		 * This method uses GMP library for the big integer representation and its arithmetic.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
 		 * @param startIndex starting index for the result polynomial
		 */
		void mulBigIntSer(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex);

	private:
		
		/**
		 * Polynomial Multiplication using Reciprocal Evaluation Points (Advanced Kronecker substitution).
		 * This trick has 3 parts (Evaluation is equal to Zero Padding)
		 *      1. Evaluate 'a' and 'b' on 2^N and 2^{-N} where N = b + e/2. (In classic KS N = 2b + e)
		 *      2. Compute c(2^N) = a(2^N).b(2^N) and 2^{2N(L-1)}c(2^{-N})=2^{N(L-1)}a(2^{-N}).=2^{N(L-1)}b(2^{-N})
		 *      3. Extract the result's coefficients.
		 * For Big Integer Multiplications, GMP Library is used.
		 * This only works for polynomials with positive coefficients.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 * @param startIndex starting index for the result polynomial
		 */
		void mulBigIntSerReciprocalUnsigned(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex);

		/**
		 * Compute the first intermediate result for KS using Reciprocal eval. points.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param digitCount maximum number of the digits for the coefficients
		 * @return the result
		 */
		mpz_class generateFirstInteger(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, int digitCount);

		/**
		 * Compute the second intermediate result for KS using Reciprocal eval. points.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param digitCount maximum number of the digits for the coefficients
		 * @return the result
		 */
                mpz_class generateSecondInteger(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, int digitCount);

};
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


