#pragma once
/**
	Implementation of polynomial multiplication using 4-Way Toom-Cook for dividing the work, 
	plus Kronecker Substitution, and GMP Integer Multiplication.

	@author Farnam Mansouri
*/

#include "MulToom.h"

class MulToom4: public MulToom{
	public:
		
		/**
		 * @return The number of workers which is 8 for 4-Way Toom-Cook
		 */
		virtual char* numberOfWorkers(){
                        return (char*) "8";
                };

		/**
		 * Multiply two polynomial using 4-Way Toom-Cook algorithm.
		 *
		 * @param a The first polynomial
		 * @param b The second polynomial
		 * @param c The result polynomial
		 */
                virtual void mulToom(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){
                        mulToom4Par(a, b, c);
                };

		/**
		 * Evaluate the polynomial with coefficients (c0,c1,c2,c3) at 7 points (0,1,-1,1/2,-1/2,2,infinity).
		 * (r0,r1,...,r6) will be the evaluated results.
		 *
		 * M = {{1,0,0,0},{1,1,1,1},{1,-1,1,-1},{8,4,2,1},{8,-4,2,-1},{1,2,4,8},{0,0,0,1}}
		 * Compute R = M * C.
		 *
		 * @param r evaluated result.
		 * @param c polynomial's coefficients.
		 */
                void evaluate(mpz_class * r, mpz_class * c);

		/**
		 * Shrink the Polynomial's coefficients with required number of zeros.
		 * (Since the 4-way Toom Cook algorithms works for polynomials with equal sizes and multiple of 4.)
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 */
                void growArrays(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b);

		/**
		 * Interpolate the result polynomial by having 7 evaluated points.
		 *
		 * Having M and W as:
		 *  M = {{1,0,0,0,0,0,0},{1,1,1,1,1,1,1},{1,-1,1,-1,1,-1,1},{64,32,16,8,4,2,1},{64,-32,16,8,-4,2,-1},{1,2,4,8,16,32,64},{0,0,0,0,0,0,1}}
		 *  W = {W0, W1, W2, W3, W4, W5, W6}
		 * We want to compute C, where:
		 *  C = M^{-1} * W = {C0,C1,C2,C3,C4,C5,C6}
		 *
		 * This will be done with a number of steps analogous to 'the Gaussian elimination' in which 
		 * we try to convert M to I (identity matrix); call it M' = I.
		 * This steps can be done on the 'W', since:
		 * We have M * C = W --> M' * C = W' --> I * C = W' --> C = W'
		 *
		 * The coefficients of the result polynomial (C) will be stored in 'w'.
		 *
		 * @param w evaluated points, and the result polynomial.
		 */
		void interpolate(mpz_class * w);

	private:

		/**
		 * Parallel Implementation of 4-Way Toom Cook.
		 * This algorithm has 5 phazes:
		 *
		 * 1. -Work Division: Divide each of the polynomials to 4 parts.
		 *    -Conversions: Convert each of the parts to Big Integers Using Kronecker Substitution.
		 *  Now we have two polynomials with very big integer coefficients with degree 3. (X = x^k)
		 *  a(x) = a0(x) + a1(x) * x^k + a2(x) * x^2k + a3(x) * x^3k --> A(X) = A0 + A1 * X + A2 * X^2 + A3 * X^3
		 *  b(x) = b0(x) + b1(x) * x^k + b2(x) * x^2k + b3(x) * x^3k --> B(X) = B0 + B1 * X + B2 * X^2 + B3 * X^3
		 *
		 * 2. Evaluation: Evaluate each of the polynomials at 7 points (0,1,-1,1/2,-1/2,2,infinity).
		 *  This step is equivalent to F = M * A and G = M * B where:
		 *  M = {{1,0,0,0},{1,1,1,1},{1,-1,1,-1},{8,4,2,1},{8,-4,2,-1},{1,2,4,8},{0,0,0,1}}
		 *  A = {A0,A1,A2,A3}, B={B0,B1,B2,B3}
		 *  F = {WA0,WA1,WA2,WA3,WA4,WA5,WA6}, G = {WB0,WB1,WB2,WB3,WB4,WB5,WB6}}
		 *
		 * 3. Multiplications: For each of the evaluated points do W(point) = A(point) * B(point):
		 *  W = {A(0)*W(0), A(1)*B(1), A(-1)*B(-1), A(1/2)*B(1/2), A(-1/2)*B(-1/2), A(2)*B(2), A(infinity)*B(infinity)}
		 *  W = {W0, W1, W2, W3, W4, W5, W6}
		 *
	 	 * 4. Interpolation: Generate a polynomial with degree 6 using the 7 evaluated points.
		 *  M = {{1,0,0,0,0,0,0},{1,1,1,1,1,1,1},{1,-1,1,-1,1,-1,1},{64,32,16,8,4,2,1},{64,-32,16,8,-4,2,-1},{1,2,4,8,16,32,64},{0,0,0,0,0,0,1}}
		 *  C = 	M^{-1} * W = {C0,C1,C2,C3,C4,C5,C6}
		 *
		 * 5. -Conversion: Convert the coefficients of the polynomial to the polynomials. (7 coefficients)
		 *    -Merge: Merge the 7 polynomials to get the final result.
		 *  C(X) = C0 + C1 * X + C2 * X^2 + C3 * X^3 + C4 * X^4 + C5 * X^5 + C6 * X^6
		 *  c(x) = c0(x) + c1(x) * x^k + c2(x) * x^2k + c3(x) * x^3k + c4(x) * x^4k + c5(x) * x^5k + c6(x) * x^6k
		 *  c(x) = c0 + c1 * x + ... + c_{8k - 1}x^{8k - 1}
		 *
		 * NOTE: This function assumes that the sizes of the polynomials are the same and divisible by 4.
		 *         A wrapper function (see 'growArrays') will handle for different sizes case.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Result polynomial
		 */
		void mulToom4Par(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c);
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


