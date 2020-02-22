#pragma once
/**
 Bivariate Polynomial class for representing Bivariate Integer Polynomial in the densed-fashion which is a
 vector of C integers.
 
 This class provide basic functions for polynomials.
 
 @author Farnam Mansouri

 */

#include "../FFT/src/modpn.h"
#include "Poly.h"

class BivariatePolynomial
{
	private:
		sfixn * coefficients;
                int d;//d = the partial degree in respect to y
                int K;//K = the partial degree in respect to x
                int M;//M = the number of bits for the coefficients
                int size;// size = d * K

        public:

		inline int getK(){
			return K;
		}

		inline sfixn * getCoefficients(){
			return coefficients;
		}

		inline int getSize(){
			return size;
		}

		inline int setSize(int s){
			size = s;
			coefficients = new sfixn[size];
		}

		inline void setCoefficients(sfixn * coeff){
			coefficients = coeff;
		}

		BivariatePolynomial(int d, int K, int M);

		BivariatePolynomial(int d, int K, int M, int n);

		/**
                 * Free the heap memory used by the coefficients.
                 */
                void freeHeap();

		/**
                 * Prints all of the coefficients of the polynomial.
                 */
                void print();

		/**
                 * Set all of the coeeficients of the polynomial to 0.
                 */
                void setToZero();

		/**
		 * Convert a very large integer to a polynomial.
		 * This will be used to converting a univariate polynomial
		 * with large coefficients to a bivariate polynomial with 
		 * relatively small coefficients.
		 *
		 * The results will be stored in the coefficients array in the 
		 * corresponding index which is starts from 'index * K' for a 
		 * specific 'index'.
		 *
		 * univariate --> a1    * y ^ 0 + a2    * y ^ 1 + ... + ad    * y ^ (d-1)
		 * bivariate  --> A1(x) * y ^ 0 + A2(x) * y ^ 1 + ... + Ad(x) * y ^ (d-1)
		 *
		 * Ai(x) = b1 * x ^ 0 + b2 * x ^ 1 + ... + bK * x ^ (K-1)
		 *
		 * @param coeff the large integer
		 * @param index index of the coefficients in the univariate polynomial
		 */
		void convertFromBigIntegerSigned(mpz_class coeff, int index);	

		/**
		 * Reconstruct the large integer corresponding to a polynomial 
		 * in respect to x. 
		 * This is a reverse function for 'convertFromBigIntegerSigned'.
		 * See its comments for more details.
		 *
		 * @param index the starting index for poly's coefficients.
		 */
		mpz_class coefficientReconstruction(int index);

		/**
		 * Adapt (moduli) the coefficients over a prime field.
		 *
		 * @param prime the prime number
		 */
		int * adapt(int prime);
		
		void writeToFile(const char* name);
		void writeToFile(const char* name, int n);
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


