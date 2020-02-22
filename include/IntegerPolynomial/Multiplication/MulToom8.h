#pragma once
/**
	Implementation of polynomial multiplication using 8-Way Toom-Cook for dividing the work, 
	plus Kronecker Substitution, and GMP Integer Multiplication.

	@author Farnam Mansouri
*/

#include "MulToom.h"

class MulToom8: public MulToom{
        public:

		/**
		 * @return The number of workers which is 16 for 8-Way Toom-Cook
		 */
		virtual char* numberOfWorkers(){
			return (char*) "16";
		};
	
		/**
		 * Multiply two polynomial using 8-Way Toom-Cook algorithm.
		 *
		 * @param a The first polynomial
		 * @param b The second polynomial
		 * @param c The result polynomial
		 */
		virtual void mulToom(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){
			mulToom8Par(a, b, c);
		};

		/**
		 * Shrink the Polynomial's coefficients with required number of zeros.
		 * (Since the 8-way Toom Cook algorithms works for polynomials with equal sizes and multiple of 8.)
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 */
		void growArrays(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b);

		/**
		 * Evaluate the polynomial with coefficients (c0,c1,c2,c3, c4, c5, c6, c7) at 
		 * 15 points (0,1,-1,1/2,-1/2,2,-2,1/4,-1/4,4,-4,1/8,-1/8,8,infinity).
		 * (r0,r1,...,r14) will be the evaluated results.
		 *
		 *  M = {{1,0,0,0,0,0,0}, {1,1,1,1,1,1,1}, {1,-1,1,-1,1,-1,1},
		 *       {128,64,32,16,8,4,2,1}, {128,-64,32,-16,8,-4,2,-1},
		 *       {1,2,4,8,16,32,64,128}, {1,-2,4,-8,16,-32,64,-128},
		 *       {16384,4096,1024,256,64,16,4,1}, {16384,-4096,1024,-256,64,-16,4,-1},
		 *       {1,4,16,64,256,1024,4096,16384}, {1,-4,16,-64,256,-1024,4096,-16384},
		 *       {2097152,262144,32768,4096,512,64,8,1}, {2097152,-262144,32768,-4096,512,-64,8,-1},
		 *       {1,8,64,512,4096,32768,262144,2097152}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}}
		 * Compute R = M * C.
		 *
		 * @param r evaluated result.
		 * @param c polynomial's coefficients.
		 */
		void evaluate(mpz_class * r, mpz_class * c);

		/**
		 * Interpolate the result polynomial by having 15 evaluated points.
		 *
		 * Having M and W as:
		 *  M = {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		 *	 {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
		 *	 {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},
		 *	 {16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1},
		 *	 {16384,-8192,4096,-2048,1024,-512,256,-128,64,-32,16,-8,4,-2,1},
		 *	 {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384},
		 *	 {1,-2,4,-8,16,-32,64,-128,256,-512,1024,-2048,4096,-8192,16384},
		 *	 {268435456,67108864,16777216,4194304,1048576,262144,65536,16384,4096,1024,256,64,16,4,1},
		 *	 {268435456,-67108864,16777216,-4194304,1048576,-262144,65536,-16384,4096,-1024,256,-64,16,-4,1},
		 *	 {1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864,268435456},
		 *	 {1,-4,16,-64,256,-1024,4096,-16384,65536,-262144,1048576,-4194304,16777216,-67108864,268435456},
		 *	 {4398046511104,549755813888,68719476736,8589934592,1073741824,134217728,16777216,2097152,262144,32768,4096,512,64,8,1},
		 *	 {4398046511104,-549755813888,68719476736,-8589934592,1073741824,-134217728,16777216,-2097152,262144,-32768,4096,-512,64,-8,1},
		 *	 {1,8,64,512,4096,32768,262144,2097152,16777216,134217728,1073741824,8589934592,68719476736,549755813888,4398046511104},
		 *	 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}}
		 *  W = {W0, W1, ... , W14}
		 * We want to compute C, where:
		 *  C = M^{-1} * W = {C0,C1, ... ,C14}
		 *
		 * We decompose the matrix multiplication and compute them in parallel.
		 *
		 * The coefficients of the result polynomial (C) will be stored in 'w'.
		 *
		 * @param w evaluated points, and the result polynomial.
		 */
		void interpolate(mpz_class *w);

        private:

		/**
		 * Parallel Implementation of 8-Way Toom Cook.
		 * This algorithm has 5 phazes:
		 *
		 * 1. -Work Division: Divide each of the polynomials to 8 parts.
		 *    -Conversions: Convert each of the parts to Big Integers Using Kronecker Substitution.
		 *  Now we have two polynomials with very big integer coefficients with degree 7.
		 *  a(x) = a0(x) + a1(x) * x^k + a2(x) * x^2k + a3(x) * x^3k + a4(x) * x^4k + a5(x) * x^5k + a6(x) * x^6k + a7(x) * x^7k.
		 *  b(x) = b0(x) + b1(x) * x^k + b2(x) * x^2k + b3(x) * x^3k + b4(x) * x^4k + b5(x) * x^5k + b6(x) * x^6k + b7(x) * x^7k.
		 *  --> A(X) = A0 + A1 * X + A2 * X^2 + A3 * X^3 + A4 * X^4 + A5 * X^5 + A6 * X^6 + A7 * X^7
		 *  --> B(X) = B0 + B1 * X + B2 * X^2 + B3 * X^3 + B4 * X^4 + B5 * X^5 + B6 * X^6 + B7 * X^7
		 *
		 * 2. Evaluation: Evaluate each of the polynomials at 15 points (0,1,-1,1/2,-1/2,2,-2,1/4,-1/4,4,-4,1/8,-1/8,8,infinity).
		 *  This step is equivalent to F = M * A and G = M * B where:
		 *  M = {{1,0,0,0,0,0,0}, {1,1,1,1,1,1,1}, {1,-1,1,-1,1,-1,1},
		 *       {128,64,32,16,8,4,2,1}, {128,-64,32,-16,8,-4,2,-1},
		 *       {1,2,4,8,16,32,64,128}, {1,-2,4,-8,16,-32,64,-128},
		 *       {16384,4096,1024,256,64,16,4,1}, {16384,-4096,1024,-256,64,-16,4,-1},
		 *       {1,4,16,64,256,1024,4096,16384}, {1,-4,16,-64,256,-1024,4096,-16384},
		 *       {2097152,262144,32768,4096,512,64,8,1}, {2097152,-262144,32768,-4096,512,-64,8,-1},
		 *       {1,8,64,512,4096,32768,262144,2097152}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}}
		 *  A = {A0,A1, ... ,A7}, B={B0,B1, ... ,B7}
		 *  F = {WA0, ... , WA15}, G = {WB0, ... , WB15}
		 *
		 * 3. Multiplications: For each of the evaluated points do W(point) = A(point) * B(point):
		 *  W = {A(0)*W(0), A(1)*B(1), A(-1)*B(-1), A(1/2)*B(1/2), A(-1/2)*B(-1/2), A(2)*B(2), A(-2)*B(-2), A(1/4)*B(1/4),
		 *       A(-1/4)*B(-1/4), A(4)*B(4), A(-4)*B(-4), A(1/8)*B(1/8), A(-1/8)*B(-1/8), A(8)*B(8), A(infinity)*B(infinity)}
		 *  W = {W0, W1, W2, W3, W4, W5, W6, W7, W8, W9, W10, W11, W12, W13, W14}
		 *
		 * 4. Interpolation: Generate a polynomial with degree 14 using the 15 evaluated points.
		 *  M = {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		 *       {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
		 *       {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},
		 *       {16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1},
		 *       {16384,-8192,4096,-2048,1024,-512,256,-128,64,-32,16,-8,4,-2,1},
		 *       {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384},
		 *       {1,-2,4,-8,16,-32,64,-128,256,-512,1024,-2048,4096,-8192,16384},
		 *       {268435456,67108864,16777216,4194304,1048576,262144,65536,16384,4096,1024,256,64,16,4,1},
		 *       {268435456,-67108864,16777216,-4194304,1048576,-262144,65536,-16384,4096,-1024,256,-64,16,-4,1},
		 *       {1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864,268435456},
		 *       {1,-4,16,-64,256,-1024,4096,-16384,65536,-262144,1048576,-4194304,16777216,-67108864,268435456},
		 *       {4398046511104,549755813888,68719476736,8589934592,1073741824,134217728,16777216,2097152,262144,32768,4096,512,64,8,1},
		 *       {4398046511104,-549755813888,68719476736,-8589934592,1073741824,-134217728,16777216,-2097152,262144,-32768,4096,-512,64,-8,1},
		 *       {1,8,64,512,4096,32768,262144,2097152,16777216,134217728,1073741824,8589934592,68719476736,549755813888,4398046511104},
		 *       {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}}
		 *  C = M^{-1} * W = {C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14}
		 *
		 * 5. -Conversion: Convert the coefficients of the polynomial to the polynomials. (15 coefficients)
		 *    -Merge: Merge the 15 polynomials to get the final result.
		 *   C(X) = C0 + C1 * X + C2 * X^2 + ... + C14 * X^14
		 *   c(x) = c0(x) + c1(x) * x^k + c2(x) * x^2k + ... c14(x) * x^14k
		 *   c(x) = c0 + c1 * x + ... + c_{16k - 1}x^{16k - 1}
		 *
		 * NOTE: This function assumes that the sizes of the polynomials are the same and divisible by 8.
		 *         A wrapper function (see 'growArrays') will handle for different sizes case.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Result polynomial
		 */
                void mulToom8Par(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c);
};

