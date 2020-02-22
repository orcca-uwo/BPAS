#pragma once
/**
	Implementation of Naive algorithm for polynomial multiplication.

	@author Farnam Mansouri
*/

#include "Mul.h"

class MulNaive: public Mul{
	public:
		virtual UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
			UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
        		mulNaiveSer(a, b, &c, 0);
                        return c;
		};
		virtual void multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){
			mulNaiveSer(a, b, c, 0);
		};
		int BASE_MUL; // Threshold for the divide and conquere's base case
                MulNaive();
	private:

		/**
		 * Naivest Implementation of polynomial multiplication C = A * B
		 * 2 Nested loops...
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
 		 * @param startIndex starting index for the result polynomial
		 */
		void mulNaiveSer(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex);

		/**
		 * A = [A0,A1], B = [B0,B1] --> 4 Recursive Calls
		 *
		 *    A0*B0 and A1*B1 will be executed concurrently in-place.
		 *   SYNC
		 *    A1*B0 amd A0*B1 will be executed concurrently with temporary arrays.
		 *   SYNC
		 *
		 * C = [C0,C1,C2,C3]
		 *      C0, from A0*B0
		 *      C1, from A1*B0, A0*B0, A0*B0
		 *      C2, from A1*B0, A0*B0, A1*B1
		 *      C3, from A1*B1
		 *
		 * This function uses Naive method for the base case's multiplication.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 * @param startIndex starting index for the result polynomial
		 */
                void mulDnRNaivePar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex);

		/**
		 * Iterative Parallel Multiplication using Naive method for the base case.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 * @param startIndex starting index for the result polynomial
		 */
		void mulIterNaivePar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex);

};
