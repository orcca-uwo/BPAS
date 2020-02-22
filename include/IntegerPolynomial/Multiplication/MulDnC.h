#pragma once
/**
	Implementation of classical Divide & Conquer for polynomial multiplication.

	@author Farnam Mansouri
*/

#include "MulKS.h"

class MulDnC: public Mul{
	public:
                virtual UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
			UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
		        BASE_MUL = MAX(a->getSize(), b->getSize()) / 4;
		        mulDnCBigIntPar(a, b, &c, 0);
		        return c;
		};

		virtual void multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){
			BASE_MUL = MAX(a->getSize(), b->getSize()) / 4;
                        mulDnCBigIntPar(a, b, c, 0);
		};


		/**
		 * Multipply two polynomial using Big Integer method for the base case in parallel.
		 * The work decomposition is like Toom-Cook implementations, so we could have a fair comparison.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @return The result polynomial
		 */
		UnivariateIntegerPolynomial mulDnC(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b);

		/**
		 * Multipply two polynomial using Big Integer method for the base case in parallel.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @return The result polynomial
		 */
		UnivariateIntegerPolynomial mulDnCBigIntParallelStatic1(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b);

		/**
		 * Multipply two polynomial using Big Integer method for the base case in parallel.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @return The result polynomial
		 */
		UnivariateIntegerPolynomial mulDnCBigIntParallelStatic2(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b);

		/**
		 * Muliplication of polynomials with signed coefficients by decomposing the multiplication
		 * into 4 multiplications.
		 * A = APlus - AMinus, B = BPlus - BMinus
		 * A * B = APlus * BPlus - APlus * BMinus - AMinus * BPlus + AMinus * BMinus
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @return The result polynomial
		 */
		UnivariateIntegerPolynomial mulSignedDecomposed(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b);

		int BASE_MUL; // Threshold for the divide and conquere's base case

                MulDnC();

	private:
		
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
		 * This method uses two auxilary arrays for computing C1 and C2 in order to get the maximum parallelism.
		 *
		 * This function uses Kronecker substitution (also known as Bruno Salvy's trick) for the base case's multiplication.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 * @param startIndex starting index for the result polynomial
		 */
		void mulDnCBigIntPar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex);

		/**
		 * Polynomial Multiplication using one auxilary array. This implementation has a sync between the
		 * four divided multiplications. As a result the parallellism is not going to be maximum; but in the
		 * sense of memory usage is working better.
		 *
		 * ATTENTION:   THIS FUNCTION HAS A BUG WHEN THE LEVEL OF RECURSION IS MORE THAN SOME CONSTANT!
		 *              TO FIX THE BUG WE HAVE TO FIX THE AUXILARY ARRAY USAGE, I GAVE UP FIXING IT SINCE
		 *              IT WAS SLOWER THAN OTHER METHODS ANYWAY.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 * @param tmp Auxilary array for intermediate result
 		 * @param startIndex starting index for the result polynomial
		 */
                void mulDnCBigIntParStatic(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, UnivariateIntegerPolynomial* tmp, int startIndex);

		/**
		 * Polynomial Multiplication using 2 auxilary arrays. This implementation has the maximum parallelism; but in the sense
		 * of memory usage, is not very good.
		 *
		 * ATTENTION:   THIS FUNCTION HAS A BUG WHEN THE LEVEL OF RECURSION IS MORE THAN SOME CONSTANT!
		 *              TO FIX THE BUG WE HAVE TO FIX THE AUXILARY ARRAYS USAGE, I GAVE UP FIXING IT SINCE
		 *              IT WAS SLOWER THAN OTHER METHODS ANYWAY.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 * @param tmp0 First auxilary array for intermediate result
		 * @param tmp1 Second auxilary array for intermediate result
		 * @param startIndex starting index for the result polynomial
		 */
                void mulDnCBigIntParStatic(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, UnivariateIntegerPolynomial* tmp0, UnivariateIntegerPolynomial* tmp1, int startIndex);

		/**
		 * Iterative Parallel multiplication using BigInt (KS) method for the base case.
		 *
		 * @param Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 * @param startIndex starting index for the result polynomial
		 */
                void mulIterBigIntPar(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int startIndex);

};
