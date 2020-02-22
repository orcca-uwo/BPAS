#pragma once
/**
	Difault class for multiplication. 
	All of the multiplication algorithms extend this class.

	@author Farnam Mansouri
*/

#include "../Poly.h"

class Mul{

	#define LOG_DEGREE(a, b) (a.getSize() == 0 || b.getSize() == 0) ? 0 : log2(MIN(a.getSize(), b.getSize()))
	#define LOG_DEGREE_P(a, b) (a->getSize() == 0 || b->getSize() == 0) ? 0 : log2(MIN(a->getSize(), b->getSize()))

	#define MAX_DIGITS(a, b) a.getCoefficientDigits() + b.getCoefficientDigits() + (LOG_DEGREE(a, b))
	#define RESULT_DEGREE(a, b) (a.getSize() == 0 || b.getSize() == 0) ? 0 : (a.getSize() + b.getSize() - 1)

	#define MAX_DIGITS_P(a, b) a->getCoefficientDigits() + b->getCoefficientDigits() + (LOG_DEGREE_P(a, b))

	#define MAX_COEFF_DIGITS(a, b) MAX(a->getCoefficientDigits(), b->getCoefficientDigits())+1


	#define RESULT_DEGREE_P(a, b) (a->getSize() == 0 || b->getSize() == 0) ? 0 : (a->getSize() + b->getSize() - 1)

	#define MAX_DIGITS_FAST_P(a, b) MAX(a->getCoefficientDigits(), b->getCoefficientDigits()) + log2(MIN(a->getSize(), b->getSize()))/2

	public:
		/**
		 * Virtual method for multiplying two univariate integer polynomials.
		 * This method will be implemented in different algorithms.
		 *
		 * @param a First polynomial
		 * @param b Second polynomial
		 * @return The result polynomial
		 */
		virtual UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b) = 0;

		/**
		 * Virtual method for multiplying two univariate integer polynomials.
		 * This method will be implemented in different algorithms.
		 *
		 * @param a First polynomial
		 * @param b Second polynomial
		 * @param c The result polynomial
		 */
		virtual void multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) = 0;

	protected:

		/**
		 * Select the representation base for the coefficients based on the number of bits.
		 * If the coefficients are small we use binary, if they are big we use base 32.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 * @param c Pointer to the result polynomial
		 */
		void setBitCounts(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c);

		/**
		 * Select the representation base for the coefficients based on the number of bits.
		 * If the coefficients are small we use binary, if they are big we use base 32.
		 *
		 * @param a Pointer to the first polynomial
		 * @param b Pointer to the second polynomial
		 */
		void setBitCounts(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b);

		/**
		 * Set the number of workers for cilk.
 		 *
		 * @param nworkers number of workers
		 */
		void setNumberOfWorkers(int nworkers);

		/**
		 * Set the number of workers for cilk.
		 *
		 * @param nworkers number of workers
		 */
		inline void setNumberOfWorkers(char* nworkers){
			__cilkrts_set_param("nworkers", nworkers);
		}

		/**
		 * Get the number of workers for cilk.
		 *
		 * @return the number of workers
		 */
		inline int getNumberOfWorkers(){
			return __cilkrts_get_nworkers();
		}

};
