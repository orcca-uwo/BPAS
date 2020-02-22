#pragma once
/**
	Representing the Toom-Cook abstract class for doing the polynomial multiplication.
	This will be extended for different 'K's (K-way Toom-Cook).

	@author Farnam Mansouri
*/

#include "Mul.h"
#include <sstream>
#include <string.h>

class MulToom: public Mul{
	public:		

		/**
		 * Shrink the Polynomial's coefficients with required number of zeros.
		 * (Since the K-way Toom Cook algorithms works for polynomials with equal sizes and multiple of K.)
		 *
		 * @param a Pointer to the first polynomial
 		 * @param b Pointer to the second polynomial
		 */
		virtual void growArrays(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b) = 0;

		/**
		 * Evaluate the polynomial with coefficients (c_0,...,c_{k-1}) at 2k-1 points (P_0,...,P_{2k-3},P_infinity).
		 * (r_0,...,r_{2k-2}) will be the evaluated results.
		 *
		 * M = {{1,...,0},{1,...,1},...,{0,...,1}}
		 * Compute R = M * C.
		 *
		 * @param r evaluated result.
		 * @param c polynomial's coefficients.
		 */
		virtual void evaluate(mpz_class * r, mpz_class * c) = 0;

		/**
		 * Interpolate the result polynomial by having 2k-1 evaluated points.
 		 *
		 * Having M and W as:
		 *  M = {{1,...,0},{1,...,1},...,{0,...,1}}
		 *  W = {W_0, ..., W_{2k-2}}
		 * We want to compute C, where:
		 *  C = M^{-1} * W = {C_0,...,C_{2k-2}}
		 *
		 * The coefficients of the result polynomial (C) will be stored in 'w'.
		 *
		 * @param w evaluated points, and the result polynomial.
		 */
		virtual void interpolate(mpz_class * w) = 0;

		/**
		 * @return The number of workers
		 */
		virtual char* numberOfWorkers() = 0;

		/**
		 * Multiply two polynomial using Toom-Cook algorithm.
		 *
		 * @param a The first polynomial
		 * @param b The second polynomial
		 * @param c The result polynomial
		 */
		virtual void mulToom(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) = 0;

		/**
		 * Multiplying two univariate integer polynomials.
		 *
		 * @param a First polynomial
		 * @param b Second polynomial
		 * @return The result polynomial
		 */
		virtual UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
			int workers = getNumberOfWorkers();
                        setNumberOfWorkers(numberOfWorkers());

                        int resultSize = RESULT_DEGREE_P(a, b);
                        int aSize = a->getSize();
                        int bSize = b->getSize();

                        growArrays(a, b);

                        UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
                        mulToom(a, b, &c);

                        c.shrink(resultSize);
                        a->shrink(aSize);
                        b->shrink(bSize);

                        setNumberOfWorkers(workers);

                        return c;

		};

		/**
		 * Multiplying two univariate integer polynomials.
		 *
		 * @param a First polynomial
		 * @param b Second polynomial
		 * @param c The result polynomial
		 **/ 
		virtual void multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c) {
			int workers = getNumberOfWorkers();
			setNumberOfWorkers(numberOfWorkers());

			int resultSize = RESULT_DEGREE_P(a, b);
			int aSize = a->getSize();
			int bSize = b->getSize();

			growArrays(a, b);

			int newSize = RESULT_DEGREE_P(a, b);
			c->grow(newSize);

			mulToom(a, b, c);

			c->shrink(resultSize);
			a->shrink(aSize);
			b->shrink(bSize);

			setNumberOfWorkers(workers);
		}

		/**
		 * Multiplying two univariate integer polynomials.
		 *
		 * @param a First polynomial
		 * @param b Second polynomial
		 * @param c The result polynomial
		 * @param k 4 or 8
		 */
		virtual void multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c, int k) {
			int workers = getNumberOfWorkers();
                        setNumberOfWorkers(numberOfWorkers());

			int max = MAX(a->getSize(), b->getSize());
			int r = max % k;
			if(r != 0) max += k - r;

			mpz_class* fcoef = new mpz_class[max];
			mpz_class* gcoef = new mpz_class[max];
			for (int i = 0; i < max; ++i) {
				if (i < a->getSize())
					fcoef[i] = a->getCoefficient(i);
				else
					fcoef[i] = 0;
				if (i < b->getSize())
					gcoef[i] = b->getCoefficient(i);
				else
					gcoef[i] = 0;
			}
			UnivariateIntegerPolynomial f(max, fcoef);
			UnivariateIntegerPolynomial g(max, gcoef);
			int size = (max << 1) -1;
			mpz_class* hcoef = new mpz_class[size];
			for (int i = 0; i < size; ++i)
				hcoef[i] = 0;
			UnivariateIntegerPolynomial h (size, hcoef);

                        mulToom(&f, &g, &h);

			for (int i = 0; i < c->getSize(); ++i)
				c->setCoefficient(i, h.getCoefficient(i));

                        setNumberOfWorkers(workers);
			delete [] fcoef;
			delete [] gcoef;
			delete [] hcoef;
		};

	protected:

		/**
		 * Multiply two GMP integers.
		 * This function wraps the GMP multiplication, so we can execute it in parallel.
		 *
		 * @result Result integer
		 * @param a First integer
		 * @param b Second integer
		 */
		void multiplyGMP(mpz_class *result, mpz_class& a, mpz_class& b){
			*result = a * b;	
		}
};
