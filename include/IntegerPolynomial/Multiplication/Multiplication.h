/**
	Top Interface for calling Univariate Integer Densed UnivariateIntegerPolynomial Multiplication.
	All of these functions will choose the algorithm based on the number of cores and size of the input polynomials.

	@author Farnam Mansouri
*/
#ifndef _MULTIPLICATION_H_
#define _MULTIPLICATION_H_

#include "Mul.h"
#include "MulKS.h"
#include "MulDnC.h"
#include "MulToom.h"
#include "MulToom4.h"
#include "MulToom8.h"
#include "MulSSA.h"
#include "MulNaive.h"

class Multiplication{

#define CONV_THRESHOLD	3000
#define KS_THRESHOLD	1024
#define NAIVE_THRESHOLD	64

private:
	MulToom4 toom4;
	MulToom8 toom8;
	MulKS ks;
	MulNaive naive;
	MulSSA ssa;
public:

	/**
	 * Multiply 2 polynomials with integer coefficients.
	 *
	 * @param r Pointer to the coefficients of the result polynomial
	 * @param a Pointer to the coefficients of the first polynomial
	 * @param b Pointer to the coefficients of the second polynomial
	 * @param aSize Size of the first polynomial
	 * @param bSize Size of the second polynomial
	 */
	void multiply(mpz_class * r, mpz_class *a, mpz_class *b, int aSize, int bSize){
		int worker = __cilkrts_get_nworkers();
	
		UnivariateIntegerPolynomial aPoly(aSize, a);
	        UnivariateIntegerPolynomial bPoly(bSize, b);

		UnivariateIntegerPolynomial rPoly;
	
		if((aSize < NAIVE_THRESHOLD || bSize < NAIVE_THRESHOLD)){
			// Too Small to use any fancy trick -> Apply NAIVE algorithm
			rPoly = naive.multiply(&aPoly, &bPoly);
		} else if((worker == 1) || (aSize <= KS_THRESHOLD && bSize <= KS_THRESHOLD)){
			// Single-Core machine, or Small enough to use any parallel code -> KRONECKER SUBSTITUTION
			rPoly = ks.multiply(&aPoly, &bPoly);
		} else if(worker <= 6){
			// Size is big enough, but the # of cores are low -> 4-WAY TOOM-COOK
			rPoly = toom4.multiply(&aPoly, &bPoly);
		} else {
			// The # of cores are hight
			if(aSize >= CONV_THRESHOLD && bSize >= CONV_THRESHOLD){
				// Size is big enough -> 2-CONVOLUTION
				rPoly = ssa.multiply(&aPoly, &bPoly);
			} else {
				// Size is not big enough -> 8-WAY TOOM-COOK
				rPoly = toom8.multiply(&aPoly, &bPoly);
			}
		}

		for(int i = 0; i < rPoly.getSize(); i++) 
			r[i] = rPoly.getCoefficient(i);

		rPoly.freeHeap();	                        
	}

	/**
	 * Multiply 2 polynomials with integer coefficients.
	 *
	 * @param a Pointer to the coefficients of the first polynomial
	 * @param b Pointer to the coefficients of the second polynomial
	 * @param aSize Size of the first polynomial
	 * @param bSize Size of the second polynomial
	 */
        mpz_class * multiply(mpz_class *a, mpz_class *b, int aSize, int bSize){

                int worker = __cilkrts_get_nworkers();

                UnivariateIntegerPolynomial aPoly(aSize, a);
                UnivariateIntegerPolynomial bPoly(bSize, b);

		UnivariateIntegerPolynomial rPoly;

		if((aSize < NAIVE_THRESHOLD || bSize < NAIVE_THRESHOLD)){
                        // Too Small to use any fancy trick -> Apply NAIVE algorithm
                        rPoly = naive.multiply(&aPoly, &bPoly);
                } else if((worker == 1) || (aSize <= KS_THRESHOLD && bSize <= KS_THRESHOLD)){ 
                        // Single-Core machine, or Small enough to use any parallel code -> KRONECKER SUBSTITUTION
                        rPoly = ks.multiply(&aPoly, &bPoly);
                } else if(worker <= 6){
                        // Size is big enough, but the # of cores are low -> 4-WAY TOOM-COOK
                        rPoly = toom4.multiply(&aPoly, &bPoly);
                } else {
                        // The # of cores are hight
                        if(aSize >= CONV_THRESHOLD && bSize >= CONV_THRESHOLD){
                                // Size is big enough -> 2-CONVOLUTION
                                rPoly = ssa.multiply(&aPoly, &bPoly);
                        } else {
                                // Size is not big enough -> 8-WAY TOOM-COOK
                                rPoly = toom8.multiply(&aPoly, &bPoly);
                        }
                }

		mpz_class * r = new mpz_class[rPoly.getSize()];

		for(int i = 0; i < rPoly.getSize(); i++)
                        r[i] = rPoly.getCoefficient(i);

		rPoly.freeHeap();

		return r;
        }


	/**
	 * Multiply 2 polynomial with integer coefficients.
	 *
	 * @param r Result polynomial
	 * @param a First polynomial
	 * @param b Second polynomial
	 */
	void multiply(UnivariateIntegerPolynomial *r, UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
	
	        int worker = __cilkrts_get_nworkers();

		if((a->getSize() < NAIVE_THRESHOLD || b->getSize() < NAIVE_THRESHOLD)){
                        // Too Small to use any fancy trick -> Apply NAIVE algorithm
                        naive.multiply(a, b, r);
                } else if((worker == 1) || (a->getSize() <= KS_THRESHOLD && b->getSize() <= KS_THRESHOLD)){ 
                        // Single-Core machine, or Small enough to use any parallel code -> KRONECKER SUBSTITUTION
                        ks.multiply(a, b, r);
                } else if(worker <= 6){
                        // Size is big enough, but the # of cores are low -> 4-WAY TOOM-COOK
                        toom4.multiply(a, b, r);
                } else {
                        // The # of cores are hight
                        if(a->getSize() >= CONV_THRESHOLD && b->getSize() >= CONV_THRESHOLD){
                                // Size is big enough -> 2-CONVOLUTION
                                ssa.multiply(a, b, r);
                        } else {
                                // Size is not big enough -> 8-WAY TOOM-COOK
                                toom8.multiply(a, b, r);
                        }
                }

	}

	/**
	 * Multiply 2 polynomial with integer coefficients.
	 *
	 * @param a First polynomial
	 * @param b Second polynomial
	 * @return The result polynomial
	 */
	UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
	
		UnivariateIntegerPolynomial r;

		int worker = __cilkrts_get_nworkers();

		if((a->getSize() < NAIVE_THRESHOLD || b->getSize() < NAIVE_THRESHOLD)){
                        // Too Small to use any fancy trick -> Apply NAIVE algorithm
                        r = naive.multiply(a, b);
                } else if((worker == 1) || (a->getSize() <= KS_THRESHOLD && b->getSize() <= KS_THRESHOLD)){
                        // Single-Core machine, or Small enough to use any parallel code -> KRONECKER SUBSTITUTION
                        r = ks.multiply(a, b);
                } else if(worker <= 6){
                        // Size is big enough, but the # of cores are low -> 4-WAY TOOM-COOK
                        r = toom4.multiply(a, b);
                } else {
                        // The # of cores are hight
                        if(a->getSize() >= CONV_THRESHOLD && b->getSize() >= CONV_THRESHOLD){
                                // Size is big enough -> 2-CONVOLUTION
                                r = ssa.multiply(a, b);
                        } else { 
                                // Size is not big enough -> 8-WAY TOOM-COOK
                                r = toom8.multiply(a, b);
                        }
                }

		r.computeExposedSize();
		return r;
	}
	
};

#endif
