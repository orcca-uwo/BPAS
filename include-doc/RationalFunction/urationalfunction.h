#ifndef _URATIONALFUNCTION_H_
#define _URATIONALFUNCTION_H_

#include <complex>
#include <mps/mpc.h>
#include "../polynomial.h"
#include "rationalfunction_euclideanmethods.h"
#include "rationalfunction_symbolicintegration.h"
#include "multiprecision_rootfinding.h"
#include "rationalfunction_integrationpostprocessing.h"
#include "rationalfunction_symbolicnumericintegration.h"
#include "rationalfunction_integrationprinting.h"
#include "rationalfunction_complexrationalnumberordering.h"

struct IntegralTerm
{
	int Sindex;
	int tindex;
	
	bool operator()(const IntegralTerm a,const IntegralTerm b) const {
		if (a.Sindex > b.Sindex)
			return 1;
		else if (a.Sindex < b.Sindex)
			return 0;
		else {
			if (a.tindex > b.tindex)
				return 1;
			else
				return 0;
		}
	}
	
	bool operator<(const IntegralTerm b) const {
		if (this->Sindex < b.Sindex)
			return 1;
		else if (this->Sindex > b.Sindex)
			return 0;
		else {
			if (this->tindex < b.tindex)
				return 1;
			else
				return 0;
		}
	}
	
	inline friend std::ostream& operator<< (std::ostream &out, IntegralTerm b) {
			out << "[" << b.Sindex << "," << b.tindex << "]"; 
			return out;
		}
};

struct complexMPF
{
	mpc_t c;
};

typedef std::map<int, IntegralTerm> RootIntegralTermMap;
typedef RootIntegralTermMap::const_iterator RITMIter;
typedef std::multimap<IntegralTerm, int> IntegralTermRootMap;
typedef IntegralTermRootMap::const_iterator ITRMIter;
typedef std::multimap<ComplexRationalNumber, int, CompareByRealThenReverseImaginary> residueRootIndexMap;
typedef residueRootIndexMap::const_iterator RRIMIter;
typedef std::multimap<int, int> TermRootMap;
typedef TermRootMap::const_iterator TRMIter;

/**
 * A univariate rational function templated by a unvariate polynomial over a field.
 * The univariate polynomial and the coefficient BPASField must be passed separately 
 * and explicitly.  
 */
template <class UnivariatePolynomialOverField, class Field>
class UnivariateRationalFunction : public BPASRationalFunction<UnivariatePolynomialOverField> {
	private:
		UnivariatePolynomialOverField den;
		UnivariatePolynomialOverField num;
		
		bool PROFILING;
		bool ERROR_ANALYSIS;
		bool PFD_LOG_PART;
		bool floatingPointPrinting;
		std::string outputFormatting;

		inline void normalize () {
			num /= den.leadingCoefficient();
			den /= den.leadingCoefficient();
		}
	public:
		mpz_class characteristic;
		static bool isPrimeField;
		static bool isSmallPrimeField;
                static bool isComplexField;
		/**
		 * Construct the zero univariate rational function
		 *
		 * @param
		 **/ 
		UnivariateRationalFunction<UnivariatePolynomialOverField,Field> () {
		  den.one();
		  num.zero();
		  PROFILING = false;
		  ERROR_ANALYSIS = false;
		  PFD_LOG_PART = false;
		  floatingPointPrinting = false;
		  outputFormatting = "MAPLE_OUTPUT";
			UnivariatePolynomialOverField e;
			characteristic = e.characteristic;
		};

		/**
		 * Copy constructor
		 *
		 * @param b: A rational function
		 **/ 
		UnivariateRationalFunction<UnivariatePolynomialOverField,Field> (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) : den(b.den), num(b.num), PROFILING(b.PROFILING), 
ERROR_ANALYSIS(b.ERROR_ANALYSIS), PFD_LOG_PART(b.PFD_LOG_PART), floatingPointPrinting(b.floatingPointPrinting), outputFormatting(b.outputFormatting) {
			characteristic = b.characteristic;
		}
 
		/**
		 *
		 * @param a: the numerator
		 * @param b: the denominator
		 **/
		UnivariateRationalFunction<UnivariatePolynomialOverField,Field> (UnivariatePolynomialOverField a, UnivariatePolynomialOverField b) {			
		  if (a.variable() != b.variable()) {
		    std::cout << "BPAS error: numerator and denominator must have the same variable." << std::endl;
                    exit(1);
		  }
		  if (b.isZero()) {
		    std::cout << "BPAS error: denominator is zero from UnivariateRationalFunction<UnivariatePolynomialOverField,Field>" << std::endl;
		    exit(1);
		  }
		  num = a;
		  den = b;
		  PROFILING = false;
		  ERROR_ANALYSIS = false;
		  PFD_LOG_PART = false;
		  floatingPointPrinting = false;
		UnivariatePolynomialOverField e;
		characteristic = e.characteristic;
		  outputFormatting = "MAPLE_OUTPUT";
		}

		/**
		 * Destroy the rational function
		 *
		 * @param
		 **/
		~UnivariateRationalFunction<UnivariatePolynomialOverField,Field> () {}
		
		inline void setVariableName(Symbol name) {
			num.setVariableName(name);
			den.setVariableName(name);
		}
		
		inline Symbol variable() {
			return num.variable();
		}
		
		inline bool isProfiling() {
			return PROFILING;
		}
		
		inline void setProfiling(bool a) {
			PROFILING = a;
		}
		
		inline bool isAnalyzingError() {
			return ERROR_ANALYSIS;
		}
		
		inline void setAnalyzingError(bool a) {
			ERROR_ANALYSIS = a;
		}
		
		inline bool isPFDLogPart() {
			return PFD_LOG_PART;
		}
		
		inline void setPFDLogPart(bool a) {
			PFD_LOG_PART = a;
		}
		
		inline bool isFloatingPointPrinting() {
			return floatingPointPrinting;
		}
		
		inline void setFloatingPointPrinting(bool a) {
			floatingPointPrinting = a;
		}
		
		inline bool isMapleOutput() {
			if (outputFormatting == "MAPLE_OUTPUT")
				return true;
			else
				return false;
		}
		
		inline void setMapleOutput() {
			outputFormatting = "MAPLE_OUTPUT";
		}
		
		inline bool isMatlabOutput() {
			if (outputFormatting == "MATLAB_OUTPUT")
				return true;
			else
				return false;
		}
		
		inline void setMatlabOutput() {
			outputFormatting = "MATLAB_OUTPUT";
		}
		
		inline void setNumerator(const UnivariatePolynomialOverField& b) {			
			if (num.variable() != b.variable()) {
				std::cout << "BPAS error: numerator and denominator must have the same variable." << std::endl;
				exit(1);
			}
			num = b;
			canonicalize();
		}
		
		inline void setDenominator(const UnivariatePolynomialOverField& b) {			
			if (num.variable() != b.variable()) {
				std::cout << "BPAS error: numerator and denominator must have the same variable." << std::endl;
				exit(1);
			}
			den = b;
			canonicalize();
		}
		
		inline void set(const UnivariatePolynomialOverField& a, const UnivariatePolynomialOverField& b) {			
			if (a.variable() != b.variable()) {
				std::cout << "BPAS error: numerator and denominator must have the same variable." << std::endl;
				exit(1);
			}
			num = a;
			den = b;
			canonicalize();
		}
		
		inline UnivariatePolynomialOverField numerator() const {
			return num;
		}
		
		inline UnivariatePolynomialOverField denominator() const {
			return den;
		}
		
		inline Field evaluate(const Field& c) {
			Field output;
			output = num.evaluate(c);
			output /= den.evaluate(c);
			return output;
		}

		inline bool operator!= (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
		  return (!(num == b.num) || !(den == b.den));
		}
		inline bool operator== (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
		  return ((num == b.num) && (den == b.den));
		}

		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> operator+ (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
		  UnivariateRationalFunction<UnivariatePolynomialOverField,Field> r(*this);
		  return (r += b);
		}
		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& operator+= (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) {
		  if (num.variable()!= b.num.variable()) {
		    std::cout << "BPAS: error, trying to add between RationalFunction[" << num.variable() << "] and RationalFunction[" << b.num.variable() << "]." << std::endl;
		    exit(1);
		  }
		  UnivariatePolynomialOverField g, r1(den), r2(b.den);
		  g = den.gcd(b.den);
		  r1 /= g;
		  r2 /= g;
		  den = r1 * r2;
		  r1 *= b.num;
		  r2 *= num;
		  r1 += r2; // r1 = ar2 + cr1;
		  r2 = r1.gcd(g);
		  r1 /= r2; // r1 = e;
		  g /= r2; // g = g';
		  den *= g;
		  num = r1;
		  normalize();
		  
		  return *this;
		}

		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> operator- (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
		  UnivariateRationalFunction<UnivariatePolynomialOverField,Field> r(*this);
		  return (r -= b);
		}
		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& operator-= (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) {
		  if (num.variable()!= b.num.variable()) {
		    std::cout << "BPAS: error, trying to subtract between RationalFunction[" << num.variable() << "] and RationalFunction[" << b.num.variable() << "]." << std::endl;
		    exit(1);
		  }
		  UnivariatePolynomialOverField g, r1(den), r2(b.den);
		  g = den.gcd(b.den);
		  r1 /= g;
		  r2 /= g;
		  den = r1 * r2;
		  r1 *= -b.num;
		  r2 *= num;
		  r1 += r2; // r1 = ar2 + cr1;
		  r2 = r1.gcd(g);
		  r1 /= r2; // r1 = e;
		  g /= r2; // g = g';
		  den *= g;
		  num = r1;
		  normalize();
 
		  return *this;
		}

		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> operator- () const {
		  UnivariateRationalFunction<UnivariatePolynomialOverField,Field> r (-num, den);
		  return r;
		}

		/**
         * Overload operator ^
         * replace xor operation by exponentiation
         *
         * @param e: The exponentiation, e > 0
         **/
        inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> operator^ (long long int e) const {
                UnivariateRationalFunction<UnivariatePolynomialOverField,Field> res;
                res.setVariableName(num.variable());
                if (isZero() || isOne() || e == 1)
                        res = *this;
                else if (e == 2)
                        res = *this * *this;
                else if (e > 2) {
                        UnivariateRationalFunction<UnivariatePolynomialOverField,Field> x (*this);
                        res.one();

                        while (e) {
                                if (e % 2) { res *= x; }
                                x = x * x;
                                e >>= 1;
                        }
                }
                else if (!e)
                        res.one();
                else {
                        res = *this;
                }
                return res;
        }

		/**
         * Overload operator ^=
         * replace xor operation by exponentiation
         *
         * @param e: The exponentiation, e > 0
         **/
        inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& operator^= (long long int e) {
                *this = *this ^ e;
                return *this;
        }

		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> inverse() const {
		  if (num.isZero()) {
		    std::cout << "BPAS error: division by zero from UnivariateRationalFunction<UnivariatePolynomialOverField,Field> inverse()" << std::endl;
		    exit(1);
		  }
		  UnivariateRationalFunction<UnivariatePolynomialOverField,Field> r(den, num);
		  r.normalize();
		  return r;
		}
		
		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> operator* (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
		  UnivariateRationalFunction<UnivariatePolynomialOverField,Field> r(*this);
		  return (r *= b);
		}
		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& operator*= (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) {
		  if (num.variable()!= b.num.variable()) {
		    std::cout << "BPAS: error, trying to multiply between RationalFunction[" << num.variable() << "] and RationalFunction[" << b.num.variable() << "]." << std::endl;
		    exit(1);
		  }

		  UnivariatePolynomialOverField g1, g2, r;
		  g1 = num.gcd(b.den);
		  g2 = den.gcd(b.num);
		  num /= g1;
		  r = b.num / g2;
		  num *= r;
		  den /= g2;
		  r = b.den / g1;
		  den *= r;
		  normalize();
		  return *this;
		}

		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> operator/ (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
		  UnivariateRationalFunction<UnivariatePolynomialOverField,Field> r(*this);
		  return (r /= b);
		}
		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& operator/= (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) {
		  if (b.isZero()) {
		    std::cout << "BPAS error: division by zero from UnivariateRationalFunction<UnivariatePolynomialOverField,Field> operator/=" << std::endl;
		    exit(1);
		  }

		  if (num.variable()!= b.num.variable()) {
		    std::cout << "BPAS: error, trying to divide between RationalFunction[" << num.variable() << "] and RationalFunction[" << b.num.variable() << "]." << std::endl;
		    exit(1);
		  }
		  UnivariateRationalFunction<UnivariatePolynomialOverField,Field> e (b.den, b.num);
		    *this *= e;
		  return *this;
		}

		inline void canonicalize() {
			UnivariatePolynomialOverField temp;
			temp = num.gcd(den);
			num /= temp;
			den /= temp;
			num /= den.leadingCoefficient();
			den /= den.leadingCoefficient();
			//UnivariateRationalFunction<UnivariatePolynomialOverField,Field> ret;
			//ret = this->unitCanonical();
			//*this = ret;
		}

		inline bool isZero() const {
			return num.isZero();
		}
		inline void zero() {
			num.zero();
			den.one();
		}
		inline bool isOne() const {
			return num.isOne() && den.isOne();
		}
		inline void one() {
			num.one();
			den.one();
		}
		inline bool isNegativeOne() const {
			return (num.isNegativeOne() && den.isOne()) || (num.isOne() && den.isNegativeOne());
		}
		inline void negativeOne() {
			num.negativeOne();
			den.one();
		}
		inline int isConstant() const {
			return num.isConstant() && den.isConstant();
		}

		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field> unitCanonical(UnivariateRationalFunction<UnivariatePolynomialOverField,Field>* u = NULL, UnivariateRationalFunction<UnivariatePolynomialOverField,Field>* v = NULL) const {
			UnivariatePolynomialOverField du,dc,dv,g,temp;
			g = num.gcd(den);
			temp = den/g;
			dc = temp.unitCanonical(&du,&dv);
			temp = du*(num/g);
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> ret;
			ret.num = temp;
			ret.den = dc;	
			if (u != NULL || v!= NULL) {
				UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp2;
				temp2.one();
				if (u != NULL) {
					*u = temp2;
				}
				if (v != NULL) {
					*v = temp2;
				}
			}

			return ret;
		}	

		/**
		 * Overload operator =
		 *
		 * @param b: A rational function
		 **/
		inline UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& operator= (const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) {
			if (this != &b) {
			  num = b.num;
			  den = b.den;
			characteristic = b.characteristic;
			}
			return *this;
		}

		ExpressionTree convertToExpressionTree() const {
			std::cerr << "UnivariateRationalFunction::convertToExpressionTree NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return ExpressionTree();
		}

		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param b: The rational function
		 **/
    	void print(std::ostream& ostream) const {
			ostream << "(" << num << ")/(" << den << ")"; 
    	}

    	/*UnivariateRationalFunction<UnivariatePolynomialOverField, Field> unitCanonical(UnivariateRationalFunction<UnivariatePolynomialOverField,Field>* u = NULL,
    																				   UnivariateRationalFunction<UnivariatePolynomialOverField,Field>* v = NULL) const {
    		UnivariatePolynomialOverField temp;
			temp = num.gcd(den);
			UnivariatePolynomialOverField tempDen = den / temp;
			Field lc = tempDen.leadingCoefficient();
			temp *= lc;
			if (u != NULL) {
				//TODO 
				*u = UnivariateRationalFunction<UnivariatePolynomialOverField,Field>();
			}
			if (v != NULL) {
				*v = UnivariateRationalFunction<UnivariatePolynomialOverField,Field>(temp, temp);
			}

			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> ret = *this;
			ret.canonicalize();
			return ret;
    	}*/

    	
    	/** BPASGCDDomain, BPASEuclideanDomain, BPASField virtual methods **/

    	/**
    	 * Get this GCD of *this and b.
    	 */
    	UnivariateRationalFunction<UnivariatePolynomialOverField, Field> gcd(const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
    		std::cerr << "UnivariateRationalFunction::gcd NOT YET IMPLEMENTED" << std::endl;
    		//TODO
    		return *this;
    	}
		
		/**
		 * Compute squarefree factorization of *this
		 */
		inline Factors<UnivariateRationalFunction> squareFree() const {
			std::cerr << "UnivariateRationalFunction::squareFree NOT YET IMPLEMENTED" << std::endl;
			//TODO
			std::vector<UnivariateRationalFunction> ret;
			ret.push_back(*this);
			return ret;
		}

    	UnivariateRationalFunction<UnivariatePolynomialOverField, Field> euclideanSize() const {
			std::cerr << "UnivariateRationalFunction::euclideanSize NOT YET IMPLEMENTED" << std::endl;
    		//TODO
    		return *this;

    	}
 
    	UnivariateRationalFunction<UnivariatePolynomialOverField, Field> euclideanDivision(const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b, 
    																							 UnivariateRationalFunction<UnivariatePolynomialOverField,Field>* q = NULL) const {
			std::cerr << "UnivariateRationalFunction::euclideanDivision NOT YET IMPLEMENTED" << std::endl;
    		//TODO
    		return *this;

    	}

    	UnivariateRationalFunction<UnivariatePolynomialOverField, Field> quotient(const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
			return (*this / b);
    	}

    	UnivariateRationalFunction<UnivariatePolynomialOverField, Field> remainder(const UnivariateRationalFunction<UnivariatePolynomialOverField,Field>& b) const {
			return UnivariateRationalFunction<UnivariatePolynomialOverField, Field>(0,1);
    	}

    	UnivariateRationalFunction<UnivariatePolynomialOverField, Field> extendedEuclidean(const UnivariateRationalFunction<UnivariatePolynomialOverField, Field>& b, 
																						     	 UnivariateRationalFunction<UnivariatePolynomialOverField, Field>* s = NULL,
																						     	 UnivariateRationalFunction<UnivariatePolynomialOverField, Field>* t = NULL) const {
    		std::cerr << "UnivariateRationalFunction::extendedEuclidean NOT YET IMPLEMENTED" << std::endl;
    		//TODO
    		return *this;
    	}

    	inline UnivariateRationalFunction<UnivariatePolynomialOverField, Field> operator%(const UnivariateRationalFunction<UnivariatePolynomialOverField, Field>& b) const {
    		UnivariateRationalFunction<UnivariatePolynomialOverField, Field> ret = this->remainder(b);
    		return ret;
    	}

    	inline UnivariateRationalFunction<UnivariatePolynomialOverField, Field>& operator%=(const UnivariateRationalFunction<UnivariatePolynomialOverField, Field>& b) {
    		*this = this->remainder(b);
    	}


		void hermiteReduce(std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > *g, UnivariateRationalFunction<UnivariatePolynomialOverField,Field> *h) {
			std::vector<UnivariatePolynomialOverField> gg,hh;
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp;
			
			temp.setVariableName(num.variable());
			h->setVariableName(num.variable());
			g->clear();
			_hermiteReduce<UnivariatePolynomialOverField,Field>(num,den,&gg,&hh);
			int i(0);
			while (i<gg.size()) {
				temp.set(gg.at(i),gg.at(i+1));
				g->push_back(temp);
				i += 2;
			}
			temp.set(hh.at(0),hh.at(1));
			*h = temp;
		}
		
		void integrateRationalFunctionLogPart(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, std::vector<UnivariatePolynomialOverField> *U) {
			U->clear();
			S->clear();
		
			_integrateRationalFunctionLogPart<UnivariatePolynomialOverField,Field>(S,U,num,den,PROFILING);
		}
		
		void differentiate() {
			/* Destructive rational function differentiation */
			// input a/d
			// output (a'*d-a*d')/d^2 = a'*e-a*f/d*e; e=d/g, f=d'/g, g=gcd(d,d')
	
			UnivariatePolynomialOverField D(den); // D = d	
			UnivariatePolynomialOverField dD(den);
			UnivariatePolynomialOverField temp;
			//std::cout << "here?" << std::endl;
			dD.differentiate(1);	// dD = d'
			//std::cout << "dD = " << dD << std::endl;
			temp = D.gcd(dD);		// temp = g
			//std::cout << "temp = " << temp << std::endl;
			D /= temp;				// D = e
			//std::cout << "D = " << D << std::endl;
			dD /= temp;				// dD = f
			//std::cout << "dD = " << dD << std::endl;
			temp = -num;			// temp = -a
			//std::cout << "temp = " << temp << std::endl;
			temp *= dD;				// temp = -a*f
			//std::cout << "temp = " << temp << std::endl;
			dD = num;				// dD = a
			//std::cout << "dD = " << dD << std::endl;
			dD.differentiate(1);	// dD = a'
			//std::cout << "dD = " << dD << std::endl;
			dD *= D;				// dD = a'*e
			//std::cout << "dD = " << dD << std::endl;
			temp += dD;				// temp = a'*e-a*f
			//std::cout << "temp = " << temp << std::endl;
			D *= den; 				// D = d*e
			//std::cout << "D = " << D << std::endl;
			
			//std::cout << "here?" << std::endl;
			num = temp;
			den = D;
			canonicalize();
		}
		
		void integrate(UnivariatePolynomialOverField *P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > *g, std::vector<UnivariatePolynomialOverField> *U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S) {
			g->clear();
			U->clear();
			S->clear();
			std::vector<UnivariatePolynomialOverField> G;
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp;
			temp.setVariableName(num.variable());	
	
			// Profiling variables
			unsigned long long start;
			float elapsed = 0;
	
			if (PROFILING){	
				std::cout << "integrate" << std::endl;
				std::cout << "--------------------------------------" << std::endl;
				startTimer(&start);
			}
					
			_integrateRationalFunction<UnivariatePolynomialOverField,Field>(num,den,P,&G,U,S,PROFILING);
	
			if (PROFILING){	
				stopTimer(&start,&elapsed);
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "integrate runtime: " << elapsed << " s" << std::endl;
			}
				
			int i(0);
			while (i<G.size()) {
				temp.set(G.at(i),G.at(i+1));
				g->push_back(temp);
				i += 2;
			}
			
		}
		
		void realSymbolicNumericIntegrate(UnivariatePolynomialOverField *P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > *g, std::vector<Field> *lg, std::vector<UnivariatePolynomialOverField> *Lg, std::vector<Field> *atn, std::vector<UnivariatePolynomialOverField> *Atn, int prec) {
			P->zero();
			g->clear();
			lg->clear();
			Lg->clear();
			atn->clear();
			Atn->clear();
			std::vector<UnivariatePolynomialOverField> G;
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp;
			temp.setVariableName(num.variable());
	
		    std::cout << "[realSymbolicNumericIntegrate (snInt): Symbolic-Numeric Integration with BPAS and MPSolve]" << std::endl;
		    std::cout << "[Integration method: Hermite reduction, LRT integration]" << std::endl;
		    std::cout << "Starting..." << std::endl;

			// Profiling variables
			unsigned long long start;
			float elapsed = 0;
	
			if (PROFILING){	
				std::cout << "--------------------------------------" << std::endl;
				startTimer(&start);
			}
		
			_realSNIntegrate<UnivariatePolynomialOverField,Field>(num,den,P,&G,lg,Lg,atn,Atn,prec,PROFILING,PFD_LOG_PART,ERROR_ANALYSIS);
	
			if (PROFILING){	
				stopTimer(&start,&elapsed);
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "realSymbolicNumericIntegrate runtime: " << elapsed << " s" << std::endl;
			}	
		
			int i(0);
			while (i<G.size()) {
				temp.set(G.at(i),G.at(i+1));
				g->push_back(temp);
				i += 2;
			}
		
		}
		
		void realSymbolicNumericIntegrate(UnivariatePolynomialOverField *P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > *g, std::vector<Field> *lg, std::vector<UnivariatePolynomialOverField> *Lg, std::vector<Field> *atn, std::vector<UnivariatePolynomialOverField> *Atn1, std::vector<UnivariatePolynomialOverField> *Atn2, int prec) {
			P->zero();
			g->clear();
			lg->clear();
			Lg->clear();
			atn->clear();
			Atn1->clear();
			Atn2->clear();
			std::vector<UnivariatePolynomialOverField> G;
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp;
			temp.setVariableName(num.variable());
	
		    std::cout << "[realSymbolicNumericIntegrate (snInt): Symbolic-Numeric Integration with BPAS and MPSolve]" << std::endl;
		    std::cout << "[Integration method: Hermite reduction, LRT integration]" << std::endl;
		    std::cout << "Starting..." << std::endl;

			// Profiling variables
			unsigned long long start;
			float elapsed = 0;
	
			if (PROFILING){	
				std::cout << "--------------------------------------" << std::endl;
				startTimer(&start);
			}
		
			_realSNIntegrate<UnivariatePolynomialOverField,Field>(num,den,P,&G,lg,Lg,atn,Atn1,Atn2,prec,PROFILING,PFD_LOG_PART,ERROR_ANALYSIS);
	
			if (PROFILING){	
				stopTimer(&start,&elapsed);
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "realSymbolicNumericIntegrate runtime: " << elapsed << " s" << std::endl;
			}	
		
			int i(0);
			while (i<G.size()) {
				temp.set(G.at(i),G.at(i+1));
				g->push_back(temp);
				i += 2;
			}
		
		}
		
		void realSymbolicNumericIntegratePFD(UnivariatePolynomialOverField *P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > *g, std::vector<Field> *lg, std::vector<UnivariatePolynomialOverField> *Lg, std::vector<Field> *atn, std::vector<UnivariatePolynomialOverField> *Atn, int prec) {
			P->zero();
			g->clear();
			lg->clear();
			Lg->clear();
			atn->clear();
			Atn->clear();
			std::vector<UnivariatePolynomialOverField> G;
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp;
			temp.setVariableName(num.variable());
	
		    std::cout << "[realSymbolicNumericIntegratePFD (snIntPFD): Symbolic-Numeric Integration with BPAS and MPSolve]" << std::endl;
		    std::cout << "[Integration method: Hermite reduction, PFD integration]" << std::endl;
		    std::cout << "Starting..." << std::endl;

			// Profiling variables
			unsigned long long start;
			float elapsed = 0;
	
			if (PROFILING){	
				std::cout << "--------------------------------------" << std::endl;
				startTimer(&start);
			}
		
			_realSNIntegratePFD<UnivariatePolynomialOverField,Field>(num,den,P,&G,lg,Lg,atn,Atn,prec,PROFILING,PFD_LOG_PART,ERROR_ANALYSIS);
	
			if (PROFILING){	
				stopTimer(&start,&elapsed);
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "realSymbolicNumericIntegratePFD runtime: " << elapsed << " s" << std::endl;
			}	
		
			int i(0);
			while (i<G.size()) {
				temp.set(G.at(i),G.at(i+1));
				g->push_back(temp);
				i += 2;
			}
		
		}
		
		void realSymbolicNumericIntegrateSimplePFD(UnivariatePolynomialOverField *P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > *g, std::vector<Field> *lg, std::vector<UnivariatePolynomialOverField> *Lg, std::vector<Field> *atn, std::vector<UnivariatePolynomialOverField> *Atn, int prec) {
			P->zero();
			g->clear();
			lg->clear();
			Lg->clear();
			atn->clear();
			Atn->clear();
			std::vector<UnivariatePolynomialOverField> G;
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp;
			temp.setVariableName(num.variable());
	
		    std::cout << "[realSymbolicNumericIntegratePFD (snIntPFD): Symbolic-Numeric Integration with BPAS and MPSolve]" << std::endl;
		    std::cout << "[Integration method: Hermite reduction, PFD integration]" << std::endl;
		    std::cout << "Starting..." << std::endl;

			// Profiling variables
			unsigned long long start;
			float elapsed = 0;
	
			if (PROFILING){	
				std::cout << "--------------------------------------" << std::endl;
				startTimer(&start);
			}
		
			_realSNIntegrateSimplePFD<UnivariatePolynomialOverField,Field>(num,den,P,&G,lg,Lg,atn,Atn,prec,PROFILING,PFD_LOG_PART,ERROR_ANALYSIS);
	
			if (PROFILING){	
				stopTimer(&start,&elapsed);
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "realSymbolicNumericIntegratePFD runtime: " << elapsed << " s" << std::endl;
			}	
		
			int i(0);
			while (i<G.size()) {
				temp.set(G.at(i),G.at(i+1));
				g->push_back(temp);
				i += 2;
			}
		
		}
		
		/*void realSymbolicNumericIntegratePFD(UnivariatePolynomialOverField *P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > *g, std::vector<Field> *lg, std::vector<UnivariatePolynomialOverField> *Lg, std::vector<Field> *atn, std::vector<UnivariatePolynomialOverField> *Atn1, std::vector<UnivariatePolynomialOverField> *Atn2, int prec) {
			P->zero();
			g->clear();
			lg->clear();
			Lg->clear();
			atn->clear();
			Atn1->clear();
			Atn2->clear();
			std::vector<UnivariatePolynomialOverField> G;
			UnivariateRationalFunction<UnivariatePolynomialOverField,Field> temp;
			temp.setVariableName(num.variable());
	
		    std::cout << "[realSymbolicNumericIntegratePFD (snIntPFD): Symbolic-Numeric Integration with BPAS and MPSolve]" << std::endl;
		    std::cout << "[Integration method: Hermite reduction, PFD integration]" << std::endl;
		    std::cout << "Starting..." << std::endl;

			// Profiling variables
			unsigned long long start;
			float elapsed = 0;
	
			if (PROFILING){	
				std::cout << "--------------------------------------" << std::endl;
				startTimer(&start);
			}
		
			_realSNIntegratePFD<UnivariatePolynomialOverField,Field>(num,den,P,&G,lg,Lg,atn,Atn1,Atn2,prec,PROFILING,PFD_LOG_PART,ERROR_ANALYSIS);
	
			if (PROFILING){	
				stopTimer(&start,&elapsed);
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "realSymbolicNumericIntegratePFD runtime: " << elapsed << " s" << std::endl;
			}	
		
			int i(0);
			while (i<G.size()) {
				temp.set(G.at(i),G.at(i+1));
				g->push_back(temp);
				i += 2;
			}
		
		}*/

		void printIntegral(UnivariatePolynomialOverField &P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > &g, std::vector<UnivariatePolynomialOverField> &U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > &S){
			std::vector<UnivariatePolynomialOverField> G;
			for (int i=0; i<g.size(); i++) {
				G.push_back(g.at(i).num);
				G.push_back(g.at(i).den);
			}
			_printFormalIntegral<UnivariatePolynomialOverField,Field>(num,den,P,G,U,S, false, floatingPointPrinting, false);
		}
		
		void printIntegral(UnivariatePolynomialOverField &P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > &g, std::vector<Field> &lg, std::vector<UnivariatePolynomialOverField> &Lg, std::vector<Field> &atn, std::vector<UnivariatePolynomialOverField> &Atn){
			std::vector<UnivariatePolynomialOverField> G;
			for (int i=0; i<g.size(); i++) {
				G.push_back(g.at(i).num);
				G.push_back(g.at(i).den);
			}
			std::vector<UnivariatePolynomialOverField> empty;
			_printIntegral<UnivariatePolynomialOverField,Field>(num,den,P,G,lg,Lg,atn,Atn,empty, false, floatingPointPrinting, false, outputFormatting);
		}
		
		void printIntegral(UnivariatePolynomialOverField &P, std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > &g, std::vector<Field> &lg, std::vector<UnivariatePolynomialOverField> &Lg, std::vector<Field> &atn, std::vector<UnivariatePolynomialOverField> &Atn1, std::vector<UnivariatePolynomialOverField> &Atn2){
			std::vector<UnivariatePolynomialOverField> G;
			for (int i=0; i<g.size(); i++) {
				G.push_back(g.at(i).num);
				G.push_back(g.at(i).den);
			}
			_printIntegral<UnivariatePolynomialOverField,Field>(num,den,P,G,lg,Lg,atn,Atn1,Atn2, false, floatingPointPrinting, false, outputFormatting);
		}

		void realSymbolicNumericIntegrate(int prec) {
                        UnivariatePolynomialOverField P;
                        std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > g;
                        std::vector<Field> lg, atn;
                        std::vector<UnivariatePolynomialOverField> Lg, Atn1, Atn2;

                        realSymbolicNumericIntegrate(&P, &g, &lg, &Lg, &atn, &Atn1, &Atn2, prec);
                        printIntegral(P, g, lg, Lg, atn, Atn1, Atn2);
                }

		void integrate() {
                        UnivariatePolynomialOverField P;
                        std::vector< UnivariateRationalFunction<UnivariatePolynomialOverField,Field> > g;
						std::vector<UnivariatePolynomialOverField> U;
						std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > S;

                        integrate(&P, &g, &U, &S);
                        printIntegral(P, g, U, S);
                }
};
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


