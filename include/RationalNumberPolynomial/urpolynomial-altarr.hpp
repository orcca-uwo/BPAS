#include "polynomial.h"
#include "Interval/interval.h"
#include "RationalNumberPolynomial/urpolynomial.h"
#include "RingPolynomial/upolynomial.h"
#include "RationalNumberPolynomial/SUQP_Support.h"
//#include "SMQP_Support_Test-AA.h"
//#include "SMQP_Support_Recursive-AA.h"
#include <gmpxx.h>
#include "ExpressionTree/ExpressionTree.hpp"
#include "DataStructures/Factors.hpp"

class SparseUnivariateRationalPolynomial : public BPASUnivariatePolynomial<RationalNumber, SparseUnivariateRationalPolynomial> {
//private:
public: 

	AltArrU_t* poly;
	Symbol name; 
  SparseUnivariateRationalPolynomial(AltArrU_t * aa,Symbol s);

public:

/* Constructors */

		/**
		 * Construct a multivariate polynomial
		 *
		 **/
     SparseUnivariateRationalPolynomial();

		

		/**
		 * Construct with a variable name such that f(x) = x;
		 *
		 * @param x: The variable name
		 **/
		SparseUnivariateRationalPolynomial (const Symbol& x);

		/**
		 * Copy Constructor.
		 * 
		 * Does not reuse underlying memory allocated by b. 
		 *
		 * @param b: A sparse multivariate polynomial
		 **/
		SparseUnivariateRationalPolynomial(const SparseUnivariateRationalPolynomial& b);

		/**
		 * Move Constructor.
		 *
		 * @params b: The r-value reference polynomial.
		 */
		SparseUnivariateRationalPolynomial(SparseUnivariateRationalPolynomial&& b);
    ///////////////////////////////////////////////////////////////////////ExpressionTree convertToExpressionTree()
    /**
 * Construct an ExprTreeNode (well, multiple) which represents a single
 * term. That is, a coefficient and a monomial.
 */

ExprTreeNode* exprTreeNodeFromAAElem( AAElemU_t* n,  Symbol& sym) const;


ExpressionTree convertToExpressionTree () const; 

		/**
		 * Output the string representation of *this to the input ostream.
		 */
		void print(std::ostream& os) const;

std::string  polyToString() const;
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////


	/**
     * Determine if *this RationalNumber element is zero, that is the additive identity.
     *
     * returns true iff *this is zero.
     */
      bool isZero() const;

    /**
     * Make *this RationalNumber element zero.
     */
      void zero();

    /**
     * Determine if *this RationalNumber element is one, that is the multiplication identity.
     *
     * returns true iff *this is one.
     */
      bool isOne() const;
    
    /**
     * Make *this RationalNumber element one.
     */
      void one();
      /**
     * Determine if *this RationalNumber element is one, that is the multiplication identity.
     *
     * returns true iff *this is one.
     */
    int isConstant() const ;
    
        /**
		 * Negate all the coefficients of *this. Note, that due to the 
		 * sharing nature of underling nodes, this may alter the Nodes of
		 * other SUQP.
		 */
		void negate();
                                      
    /**
     * Obtain the unit normal (a.k.a canonical associate) of an element. 
     * If either parameters u, v, are non-NULL then the units are returned such that 
     * b = ua, v = u^-1. Where b is the unit normal of a, and is the returned value.
     */
    
      SparseUnivariateRationalPolynomial unitCanonical(SparseUnivariateRationalPolynomial* u , SparseUnivariateRationalPolynomial* v ) const;

    /**
     * Copy assignment.
     */
      SparseUnivariateRationalPolynomial& operator= (const SparseUnivariateRationalPolynomial&b);
    
    /**
     * Addition.
     */
      SparseUnivariateRationalPolynomial operator+ (const SparseUnivariateRationalPolynomial&b) const;

    /**
     * Addition assignment.
     */
      SparseUnivariateRationalPolynomial& operator+= (const SparseUnivariateRationalPolynomial&b);
    /**
     * Negation.
     */
      SparseUnivariateRationalPolynomial operator- () const;
    /**
     * Subtraction.
     */
      SparseUnivariateRationalPolynomial operator- (const SparseUnivariateRationalPolynomial&b) const;
    
    /**
     * Subtraction assignment.
     */
      SparseUnivariateRationalPolynomial& operator-= (const SparseUnivariateRationalPolynomial&b);

    

    /**
     * Multiplication.
     */
      SparseUnivariateRationalPolynomial operator* (const SparseUnivariateRationalPolynomial&b) const;

    /**
     * Multiplication assignment.
     */
      SparseUnivariateRationalPolynomial& operator*= (const SparseUnivariateRationalPolynomial&b);
    
    /**
     * Exponentiation.
     */
      SparseUnivariateRationalPolynomial operator^ (long long int e) const;
    
    /**
     * Exponentiation assignment.
     */
      SparseUnivariateRationalPolynomial& operator^= (long long int e);

    /**
     * Equality test,
     *
     * returns true iff equal
     */
      bool operator== (const SparseUnivariateRationalPolynomial&b) const;

    /**
     * Inequality test,
     *
     * returns true iff not equal.
     */
      bool operator!= (const SparseUnivariateRationalPolynomial&b) const;



             /**
             * Add *this and the RationalNumber r.
             */

           SparseUnivariateRationalPolynomial& operator= (const RationalNumber&r);

          /**
           * Add *this and the RationalNumber r.
           */
		   SparseUnivariateRationalPolynomial operator+ (const RationalNumber&c) const;

		   SparseUnivariateRationalPolynomial& operator+= (const RationalNumber&c);
         /**
           * Subtract *this and the RationalNumber r.
           */
		   SparseUnivariateRationalPolynomial operator- (const RationalNumber&c) const;
		   SparseUnivariateRationalPolynomial& operator-= (const RationalNumber&c);
  
       /**
           * Multiply *this and the RationalNumber r.
           */
		   SparseUnivariateRationalPolynomial operator* (const RationalNumber&c) const;
		   SparseUnivariateRationalPolynomial& operator*= (const RationalNumber&c);
        /**
           * devide *this and the RationalNumber r.
           */
            SparseUnivariateRationalPolynomial& operator/= (const RationalNumber&c);
           		   SparseUnivariateRationalPolynomial operator/ (const RationalNumber&c) const;
                 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		   
         

          /**
           *    total degree
           */
       Integer degree() const;

          /**
           * Leading Coefficient.
           */

		   RationalNumber leadingCoefficient() const; 
       /**
           * the least degree coeficient.
           */
		   RationalNumber trailingCoefficient() const; 
       /**
           * if the content term is zero.
           */
		   bool isConstantTermZero() const;

       /**
           * return the number of terms.
           */
		   Integer numberOfTerms() const;
       /**
           * return the content part of the polynomials.
           */
		   RationalNumber content() const;
      /**
           * return the content part of the polynomials.
           */
		   SparseUnivariateRationalPolynomial primitivePart() const;

          /**
           *  doing differentiate.
           */
         void differentiate(); // p = dp/dx

          /**
           *  doing differentiate.
           */
		   void differentiate(int k); 

        /**
           *  return the derivative .
           */
		   SparseUnivariateRationalPolynomial derivative() const; // q = dp/dx
         /**
           *  return the derivative by k times.
           */
		   SparseUnivariateRationalPolynomial derivative(int k) const;
//////////////////////////////////////////////////////////////////////////////////// integral /////////////
       SparseUnivariateRationalPolynomial Integral(int k) const; 

       
 	 void integrate( int k) ;

/////////////////////////////////////////////////////////////////////            Evaluation            ////////////////////////////////////////////////////////////

		    RationalNumber evaluate(const RationalNumber& r) const ; 
///////////////////////////////////////////////////////////////////       
		   SparseUnivariateRationalPolynomial monicDivide(const SparseUnivariateRationalPolynomial&b);
		   SparseUnivariateRationalPolynomial monicDivide(const SparseUnivariateRationalPolynomial&b, SparseUnivariateRationalPolynomial*rem) const ;
		   SparseUnivariateRationalPolynomial lazyPseudoDivide(const SparseUnivariateRationalPolynomial&b, RationalNumber*k, RationalNumber*mult);
      
	/**
	 * Lazy pseudo dividsion
	 * Return the quotient and itself becomes remainder
	 * e is the exact number of division steps
	 *
	 * @param b: The dividend polynomial
	 * @param c: The leading coefficient of b to the power e
	 * @param d: That to the power deg(a) - deg(b) + 1 - e
	 **/
/* 	inline SparseUnivariateRationalPolynomial lazyPseudoDivide(const SparseUnivariateRationalPolynomial& b, RationalNumber* c, RationalNumber* d) {
		if (d == NULL)
			d = new RationalNumber;
		Integer da = degree(), db = b.degree();
		if (b.isZero() || db == 0) {
			std::cout << "BPAS: error, dividend is zero or constant from SparseUnivariatePolynomial<Ring>." << std::endl;
			exit(1);
		}
		c->one(), d->one();
		if (isConstant()) {
			SparseUnivariateRationalPolynomial r;
			r.zero();
			return r;
		}
		if (name != b.name) {
			std::cout << "BPAS: error, trying to pseudo divide between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}

		if (da < db) {
			SparseUnivariateRationalPolynomial  r;
			r.name = name;
			return r;
		}

		SparseUnivariateRationalPolynomial  quo;
		quo.name = name;
		RationalNumber blc = b.poly->elems[0].coef;

		int e = 0;
		Integer diff = da - db;
		while (degree() >= db) {
			UnivariateTerm <RationalNumber> at = b.poly->elems[0].coef;

			UnivariateTerm<RationalNumber> nlc;
			nlc.coef = -at.coef;
			nlc.exp = at.exp - db.get_si();

			*c *= blc;
			e++;
			pomopo(blc, nlc, b);
			at.exp = nlc.exp;
			quo.terms.insert(quo.terms.begin(), at);
		}
		for (int i = e; diff >= i; ++i)
			*d *= blc;
		return quo;
	}  */

		   SparseUnivariateRationalPolynomial lazyPseudoDivide(const SparseUnivariateRationalPolynomial&b, SparseUnivariateRationalPolynomial*rem, RationalNumber*k, RationalNumber*mult) const;
		   SparseUnivariateRationalPolynomial pseudoDivide(const SparseUnivariateRationalPolynomial&b, RationalNumber*mult);
		   SparseUnivariateRationalPolynomial pseudoDivide(const SparseUnivariateRationalPolynomial&b, SparseUnivariateRationalPolynomial*rem, RationalNumber*mult) const;
		  
      //////////////////////////////////////////////////////////////////////////     find the coeficient of deg d    ////////////////////////////////////////////////////////
      
       RationalNumber coefficient(int d) const;
///////////////////////////////////////////////////////////////////////// set coef //////////////////////////////////////////////////////////////////////////////

			   void setCoefficient(int d , const RationalNumber&r) ;


       /////////////this->degree()/////////////////////////////////////                  set variable //////////////////////////////////////////////////////////////////////////////
		   void setVariableName (const Symbol&sym) ;
       /////////////////////////////////////////////////////////////////////////////////// varieble ///////////////////////////////////////////////////////

                   Symbol variable() const;  
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       ///////////////////////////////////////////////////////////////////////////////////
 SparseUnivariateRationalPolynomial operator<< (int i) const ; // q = p * (x^i);
	 SparseUnivariateRationalPolynomial& operator<<= (int i ) ;
		   SparseUnivariateRationalPolynomial operator>> (int i) const; // q = p / (x^i);
		   SparseUnivariateRationalPolynomial& operator>>= (int i);

////////////////////////////////////////////////////////////////////////////////////////////  devide apoly by aother poly


   SparseUnivariateRationalPolynomial operator/ (const SparseUnivariateRationalPolynomial&b) const;
	 SparseUnivariateRationalPolynomial& operator/= (const SparseUnivariateRationalPolynomial&b);

   
	 Factors<SparseUnivariateRationalPolynomial> squareFree() const; //TODO




   SparseUnivariateRationalPolynomial gcd(const SparseUnivariateRationalPolynomial&b) const;

SparseUnivariateRationalPolynomial ModularGCD_AAU(const SparseUnivariateRationalPolynomial&b) const;


   /////////////////////////////////////////////////// SUQP-Specific //////////////////////////////////////////////////

  /**
		 * Determine if *this is equal to b.
		 * This takes into account the variable ordering on both poylnomials 
		 * in such a way that the same polynomial under different variable orderings
		 * are NOT equal.
		 */

	bool isEqual(const SparseUnivariateRationalPolynomial &b) const ;
  ///////////////////////////////////////         Random poly             ///////////////////////////////////////////////////

 void randomPolynomial( time_t seed,int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) ;
 
 void randomPolynomialU( int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) ;
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   Random poly Max
 void buildRandomPolyFromMax_seeded(time_t seed,  const int maxDegs, unsigned long int coefBound, float sparsity, int includeNeg) ;
void buildRandomPolyFromMax_seededU( const int maxDegs, unsigned long int coefBound, float sparsity, int includeNeg);
};




