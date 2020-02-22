#ifndef _TRIANGULARSET_H_
#define _TRIANGULARSET_H_

//#include "Polynomial/mpolynomial.h"
#include "../polynomial.h"
#include "../ring.h"
#include "../Utils/TemplateHelpers.hpp"
#include "../Utils/SymbolHelpers.hpp"
//#include "../RegularChain/regularchain.hpp"
//#include "../RegularChain/zerodimensionalregularchain.hpp"
#include <iostream>

using std::endl;
using std::cerr;

//template <class Field, class RecursivePoly> RegularChain;

/** 
 * An enumeration which describes the mode of the triangular set.
 * TS_FIXED is the ALGEBRAIC case, forces a fixed list of variables, determinining the maximum size of the TS
 * TS_VARIABLE is the DIFFERENTIAL case, allows the list of variables, hence size of the TS, to change  
 */
typedef enum TriangularSetMode {
	TS_FIXED = 0x0,
	TS_VARIABLE = 0x1
} TriangularSetMode;

/** 
 * A triangular set templated by a multivariate polynomial over a field.
 * The field should be a BPASField and the multivariate polynomial should
 * be recursively viewed, as in BPASRecursivelyViewedPolynomial.
 */
template <class Field, class RecursivePoly>
class TriangularSet : public virtual BPASTriangularSet<Field,RecursivePoly>
{

	protected:
		TriangularSetMode mode;				// fixed or variable TS structure
		mpz_class characteristic;			// characteristic of the TS (inherited from base Field)
		std::vector<RecursivePoly> set;		// the set of polynomials
		std::vector<Symbol> vars;			// list of active and inactive algebraic variables
		std::vector<Symbol> algVars;		// list of active algebraic variables
		std::vector<Symbol> trcVars;		// list of transcendental variables (cannot appear in polynomials)
		Symbol sNMaxVar;					// max variable of largest strongly normalized subset of the TS
		// TODO: remove stronglyNormalized attribute and compute it when needed using sNMaxVar
		bool stronglyNormalized;			// flag whether the TS is strongly normalized
		
//		/**
//		 * clear: empties the triangular set, preserving the list of variables and transcendentals
//		 **/
//		void clear();
//		
//		/**
//		 * erase: converts the triangular set into the result of calling the default constructor
//		 **/
//		void erase();
		
		void updateTriangularSetStates();
		void updateTriangularSetStates(const RecursivePoly& p);
		int variableIndex(const Symbol& s) const;


	public:
		/**
		 * Default constructor: creates an empty triangular set of variable size
		 * with empty list of transcendentals
		 *
		 * @param
		 **/
		TriangularSet<Field,RecursivePoly> ();
		
		/**
		 * Construct an empty triangular set of fixed size in the s decreasingly ordered 
		 * variables given by xs with empty list of transcendentals
		 *
		 * @param s: The number of polynomials
		 * @param xs: The variable names
		 **/
		TriangularSet<Field,RecursivePoly> (const std::vector<Symbol>& xs);
		
		/**
		 * Construct an empty triangular set of fixed size in the s decreasingly ordered 
		 * variables given by xs and list of transcendentals given by ts
		 *
		 * @param s: The number of polynomials
		 * @param xs: The variable names
		 * @param ts: The transcendental variable names
		 **/
		TriangularSet<Field,RecursivePoly> (const std::vector<Symbol>& xs, const std::vector<Symbol>& ts);
		
		/**
		 * Construct a variable triangular set containing p, such that the variables of p are treated as algebraic,
		 * with empty list of transcendentals
		 *
		 * @param p: The polynomial to add
		 **/
		TriangularSet<Field,RecursivePoly> (const RecursivePoly& p);
		
		/**
		 * Construct a variable triangular set containing p, such that the variables in ts are
		 * treated as transcendental, while any remaining variables of p are treated as algebraic
		 *
		 * @param p: The polynomial to add
		 * @param ts: The transcendental variable names
		 **/
		TriangularSet<Field,RecursivePoly> (const RecursivePoly& p, const std::vector<Symbol>& ts);
		
		/**
		 * Copy constructor
		 *
		 * @param a: A triangular set
		 **/
		TriangularSet<Field,RecursivePoly> (const TriangularSet<Field,RecursivePoly>& a);
		
//		/**
//		 * Copy constructor
//		 *
//		 * @param a: A regular chain
//		 **/
//		TriangularSet<Field,RecursivePoly> (const RegularChain<Field,RecursivePoly>& a);
		
//		/**
//		 * Copy constructor
//		 *
//		 * @param a: A zero dimensional regular chain
//		 **/
//		TriangularSet<Field,RecursivePoly> (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Move constructor
		 *
		 * @param a: An r-value reference triangular set
		 **/
		TriangularSet<Field,RecursivePoly> (TriangularSet<Field,RecursivePoly>&& a);
		
//		/**
//		 * Move constructor
//		 *
//		 * @param a: An r-value reference regular chain
//		 **/
//		TriangularSet<Field,RecursivePoly> (RegularChain<Field,RecursivePoly>&& a);
		
//		/**
//		 * Move constructor
//		 *
//		 * @param a: An r-value reference zero dimensional regular chain
//		 **/
//		TriangularSet<Field,RecursivePoly> (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);
		
		/**
		 * Computational constructor: creates a triangular set given all the data
		 *
		 * @param vs: variables of the triangular set
		 * @param avs: algebraic variables of the triangular set
		 * @param tvs: transcendental variables of the triangular set
		 * @param polys: polynomials of the triangular set
		 * @param tsm: whether the triangular set is variable or fixed
		 * @param c: characteristic of the triangular set
		 **/
		TriangularSet<Field,RecursivePoly> (const std::vector<Symbol>& vs, const std::vector<Symbol>& avs, const std::vector<Symbol>& tvs, const std::vector<RecursivePoly>& ts, TriangularSetMode tsm, const mpz_class& c);
		
		/**
		 * Deconstructor
		 *
		 * @param
		 **/
		~TriangularSet<Field,RecursivePoly>();

//		/**
//		 * Copy an object derived from abstract BPASTriangularSet class to type of current object
//		 *
//		 * @param ts: triangular set to copy
//		 **/
//		inline void copy(const BPASTriangularSet<Field,RecursivePoly>& ts) override {
//			if (dynamic_cast<const TriangularSet<Field,RecursivePoly>*>(&ts))
//				*this = *dynamic_cast<const TriangularSet<Field,RecursivePoly>*>(&ts);
//			else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to TriangularSet."));
//		}
		
		/**
		 * Tests if the TriangularSet is empty
		 *
		 * @param
		 **/
		bool isEmpty() const;
		
		/**
		 * Tests if the polynomial is constant relative to the TriangularSet
		 *
		 * @param
		 **/
		bool isConstantPolynomial(const RecursivePoly& p) const;
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A triangular set
		 **/
		TriangularSet<Field,RecursivePoly>& operator= (const TriangularSet<Field,RecursivePoly>& a);
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (const BPASTriangularSet<Field,RecursivePoly>& a) override;
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A triangular set
		 **/
		TriangularSet<Field,RecursivePoly>& operator= (TriangularSet<Field,RecursivePoly>&& a);
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (BPASTriangularSet<Field,RecursivePoly>&& a) override;
		
		/**
		 * Add operator +
		 * Adds a polynomial to a triangular set and returns a new triangular set
		 *
		 * @param p: A sparse multivariate polynomial
		 **/
		TriangularSet<Field,RecursivePoly> operator+ (const RecursivePoly& p);
		
		/**
		 * Add assignment operator +=
		 * Adds a polynomial to a triangular set
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		TriangularSet<Field,RecursivePoly>& operator+= (const RecursivePoly& p);
		
		/**
		 * Identity operator ==
		 *
		 *
		 * @param a: A triangular set
		 **/
		bool operator== (TriangularSet<Field,RecursivePoly>& a);

		/**
		 * Negated identity operator !=
		 *
		 *
		 * @param a: A triangular set
		 **/
		bool operator!= (TriangularSet<Field,RecursivePoly>& a);
		
		/**
		 * Get the number of variables
		 *
		 **/
		inline int numberOfVariables() const {
			return vars.size();
		}
		
		/**
		 * Get the size of the triangular set
		 *
		 **/
		inline int size() const {
			return algVars.size();
		}
		
		/**
		 * Get the number of algebraic variables
		 *
		 **/
		inline int numberOfAlgebraicVariables() const {
			return algVars.size();
		}
		
		/**
		 * Get the number of transcendental variables
		 *
		 * @param
		 **/
		inline int numberOfTranscendentalVariables() const {
			return trcVars.size();
		}
		
		/**
		 * Get the variable names in decreasing order
		 *
		 **/
		inline std::vector<Symbol> variables() const {
			return vars;
		}
		
		/**
		 * Get algebraic variables
		 *
		 **/
		inline std::vector<Symbol> mainVariables() const {
			return algVars;
		}
		
		/**
		 * Get transcendentalVariables variables
		 *
		 **/
		inline std::vector<Symbol> transcendentalVariables() const {
			return trcVars;
		}
		
		/**
		 * Get list of variables followed by the transcendental variables
		 *
		 **/
		std::vector<Symbol> allVariables() const;
		
		inline bool isAlgebraic(const Symbol& s) const {
			if (std::find(algVars.begin(),algVars.end(),s) != algVars.end())
				return true;
			else
				return false;
		}
		
		inline bool isStronglyNormalized() const {
			return stronglyNormalized;
		}
		
		inline std::vector<RecursivePoly> polynomials() const {
			return set;
		}
		
		/**
		 * The dimension of the triangular set (understood in terms of the space of potentially algebraic variables).
		 **/
		inline int dimension() const {
			return (vars.size() - algVars.size());
		}
		
		/**
		 * The dimension of the triangular set lower(v) (understood in terms of the space of potentially algebraic variables).
		 **/
		inline int dimensionLower(Symbol v) const {
			if (!isAMemberOf(v,vars)) {
				return 0;
			}
			int vi(variableIndex(v)),varsSize(vi),algVarsSize;
			varsSize = vars.size() - 1 - varsSize;
			for (int i=0; i<algVars.size(); ++i) {
				if (variableIndex(algVars[i])>vi)
					algVarsSize++;
			}
			return (varsSize - algVarsSize);
		}
		
		/**
		 * The codimension of the triangular set (understood in terms of the space of potentially algebraic variables).
		 **/
		inline int codimension() const {
			return algVars.size();
		}
		
		/**
		 * Test to determine whether the triangular set can be treated as zero dimensional.
		 **/
		bool canComputeInDimensionZero(const RecursivePoly& p, bool excludeMainVariable = false) const;
		
		/**
		 * Test to determine if only algebraic variables (aside from transcendentals) appear in the polynomials of the set.
		 **/
		bool isZeroDimensionalMathematically() const;
		
		/**
		 * Select a polynomial given the leading variable;
		 * if no such polynomial, 0 is returned
		 * @param x: The leading variable name
		 **/
		RecursivePoly select(const Symbol& s) const;
		
		/**
		 * Replace each polynomial of the regular chain with its primitive part
		 *
		 **/
		void makePrimitive();
		
		/**
		 * Returns the ts consisting of polynomials with 
		 * main variable strictly less than s
		 *
		 * @param s: Symbol of the main variable of specified element of the triangular set
		 * @param ts: The returned triangular set
		 **/
		void lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;
		
		/**
		 * Returns the ts consisting of polynomials with
		 * main variable strictly greater than s
		 *
		 * @param s: Symbol of the main variable of specified element of the triangular set
		 * @param ts: The returned triangular set
		 **/
		void upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;
		
		/**
		 * Pseudo division
		 * Return the pseudo-remainder, the pseudo-quotients and
		 * c such that c*p = ∑(q_i T_i) + r 
		 * @param p: An input polynomial
		 * @param quo: The array of quotients
		 * @param c: The constant multiplied to the input polynomial
		 **/
		RecursivePoly pseudoDivide (const RecursivePoly& p, std::vector<RecursivePoly>* quo=NULL, RecursivePoly* c=NULL) const;
		
		/**
		 * Monic division
		 * Returns the remainder; requires the triangular set be monic
		 *
		 * @param p: An input polynomial
		 * @param quo: The quotients
		 **/
//		RecursivePoly monicDivide (const RecursivePoly& p, std::vector<RecursivePoly>* quo=NULL) const;
		
		/**
		 * normalForm in the sense of Groebner basis
		 *
		 * @param p: Input polynomial
		 * 
		 **/
		RecursivePoly normalForm(const RecursivePoly& p, std::vector<RecursivePoly>* Q=NULL) const;
		
		/**
		 * reduce the polynomial p modulo the triangular set
		 *
		 * @param p: Input polynomial
		 * 
		 **/
		RecursivePoly reduce(const RecursivePoly& p) const;
		
		/**
		 * returns (c, r) such that c*r = p modulo sat(T) such that
		 * c has no algebraic variables
		 *
		 * @param p: Input polynomial
		 * 
		 **/
		RecursivePoly reduce(const RecursivePoly& p, RecursivePoly& c, bool usePrimitiveFactorization = true, bool onlyInDimZero = false) const;
		
		/**
		 * Make monic
		 *
		 * @param
		 **/
//		inline void makeMonic();

		/**
		 * Generate a random triangular set polynomial
		 *
		 * @param variables: variables of the triangular set
		 * @param algVar: index of the algebraic variable
		 * @param transcendentalVariables: transcendental variables of the triangular set
		 * @param nTerms: number of terms in the polynomial
		 * @param coefBound: maximum size of the coefficients
		 * @param pSparsity: sparsity of the polynomial
		 * @param includeNeg: whether to include negative coefficients
		 **/
		RecursivePoly randomTriangularSetPolynomial(std::vector<Symbol> variables, int algVar, std::vector<Symbol> transcendentalVariables, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg);

		/**
		 * Generate a random triangular set polynomial
		 *
		 * @param variables: variables of the triangular set
		 * @param algVar: index of the algebraic variable
		 * @param transcendentalVariables: transcendental variables of the triangular set
		 * @param maxDegs: maximum degrees among the set of variables
		 * @param coefBound: maximum size of the coefficients
		 * @param pSparsity: proportional sparsity of the polynomial
		 * @param includeNeg: whether to include negative coefficients
		 **/
		RecursivePoly randomTriangularSetPolynomial(std::vector<Symbol> variables, int algVar, std::vector<Symbol> transcendentalVariables, std::vector<int> maxDegs, unsigned long int coefBound, double pSparsity, bool includeNeg);

		/**
		 * Generate a random triangular set
		 *
		 * @param nVars: number of variables
		 * @param nAlgVars: number of algebraic variables
		 * @param nTrcVars: number of transcendental variables
		 * @param nTerms: number of terms in the polynomial
		 * @param coefBound: maximum size of the coefficients
		 * @param pSparsity: sparsity of the polynomial
		 * @param includeNeg: whether to include negative coefficients
		 **/
		void randomTriangularSet(int nVars, int nAlgVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg);

		/**
		 * Generate a random triangular set
		 *
		 * @param nVars: number of variables
		 * @param nAlgVars: number of algebraic variables
		 * @param nTrcVars: number of transcendental variables
		 * @param nTerms: number of terms in the polynomial
		 * @param coefBound: maximum size of the coefficients
		 * @param pSparsity: sparsity of the polynomial
		 * @param includeNeg: whether to include negative coefficients
		 **/
		void randomStronglyNormalizedTriangularSet(int nVars, int nAlgVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg);
		
		/**
		 * Display the triangular set
		 *
		 **/
		void display();

		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param a: A triangular set
		 **/
		inline friend std::ostream& operator<< (std::ostream& out, TriangularSet<Field,RecursivePoly>& a) {
			bool isNotFirst = 0;
			out << "[";
			for (int i = 0; i < a.set.size(); ++i) {
				if (!a.set[i].isZero()) {
					if (isNotFirst) { out << ", "; }
					out << a.set[i];
					isNotFirst = 1;
				}
			}
			out << "]";
			return out;
		}
		
		inline ExpressionTree convertToExpressionTree() const {
			if (!set.size()) {
				ExprTreeNode* node = new ExprTreeNode(EXPR_ARRAY, NULL, NULL, NULL);
				return ExpressionTree(node);
			}
			else {
				std::vector<RecursivePoly> tsp;
				for (int i=0; i<set.size(); ++i) {
					if (!set[i].isZero())
						tsp.push_back(set[i]);
				}
				ExpressionTree t;
				t.fromVector<RecursivePoly>(tsp);
				return(t);
			}
		}
};

#endif
