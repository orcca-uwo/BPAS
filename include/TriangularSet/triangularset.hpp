#ifndef _TRIANGULARSET_H_
#define _TRIANGULARSET_H_

//#include "Polynomial/mpolynomial.h"
#include "../ring.h"
#include "../Utils/TemplateHelpers.hpp"
#include "../Utils/SymbolHelpers.hpp"
#include "BPASTriangularSet.hpp"
//#include "../RegularChain/regularchain.hpp"
//#include "../RegularChain/zerodimensionalregularchain.hpp"
#include <iostream>

using std::endl;
using std::cerr;

extern long long unsigned int rcProfilingStart;
extern float primitivePartTime;
extern float squareFreePartTime;
extern float subresultantChainTime;
extern float zerodimensionalregularchainTime;
extern float pseudoDivideTime;
extern float normalFormTime;

extern float tsCopyTime;

//template <class Field, class RecursivePoly> RegularChain;

/**
 * An enumeration which describes the mode of the triangular set.
 * TS_FIXED is the ALGEBRAIC case, forces a fixed list of variables, determinining the maximum size of the TS.
 * TS_VARIABLE is the DIFFERENTIAL case, allows the list of variables, hence size of the TS, to change.
 **/
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

		/**
		 * Determine whether the triangular set is strongly normalized and the maximum variable
		 * of the largest strongly normalized subset of the triangular set.
		 *
		 * @param
		 **/
		void updateTriangularSetStates();

		/**
		 * Determine whether after adding p the triangular set is strongly normalized and update the maximum variable
		 * of the largest strongly normalized subset of the triangular set.
		 *
		 * @param p: a recursively viewed polynomial that has just been added to the set
		 **/
		void updateTriangularSetStates(const RecursivePoly& p);

		/**
		 * Return the index of the input symbol in the array of (potentially algebraic) variables of the triangular set.
		 *
		 * @param s: a symbol
		 **/
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
		 * variables given by xs with empty list of transcendentals.
		 *
		 * @param xs: The variable names
		 **/
		TriangularSet<Field,RecursivePoly> (const std::vector<Symbol>& xs);

		/**
		 * Construct an empty triangular set of fixed size in the s decreasingly ordered
		 * variables given by xs and list of transcendentals given by ts.
		 *
		 * @param xs: The variable names
		 * @param ts: The transcendental variable names
		 **/
		TriangularSet<Field,RecursivePoly> (const std::vector<Symbol>& xs, const std::vector<Symbol>& ts);

		/**
		 * Construct a variable triangular set containing p, such that the variables of p are treated as (potentially algebraic) variables,
		 * with empty list of transcendentals.
		 *
		 * @param p: The recursively viewed polynomial to add
		 **/
		TriangularSet<Field,RecursivePoly> (const RecursivePoly& p);

		/**
		 * Construct a variable triangular set containing p, such that the variables in ts are
		 * treated as transcendental, while any remaining variables of p are treated as (potentially algebraic) variables.
		 *
		 * @param p: The recursively viewed polynomial to add
		 * @param ts: The transcendental variable names
		 **/
		TriangularSet<Field,RecursivePoly> (const RecursivePoly& p, const std::vector<Symbol>& ts);

		/**
		 * Copy constructor.
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
		 * Move constructor.
		 *
		 * @param a: An r-value triangular set
		 **/
		TriangularSet<Field,RecursivePoly> (TriangularSet<Field,RecursivePoly>&& a);

//		/**
//		 * Move constructor
//		 *
//		 * @param a: An r-value regular chain
//		 **/
//		TriangularSet<Field,RecursivePoly> (RegularChain<Field,RecursivePoly>&& a);

//		/**
//		 * Move constructor
//		 *
//		 * @param a: An r-value zero dimensional regular chain
//		 **/
//		TriangularSet<Field,RecursivePoly> (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);

		/**
		 * Computational constructor: creates a triangular set given all the data.
		 *
		 * @param vs: rvalue reference to variables of the triangular set
		 * @param avs: rvalue reference to algebraic variables of the triangular set
		 * @param tvs: rvalue reference to transcendental variables of the triangular set
		 * @param polys: rvalue reference to polynomials of the triangular set
		 * @param tsm: whether the triangular set is variable or fixed
		 * @param c: characteristic of the triangular set
		 **/
		TriangularSet<Field,RecursivePoly> (const std::vector<Symbol>&& vs, const std::vector<Symbol>&& avs, const std::vector<Symbol>&& tvs, const std::vector<RecursivePoly>&& ts, TriangularSetMode tsm, const mpz_class& c);

		/**
		 * Deconstructor.
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
		 * Tests if the TriangularSet is empty.
		 *
		 * @param
		 **/
		bool isEmpty() const;

		/**
		 * Tests if the polynomial is constant relative to the TriangularSet, i.e., whether it is and element of
		 * the Field or its only variables are transcendental.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		bool isConstantPolynomial(const RecursivePoly& p) const;

		/**
		 * Assignment operator =.
		 *
		 * @param a: A triangular set
		 **/
		TriangularSet<Field,RecursivePoly>& operator= (const TriangularSet<Field,RecursivePoly>& a);

		/**
		 * Assignment operator =.
		 *
		 * @param a: A triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (const BPASTriangularSet<Field,RecursivePoly>& a) override;

		/**
		 * Move assignment operator =.
		 *
		 * @param a: A triangular set
		 **/
		TriangularSet<Field,RecursivePoly>& operator= (TriangularSet<Field,RecursivePoly>&& a);

		/**
		 * Move assignment operator =.
		 *
		 * @param a: A triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (BPASTriangularSet<Field,RecursivePoly>&& a) override;

		/**
		 * Add operator +.
		 * Adds a polynomial to a triangular set and returns a new triangular set.
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		TriangularSet<Field,RecursivePoly> operator+ (const RecursivePoly& p);

		/**
		 * Add assignment operator +=.
		 * Adds a polynomial to a triangular set.
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		TriangularSet<Field,RecursivePoly>& operator+= (const RecursivePoly& p);

		/**
		 * Identity operator ==.
		 *
		 * @param a: A triangular set
		 **/
		bool operator== (const TriangularSet<Field,RecursivePoly>& a) const;

		/**
		 * Negated identity operator !=.
		 *
		 * @param a: A triangular set
		 **/
		bool operator!= (const TriangularSet<Field,RecursivePoly>& a) const;

		/**
		 * Get the number of variables.
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return vars.size();
		}

		/**
		 * Get the size of the triangular set.
		 *
		 * @param
		 **/
		inline int size() const {
			return algVars.size();
		}

		/**
		 * Get the number of algebraic variables.
		 *
		 * @param
		 **/
		inline int numberOfAlgebraicVariables() const {
			return algVars.size();
		}

		/**
		 * Get the number of transcendental variables.
		 *
		 * @param
		 **/
		inline int numberOfTranscendentalVariables() const {
			return trcVars.size();
		}

		/**
		 * Get the variable names in decreasing order.
		 *
		 * @param
		 **/
		inline std::vector<Symbol> variables() const {
			return vars;
		}

		/**
		 * Get the algebraic variables.
		 *
		 * @param
		 **/
		inline std::vector<Symbol> mainVariables() const {
			return algVars;
		}

		/**
		 * Get the transcendentalVariables variables.
		 *
		 * @param
		 **/
		inline std::vector<Symbol> transcendentalVariables() const {
			return trcVars;
		}

		/**
		 * Get the list of variables followed by the transcendental variables.
		 *
		 * @param
		 **/
		std::vector<Symbol> allVariables() const;


		/**
		 * Determine if the input variable s is algebraic, i.e., if the triangular set
		 * contains a polynomial with s its as leading variable.
		 *
		 * @param s: the input variable
		 **/
		inline bool isAlgebraic(const Symbol& s) const {
			if (std::find(algVars.begin(),algVars.end(),s) != algVars.end())
				return true;
			else
				return false;
		}

		/**
		 * Return true if the triangular set is strongly normalized, i.e., the initals of all polynomials
		 * are in the Field; return false otherwise.
		 *
		 * @param
		 **/
		inline bool isStronglyNormalized() const {
			return stronglyNormalized;
		}

		/**
		 * Get the vector of polynoials in the triangular set.
		 *
		 * @param
		 **/
		inline std::vector<RecursivePoly> polynomials() const {
			return set;
		}

		/**
		 * Return the dimension of the triangular set (understood in terms of the space of (potentially algebraic) variables).
		 *
		 * @param
		 **/
		inline int dimension() const {
			return (vars.size() - algVars.size());
		}

		/**
		 * Return the dimension of the triangular set lower(v) (understood in terms of the space of (potentially algebraic) variables).
		 *
		 * @param
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
		 * Return the codimension of the triangular set (understood in terms of the space of (potentially algebraic) variables).
		 *
		 * @param
		 **/
		inline int codimension() const {
			return algVars.size();
		}

		/**
		 * Test to determine whether the triangular set can be treated as zero dimensional, i.e., whether the triangular set
		 * becomes zero dimensional when all non-algebraic variables are removed and whether the polynomial p contains only
		 * algebraic variables.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param p: (optional) flag to exclude the main variable of p when determining whether it contains only algebraic variables (default false)
		 **/
		bool canComputeInDimensionZero(const RecursivePoly& p, bool excludeMainVariable = false) const;

		/**
		 * Test to determine if only algebraic variables (aside from transcendentals) appear in the polynomials of the triangular set.
		 *
		 * @param
		 **/
		bool isZeroDimensionalMathematically() const;

		/**
		 * Select a polynomial given the leading variable;
		 * if no such polynomial, 0 is returned.
		 *
		 * @param x: The leading variable name
		 **/
		RecursivePoly select(const Symbol& s) const;

		/**
		 * Replace each polynomial of the triangular set with its primitive part.
		 *
		 * @param
		 **/
		void makePrimitive();

		/**
		 * Returns the triangular set consisting of polynomials with
		 * main variable strictly less than s.
		 *
		 * @param s: symbol of the main variable of specified element of the triangular set
		 * @param ts: The returned triangular set
		 **/
		void lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;

		/**
		 * Returns the triangular set consisting of polynomials with
		 * main variable strictly greater than s.
		 *
		 * @param s: symbol of the main variable of specified element of the triangular set
		 * @param ts: The returned triangular set
		 **/
		void upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;

		/**
		 * Cut an input triangular set at the symbol v, returning the subchain below v, the polynomial with main variable v and the subchain above v.
		 *
		 * @param T: a triangular set
		 * @param v: a symbol
		 * @param Tlv: the subchain of T below v
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 * @param Tgv: the subchain of T above v
		 **/
		void cutChain(const TriangularSet<Field,RecursivePoly>& T, const Symbol& v, TriangularSet<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, TriangularSet<Field,RecursivePoly>& Tgv) const;

		/**
		 * Cut the current object at the symbol v, returning the polynomial with main variable v and the subchain above v.
		 *
		 * @param v: a symbol
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 * @param Tgv: the subchain of the current object above v
		 **/
		void cutChain(const Symbol& v, RecursivePoly& Tv, TriangularSet<Field,RecursivePoly>& Tgv) const;

		/**
		 * Cut the current object at the symbol v, returning the subchain below v and the polynomial with main variable v.
		 *
		 * @param v: a symbol
		 * @param Tlv: the subchain of the current object below v
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 **/
		void cutChain(const Symbol& v, TriangularSet<Field,RecursivePoly>& Tlv, RecursivePoly& Tv) const;

		/**
		 * Pseudo division: return the pseudo-remainder, the pseudo-quotients and
		 * c such that c*p = ∑(q_i T_i) + r.
		 *
		 * @param p: an input recursively viewed polynomial
		 * @param quo: (optional) the array of quotients
		 * @param c: (optional) the constant multiplied to the input polynomial
		 **/
		RecursivePoly pseudoDivide (const RecursivePoly& p, std::vector<RecursivePoly>* quo=NULL, RecursivePoly* c=NULL) const;

		/**
		 * Return the normalForm of the input polynomial modulo the triangular set in the sense of Groebner basis
		 *
		 * @param p: innput recursively viewed polynomial
		 * @param Q: (optional) the array of quotient
		 **/
		RecursivePoly normalForm(const RecursivePoly& p, std::vector<RecursivePoly>* Q=NULL) const;

		/**
		 * Reduce the input polynomial modulo the triangular set.
		 *
		 * @param p: input recursively viewed polynomial
		 **/
		RecursivePoly reduce(const RecursivePoly& p) const;

		/**
		 * returns r such that c*r = p modulo sat(T) such that
		 * c has no algebraic variables, and c is returned as an input parameter.
		 *
		 * @param p: input recursively viewed polynomial
		 * @param c: returned value of the content of p modulo sat(T)
		 * @param usePrimitiveFactorization: (optional) whether to use primitive factorization to compute c (default true)
		 * @param onlyInDimZero: (optional) only perform the reduction if the canComputeInDimensionZero(p) is true (default false)
		 **/
		RecursivePoly reduce(const RecursivePoly& p, RecursivePoly& c, bool takeMainPrimitivePart = false, bool onlyInDimZero = false) const;

		/**
		 * Generate a random triangular set polynomial based on its number of terms.
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
		 * Generate a random triangular set polynomial based on its maximum degrees in its variables.
		 *
		 * @param variables: variables of the triangular set
		 * @param algVar: index of the algebraic variable
		 * @param transcendentalVariables: transcendental variables of the triangular set
		 * @param maxDegs: vector of maximum degrees among the set of variables
		 * @param coefBound: maximum size of the coefficients
		 * @param pSparsity: proportional sparsity of the polynomial
		 * @param includeNeg: whether to include negative coefficients
		 **/
		RecursivePoly randomTriangularSetPolynomial(std::vector<Symbol> variables, int algVar, std::vector<Symbol> transcendentalVariables, std::vector<int> maxDegs, unsigned long int coefBound, double pSparsity, bool includeNeg);

		/**
		 * Generate a random triangular set based on the number of terms of its polynomials.
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
		 * Generate a random strongly normalized triangular set based on the number of terms in its polynomials.
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
		 * Display the triangular set.
		 *
		 * @param
		 **/
		void display();

		/**
		 * Overload stream operator << for triangular sets.
		 *
		 * @param out: Stream object
		 * @param a: A triangular set
		 **/
		inline friend std::ostream& operator<< (std::ostream& out, const TriangularSet<Field,RecursivePoly>& a) {
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

	   /**
	    * Convert a triangular set to an expression tree (array of its polynomials).
	    *
		* @param
		**/
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
