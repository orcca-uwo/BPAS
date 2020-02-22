#ifndef _ZERODIMENSIONALREGULARCHAIN_H_
#define _ZERODIMENSIONALREGULARCHAIN_H_


#include "../polynomial.h"
#include "../Ring/BPASField.hpp"
#include "regularchain.hpp"
#include "chainstructures.hpp"
#include "../TriangularSet/triangularset.hpp"
#include "../SubResultantChain/subresultantchain.hpp"

/** 
 * A zero-dimensional RegularChain.
 */
template <class Field, class RecursivePoly>
class ZeroDimensionalRegularChain : public RegularChain<Field,RecursivePoly>,
									public virtual BPASZeroDimensionalRegularChain<Field,RecursivePoly>
{
	protected:
	
		using RegularChain<Field,RecursivePoly>::mode;
		using RegularChain<Field,RecursivePoly>::set;
		using RegularChain<Field,RecursivePoly>::vars;
		using RegularChain<Field,RecursivePoly>::algVars;
		using RegularChain<Field,RecursivePoly>::trcVars;
		using RegularChain<Field,RecursivePoly>::stronglyNormalized;
		using RegularChain<Field,RecursivePoly>::characteristic;
		using RegularChain<Field,RecursivePoly>::regularChainOptions;
		
		using RegularChain<Field,RecursivePoly>::updateTriangularSetStates;
		using RegularChain<Field,RecursivePoly>::variableIndex;
		
		void cutRegularChain(const ZeroDimensionalRegularChain<Field,RecursivePoly>& T, const Symbol& v, ZeroDimensionalRegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const;
		void cutRegularChain(const Symbol& v, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv);
		void cutRegularChain(const Symbol& v, ZeroDimensionalRegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv);

		/**
		 * Merge two compatible zero-dimensional regular chains, where the polynomials in the current object are above all the elements in rc.
		 *
		 * @param rc: input regular chain
		 **/
		void mergeZeroDimensionalRegularChains(const ZeroDimensionalRegularChain<Field,RecursivePoly>& rc);
		
		/**
		 * Add an input polynomial to the current object
		 *
		 * @param p: A polynomial
		 **/
		void constructChain(const RecursivePoly& p, int options=ASSUME_REGULAR);
		
		/**
		 * Add an an upper regular chain to the current object
		 *
		 * @param T: A polynomial
		 **/
		void constructChain(const RegularChain<Field,RecursivePoly>& T, int options=ASSUME_REGULAR);
		
		/**
		 * Construct a set of regular chains from the current object and an input polynomial
		 *
		 * @param p: A polynomial
		 **/
		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> constructChains(const RecursivePoly& p, int options=ASSUME_REGULAR) const;
		
		/**
		 * Construct a set of regular chains from the current object and an input polynomial
		 *
		 * @param p: A polynomial
		 **/
		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> constructChains(const RegularChain<Field,RecursivePoly>& T, int options=ASSUME_REGULAR) const;
		
	public:
	
		using RegularChain<Field,RecursivePoly>::allVariables;
		
		/**
		 * Default constructor: creates an empty triangular set of variable size
		 * with empty list of transcendentals
		 *
		 * @param
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>();
		
		/**
		 * Construct an empty triangular set of variable size with
		 * variables given by xs with empty list of transcendentals
		 *
		 * @param ps: The transcendental variable names
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const std::vector<Symbol>& ps);
		
//		/**
//		 * Construct an empty triangular set of fixed size with
//		 * variables given by xs and list of transcendentals given by ts
//		 *
//		 * @param xs: The variable names
//		 * @param ts: The transcendental variable names
//		 **/
//		ZeroDimensionalRegularChain<Field,RecursivePoly> (const std::vector<Symbol>& xs, const std::vector<Symbol>& ts);
		
		/**
		 * Construct a variable triangular set containing p, such that the variables of p are treated as algebraic,
		 * with empty list of transcendentals
		 *
		 * @param p: The polynomial to add
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const RecursivePoly& p);
		
		/**
		 * Construct a variable triangular set containing p, such that the variables in ts are
		 * treated as transcendental, while any remaining variables of p are treated as algebraic
		 *
		 * @param p: The polynomial to add
		 * @param ts: The transcendental variable names
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const RecursivePoly& p, const std::vector<Symbol>& ts);
		
		/**
		 * Copy constructor
		 *
		 * @param a: A triangular set
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Copy constructor
		 *
		 * @param a: A triangular set
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const RegularChain<Field,RecursivePoly>& a, int options=0);
		
		/**
		 * Move constructor
		 *
		 * @param a: An r-value reference triangular set
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);
		
		/**
		 * Move constructor
		 *
		 * @param a: An r-value reference triangular set
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (RegularChain<Field,RecursivePoly>&& a, int options=0);

		/**
		 * Computational constructor: creates a triangular set given all the data
		 *
		 * @param vs: variables of the triangular set
		 * @param avs: algebraic variables of the triangular set
		 * @param tvs: transcendental variables of the triangular set
		 * @param polys: polynomials of the triangular set
		 * @param tsm: whether the triangular set is variable or fixed
		 * @param c: characteristic of the triangular set
		 * @param normal: whether the triangular set is strongly normalized
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const std::vector<Symbol>& vs, const std::vector<Symbol>& avs, const std::vector<Symbol>& tvs, const std::vector<RecursivePoly>& ts, TriangularSetMode tsm, const mpz_class& c);
		

		/**
		 * Copy an object derived from abstract BPASTriangularSet class to type of current object
		 *
		 * @param ts: triangular set to copy
		 **/
//		inline void copy(const BPASRegularChain<Field,RecursivePoly>& ts) override {
//			if (dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&ts))
//				*this = *dynamic_cast<const ZeroDimensionalRegularChain<Field,RecursivePoly>*>(&ts);
//			else throw (std::invalid_argument("BPAS: Cannot cast BPASTriangularSet to RegularChain."));
//		}
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A triangular set
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A BPASTriangularSet
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (const BPASTriangularSet<Field,RecursivePoly>& a) override;
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A BPASRegularChain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (const BPASRegularChain<Field,RecursivePoly>& a) override;
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A BPASZeroDimensionalRegularChain
		 **/
		BPASZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (const BPASZeroDimensionalRegularChain<Field,RecursivePoly>& a) override;
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A triangular set
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A BPASTriangularSet
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (BPASTriangularSet<Field,RecursivePoly>&& a) override;
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A BPASRegularChain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (BPASRegularChain<Field,RecursivePoly>&& a) override;
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A BPASZeroDimensionalRegularChain
		 **/
		BPASZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (BPASZeroDimensionalRegularChain<Field,RecursivePoly>&& a) override;
		
		/**
		 * Add operator +
		 * Adds a polynomial to a triangular set and returns a new triangular set
		 *
		 * @param p: A sparse multivariate polynomial
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> operator+ (const RecursivePoly& p);
		
		/**
		 * Add assignment operator +=
		 * Adds a polynomial to a triangular set
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator+= (const RecursivePoly& p);
		
		/**
		 * Add operator +
		 * Adds a polynomial to a regular chain and returns a new regular chain
		 *
		 * @param p: A sparse multivariate polynomial
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> operator+ (const RegularChain<Field,RecursivePoly>& T) const;

		/**
		 * Add assignment operator +=
		 * Adds a polynomial to a regular chain
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator+= (const RegularChain<Field,RecursivePoly>& T);
		
		/**
		 * Identity operator ==
		 *
		 *
		 * @param a: A triangular set
		 **/
		bool operator== (ZeroDimensionalRegularChain<Field,RecursivePoly>& a);

		/**
		 * Negated identity operator !=
		 *
		 *
		 * @param a: A triangular set
		 **/
		bool operator!= (ZeroDimensionalRegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Get the number of variables
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return RegularChain<Field,RecursivePoly>::numberOfVariables();
		}
//		using RegularChain<Field,RecursivePoly>::numberOfVariables;
		
		/**
		 * Get the size of the triangular set
		 *
		 * @param
		 **/
		inline int size() const {
			return RegularChain<Field,RecursivePoly>::size();
		}
		
		/**
		 * Get the number of algebraic variables
		 *
		 * @param
		 **/
		inline int numberOfAlgebraicVariables() const {
			return RegularChain<Field,RecursivePoly>::numberOfAlgebraicVariables();
		}
		
		/**
		 * Get the number of transcendental variables
		 *
		 * @param
		 **/
		inline int numberOfTranscendentalVariables() const {
			return RegularChain<Field,RecursivePoly>::numberOfTranscendentalVariables();
		}
		
		/**
		 * Get the variable names in decreasing order
		 *
		 * @param
		 **/
		inline std::vector<Symbol> variables() const {
			return RegularChain<Field,RecursivePoly>::variables();
		}
		
		/**
		 * Get algebraic variables
		 *
		 * @param
		 **/
		inline std::vector<Symbol> mainVariables() const {
			return RegularChain<Field,RecursivePoly>::mainVariables();
		}
		
		/**
		 * Get transcendentalVariables variables
		 *
		 * @param
		 **/
		inline std::vector<Symbol> transcendentalVariables() const {
			return RegularChain<Field,RecursivePoly>::transcendentalVariables();
		}
		
		inline bool isAlgebraic(const Symbol& s) const {
			return RegularChain<Field,RecursivePoly>::isAlgebraic(s);
		}
		
		inline bool isEmpty() const {
			return RegularChain<Field,RecursivePoly>::isEmpty();
		}
		
		inline std::vector<RecursivePoly> polynomials() const {
			return RegularChain<Field,RecursivePoly>::polynomials();
		}
		
		inline RecursivePoly select(const Symbol& s) const {
			return RegularChain<Field,RecursivePoly>::select(s);
		}
		
		void lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;
		
		void upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;
		
		/**
		 * Compute the intersection of the varieties of the input polynomial and the 
		 * current regular chain, i.e., compute a splitting of (T_i) such that all 
		 * zeros of T_i are common to p and all common zeros of p and (T_i) are included
		 * among the T_i.
		 *
		 * @param p: input polynomial
		 **/
		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> intersect(const RecursivePoly& p) const;
		
		/**
		 * Regularize the input polynomial with respect to the current regular chain,
		 * i.e., compute a splitting of pairs (f_i,T_i) such that h_i has invertible
		 * initial with respect to T_i.
		 *
		 * @param p: input polynomial
		 **/
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> regularize(const RecursivePoly& p, bool regularizeInitial=true) const;
		
		/**
		 * Determine whether a polynomial is invertible with respect to the current regular chain,
		 * i.e. compute a splitting of pairs (b_i,T_i) such that b_i is true if f is invertible
		 * modulo T_i and b_i is false if f is zero modulo T_i.
		 *
		 * @param p: input polynomial
		 **/
		std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> isInvertible(const RecursivePoly& p) const;
		
		/**
		 * Compute the last nonzero subresultant of the subresultant chain src of f and g modulo the current regular chain,
		 * i.e., compute a splitting of pairs (R_i,T_i) such that R_i is the last nonzero subresultant modulo T_i.
		 *
		 * @param f: input polynomial
		 * @param g: input polynomial
		 * @param v: common main variable of f and g
		 * @param isGCD: boolean flag specifying whether the function is called for GCD computation
		 **/
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> lastNonZeroSubResultant(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, bool isGCD) const;
		
		/**
		 * Compute the last nonzero subresultant of the subresultant chain src of f and g modulo the current regular chain,
		 * i.e., compute a splitting of pairs (R_i,T_i) such that R_i is the last nonzero subresultant modulo T_i.
		 *
		 * @param f: input polynomial
		 * @param g: input polynomial
		 * @param v: common main variable of f and g
		 * @param src: subresultant chain of f and g
		 * @param isGCD: boolean flag specifying whether the function is called for GCD computation
		 **/
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> lastNonZeroSubResultant_inner(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool isGCD) const;
//		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> lastNonZeroSubResultant_inner(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, const SubResultantChain<Field,SparseUnivariatePolynomial<RecursivePoly>>& src, bool isGCD) const;
		
		/**
		 * Compute the gcd of polynomials f and g with respect to v modulo the current object
		 * regular chain
		 *
		 * @param f: input polynomial
		 * @param g: input polynomial
		 * @param v: common main variable of f and g
		 **/
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> regularGCD(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v);
		
		/**
		 * Get the number of variables
		 *
		 * @param
		 **/
//		using RegularChain<Field,RecursivePoly>::numberOfVariables;

//		inline void updateRegularChainsStates(RecursivePoly smp, int k) {}
		
//		int regularizeDim0(RecursivePoly poly, std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> &newRcOut, std::vector<RecursivePoly> &newPolys );
		
//		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> intersect(RecursivePoly poly);
		
//		int regularize(RecursivePoly p, std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> &newRcOut, std::vector<RecursivePoly> &newPolys );
		
		/**
		 * Generate a random regular chain
		 *
		 * @param nVars: number of variables = number of algebraic variables
		 * @param nTrcVars: number of transcendental variables
		 * @param nTerms: maximum number of terms in the polynomials
		 * @param coefBound: maximum coefficient size
		 * @param pSparsity: sparsity of the polynomials
		 * @param includeNeg: whether to include negative coefficients
		 **/
		void randomZeroDimensionalRegularChain(int nVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg);

		/**
		 * Generate a random zero dimensional regular chain
		 *
		 * @param nVars: number of variables = number of algebraic variables
		 * @param nTrcVars: number of transcendental variables
		 * @param maxDegs: maximum degrees among the full list of variables
		 * @param coefBound: maximum coefficient size
		 * @param pSparsity: sparsity of the polynomials
		 * @param includeNeg: whether to include negative coefficients
		 **/
		void randomZeroDimensionalRegularChain(int nVars, int nTrcVars, std::vector<int> maxDegs, unsigned long int coefBound, double pSparsity, bool includeNeg);
		
};

#endif
