#ifndef _REGULARCHAIN_H_
#define _REGULARCHAIN_H_

#include "../TriangularSet/triangularset.hpp"
#include "../SubResultantChain/subresultantchain.hpp"
#include "../polynomial.h"
#include "chainstructures.hpp"

template <class Field, class RecursivePoly>
class ZeroDimensionalRegularChain;

enum RegularChainOption {
	ASSUME_REGULAR         = 0x001,
	ASSUME_REDUCED         = 0x002,
	ASSUME_PRIMITIVE       = 0x004,
	ASSUME_SQUAREFREE      = 0x008,
	ASSUME_ZERODIMENSIONAL = 0x010,
	ASSUME_MAKESCHAIN      = 0x020,
	MAINTAIN_SQUAREFREE    = 0x040,
	MAINTAIN_NORMALIZED    = 0x080,
	MAINTAIN_PRIME         = 0x100,
	CHECK_REGULAR          = 0xffe
};

/**
 * A RegularChain over a BPASRecursivelyViewedPolynomial with coefficients in a field.
 */
template <class Field, class RecursivePoly>
class RegularChain : public TriangularSet<Field,RecursivePoly>,
					 public virtual BPASRegularChain<Field,RecursivePoly>
{
	protected:
		using TriangularSet<Field,RecursivePoly>::mode;
		using TriangularSet<Field,RecursivePoly>::set;
		using TriangularSet<Field,RecursivePoly>::vars;
		using TriangularSet<Field,RecursivePoly>::algVars;
		using TriangularSet<Field,RecursivePoly>::trcVars;
		using TriangularSet<Field,RecursivePoly>::stronglyNormalized;
		using TriangularSet<Field,RecursivePoly>::characteristic;
		using TriangularSet<Field,RecursivePoly>::updateTriangularSetStates;
		using TriangularSet<Field,RecursivePoly>::variableIndex;
		using TriangularSet<Field,RecursivePoly>::canComputeInDimensionZero;
		using TriangularSet<Field,RecursivePoly>::isZeroDimensionalMathematically;
		
		int regularChainOptions = MAINTAIN_SQUAREFREE;
		bool isSquareFree;

		void updateRegularChainStates();
		void updateRegularChainStates(const RecursivePoly& p);
		
		std::vector<RegularChain<Field,RecursivePoly>> triangularize_inner(std::vector<RecursivePoly>& F) const;
		bool cleanSet(std::vector<RecursivePoly>& polys) const;
		
		void cutRegularChain(const RegularChain<Field,RecursivePoly>& T, const Symbol& v, RegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const;
		void cutRegularChain(const Symbol& v, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv);
		void cutRegularChain(const Symbol& v, RegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv);
		
		void regularizeList(const std::vector<RecursivePoly>& knownRegular, const std::vector<RecursivePoly>& unknownIfRegular, std::vector<RegularChain<Field,RecursivePoly>>& singularComponents, std::vector<RegularChain<Field,RecursivePoly>>& regularComponents) const;
//		std::vector<BoolChainPair<RegularChain<Field,RecursivePoly>>> areFactorsRegular(const std::vector<RecursivePoly>& factors) const;
		
		/**
		 * Add an input polynomial to the current object
		 *
		 * @param p: A polynomial
		 **/
		void constructChain(const RecursivePoly& p, int options=ASSUME_REGULAR);
		
		/**
		 * Construct a regular chain from the current object and an input upper regular chain
		 *
		 * @param T: An upper regular chain
		 **/
		void constructChain(const RegularChain<Field,RecursivePoly>& T, int options);
		
		/**
		 * Construct a set of regular chains from the current object and an input polynomial
		 *
		 * @param p: A polynomial
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> constructChains(const RecursivePoly& p, int options=ASSUME_REGULAR) const;
		
		/**
		 * Construct a set of regular chains from the current object and an input polynomial
		 *
		 * @param p: A polynomial
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> constructChains(const RegularChain<Field,RecursivePoly>& T, int options=ASSUME_REGULAR) const;


	public:
	
		using TriangularSet<Field,RecursivePoly>::allVariables;
	
		/**
		 * Default constructor: creates an empty regular chain of variable size
		 * with empty list of transcendentals
		 *
		 * @param
		 **/
		RegularChain<Field,RecursivePoly>();
		
		/**
		 * Construct an empty regular chain of fixed size in the s decreasingly ordered 
		 * variables given by xs with empty list of transcendentals
		 *
		 * @param s: The number of polynomials
		 * @param xs: The variable names
		 **/
		RegularChain<Field,RecursivePoly> (const std::vector<Symbol>& xs);
		
		/**
		 * Construct an empty regular chain of fixed size in the s decreasingly ordered 
		 * variables given by xs and list of transcendentals given by ts
		 *
		 * @param s: The number of polynomials
		 * @param xs: The variable names
		 * @param ts: The transcendental variable names
		 **/
		RegularChain<Field,RecursivePoly> (const std::vector<Symbol>& xs, const std::vector<Symbol>& ts);
		
		/**
		 * Construct a variable regular chain containing p, such that the variables of p are treated as algebraic,
		 * with empty list of transcendentals
		 *
		 * @param p: The polynomial to add
		 **/
		RegularChain<Field,RecursivePoly> (const RecursivePoly& p);
		
		/**
		 * Construct a variable regular chain containing p, such that the variables in ts are
		 * treated as transcendental, while any remaining variables of p are treated as algebraic
		 *
		 * @param p: The polynomial to add
		 * @param ts: The transcendental variable names
		 **/
		RegularChain<Field,RecursivePoly> (const RecursivePoly& p, const std::vector<Symbol>& ts);

//		/**
//		 * Merge constructor
//		 *
//		 * @param a: A regular chain
//		 **/
//		RegularChain<Field,RecursivePoly> (const RegularChain<Field,RecursivePoly>& a, const RegularChain<Field,RecursivePoly>& b);
		
		/**
		 * Copy constructor
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly> (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Copy constructor
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly> (const RegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Copy constructor
		 *
		 * @param a: A triangular set
		 **/
		RegularChain<Field,RecursivePoly> (const TriangularSet<Field,RecursivePoly>& a);
		
		/**
		 * Move constructor
		 *
		 * @param a: An r-value reference regular chain
		 **/
		RegularChain<Field,RecursivePoly> (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);
		
		/**
		 * Move constructor
		 *
		 * @param a: An r-value reference regular chain
		 **/
		RegularChain<Field,RecursivePoly> (RegularChain<Field,RecursivePoly>&& a);
		
		/**
		 * Move constructor
		 *
		 * @param a: An r-value reference triangular set
		 **/
		RegularChain<Field,RecursivePoly> (TriangularSet<Field,RecursivePoly>&& a);

		/**
		 * Computational constructor: creates a regular chain given all the data
		 *
		 * @param vs: variables of the regular chain
		 * @param avs: algebraic variables of the regular chain
		 * @param tvs: transcendental variables of the regular chain
		 * @param polys: polynomials of the regular chain
		 * @param tsm: whether the regular chain is variable or fixed
		 * @param c: characteristic of the regular chain
		 * @param normal: whether the regular chain is strongly normalized
		 **/
		RegularChain<Field,RecursivePoly> (const std::vector<Symbol>& vs, const std::vector<Symbol>& avs, const std::vector<Symbol>& tvs, const std::vector<RecursivePoly>& ts, TriangularSetMode tsm, const mpz_class& c);
		
		/**
		 * Construct a set of regular chains from an input triangular set
		 *
		 * @param T: A triangular set
		 **/
		static std::vector<RegularChain<Field,RecursivePoly>> constructChains(const TriangularSet<Field,RecursivePoly>& T);
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (const RegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (const BPASTriangularSet<Field,RecursivePoly>& a) override;
		
		/**
		 * Assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (const BPASRegularChain<Field,RecursivePoly>& a) override;
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (RegularChain<Field,RecursivePoly>&& a);
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (BPASTriangularSet<Field,RecursivePoly>&& a) override;
		
		/**
		 * Move assignment operator =
		 *
		 * @param a: A regular chain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (BPASRegularChain<Field,RecursivePoly>&& a) override;
		
		/**
		 * Add operator +
		 * Adds a polynomial to a regular chain and returns a new regular chain
		 *
		 * @param p: A sparse multivariate polynomial
		 **/
		RegularChain<Field,RecursivePoly> operator+ (const RecursivePoly& p) const;
		
		/**
		 * Add assignment operator +=
		 * Adds a polynomial to a regular chain
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		RegularChain<Field,RecursivePoly>& operator+= (const RecursivePoly& p);
		
		/**
		 * Add operator +
		 * Adds a regular chain to a regular chain and returns a new regular chain
		 *
		 * @param p: A regular chain
		 **/
		RegularChain<Field,RecursivePoly> operator+ (const RegularChain<Field,RecursivePoly>& T) const;
		
		/**
		 * Add assignment operator +=
		 * Adds a regular chain to a regular chain
		 *
		 * @param p: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator+= (const RegularChain<Field,RecursivePoly>& T);
		
		/**
		 * Identity operator ==
		 *
		 *
		 * @param a: A regular chain
		 **/
		bool operator== (RegularChain<Field,RecursivePoly>& a);

		/**
		 * Negated identity operator !=
		 *
		 *
		 * @param a: A regular chain
		 **/
		bool operator!= (RegularChain<Field,RecursivePoly>& a);
		
		/**
		 * Get the number of variables
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return TriangularSet<Field,RecursivePoly>::numberOfVariables();
		}
		
		/**
		 * Get the encoded options of the regular chain
		 *
		 **/
		int options() const;
		
		/**
		 * Set the encoded options of the regular chain
		 *
 		* @param opts: bitwise or of RegularChainOption values
		 **/
		void setOptions(int opts);

		
		/**
		 * Find out if the regular chain is known to be squarefree
		 *
		 **/
		bool isKnownToBeSquareFree() const;

		
//		/**
//		 * Find out if the regular chain is known to be irreducible
//		 *
//		 **/
//		bool isKnownToBeIrreducible() const;
		
		/**
		 * Find out if the input polynomial is in the saturated ideal of the current regular chain
		 *
		 * @param p: polynomial
		 **/
		bool isInSaturatedIdeal(const RecursivePoly& p) const;
		
		/**
		 * Find out if the input polynomial is in the saturated ideal of the current regular chain
		 * and return the reduced input polynomial
		 *
		 * @param p: polynomial
		 **/
		bool isInSaturatedIdeal(const RecursivePoly& p, RecursivePoly& redp) const;
		
		/**
		 * Find out if the input polynomial is regular modulo the current regular chain
		 *
		 * @param p: polynomial
		 **/
		 bool isRegular(const RecursivePoly& p) const;
		
		
		
//		using TriangularSet<Field,RecursivePoly>::numberOfVariables;
		
		/**
		 * Get the size of the regular chain
		 *
		 * @param
		 **/
//		inline int size() const {
//			return TriangularSet<Field,RecursivePoly>::size();
//		}
		
		/**
		 * Get the number of algebraic variables
		 *
		 * @param
		 **/
//		inline int numberOfAlgebraicVariables() const {
//			return TriangularSet<Field,RecursivePoly>::numberOfAlgebraicVariables();
//		}
		
		/**
		 * Get the number of transcendental variables
		 *
		 * @param
		 **/
//		inline int numberOfTranscendentalVariables() const {
//			return TriangularSet<Field,RecursivePoly>::numberOfAlgebraicVariables();
//		}
		
		/**
		 * Get the variable names in decreasing order
		 *
		 * @param
		 **/
		inline std::vector<Symbol> variables() const {
			return TriangularSet<Field,RecursivePoly>::variables();
		}
		
		/**
		 * Get algebraic variables
		 *
		 * @param
		 **/
//		inline std::vector<Symbol> mainVariables() const {
//			return TriangularSet<Field,RecursivePoly>::mainVariables();
//		}
		
		/**
		 * Get transcendentalVariables variables
		 *
		 * @param
		 **/
//		inline std::vector<Symbol> transcendentalVariables() const {
//			return TriangularSet<Field,RecursivePoly>::transcendentalVariables();
//		}
		
//		inline bool isAlgebraic(const Symbol& s) const {
//			return TriangularSet<Field,RecursivePoly>::isAlgebraic(s);
//		}
		
//		inline bool isEmpty() const {
//			return TriangularSet<Field,RecursivePoly>::isEmpty();
//		}
		
//		inline std::vector<RecursivePoly> polynomials() const {
//			return TriangularSet<Field,RecursivePoly>::polynomials();
//		}
		
		inline RecursivePoly select(const Symbol& s) const {
			return TriangularSet<Field,RecursivePoly>::select(s);
		}
		
		void lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;
		
		void upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;
		
		void lowerSlice(const Symbol& s);
		
		std::vector<RecursivePoly> GCDFreeFactorization(const RecursivePoly& p, int type = 0) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> triangularize(const std::vector<RecursivePoly>& F);
		
		std::vector<RegularChain<Field,RecursivePoly>> intersect(const RecursivePoly& p) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> intersectFree(const RecursivePoly& p, const Symbol& v) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> intersectAlgebraic(const RecursivePoly& p, const RegularChain<Field,RecursivePoly>& T, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src) const;
//		std::vector<RegularChain<Field,RecursivePoly>> intersectAlgebraic(const RecursivePoly& p, const RegularChain<Field,RecursivePoly>& T, const Symbol& v, const SubResultantChain<Field,SparseUnivariatePolynomial<RecursivePoly>>& src) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> cleanChain(const RegularChain<Field,RecursivePoly>& C, const Symbol& v) const;
		
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularizeSingle(const RecursivePoly& p) const;
		
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularize(const RecursivePoly& p) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> extend(const RecursivePoly& p, const Symbol& v) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> extend(const std::vector<RecursivePoly>& T, const Symbol& v) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> extend(const RegularChain<Field,RecursivePoly>& T, const Symbol& v) const;
		
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src) const;
//		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<Field,SparseUnivariatePolynomial<RecursivePoly>>& src) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> squareFreePart(const RecursivePoly& p, const Symbol& v, int options=ASSUME_REGULAR) const;
		
		std::vector<RegularChain<Field,RecursivePoly>> squareFreePart(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, int options=ASSUME_REGULAR) const;
//		std::vector<RegularChain<Field,RecursivePoly>> squareFreePart(const RecursivePoly& p, const Symbol& v, const SubResultantChain<Field,SparseUnivariatePolynomial<RecursivePoly>>& src, int options=ASSUME_REGULAR) const;
		
//		std::vector<RegularChain<Field,RecursivePoly>> squareFreePart() const;

		int options();
		

		/**
		 * Generate a random regular chain
		 *
		 * @param nVars: number of variables
		 * @param nAlgVars: number of algebraic variables
		 * @param nTrcVars: number of transcendental variables
		 * @param nTerms: number of terms in each polynomial
		 * @param coefBound: bound on the coefficient size
		 * @param pSparsity: sparsity of the polynomial in the term separation sense
		 * @param includeNeg: whether to include negative coefficients
		 *
		 **/
		void randomRegularChain(int nVars, int nAlgVars, int nTrcVars, int nTerms, unsigned long int coefBound, int pSparsity, bool includeNeg);
		
		/**
		 * Generate a random regular chain
		 *
		 * @param nVars: number of variables
		 * @param nAlgVars: number of algebraic variables
		 * @param nTrcVars: number of transcendental variables
		 * @param maxDegs: maximum degrees among the full list of variables
		 * @param coefBound: bound on the coefficient size
		 * @param pSparsity: sparsity of the polynomial in the term separation sense
		 * @param includeNeg: whether to include negative coefficients
		 *
		 **/
		void randomRegularChain(int nVars, int nAlgVars, int nTrcVars, std::vector<int> maxDegs, unsigned long int coefBound, double pSparsity, bool includeNeg);
		
};

#endif
