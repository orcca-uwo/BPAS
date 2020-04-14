#ifndef _ZERODIMENSIONALREGULARCHAIN_H_
#define _ZERODIMENSIONALREGULARCHAIN_H_

#include "regularchain_macros.hpp"

#include "BPASZeroDimRegularChain.hpp"
#include "../Ring/BPASField.hpp"
#include "regularchain.hpp"
#include "chainstructures.hpp"
#include "../TriangularSet/triangularset.hpp"
#include "../SubResultantChain/subresultantchain.hpp"



extern long long unsigned int rcProfilingStart;
extern float primitivePartTime;
extern float squareFreePartTime;
extern float subresultantChainTime;
extern float zerodimensionalregularchainTime;
extern float pseudoDivideTime;
extern float normalFormTime;

extern float zdrcCopyTime;

/**
 * A class for handling regular chains in dimension zero.
 * A ZeroDimensionalRegularChain contains polynomials of type BPASRecursivelyViewedPolynomial, which have coefficients in a BPASField.
 **/
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
		using RegularChain<Field,RecursivePoly>::reduceMinimal;

		/**
		 * Cut an input zero-dimensional regular chain at the symbol v, returning the zero-dimensional subchain below v,
		 * the polynomial with main variable v and the potentially positive-dimensional subchain above v.
		 *
		 * @param T: a regular chain
		 * @param v: a symbol
		 * @param Tlv: the zero-dimensional subchain of T below v
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 * @param Tgv: the (potentially positive-dimensional) subchain of T above v
		 **/
		void cutChain(const ZeroDimensionalRegularChain<Field,RecursivePoly>& T, const Symbol& v, ZeroDimensionalRegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const;

		/**
		 * Cut the current object at the symbol v, returning the polynomial with main variable v and the subchain above v, which may no longer be zero-dimensional.
		 *
		 * @param v: a symbol
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 * @param Tgv: the regular chain obtained from the polynomials of the current object above v
		 **/
		void cutChain(const Symbol& v, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const;

		/**
		 * Cut the current object at the symbol v, returning the zero-dimensional subchain below v and the polynomial with main variable v.
		 *
		 * @param v: a symbol
		 * @param Tlv: the zero-dimensional subchain of the current object below v
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 **/
		void cutChain(const Symbol& v, ZeroDimensionalRegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv) const;

		/**
		 * Add an input polynomial to the current object
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		void constructChain(const RecursivePoly& p, int options=ASSUME_REGULAR);

		/**
		 * Add an an upper regular chain to the current object
		 *
		 * @param T: a recursively viewed polynomial
		 **/
		void constructChain(const RegularChain<Field,RecursivePoly>& T, int options=ASSUME_REGULAR);

		/**
		 * Construct a set of regular chains from the current object and an input polynomial
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> constructChains(const RecursivePoly& p, int options=ASSUME_REGULAR) const;

		/**
		 * Construct a set of regular chains from the current object and an input polynomial
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> constructChains(const RegularChain<Field,RecursivePoly>& T, int options=ASSUME_REGULAR) const;

		/*******************************
		 * Generator methods
		 *
		 *******************************/

		/**
		 * Compute the last nonzero subresultant of the subresultant chain src of f and g modulo the current regular chain,
		 * i.e., compute a splitting of pairs (R_i,T_i) such that R_i is the last nonzero subresultant modulo T_i.
		 *
		 * @param f: input polynomial
		 * @param g: input polynomial
		 * @param v: common main variable of f and g
		 * @param assumeRegular: boolean flag specifying whether the input polynomials are assumed to be regular modulo the current chain.
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void lastNonZeroSubResultant(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, bool assumeRegular, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> lastNonZeroSubResultant(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, bool assumeRegular) const;
		#endif

		/**
		 * Compute the last nonzero subresultant of the subresultant chain src of f and g modulo the current regular chain,
		 * i.e., compute a splitting of pairs (R_i,T_i) such that R_i is the last nonzero subresultant modulo T_i.
		 *
		 * @param f: input polynomial
		 * @param g: input polynomial
		 * @param v: common main variable of f and g
		 * @param src: subresultant chain of f and g
		 * @param assumeRegular: boolean flag specifying whether the input polynomials are assumed to be regular modulo the current chain.
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void lastNonZeroSubResultant_inner(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool assumeRegular, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> lastNonZeroSubResultant_inner(const RecursivePoly& f, const RecursivePoly& g, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool assumeRegular) const;
		#endif

		/**
		 * Compute a decomposition of the current object such that on each component the initial of the input polynomial is either
		 * zero or regular modulo the saturated ideal of that component. If the initial is zero, (0, T_i) is returned;
		 * if the initial is regular, (p_i, T_i), is returned, where p_i is equivalent to p modulo the saturated ideal of T_i and
		 * p_i is reduced modulo Sat(T_i).
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _regularizeInitial(const RecursivePoly& f, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>&) const;
		void _isInvertible(const RecursivePoly& f, AsyncGenerator<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> _regularizeInitial(const RecursivePoly& p) const;
		std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> _isInvertible(const RecursivePoly& f) const;
		#endif

	public:

		// TODO: the underscore routines should be protected+friends

		/**
		 * Compute the common solutions of the input polynomial and the current regular chain.
		 * More precisely, this routine computes the intersection of the varieties of the input polynomial and the
		 * current zero-dimensional regular chain, expressed as a set of zero-dimensional regular chains.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _intersect(const RecursivePoly& p, AsyncGenerator<ZeroDimensionalRegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> _intersect(const RecursivePoly& p) const;
		#endif

		/**
		 * Compute the gcd of two input polynomials p and q modulo the saturated ideal of the current object. The result is a list of pairs (g_i,T_i) of
		 * polynomials g_i and regular chains T_i, where g_i is gcd(p,q) modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param q: a recursively viewed polynomial
		 * @param v: common main variable of f and g
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results);
		#else
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> _regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v);
		#endif

		/**
		 * Compute a decomposition of the current object such that on each component the input polynomial is either
		 * zero or regular modulo the saturated ideal of that component. If the polynomial is zero, (0, T_i) is returned;
		 * if the polynomial is regular, (p_i, T_i), is returned, where p_i is equivalent to p modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _regularize(const RecursivePoly& p, AsyncGenerator<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> _regularize(const RecursivePoly& p) const;
		#endif

		using RegularChain<Field,RecursivePoly>::allVariables;

		/**
		 * Default constructor: creates an empty "zero-dimensional" regular chain with variable size
		 * with empty list of transcendentals.
		 *
		 * @param
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>();

		/**
		 * Construct an empty zero-dimensional regular chain of variable size with
		 * variables given by xs with empty list of transcendentals.
		 *
		 * @param ps: the transcendental variable names
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const std::vector<Symbol>& ps);

		/**
		 * Construct a zero-dimensional regular chain of variable size containing a univariate polynomial p
		 * with empty list of transcendentals.
		 *
		 * @param p: a recursively viewed polynomial with only one variable
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const RecursivePoly& p);

		/**
		 * Construct a variable zero-dimensional regular chain containing p, such that the supplied list of transcendental variable names
		 * includes all of the variables in p except its main variable, which becomes the only algebraic variable of the chain.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param ts: the transcendental variable names
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const RecursivePoly& p, const std::vector<Symbol>& ts);

		/**
		 * Construct a fixed zero-dimensional regular chain containing the polynomials in polys.
		 * @note: It is assumed that the polynomials in polys form a valid zero-dimensional regular chain.
		 *
		 * @param polys: a list of recursively viewed polynomials.
		 *
		 */
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const std::vector<RecursivePoly> polys);


		/**
		 * Copy constructor.
		 *
		 * @param a: a zero-dimensional regular chain
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);

		/**
		 * Copy constructor taking a regular chain as input. Note that this routine allow the construction of fixed zero-dimensional regular
		 * chains that are technically positive-dimensional because the polynomials above the maximum algebraic variable are zero.
		 * This is allowed to optimize the performance of routines in the RegularChain class.
		 *
		 * @param a: a regular chain
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const RegularChain<Field,RecursivePoly>& a, int options=0);
		// TODO: make this protected

		/**
		 * Move constructor.
		 *
		 * @param a: an r-value reference zero-dimensional regular chain
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);

		/**
		 * Move constructor taking a regular chain as input.  Note that this routine allow the construction of fixed zero-dimensional regular
		 * chains that are technically positive-dimensional because the polynomials above the maximum algebraic variable are zero.
		 * This is allowed to optimize the performance of routines in the RegularChain class.
		 *
		 * @param a: an r-value reference regular chain
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (RegularChain<Field,RecursivePoly>&& a, int options=0);
		// TODO: make this protected

		/**
		 * Computational constructor: creates a triangular set given all the data
		 *
		 * @param vs: rvalue reference to variables of the zero-dimensional regular chain
		 * @param avs: rvalue reference to algebraic variables of the zero-dimensional regular chain
		 * @param tvs: rvalue reference to transcendental variables of the zero-dimensional regular chain
		 * @param polys: rvalue reference to polynomials of the zero-dimensional regular chain
		 * @param tsm: whether the zero-dimensional regular chain is variable or fixed
		 * @param c: characteristic of the zero-dimensional regular chain
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> (const std::vector<Symbol>&& vs, const std::vector<Symbol>&& avs, const std::vector<Symbol>&& tvs, const std::vector<RecursivePoly>&& ts, TriangularSetMode tsm, const mpz_class& c);
		// TODO: do we really need this?

		/**
		 * Assignment operator =.
		 *
		 * @param a: a triangular set
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);

		/**
		 * Assignment operator = imposed by abstract class BPASTriangularSet.
		 *
		 * @param a: a triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (const BPASTriangularSet<Field,RecursivePoly>& a) override;

		/**
		 * Assignment operator = imposed by abstract class BPASRegularChain.
		 *
		 * @param a: a regular chain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (const BPASRegularChain<Field,RecursivePoly>& a) override;

		/**
		 * Assignment operator = imposed by abstract class BPASZeroDimensionalRegularChain.
		 *
		 * @param a: a zero-dimensional regular chain
		 **/
		BPASZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (const BPASZeroDimensionalRegularChain<Field,RecursivePoly>& a) override;

		/**
		 * Move assignment operator = imposed by class ZeroDimensionalRegularChain.
		 *
		 * @param a: a zero-dimensional regular chain
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);

		/**
		 * Move assignment operator = imposed by abstract class BPASTriangularSet.
		 *
		 * @param a: a triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (BPASTriangularSet<Field,RecursivePoly>&& a) override;

		/**
		 * Move assignment operator =
		 *
		 * @param a: a BPASRegularChain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (BPASRegularChain<Field,RecursivePoly>&& a) override;

		/**
		 * Move assignment operator = imposed by abstract class BPASZeroDimensionalRegularChain.
		 *
		 * @param a: a zero-dimensional regular chain
		 **/
		BPASZeroDimensionalRegularChain<Field,RecursivePoly>& operator= (BPASZeroDimensionalRegularChain<Field,RecursivePoly>&& a) override;

		/**
		 * Add operator +:
		 * Adds a polynomial p to a zero-dimensional regular chain and returns a new zero-dimensional regular chain,
		 * assuming that the main variable of p is neither algebraic nor transcendental, p contains no other non-transcendental
		 * variables, and that init(p) is regular modulo the saturated ideal of the current object.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> operator+ (const RecursivePoly& p) const;

		/**
		 * Add assignment operator +=:
		 * Adds a polynomial p to a zero-dimensional regular chain, assuming that the main variable of p
		 * is neither algebraic nor transcendental, p contains no other non-transcendental
		 * variables, and that init(p) is regular modulo the saturated ideal of the current object.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator+= (const RecursivePoly& p);

		/**
		 * Add operator +:
		 * Adds the polynomials of an input regular chain to the current object and returns a new zero-dimensional
		 * regular chain, assuming that the result of the addition is both a regular chain and zero-dimensional.
		 *
		 * @param p: a regular chain
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly> operator+ (const RegularChain<Field,RecursivePoly>& T) const;

		/**
		 * Add assignment operator +=:
		 * Adds the polynomials of an input regular chain to the current object, assuming that the result of
		 * the addition is both a regular chain and zero-dimensional.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		ZeroDimensionalRegularChain<Field,RecursivePoly>& operator+= (const RegularChain<Field,RecursivePoly>& T);

		/**
		 * Identity operator ==.
		 *
		 *
		 * @param a: A triangular set
		 **/
		bool operator== (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) const;

		/**
		 * Negated identity operator !=.
		 *
		 *
		 * @param a: A triangular set
		 **/
		bool operator!= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a) const;

		/**
		 * Get the number of (potentially algebraic) variables. This can only be different than the number
		 * of algebraic variables for fixed zero-dimensional regular chains.
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return RegularChain<Field,RecursivePoly>::numberOfVariables();
		}

		/**
		 * Get the number of algebraic variables in the current object.
		 *
		 * @param
		 **/
		inline int numberOfAlgebraicVariables() const {
			return RegularChain<Field,RecursivePoly>::numberOfAlgebraicVariables();
		}

		/**
		 * Get the number of transcendental variables in the current object.
		 *
		 * @param
		 **/
		inline int numberOfTranscendentalVariables() const {
			return RegularChain<Field,RecursivePoly>::numberOfTranscendentalVariables();
		}

		/**
		 * Get the (potentially algebraic) variable names for the current object in decreasing order.
		 * This can only be different from the list of algebraic variable names for fixed zero-dimensional
		 * regular chains.
		 *
		 * @param
		 **/
		inline std::vector<Symbol> variables() const {
			return RegularChain<Field,RecursivePoly>::variables();
		}

		/**
		 * Get algebraic variables in the current object.
		 *
		 * @param
		 **/
		inline std::vector<Symbol> mainVariables() const {
			return RegularChain<Field,RecursivePoly>::mainVariables();
		}

		/**
		 * Get transcendental variables in the current object.
		 *
		 * @param
		 **/
		inline std::vector<Symbol> transcendentalVariables() const {
			return RegularChain<Field,RecursivePoly>::transcendentalVariables();
		}

		/**
		 * Find out if the input symbol is an algebraic variable of the current object.
		 *
		 * @param s: a symbol
		 **/
		inline bool isAlgebraic(const Symbol& s) const {
			return RegularChain<Field,RecursivePoly>::isAlgebraic(s);
		}

		/**
		 * Find out if the current object is the empty chain.
		 *
		 * @param
		 **/
		inline bool isEmpty() const {
			return RegularChain<Field,RecursivePoly>::isEmpty();
		}

		/**
		 * Get the list of polynomials in the current object.
		 *
		 * @param
		 **/
		inline std::vector<RecursivePoly> polynomials() const {
			return RegularChain<Field,RecursivePoly>::polynomials();
		}

		/**
		 * Select the polynomial in the current object with main variable s, if it exists.
		 *
		 * @param s: a symbol
		 **/
		inline RecursivePoly select(const Symbol& s) const {
			return RegularChain<Field,RecursivePoly>::select(s);
		}

		/**
		 * Returns the zero-dimensional regular chain consisting of polynomials with
		 * main variable strictly less than s. NB: the type of the returned object is always
		 * a ZeroDimensionalRegularChain, since the lower chain is always genuinely zero-dimensional;
		 * however, if the current object is a variable zero-dimensional regular chain, the returned
		 * object really is zero-dimensional, in the sense that all variables are algebraic, but if
		 * the current object is a fixed type, then the chain may only be morally zero-dimensional
		 * in the sense that all variables below some variable are algebraic but the chain has
		 * non-algebraic variables.
		 *
		 * @param s: symbol of the main variable of specified element of the regular chain
		 * @param ts: The returned regular chain
		 **/
		void lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;

		/**
		 * Returns the regular chain consisting of polynomials with
		 * main variable strictly greater than s. The type of the returned object is always
		 * a RegularChain, since the upper chain can be genuinely positive dimensional.
		 *
		 * @param s: symbol of the main variable of specified element of the regular chain
		 * @param ts: The returned regular chain
		 **/
		void upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;

		/**
		 * Compute the common solutions of the input polynomial and the current regular chain.
		 * More precisely, this routine computes the intersection of the varieties of the input polynomial and the
		 * current zero-dimensional regular chain, expressed as a set of zero-dimensional regular chains.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<ZeroDimensionalRegularChain<Field,RecursivePoly>> intersect(const RecursivePoly& p) const;

		/**
		 * Compute a decomposition of the current object such that on each component the input polynomial is either
		 * zero or regular modulo the saturated ideal of that component. If the polynomial is zero, (0, T_i) is returned;
		 * if the polynomial is regular, (p_i, T_i), is returned, where p_i is equivalent to p modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> regularize(const RecursivePoly& p) const;

		/**
		 * Compute a decomposition of the current object such that on each component the initial of the input polynomial is either
		 * zero or regular modulo the saturated ideal of that component. If the initial is zero, (0, T_i) is returned;
		 * if the initial is regular, (p_i, T_i), is returned, where p_i is equivalent to p modulo the saturated ideal of T_i and
		 * p_i is reduced modulo Sat(T_i).
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> regularizeInitial(const RecursivePoly& p) const;

		/**
		 * Determine whether a recursively viewed polynomial is invertible with respect to the current regular chain.
		 * More precisely, compute a splitting consisting of pairs (b_i,T_i) such that b_i is true if f is invertible
		 * modulo T_i and b_i is false if f is zero modulo T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<BoolChainPair<ZeroDimensionalRegularChain<Field,RecursivePoly>>> isInvertible(const RecursivePoly& p) const;

		/**
		 * Compute the gcd of two input polynomials p and q modulo the saturated ideal of the current object. The result is a list of pairs (g_i,T_i) of
		 * polynomials g_i and regular chains T_i, where g_i is gcd(p,q) modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param q: a recursively viewed polynomial
		 * @param v: common main variable of f and g
		 **/
		std::vector<PolyChainPair<RecursivePoly,ZeroDimensionalRegularChain<Field,RecursivePoly>>> regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v);

		/**
		 * Generate a random zero-dimensional regular chain based on the number of terms of the polynomials of the chain.
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
		 * Generate a random zero-dimensional regular chain based on a list of maximum degrees of variables in the polynomials in the chain.
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
