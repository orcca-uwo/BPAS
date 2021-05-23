#ifndef _REGULARCHAIN_H_
#define _REGULARCHAIN_H_

#include "regularchain_macros.hpp"

#include "../TriangularSet/triangularset.hpp"
#include "../SubResultantChain/subresultantchain.hpp"
#include "chainstructures.hpp"

#include "BPASRegularChain.hpp"

#include <memory>

/* Debugging utils  */

extern int intersectDepth;
extern int regularGCDDepth;
extern int intersectFreeDepth;
extern int intersectAlgebraicDepth;
extern int regularizeDepth;
extern int regularizeSingleDepth;
extern int extendDepth;
extern int cleanChainDepth;
extern int depth;

extern float primitivePartTime;
extern float squareFreePartTime;
extern float subresultantChainTime;
extern float zerodimensionalregularchainTime;
extern float pseudoDivideTime;
extern float normalFormTime;
extern float removeRedundantChainsTime;
extern float factorTime;

extern float intersectTime;
extern float regularGCDTime;
extern float intersectFreeTime;
extern float intersectAlgebraicTime;
extern float regularizeTime;
extern float extendTime;
extern float cleanChainTime;
extern float constructChainTime;
extern float constructChainsTime;
extern float GCDFreeFactorizationTime;
extern float ZDIntersectTime;
extern float ZDRegularizeTime;
extern float isRegularTime;
extern float isInSaturatedIdealTime;
extern float cleanSetTime;

extern float tsCopyTime;
extern float rcCopyTime;
extern float zdrcCopyTime;

#define RC_TRIANGULARIZE_TASKTREEDATA 0


/* forward declares */

template <class Field, class RecursivePoly>
class ZeroDimensionalRegularChain;

template <class Value>
class SynchronizedWriteVector;

class TaskScheduler;

template <class Field, class RecursivePoly>
class RegularChain;

/* internal helper functions */

template<class Field, class RecursivePoly>
void triangularizeTask(const RegularChain<Field, RecursivePoly>& rc, std::vector<RecursivePoly>& polys, bool lazardDecompose, int heightBound, TaskScheduler* tasks, std::shared_ptr<SynchronizedWriteVector<RegularChain<Field,RecursivePoly>>> results);

template<class Field, class RecursivePoly>
void intersectOne(int j, const RegularChain<Field, RecursivePoly>& T, const RecursivePoly& p, int lazardDecompose, int heightBound, SynchronizedWriteVector<RegularChain<Field,RecursivePoly>>&);


/**
 * An enumeration for various options passed to certain routines of the regular chain classes.
 **/
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
	CONSTRUCT_FACTORIZE    = 0x200,
};

/**
 * A class for handling regular chains of arbitrary dimension.
 * A RegularChain contains polynomials of type BPASRecursivelyViewedPolynomial, which have coefficients in a BPASField.
 **/
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

#if defined(RC_TRIANGULARIZE_TASKTREEDATA) && RC_TRIANGULARIZE_TASKTREEDATA
		mutable int RegChain_UniqueID;
#endif

		int regularChainOptions = CONSTRUCT_FACTORIZE;
		// int regularChainOptions = MAINTAIN_SQUAREFREE | CONSTRUCT_FACTORIZE;

		//TODO: These don't really work. They only count up from the empty chain.
		//Ex: RC over [x>y>z] with polys [x,y]. Adding a non-square free poly in z would make squarefree-ness = 0
		int squareFreeLevel;
		int saturatedIdealPrimeLevel;

		/**
		 * Determine whether the regular chain is strongly normalized and the maximum variable
		 * of the largest strongly normalized subset of the regular chain.
		 *
		 * @param
		 **/
		void updateRegularChainStates();

		/**
		 * Determine whether after adding p the regular chain is strongly normalized and update the maximum variable
		 * of the largest strongly normalized subset of the triangular set.
		 *
		 * @param p: a recursively viewed polynomial that has just been added to the chain
		 **/
		void updateRegularChainStates(const RecursivePoly& p);

//		/**
//		 * A protected internal routine used by triangularize to
//		 *
//		 * @param p: a recursively viewed polynomial that has just been added to the set
//		 **/
//		std::vector<RegularChain<Field,RecursivePoly>> triangularize_inner(std::vector<RecursivePoly>& F) const;

		/**
		 * Compute a canonical representation of the input polynomial.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		RecursivePoly makeCanonical(const RecursivePoly& p_in) const;

		/**
		 * Reduce the input polynomial with respect to the polynomials in the the current object with constant
		 * initial and remove parts of the result that reduce to zero modulo the current object. This is a less
		 * involved computation than reduce from the TriangularSet class.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		RecursivePoly reduceMinimal(const RecursivePoly& p) const;

		/**
		 * Do a minimal reduction (calling reduceMinimal) and take the mainPrimitivePart.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		RecursivePoly reduceMinimalPrimitivePart(const RecursivePoly& p) const;

		/**
		 * Do a minimal reduction (calling reduceMinimal) and take the mainPrimitivePart and squareFreePart.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		RecursivePoly reduceMinimalPrimitiveSquareFreePart(const RecursivePoly& p_in) const;

		/**
		 * Reduce the input polynomial such that p = 0 iff p is in the saturated ideal of the current object and such that the p - ret is in the same saturated ideal, where ret is the returned value.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		RecursivePoly removeZero(const RecursivePoly& p) const;
		// TODO: move this to the TriangularSet class.

		/**
		 * Cut an input regular chain at the symbol v, returning the subchain below v, the polynomial with main variable v and the subchain above v.
		 *
		 * @param T: a regular chain
		 * @param v: a symbol
		 * @param Tlv: the subchain of T below v
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 * @param Tgv: the subchain of T above v
		 **/
		void cutChain(const RegularChain<Field,RecursivePoly>& T, const Symbol& v, RegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const;

		/**
		 * Cut the current object at the symbol v, returning the polynomial with main variable v and the subchain above v.
		 *
		 * @param v: a symbol
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 * @param Tgv: the subchain of the current object above v
		 **/
		void cutChain(const Symbol& v, RecursivePoly& Tv, RegularChain<Field,RecursivePoly>& Tgv) const;

		/**
		 * Cut the current object at the symbol v, returning the subchain below v and the polynomial with main variable v.
		 *
		 * @param v: a symbol
		 * @param Tlv: the subchain of the current object below v
		 * @param Tv: the recursively viewed polynomial with main variable v, if it exists
		 **/
		void cutChain(const Symbol& v, RegularChain<Field,RecursivePoly>& Tlv, RecursivePoly& Tv) const;

		/**
		 * Cut the current object at the symbol v and return the subchain above v.
		 *
		 * @param v: a symbol
		 * @param[out] Tgv: the upper subchain
		 */
		void upper(const Symbol& v, RegularChain<Field, RecursivePoly>& Tgv) const;

		/**
		 * Cut the current object at the symbol v and return the subchain below v.
		 *
		 * @param v: a symbol
		 * @param[out] Tlv: the lower subchain
		 */
		void lower(const Symbol& v, RegularChain<Field, RecursivePoly>& Tlv) const;

		/**
		 * A protected internal routine to regularize a list of polynomials with respect to the current object.
		 *
		 * @param knownRegular: a list of polynomials known to be regular modulo the saturated ideal of the current object
		 * @param unknownIfRegular: a list of polynomials for which it is not known if they are regular or not wrt the current object
		 * @param singularComponents: a list of regular chains consisting of those components of the current object over which some element
		 * of both knownRegular and unknownIfRegular is zero
		 * @param regularComponents: a list of regular chains consisting of those components of the current object over which each element
		 * of both knownRegular and unknownIfRegular is regular
		 **/
		void regularizeList(const std::vector<RecursivePoly>& knownRegular, const std::vector<RecursivePoly>& unknownIfRegular, std::vector<RegularChain<Field,RecursivePoly>>& singularComponents, std::vector<RegularChain<Field,RecursivePoly>>& regularComponents, bool lazardDecompose, int heightBound) const;

		/**
		 * A protected internal routine to regularize a list of polynomials with respect to the current object.
		 *
		 * @param knownRegular: a list of polynomials known to be regular modulo the saturated ideal of the current object
		 * @param unknownIfRegular: a list of polynomials for which it is not known if they are regular or not wrt the current object
		 * @param results: a list of pairs of polynomials and regular chains
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void regularizeListStream(const RecursivePoly& p, const std::vector<RecursivePoly>& knownRegular, const std::vector<RecursivePoly>& unknownIfRegular, bool lazardDecompose, int heightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularizeListStream(const RecursivePoly& p, const std::vector<RecursivePoly>& knownRegular, const std::vector<RecursivePoly>& unknownIfRegular, bool lazardDecompose, int heightBound) const;
		#endif

		/**
		 * Add an input polynomial to the current object
		 *
		 * @param p: A recursively viewed polynomial
		 * @param options: a bitwise or of RegularChainOption values (default assume that init(p) is regular modulo the saturated ideal of the current object)
		 **/
		void constructChain(const RecursivePoly& p, int options=ASSUME_REGULAR);

		/**
		 * Construct a regular chain from the current object and an input regular chain known to have polynomials
		 * with main variable greater than any in the current object.
		 *
		 * @param T: An upper regular chain
		 * @param options: a bitwise or of RegularChainOption values
		 **/
		void constructChain(const RegularChain<Field,RecursivePoly>& T, int options);

		/**
		 * Construct a set of regular chains from the current object and an input polynomial known to have main
		 * variable greater than any in the current object.
		 * Here, the parameter options is information about the p being added.
		 *
		 * @param p: A recursively viewed polynomial
		 * @param options: a bitwise or of RegularChainOption values (default assume init(p) is regular modulo the saturated ideal of current object)
		 **/
		 #if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void constructChainsFromPoly(const RecursivePoly& p, bool lazardDecompose, int heightBound, int options, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> constructChainsFromPoly(const RecursivePoly& p, bool lazardDecompose, int heightBound, int options=ASSUME_REGULAR) const;
		#endif

		/**
		 * Construct a set of regular chains from the current object and an input regular chain above the current object, when splitting can occur to impose
		 * the condition that the regular chains are squarefree.
		 * Here, the parameter options is information about the polys in the chain being added on top.
		 *
		 * @param T: An upper regular chain
		 * @param options: a bitwise or of RegularChainOption values (default assume that the initials of T are regular modulo the saturated ideal of the current object)
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> constructChainsFromChain(const RegularChain<Field,RecursivePoly>& T, bool lazardDecompose, int heightBound, int options) const;

		/**
		 * Return the intersection of the zero sets of an input polynomial and the current regular chain, when the main variable
		 * of the input is not algebraic.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param v: the main variable of p
		 **/
		 #if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void intersectFree(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> intersectFree(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound) const;
		#endif

		/**
		 * Return the intersection of the zero sets of an input polynomial, an input regular chain and the current regular chain, when the main variable
		 * of the input is algebraic.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param T: a regular chain acting as an additional constraint
		 * @param v: the main variable of p
		 * @param src: the subresultant chain of p and the polynomial of the current object with main variable v
		 **/
		 #if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void intersectAlgebraic(const RecursivePoly& p, const RegularChain<Field,RecursivePoly>& T, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> intersectAlgebraic(const RecursivePoly& p, const RegularChain<Field,RecursivePoly>& T, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound) const;
		#endif

		/**
		 * A routine to clean an input regular chain wrt the current object.
		 * More precisely, if the current object has a polynomial with main variable v, the routine returns a set of regular
		 * chains T_j with larger radical saturated ideal than that of in the input chain C (specializations of C) such  that T_j+Tv is
		 * a regular chain, and the quasi-component of C minus the variety of init(Tv) is contained within the
		 * union of the quasi-components of the T_j; otherwise it returns C.
		 *
		 * @param C: a regular chain with larger radical saturated ideal than that of T<v (specialization of T<v), where T is the current object
		 * @param v: a (potentially algebriac) variable of the current object
		 **/
		 #if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void cleanChain(const RegularChain<Field,RecursivePoly>& C, const Symbol& v, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> cleanChain(const RegularChain<Field,RecursivePoly>& C, const Symbol& v, bool lazardDecompose, int heightBound) const;
		#endif

		/**
		 * A routine to regularize a polynomial modulo the saturated ideal of the current object.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void regularizeSingle(const RecursivePoly& p, bool lazardDecompose, int heightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularizeSingle(const RecursivePoly& p, bool lazardDecompose, int heightBound) const;
		#endif

		/**
		 * Compute a factorization of the input such that the returned vector of polynomials
		 * pairwise gcd 1.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param type: a flag to determine the type of factorization [0: none, 1: primitive, 2: irreducible, 3: squareFree] (default is 1)
		 **/
		std::vector<RecursivePoly> GCDFreeFactorization(const RecursivePoly& p, int type = 2) const;

		/**
		 * A routine to extend the current object to a list of regular chains corresponding to adding the polynomial p and
		 * requiring that p is regular modulo the saturated ideal of the current object.
		 *
		 * @param p: a recursively viewed polynomial with main variable greater than any algebraic variable of the current object
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> extend(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound) const;
//		std::vector<RegularChain<Field,RecursivePoly>> extend(const RecursivePoly& p, const Symbol& v) const;

		/**
		 * A routine to extend the current object with a list of polynomials
		 *
		 * @param T: a vector of recursively viewed polynomials with main variables greater than any algebraic variable of the current object
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> extend(const std::vector<RecursivePoly>& T, const Symbol& v, bool lazardDecompose, int heightBound) const;
//		std::vector<RegularChain<Field,RecursivePoly>> extend(const std::vector<RecursivePoly>& T, const Symbol& v) const;

		/**
		 * A routine to extend the current object with an upper regular chain.
		 *
		 * @param T: a regular chain with algebraic variables greater than any algebraic variable of the current object.
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> extend(const RegularChain<Field,RecursivePoly>& T, const Symbol& v, bool lazardDecompose, int heightBound) const;
//		std::vector<RegularChain<Field,RecursivePoly>> extend(const RegularChain<Field,RecursivePoly>& T, const Symbol& v) const;

		/**
		 * A routine that decomposes the regular chain formed from the current object and an input polynomial into a set of squarefree regular chains.
		 *
		 * @param p: a squarefree recursively viewed polynomial such that init(p) is regular modulo the saturated ideal of the current object
		 * @param v: the main variable of p
		 * @param src: the subresultant chain of p and p'
		 * @param options: a bitwise or of RegularChainOption values (default assume that T+p is a regular chain, where T is the current object)
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _squareFreePartInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> _squareFreePartInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options) const;
		#endif

		/**
		 * A routine that decomposes the regular chain formed from the current object and an input polynomial into a set of squarefree regular chains.
		 *
		 * @param p: a recursively viewed polynomial such that init(p) is regular modulo the saturated ideal of the current object
		 * @param v: the main variable of p
		 * @param options: a bitwise or of RegularChainOption values (default assume that T+p is a regular chain, where T is the current object)
		 **/
		 #if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _squareFreePart(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound, int options, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> _squareFreePart(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound, int options) const;
		#endif

		/**
		 * A routine that decomposes the current object into a set of pairs of polynomials and regular chains, where the
		 * polynomial is the squarefree part of the input modulo the corresponding regular chain.
		 *
		 * @param p: a squarefree recursively viewed polynomial such that init(p) is regular modulo the saturated ideal of the current object
		 * @param v: the main variable of p
		 * @param src: the subresultant chain of p and p'
		 * @param options: a bitwise or of RegularChainOption values (default assume that T+p is a regular chain, where T is the current object)
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _squareFreePartPolynomialInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> _squareFreePartPolynomialInner(const RecursivePoly& p, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose, int heightBound, int options) const;
		#endif

		/**
		 * A routine that decomposes the current object into a set of pairs of polynomials and regular chains, where the
		 * polynomial is the squarefree part of the input modulo the corresponding regular chain.
		 *
		 * @param p: a recursively viewed polynomial such that init(p) is regular modulo the saturated ideal of the current object
		 * @param v: the main variable of p
		 * @param options: a bitwise or of RegularChainOption values (default assume that T+p is a regular chain, where T is the current object)
		 **/
		 #if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _squareFreePartPolynomial(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound, int options, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> _squareFreePartPolynomial(const RecursivePoly& p, const Symbol& v, bool lazardDecompose, int heightBound, int options) const;
		#endif

		/**
//		 * Merge two lists of regular chains that are each pairwise irredundant to form a single list of pairwise irredundant chains.
//		 *
//		 * @param lrc1: a list of pairwise irredundant regular chains
//		 * @param lrc2: a list of pairwise irredundant regular chains
//		 **/
//		static std::vector<RegularChain<Field,RecursivePoly>> mergeIrredundantLists(const std::vector<RegularChain<Field,RecursivePoly>>& lrc1, const std::vector<RegularChain<Field,RecursivePoly>>& lrc2);

//		/**
//		 * Compare two lists of regular chains to find those components of the first list are not contained in a component of an element in the second list.
//		 *
//		 * @param lrc1: the first list of regular chains
//		 * @param lrc2: the second list of regular chains
//		 **/
//		static std::vector<RegularChain<Field,RecursivePoly>> oneSideCompare(const std::vector<RegularChain<Field,RecursivePoly>>& lrc1, const std::vector<RegularChain<Field, RecursivePoly>>& lrc2);
//
//		/**
//		 * Determine which components of a list of regular chains are not contained in a given regular chain
//		 *
//		 * @param lrc: a list of regular chains
//		 * @param rc: a regular chain
//		 **/
//		static std::vector<RegularChain<Field,RecursivePoly>> oneSideCompare(const std::vector<RegularChain<Field,RecursivePoly>>& lrc, const RegularChain<Field,RecursivePoly>& rc);
//
//		/**
//		 * Find any components of the first regular chain that are not contained in a component of the second regular chain.
//		 *
//		 * @param rc1: the first regular chain
//		 * @param rc2: the second regular chain
//		 **/
//		static std::vector<RegularChain<Field,RecursivePoly>> oneSideCompare(const RegularChain<Field,RecursivePoly>& rc1, const RegularChain<Field,RecursivePoly>& rc2);


		/**
		 * Compute the gcd of two input polynomials p and q modulo the saturated ideal of the current object. The result is a list of pairs (g_i,T_i) of
		 * polynomials g_i and regular chains T_i, where g_i is gcd(p,q) modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param q: a recursively viewed polynomial
		 * @param v: the common main variable of p and q
		 * @param src: the subresultant chain of p and q
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& S, bool lazardDecompose, int inputHeightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> _regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose = false, int heightBound = 0) const;
		#endif

		/**
		 * Compute a decomposition of the current object such that on each component the input polynomial is either
		 * zero or regular modulo the saturated ideal of that component. If the polynomial is zero, (0, T_i) is returned;
		 * if the polynomial is regular, (p_i, T_i), is returned, where p_i is equivalent to p modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _regularize(const RecursivePoly& p, bool lazardDecompose, int heightBound, AsyncGenerator<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>>& results) const;
		#else
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> _regularize(const RecursivePoly& p, bool lazardDecompose, int heightBound) const;
		#endif

		/**
		 * Compute the common solutions of the input polynomial and the current regular chain.
		 * More precisely, this routine computes the intersection of the varieties of the input polynomial and the
		 * current regular chain, expressed as a set of regular chains.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _intersect(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> _intersect(const RecursivePoly& p, bool lazardDecompose, int inputHeightBound) const;
		#endif

		/**
		 * Compute a triangular decomposition of the list of input polynomials.
		 *
		 * @param polys: a vector of recursively viewed polynomials
		 * @param chains: a vector of regular chains
		 * @param lazardDecompose: whether to compute a Lazard decomposition
		 **/
		#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS
		void _triangularize(const std::vector<RecursivePoly>& polys, const std::vector<RegularChain<Field,RecursivePoly>>& chains, bool lazardDecompose, int heightBound, AsyncGenerator<RegularChain<Field,RecursivePoly>>& results) const;
		#else
		std::vector<RegularChain<Field,RecursivePoly>> _triangularize(const std::vector<RecursivePoly>& polys, std::vector<RegularChain<Field,RecursivePoly>>& chains, bool lazardDecompose, int heightBound) const;
		#endif

		std::vector<RegularChain<Field,RecursivePoly>> _triangularizeByTasks(std::vector<RecursivePoly>& polys, bool lazardDecompose, int heightBound) const;

		friend void triangularizeTask<Field,RecursivePoly>(const RegularChain<Field, RecursivePoly>& rc, std::vector<RecursivePoly>& polys, bool lazardDecompose, int heightBound, TaskScheduler* tasks, std::shared_ptr<SynchronizedWriteVector<RegularChain<Field,RecursivePoly>>> results);

		friend void intersectOne<Field,RecursivePoly>(int j, const RegularChain<Field, RecursivePoly>& T, const RecursivePoly& p, int lazardDecompose, int heightBound, SynchronizedWriteVector<RegularChain<Field,RecursivePoly>>&);

		std::vector<RegularChain<Field,RecursivePoly>> intersectTrivial(const RecursivePoly& p) const;

		bool isIntersectionTrivial(const RecursivePoly& p, RecursivePoly& pReduced) const;

		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularizeTrivial(const RecursivePoly& p, const RecursivePoly& pReduced) const;

		bool isRegularizationTrivial(const RecursivePoly& p, RecursivePoly& pReduced) const;

		RecursivePoly moduloPolysWithConstantInitial(const RecursivePoly& p_in) const;

		std::vector<RecursivePoly> factorPolynomial(const RecursivePoly& p) const;

	public:

		using TriangularSet<Field,RecursivePoly>::allVariables;

		/**
		 * Default constructor: creates an empty regular chain of variable size
		 * with empty list of transcendentals.
		 *
		 * @param
		 **/
		RegularChain<Field,RecursivePoly>();

		/**
		 * Construct an empty fixed variable list regular chain in the decreasingly ordered
		 * variables given by xs with empty list of transcendentals.
		 *
		 * @param xs: The variable names
		 **/
		RegularChain<Field,RecursivePoly> (const std::vector<Symbol>& xs);

		/**
		 * Construct an empty fixed variable list regular chain in the decreasingly ordered
		 * variables given by xs and list of transcendentals given by ts.
		 *
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

		/**
		 * Construct a fixed regular chain containing the polynomials in polys.
		 * It is assumed that the polynomials in polys form a valid regular chain.
		 *
		 * @param polys: a list of recursively viewed polynomials.
		 *
		 */
		RegularChain<Field,RecursivePoly> (const std::vector<RecursivePoly> polys);

		/**
		 * Copy constructor taking a zero-dimensional regular chain as input.
		 *
		 * @param a: A zero-dimensional regular chain
		 **/
		RegularChain<Field,RecursivePoly> (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);

		/**
		 * Copy constructor.
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly> (const RegularChain<Field,RecursivePoly>& a);

		/**
		 * Copy constructor taking a triangular set as input, assuming that the triangular set is a regular chain.
		 *
		 * @param a: A triangular set
		 **/
		RegularChain<Field,RecursivePoly> (const TriangularSet<Field,RecursivePoly>& a);

		/**
		 * Move constructor taking an r-value zero-dimensional regular chain as input.
		 *
		 * @param a: An r-value reference zero-dimensional regular chain
		 **/
		RegularChain<Field,RecursivePoly> (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);

		/**
		 * Move constructor.
		 *
		 * @param a: An r-value reference regular chain
		 **/
		RegularChain<Field,RecursivePoly> (RegularChain<Field,RecursivePoly>&& a);

		/**
		 * Move constructor taking an r-value triangular set as input, assuming that the triangular set is a regular chain.
		 *
		 * @param a: An r-value reference triangular set
		 **/
		RegularChain<Field,RecursivePoly> (TriangularSet<Field,RecursivePoly>&& a);

		/**
		 * Computational constructor: creates a regular chain given all the data.
		 *
		 * @param vs: rvalue reference to variables of the regular chain
		 * @param avs: rvalue reference to algebraic variables of the regular chain
		 * @param tvs: rvalue reference to transcendental variables of the regular chain
		 * @param polys: rvalue reference to polynomials of the regular chain
		 * @param tsm: whether the regular chain is variable or fixed
		 * @param c: characteristic of the regular chain
		 **/
		RegularChain<Field,RecursivePoly> (const std::vector<Symbol>&& vs, const std::vector<Symbol>&& avs, const std::vector<Symbol>&& tvs, const std::vector<RecursivePoly>&& ts, TriangularSetMode tsm, const mpz_class& c);

		/**
		 * Construct a set of regular chains from an input triangular set by triangularizing the elements of the input.
		 *
		 * @param T: A triangular set
		 **/
		static std::vector<RegularChain<Field,RecursivePoly>> constructChains(const TriangularSet<Field,RecursivePoly>& T);

		/**
		 * Assignment operator = for a zero-dimensional regular chain.
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (const ZeroDimensionalRegularChain<Field,RecursivePoly>& a);

		/**
		 * Assignment operator =.
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (const RegularChain<Field,RecursivePoly>& a);


		RegularChain<Field, RecursivePoly>& operator= (const TriangularSet<Field,RecursivePoly>& a);

		/**
		 * Assignment operator = imposed by abstract class BPASTriangularSet.
		 *
		 * @param a: A triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (const BPASTriangularSet<Field,RecursivePoly>& a) override;

		/**
		 * Assignment operator = imposed by abstract class BPASRegularChain.
		 *
		 * @param a: A regular chain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (const BPASRegularChain<Field,RecursivePoly>& a) override;

		/**
		 * Move assignment operator = taking an r-value zero-dimensional regular chain as input.
		 *
		 * @param a: An r-value reference zero-dimensional regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (ZeroDimensionalRegularChain<Field,RecursivePoly>&& a);

		/**
		 * Move assignment operator =.
		 *
		 * @param a: An r-value reference regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator= (RegularChain<Field,RecursivePoly>&& a);

		RegularChain<Field, RecursivePoly>& operator= (TriangularSet<Field,RecursivePoly>&& a);

		/**
		 * Move assignment operator = imposed by abstract class BPASTriangularSet.
		 *
		 * @param a: An r-value reference triangular set
		 **/
		BPASTriangularSet<Field,RecursivePoly>& operator= (BPASTriangularSet<Field,RecursivePoly>&& a) override;

		/**
		 * Move assignment operator = imposed by abstract class BPASRegularChain.
		 *
		 * @param a: An r-value reference regular chain
		 **/
		BPASRegularChain<Field,RecursivePoly>& operator= (BPASRegularChain<Field,RecursivePoly>&& a) override;

		/**
		 * Add operator +:
		 * Adds a polynomial to a regular chain and returns a new regular chain, assuming that the
		 * main variable of p is above any in the current object and that init(p) is regular modulo
		 * the saturated ideal of the current object.
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		RegularChain<Field,RecursivePoly> operator+ (const RecursivePoly& p) const;

		/**
		 * Add assignment operator +=:
		 * Adds a polynomial to a regular chain, assuming that the main variable of p is above any in
		 * the current object and that init(p) is regular modulo the saturated ideal of the current object.
		 *
		 * @param p: A recursively viewed polynomial
		 **/
		RegularChain<Field,RecursivePoly>& operator+= (const RecursivePoly& p);

		/**
		 * Add operator +:
		 * Adds a regular chain to a regular chain and returns a new regular chain, assuming that the
		 * input is above the current object and the result of adding the chains is a regular chain.
		 *
		 * @param p: A regular chain
		 **/
		RegularChain<Field,RecursivePoly> operator+ (const RegularChain<Field,RecursivePoly>& T) const;

		/**
		 * Add assignment operator +=:
		 * Adds a regular chain to a regular chain, assuming that the
		 * input is above the current object and the result adding the chains is a regular chain.
		 *
		 * @param p: A regular chain
		 **/
		RegularChain<Field,RecursivePoly>& operator+= (const RegularChain<Field,RecursivePoly>& T);

		/**
		 * Identity operator ==.
		 *
		 * @param a: A regular chain
		 **/
		bool operator== (const RegularChain<Field,RecursivePoly>& a) const;

		/**
		 * Negated identity operator !=.
		 *
		 * @param a: A regular chain
		 **/
		bool operator!= (const RegularChain<Field,RecursivePoly>& a) const;

		/**
		 * Get the number of (potentially algebraic) variables in the current object.
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return TriangularSet<Field,RecursivePoly>::numberOfVariables();
		}

		/**
		 * Get the encoded options of the regular chain, a bitwise or of RegularChainOption values.
		 *
		 * @param
		 **/
		int options() const;

		/**
		 * Set the encoded options of the regular chain.
		 *
 		 * @param opts: bitwise or of RegularChainOption values
		 **/
		void setOptions(int opts);


		/**
		 * Find out if the regular chain is known to be squarefree.
		 *
		 * @param
		 **/
		bool isSquareFree() const {
			return squareFreeLevel >= algVars.size();
		}

		bool isSaturatedIdealPrime() const {
			return saturatedIdealPrimeLevel >= algVars.size();
		}

		/**
		 * Efficiently find out if the input polynomial is in the saturated ideal of the current regular chain.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		bool isInSaturatedIdealMinimal(const RecursivePoly& p) const;

		bool isInSaturatedIdealMinimal_inner(const RecursivePoly& p) const;


		/**
		 * Find out if the input polynomial is in the saturated ideal of the current regular chain.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		bool isInSaturatedIdeal(const RecursivePoly& p) const;

		/**
		 * Find out if the input polynomial is in the saturated ideal of the current regular chain.
		 * and return the reduced input polynomial.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		bool isInSaturatedIdeal(const RecursivePoly& p, RecursivePoly& redp) const;

		/**
		 * Find out if the input polynomial is in the radical saturated ideal of the current regular chain.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		bool isInRadicalSaturatedIdeal(const RecursivePoly& p) const;

		/**
		 * Find out if the input polynomial is regular modulo the saturated ideal of the current regular chain.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		bool isRegular(const RecursivePoly& p) const;

		/**
		 * Using a modular method determine if the iterated resultant of p and this regular chain is zero or not.
		 *
		 * @param p: a recursively viewed polynomial
		 * @return true iff the iterated resultant is zero.
		 */
		bool isIteratedResultantZeroModular(const RecursivePoly& p) const;


		/**
		 * Determine whether or not the quasicomponent of the first regular chain is contained in the quasicomponent of the second
		 * using a certified method that returns true if the first is contained in the second and false if not.
		 *
		 * @param rc1: the first regular chain
		 * @param rc1: the second regular chain
		 **/
		static bool compareCertifiedNoSplit(const RegularChain<Field,RecursivePoly>& rc1, const RegularChain<Field,RecursivePoly>& rc2);

		/**
		 * Determine whether or not the quasicomponent of the first regular chain is contained in the quasicomponent of the second
		 * using a heuristic method that returns true when the first is contained in the second and false when no conclusion is possible.
		 *
		 * @param rc1: the first regular chain
		 * @param rc1: the second regular chain
		 **/
		static bool compareHeuristicNoSplit(const RegularChain<Field,RecursivePoly>& rc1, const RegularChain<Field,RecursivePoly>& rc2);

		/**
		 * Remove redundancy from the input list of regular chains.
		 *
		 * @param lrc: a list of regular chains
		 **/
		static void removeRedundantChains(const std::vector<RegularChain<Field,RecursivePoly>>& lrc, std::vector<RegularChain<Field,RecursivePoly>>& results);
//		static std::vector<RegularChain<Field,RecursivePoly>> removeRedundantChains(const std::vector<RegularChain<Field,RecursivePoly>>& lrc);

		/**
		 * Get the (potentially algebriac) variable names in decreasing order.
		 *
		 * @param
		 **/
		inline std::vector<Symbol> variables() const {
			return TriangularSet<Field,RecursivePoly>::variables();
		}

		/**
		 * Select the polynomial in the current object with main variable s, if it exists.
		 *
		 * @param s: a symbol
		 **/
		inline RecursivePoly select(const Symbol& s) const {
			return TriangularSet<Field,RecursivePoly>::select(s);
		}

		/**
		 * Returns the regular chain consisting of polynomials with
		 * main variable strictly less than s.
		 *
		 * @param s: symbol of the main variable of specified element of the regular chain
		 * @param ts: The returned regular chain
		 **/
		void lower(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;

		/**
		 * Returns the regular chain consisting of polynomials with
		 * main variable strictly greater than s.
		 *
		 * @param s: symbol of the main variable of specified element of the regular chain
		 * @param ts: The returned regular chain
		 **/
		void upper(const Symbol& s, BPASTriangularSet<Field,RecursivePoly>& ts) const;

		/**
		 * Destructively converts the current object into lower(s) by changing the set of
		 * (potentially algebraic) variables to be only those below s in the variable order.
		 *
		 * @param s: symbol of the main variable of specified element of the regular chain
		 **/
		void lowerSlice(const Symbol& s);

		/**
		 * Compute a triangular decomposition of the list of input polynomials.
		 *
		 * @param F: a vector of recursively viewed polynomials
		 * @param type: a flag to determine what strategy to use (type = 0 processes tasks as they are generated, type = 1 process by level of the
		 *              tree of chains generted in the decomposition process)
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> triangularize(const std::vector<RecursivePoly>& F, bool lazardDecompose = false, int type = 0);

		/**
		 * Compute the common solutions of the input polynomial and the current regular chain.
		 * More precisely, this routine computes the intersection of the varieties of the input polynomial and the
		 * current regular chain, expressed as a set of regular chains.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> intersect(const RecursivePoly& p, bool lazardDecompose = false, int heightBound = 0) const;
//		std::vector<RegularChain<Field,RecursivePoly>> intersect(const RecursivePoly& p) const;

		/**
		 * Compute a decomposition of the current object such that on each component the input polynomial is either
		 * zero or regular modulo the saturated ideal of that component. If the polynomial is zero, (0, T_i) is returned;
		 * if the polynomial is regular, (p_i, T_i), is returned, where p_i is equivalent to p modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 **/
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularize(const RecursivePoly& p, bool lazardDecompose = false, int heightBound = 0) const;
//		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularize(const RecursivePoly& p) const;

		/**
		 * Compute the gcd of two input polynomials p and q modulo the saturated ideal of the current object. The result is a list of pairs (g_i,T_i) of
		 * polynomials g_i and regular chains T_i, where g_i is gcd(p,q) modulo the saturated ideal of T_i.
		 *
		 * @param p: a recursively viewed polynomial
		 * @param q: a recursively viewed polynomial
		 * @param v: the common main variable of p and q
		 * @param src: the subresultant chain of p and q
		 **/
		std::vector<PolyChainPair<RecursivePoly,RegularChain<Field,RecursivePoly>>> regularGCD(const RecursivePoly& p, const RecursivePoly& q, const Symbol& v, const SubResultantChain<RecursivePoly,RecursivePoly>& src, bool lazardDecompose = false, int heightBound = 0) const;

		/**
		 * A routine that decomposes the regular chain formed from the current object and an input polynomial into a set of squarefree regular chains.
		 *
		 * @param p: a recursively viewed polynomial such that init(p) is regular modulo the saturated ideal of the current object
		 * @param v: the main variable of p
		 * @param options: a bitwise or of RegularChainOption values (default assume that T+p is a regular chain, where T is the current object)
		 **/
		std::vector<RegularChain<Field,RecursivePoly>> squareFreePart(const RecursivePoly& p, const Symbol& v, bool lazardDecompose = false, int heightBound = 0, int options=ASSUME_REGULAR) const;


		/**
		 * Generate a random regular chain based on the number of terms of the polynomials in the chain.
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
		 * Generate a random regular chain based on a list of maximum degrees of variables in the polynomials in the chain.
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

		/**
		 * Reduce the polynomials of the input vector modulo the saturated ideal of the current object and detect any obvious inconsistency among the set.
		 *
		 * @param polys: a vector of recursively viewed polynomials
		 **/
		bool cleanSet(std::vector<RecursivePoly>& polys) const;
//		bool cleanSet(const std::vector<RecursivePoly>& inPolys, std::vector<RecursivePoly>& outPolys) const;

};

#endif
