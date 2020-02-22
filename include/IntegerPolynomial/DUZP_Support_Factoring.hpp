
/*****
 * Supporting methods for univariate polynomial factoring.
 *
 * NOTE: functions herein are still experimental, may not be bug-free, 
 *       and may change in future releases.
 *****/

#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#ifdef __cplusplus
/*extern "C" { */
#include <gmpxx.h>
#endif

#include <stdlib.h>
#include <gmp.h>
#include "IntegerPolynomial/DUZP_Hensel.h"
#include "IntegerPolynomial/DUZP_Support.h"
#include "ModularPolynomial/DUSP_Support.h"

//#include "NTL/ZZ.h"
//#include "NTL/vec_long.h"
//#include "NTL/vec_vec_long.h"
//#include "NTL/vec_ZZ.h"
//#include "NTL/ZZX.h"
//#include "NTL/mat_ZZ.h"
//#include "NTL/LLL.h"

//#define CARDINALITY_THRESHOLD 1
//#define SIZE_THRESHOLD 1
#define CARDINALITY_THRESHOLD 3			// NTL uses 3
#define SIZE_THRESHOLD 12				// NTL uses 12

#define DUZP_OVERLIFT 64
#define BPAS_SP_NBITS 63
#define BPAS_SPFACT_INIT_NUM_PRIMES 2	// NTL uses 7

namespace DUZP {

    namespace Factoring {


        typedef long long int long_t;

        /**
         * A vector container for duspoly_t elements.
         *
         * The container can be used in two ways:
         *
         *  - Automatic memory control:
         *
         *    after initialization with the init 
         *    method (either the default allocation 
         *    of 1 and size of 0 [init], or a 
         *    supplied allocation size with a size 
         *    of 0 [init2]), elements are pushed and 
         *    popped using corresponding methods, 
         *    with the allocated size of the 
         *    container automatically adjusted.
         *
         *  - Manual memory control:
         *
         *    after initialization with the init 
         *    method to set the allocation size, or 
         *    the init_set method to set the size to 
         *    be the allocation size, the underlying 
         *    array of duspoly_t can be managed 
         *    directly for finer control of memory; 
         *    the size of the allocation can be 
         *    controlled by the set method (which 
         *    preserves information if allocation 
         *    size increases).
         */
        typedef struct vec_duspoly {
            duspoly_t** polys;
            Prime_ptr* Pptr;
            long alloc;
            long size;
        } vec_duspoly_t;

        /**
         * A vector container for DUZP_t elements.
         *
         * The container can be used in two ways:
         *
         *  - Automatic memory control:
         *
         *    after initialization with the init 
         *    method (either the default allocation 
         *    of 1 and size of 0 [init], or a 
         *    supplied allocation size with a size 
         *    of 0 [init2]), elements are pushed and 
         *    popped using corresponding methods, 
         *    with the allocated size of the 
         *    container automatically adjusted.
         *
         *  - Manual memory control:
         *
         *    after initialization with the init 
         *    method to set the allocation size, or 
         *    the init_set method to set the size to 
         *    be the allocation size, the underlying 
         *    array of DUZP_t can be managed directly 
         *    for finer control of memory; the size 
         *    of the allocation can be controlled by 
         *    the set method (which preserves 
         *    information if allocation size 
         *    increases).
         */
        typedef struct vec_DUZP {
            DUZP_t** polys;
            long alloc;
            long size;
        } vec_DUZP_t;

        /**
         * A pair struct for use by the 
         * vec_DUZP_long_t container.
         */
        typedef struct DUZP_long {
        	DUZP_t* a;
        	long b;
        } DUZP_long_t;

        /**
         * A vector container for DUZP_t/long pair 
         * [DUZP_long_t] elements.
         *
         * The container can be used in two ways:
         *
         *  - Automatic memory control:
         *
         *    after initialization with the init 
         *    method (either the default allocation 
         *    of 1 and size of 0 [init], or a 
         *    supplied allocation size with a size 
         *    of 0 [init2]), elements are pushed and 
         *    popped using corresponding methods, 
         *    with the allocated size of the 
         *    container automatically adjusted.
         *
         *  - Manual memory control:
         *
         *    after initialization with the init 
         *    method to set the allocation size, or 
         *    the init_set method to set the size to 
         *    be the allocation size, the underlying 
         *    array of pairs can be managed directly 
         *    for finer control of memory; the size 
         *    of the allocation can be controlled by 
         *    the set method (which preserves 
         *    information if allocation size 
         *    increases).
         */
        typedef struct vec_DUZP_long {
            DUZP_long_t* pairs;
            long alloc;
            long size;
        } vec_DUZP_long_t;

        /**
         * A vector container for vectors of mpz_t of 
         * the same depth (similar to a matrix where 
         * each row is stored in a separately allocated 
         * array). The length gives the outer vector 
         * size (number of rows) and the depth gives 
         * the innner vector size (number of columns).
         *
         * The container is designed to first be 
         * initialized with a given length and depth 
         * of 1. The length and depth can then be 
         * changed with the set method. The set method 
         * is information-preserving if the length and 
         * depth are non-decreasing. The container must
         * be freed when it is no longer needed.
         */
        typedef struct vec_vec_mpz {
        	mpz_t** array;
        	long length;
        	long depth;
        } vec_vec_mpz_t;

        /**
         * A vector container for GMP mpz_t elements.
         *
         * The container can be used in two ways:
         *
         *  - Automatic memory control:
         *
         *    after initialization with the init 
         *    method (either the default allocation 
         *    of 1 and size of 0 [init], or a 
         *    supplied allocation size with a size 
         *    of 0 [init2]), elements are pushed and 
         *    popped using corresponding methods, 
         *    with the allocated size of the 
         *    container automatically adjusted.
         *
         *  - Manual memory control:
         *
         *    after initialization with the init 
         *    method to set the allocation size, or 
         *    the init_set method to set the size to 
         *    be the allocation size, the underlying 
         *    array of mpz_t can be managed directly 
         *    for finer control of memory; the size 
         *    of the allocation can be controlled by 
         *    the set method (which preserves 
         *    information if allocation size 
         *    increases).
         */
        typedef struct vec_mpz {
        	mpz_t* array;
        	long alloc;
        	long size;
        } vec_mpz_t;

        /**
         * Print a vector of mpz_t.
         **/
        void print_vec_mpz(const vec_mpz_t a);

        /**
         * Print a vector of DUZP_t.
         **/ 
        void print_vec_DUZP(const vec_DUZP_t a);

        /**
         * Print a vector of duspoly_t.
         **/ 
        void print_vec_duspoly(const vec_duspoly_t a);

        /**
         * Print a vector of DUZP, long pairs.
         **/ 
        void print_vec_DUZP_long(const vec_DUZP_long_t a);

        /**
         * print an array of mpz_t.
         **/ 
        void print_mpz_array(const mpz_t* a, long len);

        /**
         * Print an array of long_t.
         **/ 
        void print_long_t_array(const long_t* a, long len);

        /**
         * Print a vector of long (currently takes an NTL 
         * type, but will eventually be a BPAS vec_long.
         **/ 
        //void print_vec_long(const NTL::vec_long& a);

        /**
         * Print a vector of vectors of mpz_t.
         **/ 
        void print_vec_vec(const vec_vec_mpz_t* vec);

        /**
         * Initialize a vector of mpz_t (allocation 1; size 0).
         **/
        void vec_mpz_t_init(vec_mpz_t* a);

        /**
         * Initialize a vector of mpz_t (allocation alloc; size 0).
         **/
        void vec_mpz_t_init2(vec_mpz_t* a, long alloc);

        /**
         * Set the size of a vector of mpz_t; if alloc < size, then 
         * the vector will be realloc'd to allow the size to be set.
         *
         * This function is intended for use only when the user is 
         * manually maintaining the size of the vector for efficiency
         * or control purposes. For automatic size handling the push 
         * and pop functions should be used instead.
         **/
        void vec_mpz_t_set(vec_mpz_t* a, long size);

        /**
         * Push a mpz_t onto the end of a vector of mpz_t.
         *
         * This function automatically expands the vector (by doubling)
         * if the allocated space is full.
         **/
        void vec_mpz_t_push(vec_mpz_t* a, const mpz_t b);

        /**
         * Push a long int onto a vector of mpz_t.
         *
         * This function automatically expands the vector (by doubling)
         * if the allocated space is full.
         **/
        void vec_mpz_t_push_l(vec_mpz_t* a, long b);

        /**
         * Pop a mpz_t from the end of a vector of DUZP_t.
         *
         * This function automatically shrinks the vector when
         * the allocated space is 3 times or more than the 
         * size of the vector.
         **/
        void vec_mpz_t_pop(mpz_t* r, vec_mpz_t* a);

        /**
         * Free the allocated space of a vector of mpz_t.
         **/
        void vec_mpz_t_free(vec_mpz_t* a);

        /**
         * Initialize a vector of vectors of mpz_t.
         * 
         * The outer vector is set to have the supplied length, with
         * each inner vector in it being initialized to size 1.
         **/
        void vec_vec_mpz_t_init(vec_vec_mpz_t* a, long length);

        /**
         * Set the size of a vector of vectors of mpz_t based on a 
         * supplied length of the outer vector and depth of the inner
         * vector.
         *
         * This routine will automatically extend the current length 
         * and depth of the vector of vectors to acommodate the required
         * new length and depth by reallocating the underlying arrays.
         * Accordingly there is no information loss when the length and 
         * depth are non-decreasing, with information only being lost 
         * when either the length or depth decreases. Note that each 
         * inner vector in the outer vector has the same depth.
         **/
        void vec_vec_mpz_t_set(vec_vec_mpz_t* a, long length, long depth);

        /**
         * Free the allocated space of a vector of vectors of mpz_t.
         **/
        void vec_vec_mpz_t_free(vec_vec_mpz_t* a);

        /**
         * Initialize a vector of DUZP_t/long pairs 
         * (allocation of 1, size 0).
         **/
        void vec_DUZP_long_t_init(vec_DUZP_long_t* a);
         
        /**
         * Initialize a vector of DUZP_t/long pairs 
         * (allocation alloc; size 0).
         **/
        void vec_DUZP_long_t_init2(vec_DUZP_long_t* a, long alloc);

        /**
         * Initialize a vector of DUZP_t/long pairs 
         * (allocation alloc; size alloc).
         *
         * This function is intended for use only when the user is 
         * manually maintaining the size of the vector for efficiency
         * or control purposes. For automatic size handling the push 
         * and pop functions should be used instead.
         **/
        void vec_DUZP_long_t_init_set(vec_DUZP_long_t* a, long alloc);

        /**
         * Set the size of a vector of DUZP_t/long pairs. If 
         * alloc < size, then the vector will be realloc'd to allow
         * the size to be set.
         *
         * This function is intended for use only when the user is 
         * manually maintaining the size of the vector for efficiency
         * or control purposes. For automatic size handling the push 
         * and pop functions should be used instead.
         **/
        void vec_DUZP_long_t_set(vec_DUZP_long_t* a, long size);

        /**
         * Clear a vector of DUZP_t/long pairs. The size will be set 
         * to zero and the alloc to 1.
         **/ 
        void vec_DUZP_long_t_clear(vec_DUZP_long_t* a);

        /**
         * Push a DUZP_t/long pair onto the end of a vector of 
         * DUZP_t/long pairs.
         *
         * This function automatically expands the vector (by doubling)
         * if the allocated space is full.
         **/
        void vec_DUZP_long_t_push(vec_DUZP_long_t* a, const DUZP_long_t b);
         
        /**
         * Pop a DUZP_t/long pair from the end of a vector of 
         * DUZP_t/long pairs.
         *
         * This function automatically shrinks the vector when
         * the allocated space is 3 times or more than the 
         * size of the vector.
         **/
        DUZP_long_t vec_DUZP_long_t_pop(vec_DUZP_long_t* a);
         
        /**
         * Free the allocated space of a vector of DUZP_t/long pairs.
         **/
        void vec_DUZP_long_t_free(vec_DUZP_long_t* a);

        /**
         * Initialize a vector of DUZP_t (allocation of 1, size 0).
         **/
        void vec_DUZP_t_init(vec_DUZP_t* a);

        /**
         * Initialize a vector of DUZP_t (allocation alloc; size 0).
         **/
        void vec_DUZP_t_init2(vec_DUZP_t* a, long alloc);

        /**
         * Initialize a vector of DUZP_t (allocation alloc; size alloc).
         *
         * This function is intended for use only when the user is 
         * manually maintaining the size of the vector for efficiency
         * or control purposes. For automatic size handling the push 
         * and pop functions should be used instead.
         **/
        void vec_DUZP_t_init_set(vec_DUZP_t* a, long alloc);

        /**
         * Set the size of a vector of DUZP_t; if alloc < size, then 
         * the vector will be realloc'd to allow the size to be set.
         *
         * This function is intended for use only when the user is 
         * manually maintaining the size of the vector for efficiency
         * or control purposes. For automatic size handling the push 
         * and pop functions should be used instead.
         **/
        void vec_DUZP_t_set(vec_DUZP_t* a, long size);

        /**
         * Clear a vector of DUZP_t. The size will be set to 
         * zero and the alloc to 1.
         **/
        void vec_DUZP_t_clear(vec_DUZP_t* a);

        /**
         * Push a DUZP_t onto the end of a vector of DUZP_t.
         *
         * This function automatically expands the vector (by doubling)
         * if the allocated space is full.
         **/
        void vec_DUZP_t_push(vec_DUZP_t* a, const DUZP_t* b);

        /**
         * Pop a DUZP_t from the end of a vector of DUZP_t.
         *
         * This function automatically shrinks the vector when
         * the allocated space is 3 times or more than the 
         * size of the vector.
         **/
        DUZP_t* vec_DUZP_t_pop(vec_DUZP_t* a);

        /**
         * Free the allocated space of a vector of DUZP_t.
         **/
        void vec_DUZP_t_free(vec_DUZP_t* a);

        /**
         * Initialize a vector of duspoly_t (allocation of 1, size 0).
         **/
        void vec_duspoly_t_init(vec_duspoly_t* a, Prime_ptr* Pptr);
         
        /**
         * Initialize a vector of duspoly_t (allocation alloc; size 0).
         **/
        void vec_duspoly_t_init2(vec_duspoly_t* a, Prime_ptr* Pptr, long alloc);

        /**
         * Initialize a vector of duspoly_t (allocation alloc; size alloc).
         *
         * This function is intended for use only when the user is 
         * manually maintaining the size of the vector for efficiency
         * or control purposes. For automatic size handling the push 
         * and pop functions should be used instead.
         **/
        void vec_duspoly_t_init_set(vec_duspoly_t* a, Prime_ptr* Pptr, long alloc);
         
        /**
         * Set the size of a vector of duspoly_t; if alloc < size, then 
         * the vector will be realloc'd to allow the size to be set.
         *
         * This function is intended for use only when the user is 
         * manually maintaining the size of the vector for efficiency
         * or control purposes. For automatic size handling the push 
         * and pop functions should be used instead.
         **/
        void vec_duspoly_t_set(vec_duspoly_t* a, long size);

        /**
         * Push a DUZP_t onto the end of a vector of duspoly_t.
         *
         * This function automatically expands the vector (by doubling)
         * if the allocated space is full.
         **/
        void vec_duspoly_t_push(vec_duspoly_t* a, const duspoly_t* b);

        /**
         * Pop a duspoly_t from the end of a vector of duspoly_t.
         *
         * This function automatically shrinks the vector when
         * the allocated space is 3 times or more than the 
         * size of the vector.
         **/
        duspoly_t* vec_duspoly_t_pop(vec_duspoly_t* a);
         
        /**
         * Free the allocated space of a vector of duspoly_t.
         **/
        void vec_duspoly_t_free(vec_duspoly_t* a);

        /**
         * Computes a cut-off value determining when the 
         * factorization algorithm should switch from 
         * sparse to dense mode (the latter using a 
         * random linear transformation A of traces of the 
         * polynomial being factored. After this cut-off,
         * the value return determines the number of 
         * random mixtures of the traces computed 
         * (i.e., the dimension of the image of A).
         **/
        long d1_value(long deltaInBits, long r, long s);

        /**
         * Compute the excess precision needed for the current 
         * iteration of the van Hoeij algorithm given as the 
         * prime p to the power of delta, where delta is given 
         * as a number of bits. Then delta and p^delta are 
         * returned.
         **/
        void Computepdelta(long* delta, mpz_t* pdelta, long_t p, long deltaInBits);

        /**
         * Compute the vector of values b for p^b such that p^{b_i} 
         * is a bound for the ith trace of a rational factor of the 
         * polynomial f being factored. This computation requires a 
         * bound on the roots of f, with n = deg(f).
         **/
        //void Computepb(NTL::vec_long& b, vec_mpz_t* pb, long_t p, long d, const mpz_t rootBound, long n);

        /**
         * Compute an effective value b_eff for the bound p^{b_eff} on 
         * traces of a rational factor of the polynomial being factored
         * when a random linear transformation of the traces is being 
         * used. This computation requires a bound on the roots of f, 
         * with n = deg(f).
         **/
        void Computepb_eff(long* b_eff, mpz_t* pb_eff, long_t p, long d, const mpz_t rootBound, long n, long ran_bits);

        /**
         * Computes the dth trace of of the kth p-adic factor f mod P, 
         * with Tr storing all of the computed traces of all of the 
         * p-adic factors mod P.
         *
         * The routine ensures that the input f is monic, and deg(f)>0.
         * The prime power P must be > 1.
         * Tr->length must be >= d and Tr->array[i-1], for i = 1..d-1, 
         * should be the Tr_i(f) mod P (in van Hoeij's notation).
         * The quantity Tr_d(f) mod P is computed, and stored in 
         * Tr->array[d-1].
         */
        void ComputeTrace(vec_vec_mpz_t* Tr, long k, const DUZP_t* f, long d, const mpz_t P);

        /**
         * Computes the two-sided cut (see van Hoeij's paper) of d
         * traces of the kth p-adic factor mod P of the polynomial 
         * being factored. This computation requires the vector of 
         * bounds p^{b_i} on the traces, the excess precision p^delta 
         * above the largest bound p^{b_i} to which the p-adic factors 
         * have been computed, as well as the leading coefficient lc 
         * of the polynomial being factored.
         *
         * The routine requires that the number of traces d is > 0.
         * The depth of the vectors cutTraces and traces must be >= d.
         * The length of the vector pb must have length >= d
         */
        void CutTraces(vec_vec_mpz_t* cutTraces, const vec_vec_mpz_t* traces, long k, long d, const vec_mpz_t pb, const mpz_t pdelta, const mpz_t P, const mpz_t lc);

        /**
         * Computes the basis matrix for the knapsack problem defined 
         * by van Hoeij (2002). This is formed from the lattice 
         * basis matrix B_L for the lattice of (powers of) rational 
         * factors of the input polynomial, the two-sided cuts of the 
         * (possibly transformed) traces cutTraces of the p-adic 
         * factors, the excess precision pdelta that determines the 
         * precision of the two-sided cuts, the number d of traces 
         * used and the number r of p-adic factors.
         *
         * The routine computes and returns a `balancing constant' C
         * to prevent significant size differences between different 
         * entries in the basis matrix, and returns the reduction 
         * matrix M.
         */
        //void BuildReductionMatrix(NTL::mat_ZZ& M, long* C, long r, long d, const mpz_t pdelta, const vec_vec_mpz_t* cutTraces, const NTL::mat_ZZ& B_L);

        /**
         * Computes the projection of the reduced basis matrix 
         * computed by LLL onto the first r coordinates, generating 
         * an updated basis for the the lattice of (powers of) 
         * rational factors of the input polynomial, after 
         * determining the number of rows that are `short' in 
         * the sense of LLL reduction, which is determined from 
         * the vector D of sizes of the Gram-Schmidt basis vectors
         * returned by the LLL_plus algorithm.
         *
         * This computation requires the `balancing constant' C, 
         * the number r of p-adic factors, the number d of 
         * (transformed) traces computed in order to compute the 
         * size bound on all vectors in the solution set of the 
         * knapsack problem must have. This determines the number 
         * of rows in the projected basis.
         */
        //void ProjectBasis(NTL::mat_ZZ& B1, NTL::vec_ZZ& D, NTL::mat_ZZ& M, long C, long r, long d);

        /**
         * Computes an integer d together with an n x m matrix R
         * so that R/d is the reduced row echelon form of the 
         * input n x m integer matrix M, which is assumed to have 
         * linearly independent rows.
         *
         * This routine is probabilistic in the following sense: the
         * result is always correct, but with a negligible 
         * probability (specifically, if NTL::GenPrime returns a 
         * composite, and the modular NTL::gauss routine can't 
         * invert a non-zero element).
         */
        //void ReducedRowEchelonForm(NTL::ZZ& d_out, NTL::mat_ZZ& R_out, const NTL::mat_ZZ& M);

        /**
         * Computes a potential rational factor of the polynomial 
         * being factored from the p-adic factors mod P according 
         * to the index vector I.
         */
        //void ComputePotentialRationalFactor(DUZP_t** factor, const vec_DUZP_t* padicFactors, const mpz_t P, const NTL::vec_long& I);

        /**
         * Checks that the success conditions A and B for van 
         * Hoeij's algorithm are met. If both conditions are 
         * met then the vector of rational factors of the 
         * polynomial f being factored is returned.
         *
         * Condition A is that each column of B_L contains
         * exactly one 1 with all other entries being zero,
         * required if the basis vectors are to determine a
         * rational factor.
         *
         * Condition B is that each row provides the indices 
         * of the p-adic factors mod P that when multiplied 
         * together and reduced mod P produce a rational 
         * factor of f.
         */
        //long ConditionsAreMet(vec_DUZP_t* factors, const NTL::mat_ZZ& B_L, const vec_DUZP_t* padicFactors, const mpz_t P, const DUZP_t* f, long bound);

        /**
         * Computes the recombination step for univariate 
         * factorization following the algorithm of van Hoeij 
         * (2002) and the implementation strategy of the NTL 
         * library.
         *
         * The routine takes as input the polynomial f_in to 
         * be factored, the factorization of f_in modularFactors
         * computed modulo p, the vector of p-adic factors 
         * lifted modulo P_in = p^{e_in}, the lifting bound, as 
         * well as the linear (sigmas2) and quadratic (sigmas) 
         * coefficients in the multi-term diophantine problem 
         * solved in the lifting process.
         *
         * Currently linear lifting is being used, so sigmas 
         * is not being used.
         */
        void Recombine(vec_DUZP_t* rationalFactors, const DUZP_t* f_in, const vec_DUZP_t* padicFactors_in, DUZP_t** sigmas, const mpz_t P_in, long_t p, int e_in, long_t bound, const factors_t* modularFactors, duspoly_t** sigmas2);

        /**
         * Computes the factorization of the input polynomial f
         * modulo a number of small primes (small in the sense 
         * of starting from 3 and working upwards) and computes
         * the factorization of modulo the prime that has the 
         * smallest number of factors of the primes tried.
         *
         * The input f is partially factored using distinct 
         * degree factoriztion if the prime does not divide 
         * lc(f) and f is squarefree modulo that prime. Of the 
         * primes for which distinct degree factorization is 
         * performed (numbering BPAS_SPFACT_INIT_NUM_PRIMES), 
         * the prime with the smallest number of modular factors 
         * is the factored fully using same degree factorization.
         * The result of this computation is returned.
         *
         * Degree pattern computations are used for early 
         * detection of irreducible input.
         */
        //long SmallPrimeFactorize(factors_t** modularFactors, factoring_info_t* factoringInfo, const DUZP_t* f);

        /**
         * Computes the squarefree decomposition of the input 
         * polynomial ff using Yun's squarefree factorization 
         * algorithm.
         */
        void SquareFreeDecompose(vec_DUZP_long_t* u, const DUZP_t* ff);

        /**
         * Factor a primitive square free univariate polynomial ff into irreducible factors.
         * A lifting bound bnd must be supplied. Unlike factorWithBound, no
         * default bound will be computed or used.
         */
        void factor_prim_sqf(vec_DUZP_t* factors, const DUZP_t* ff, long bnd);

        /**
         * Exported univariate factorization routine with an 
         * supplied lifting bound bnd (if bnd is set to zero,
         * the Landau-Mignotte bound will be used).
         *
         * The routine returns the content c along with the 
         * vector of factors of ff over Z.
         */
        void factorWithBound(mpz_t* c, vec_DUZP_long_t* factors, const DUZP_t* ff, long bnd);

        /**
         * Exported univariate factorization routine using
         * the Landau-Mignotte bound.
         *
         * The routine returns the content c along with the 
         * vector of factors of ff over Z.
         */
        static inline void factor(mpz_t* c, vec_DUZP_long_t* factors, const DUZP_t* ff) {

        	factorWithBound(c,factors,ff,0);
        }

        /**
         * Possible degree computation: rational version.
         *
         * Computes all the possible degrees of products of 
         * p-adic factors in polys in the following way:
         *
         *   - degs[i] encodes as a bit vector the possible; 
         *   - degrees of products of size m in the set 
         *     {polys[i],...,polys[len-1]}.
         */
        void ComputePossibleDegreesRational(mpz_t* degs, const vec_DUZP_t* factors, long m);

        /**
         * Possible degree computation: rational version.
         *
         * Computes all the possible degrees of products of 
         * modular factors in polys in the following way:
         *
         *   - degs[i] encodes as a bit vector the possible; 
         *   - degrees of products of size m in the set 
         *     {polys[i],...,polys[len-1]}.
         */
        void ComputePossibleDegreesModular(mpz_t* degs, const factors_t* factors, long m, const Prime_ptr* pp);

        /**
         * Computes a potential factor g of the polynomial being
         * factored from m of the p-adic factors (polys) mod P
         *  according to the index array idxs.
         **/
        void ComputePotentialFactorRational(DUZP_t** g, DUZP_t** polys, const long* idxs, long m, const mpz_t P);

        /**
         * Remove m modular factors from a vector according 
         * to an index array.
         **/
        void RemoveFactorsModular(factors_t* modularFactors, const long* idxs, long m);

        /**
         * Convert a bit vector v to an array x of integers 
         * of size n.
         **/
        void UnpackBitVector(int* x, const mpz_t v, long n);

        /**
         * Returns the maximum number of bits needed to store the 
         * coefficients of a polynomial f.
         **/
        long MaxBits(const DUZP_t* f);

        /**
         * Computes the number of bits needed to store a long int.
         */
        static long NumBits(long n) {

            return (long) ceil(log2(abs(n)));

        }

        /**
         * Computes the number of bits needed to store an mpz_t.
         */
        static long NumBits(const mpz_t n) {
            
            return (long) mpz_sizeinbase(n,2);

        }

        /**
         * Computes the log2 of the Landau-Mignotte bound 
         *      ∥h∥∞ ≤ (n + 1)^0.5 * 2^k * ∥f∥∞
         * where n = deg(f) and k = deg(h), h a rational factor of f.
         **/
        long LandauMignotteBoundInBits(const DUZP_t* f, long k);

        /**
         * Exponential complexity search for the rational factors of 
         * the polynomial f_in from the p-adic factors of f_in mod P,
         * P = p^e. The algorithm will search for rational factors 
         * formed from m p-adic factors, based on whether the product 
         * of p-adic factors has coefficients less than a supplied 
         * bound.
         **/
        void NaiveFactorSearch(DUZP_t** f_in, vec_DUZP_t* rationalFactors, vec_DUZP_t* padicFactors, mpz_t P, long_t p, int e, long m, long_t bound);

        //void Recombine(vec_DUZP_t* rationalFactors, const DUZP_t* f_in, const vec_DUZP_t* padicFactors_in, const mpz_t P_in, long_t p, int e, long_t bound);

        /**
         * Compute a bound bnd on the size of the roots of the 
         * polynomial f.
         **/
        void computeRootBound(mpz_t* bnd, const DUZP_t* f);


    } //namespace Factoring
} //namespace DUZP
#endif
