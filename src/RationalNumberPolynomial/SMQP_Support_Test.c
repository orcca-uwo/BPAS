
#include "RationalNumberPolynomial/SMQP_Support_Test.h"
#include <time.h>

/** 
 * Get the next degree_t for a polynomial given the previous term's degree_t
 * and a "step" to take. In the univariate case this is prev+step. 
 * In the multivatirate case we consider an integer of radix maxUniDeg with 
 * coefficients described our degrees_t. We step such that the returned value
 * is equal to prev + step in this radix maxUniDeg representation.
 * e.g: prev = [1,2,7], step = 5, maxUniDeg = 10. Then next is [1,3,2];
 */
degrees_t getNextDegrees (degrees_t prev, degree_t step, degree_t maxUniDeg, int nvar) {
	degree_t weightedDegree = 0;
	degree_t weight;
	for (int i = 0; i < nvar; ++i) {
		weight = (degree_t) pow(maxUniDeg, (nvar-1 - i));
		weightedDegree += (weight*prev[i]);
	}
	weightedDegree += step;

	degrees_t nextDegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
	for (int i = 0; i < nvar; ++i) {
		weight = (degree_t) pow(maxUniDeg,(nvar-1 - i));
		if (weight != 0) {
			nextDegs[i] = weightedDegree / weight;
		} else {
			nextDegs[i] = weightedDegree;
		}
		weightedDegree -= weight * nextDegs[i];
	}
	return nextDegs;
}


/** 
 * Build a random polynomial given the the number of variables, nvar,
 * the number of terms, nterms, an (exclusive) upper bound on the absolute value
 * of the cofficients and a sparsity factor. 
 *
 * The sparsity factor is such that the difference in degree_t between sucsessive terms 
 * in the generated polynomial is 1 <= diff < sparsity;
 *
 * NOTE it is assumed nvar >= 1.
 *
 */
Node* buildRandomPoly(int nvar, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
	if (nterms <= 0) {
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = time(NULL);
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		initRand = 1;
	}

	degree_t maxTotalDeg = sparsity * nterms;
	degree_t maxUniDeg = (degree_t) ceil(pow(maxTotalDeg, 1.0 / nvar));

	Node* head;
	Node* tail;

	Node* n = (Node*) malloc(sizeof(Node));
	n->next = NULL;
	degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
	degrees_t lastDegs = degs;
	int workingIndex = nvar-1;
	if (sparsity == 2) {
		//if dense always include the constant
		degs[workingIndex] = 0;
	} else {
		degs[workingIndex] = rand() % 2; 		//50/50 chance of first term being a constant;
	}

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);

	while(mpz_sgn(mpzNum) == 0) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}
	while(mpz_sgn(mpzDen) == 0) {
		mpz_urandomb(mpzDen, R_STATE, coefBound);
	}

	// long int coef_l = rand() % coefBound;
	// if (coef_l == 0) {
	// 	coef_l = 1l;
	// }
	if (includeNeg && rand() % 2) {
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
		// coef_l *= -1;
	}
	n->degs = degs;
	mpq_init(n->coef);
	mpz_set(mpq_numref(n->coef), mpzNum);
	mpz_set(mpq_denref(n->coef), mpzDen);
	mpq_canonicalize(n->coef);
	// mpq_set_si(n->coef, coef_l, 1ul);
		// n->coef = coef_ul;
	head = n;
	tail = n;
	--nterms;

	degree_t step = 0;
	for(int i = 0; i < nterms; ++i) {
		n = (Node*) malloc(sizeof(Node));
		n->next = NULL;

		step = 0;
		while(step == 0) {
			step = rand() % sparsity;
			step = step < sparsity / 2 ? step + (sparsity / 2) : step;
		}
		degs = getNextDegrees(lastDegs, step, maxUniDeg,nvar);

	
		mpz_urandomb(mpzNum, R_STATE, coefBound);
		while(mpz_sgn(mpzNum) == 0) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
		}
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while(mpz_sgn(mpzDen) == 0) {
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}

		// coef_l = 0l; 
		// while (coef_l == 0) {
		// 	coef_l = rand() % coefBound;
		// }
		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
			// coef_l *= -1;
		}

		n->degs = degs;
		lastDegs = degs;
		mpq_init(n->coef);
		mpz_set(mpq_numref(n->coef), mpzNum);
		mpz_set(mpq_denref(n->coef), mpzDen);
		mpq_canonicalize(n->coef);
		// mpq_set_si(n->coef, coef_l, 1l);
		// n->coef = coef_ul;

		tail->next = n;
		tail = n;
	}

		//Now reverse the poly so it is in decreasing order;
	Node* cur = head;
	Node* prev = head;
	Node* next = head->next;
	cur->next = NULL;
	while(next != NULL) {
		cur = next;
		next = cur->next;
		cur->next = prev;
		prev = cur;
	}
	head = cur;

	return head;
}

/****************
* Random Polynomial Algorithm with Maximum Degree
*****************/

/**
 * @param[in] lastDegs Previous degree
 * @param[in] step The distance between the previous and returned degree
 * @param[in] nvar The number of variables
 * \example if the lastDegs = (2 2 4) and step = 3 then output is (1 2 4), and
 *  the sequence of in/out-puts is -> (0 2 4) -> (0 0 4) -> (0 0 3) -> (0 0 2).
 */
degrees_t getNextMaxDegs(degrees_t lastDegs, degree_t step, int nvar)
{
	degrees_t nextDegs = (degrees_t) calloc (nvar, sizeof(degree_t));
	degree_t weight = step;

	for (int i = 0; i < nvar ; ++i)
	{
		if (weight > 0)
		{
			if (lastDegs[i] >= weight){
				nextDegs[i] = lastDegs[i] - weight;
				weight = 0;
			} else {
				weight -= lastDegs[i];
				nextDegs[i] = 0;
			}

		} else {
			nextDegs[i] = lastDegs[i];
		}
	}
	return nextDegs;
}

/**
 * Random Polynomial Algorithm with MaxDegs
 * @param[in] nvar The number of variables
 * @param[in] maxDegs The Maximum degrees of the output polynomial
 * @param[in] sparsity
 * @param[in] coefBound
 * @param[in] includeNeg
 * \brief This algorithm generates a random polynomial where the maximum degree is maxDegs.
 */
Node* randomMaxDegsPoly(int nvar, degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg)
{

	if (nvar < 1)
		fprintf(stderr,"nvar is %d < 1", nvar);

	static int initRand = 0; // initialize seed!
	if (!initRand){
		time_t t = time(NULL);
		srand(t);
		initRand = 1;
	}

	// define head and tail:
	Node* head;
	Node* tail;

	Node* n = (Node*) malloc (sizeof(Node));
	degrees_t degs = (degrees_t) calloc (nvar, sizeof(degree_t));
	degrees_t lastDegs = maxDegs;
	long int coef_l = rand() % coefBound;
	if (coef_l == 0)
		coef_l = 1l;
	if (includeNeg && rand()%2)
		coef_l *= -1;


	// set the first term's degree maxDegs:
	mpq_init (n->coef);
	mpq_set_si(n->coef, coef_l, 1ul);
	n->degs = lastDegs;

	head = n;
	tail = n;

	degree_t step = 0;
	while (lastDegs[nvar-1] != 0){
		n = (Node*) malloc (sizeof(Node));
		n->next = NULL;

		step = 0;
		do{
			step = rand() % sparsity;
		} while(step == 0);
		degs = getNextMaxDegs(lastDegs, step, nvar);
		coef_l = 0l;
		do{
			coef_l = rand() % coefBound;
		} while (coef_l == 0);
		if (includeNeg && rand() % 2)
			coef_l *= -1;

		n->degs = degs;
		lastDegs = degs;
		mpq_init(n->coef);
		mpq_set_si(n->coef, coef_l, 1l);

		tail->next = n;
		tail = n;
	}

	return head;
}

/***************
* Random Triangular Set
***************/

/**
 * Convert the exponents of a polynomial
 * @param[in] inPoly The polynomial
 * @param[in] inNvar The number of variables of inPoly
 * @param[in] outNvar The number of variables of output polynomial.
 */
void convertExponentPoly(Node* inPoly, int inNvar, int outNvar)
{
	if (inNvar == outNvar)
		return;

	degrees_t outDegs = (degrees_t) calloc (outNvar, sizeof(degree_t));
	Node* cur = inPoly;

	while(cur != NULL)
	{
		// Ex: (w,v) ---> (x=0,y=0,z=0,w,v)
		for(int i = 0; i < inNvar; ++i)
			outDegs[outNvar - inNvar + i] = cur->degs[i];

		free(cur->degs);
		cur->degs = outDegs;
		outDegs = (degrees_t) calloc (outNvar, sizeof(degree_t));
		cur = cur->next;
	}
	return;
}

/** Build a Random Triangular Set based on maxDegs
 * @param[in] T A pointer to the triangular set
 * @param[in] outNvar The number of variables of output polynomial
 * @param[in] maxDegs The Maximum degree of polynomials
 * @param[in] sparsity
 * @param[in] coefBound
 * @param[in] includeNeg
 * \brief Output is a Lazard or general triangular set (w.r.t the lazard_flag)
 * such that if lazard_flag != 0, then the output set is Lazard triangular set.
 */
void randomTriangularSet (Node** T, int outNvar,degrees_t maxDegs, degree_t sparsity, unsigned long coefBound, int includeNeg)
{
	const int lazard_flag = 1;
	time_t sed = time(NULL);
    srand(sed);

	// generate the members
	for (int i = 0; i < outNvar; i++)
	{

		if (i != outNvar-1){
			degrees_t tmpDegs = (degrees_t) calloc (i+1, sizeof(degree_t));
			for (int j = 0; j <= i; ++j)
				tmpDegs[j] = maxDegs[outNvar - 1 - i + j];

			T[i] = randomMaxDegsPoly(i+1, tmpDegs, sparsity, coefBound, includeNeg);
			convertExponentPoly(T[i], i+1, outNvar);
		} else {
			T[i] = randomMaxDegsPoly(i+1, maxDegs, sparsity, coefBound, includeNeg);
		}

		//Lazard Triangular Set:
		if (lazard_flag != 0){
			mpq_set_si(T[i]->coef, 1l, 1ul);
			for (int j = 0; j < i; ++j){
				T[i]->degs[outNvar - 1 - j] = 0;
			}
		}
	}
	return;
}
