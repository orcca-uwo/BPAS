
#include "../../include/RationalNumberPolynomial/SUQP_Support.h"
//#include "../../include/IntegerPolynomial/primes64.h"
//#include "../../include/IntegerPolynomial/PP-primes.h"
#include "../../include/FFT/src/modpn_hfiles/inlineFuncs.h"
#include "../../include/Utils/Unix_Timer.h"
//#include"../../include/FiniteFields/SmallPrimeField_Support.h"
#include "../../include/ModularPolynomial/DUSP_Support.h"
//#include "../../include/ModularPolynomial/DUSP_Support_Test.h"

//#include "FiniteFields/prime64_constants.h"

//#include "../../include/FiniteFields/SmallPrimeField_Support.h"

//#include "emmintrin.h"
//#include "../../include/IntegerPolynomial/modpoly.h"

/*****************
 * Exponent Vector functions.
 *****************/ 

/**
 * Compare two exponent vectors for lexicographical term order.
 * a: int array representing exponent vector
 * b: int array representing exponent vector
 * v: size of the arrays
 * returns -1 if a < b, 0 if equal, 1 if a > b
 */

void printAAU(AltArrU_t *aa)
{
	if (aa == NULL || aa->size == 0)
	{
		fprintf(stderr, "%d", 0);
		return;
	}

	for (int i = 0; i < AA_SIZE(aa); ++i)
	{
		gmp_fprintf(stderr, "%Qd*x^%lld + ", aa->elems[i].coef, aa->elems[i].deg);
	}

	fprintf(stderr, "\n");
}

AltArrU_t *deepCopyPolynomial_AAU(AltArrU_t *aa)
{

	if (aa == NULL)
	{
		return NULL;
	}
	AltArrU_t *newAA = (AltArrU_t *)malloc(sizeof(AltArrU_t));
	AA_SIZE(newAA) = AA_SIZE(aa);
	AA_ALLOC(newAA) = AA_SIZE(aa);
	//	newAA->alloc = aa->alloc;

	if (newAA->alloc > 0)
	{

		newAA->elems = (AAElemU_t *)malloc(sizeof(AAElemU_t) * newAA->alloc);
		
		int size = AA_SIZE(newAA);
		AAElemU_t *elems = newAA->elems;
		AAElemU_t *oldelems = aa->elems;
		for (int i = 0; i < size; ++i)
		{
			mpq_init(elems[i].coef);

			mpq_set(elems[i].coef, oldelems[i].coef);
			elems[i].deg = oldelems[i].deg;
		}
	}
	else
	{
		newAA->elems = NULL;
	}

	return newAA;
}

AltArrU_t *sortPolynomial_AAU(AltArrU_t *aa)
{
	//TODO not insertion sort.
	AAElemU_t *elems = aa->elems;
	int size = AA_SIZE(aa);

	degree_t swapDegs;
	for (int i = 1; i < size; ++i)
	{
		for (int j = i; j > 0 && (elems[j - 1].deg) <(elems[j].deg); --j)
		{
			mpq_swap(elems[j - 1].coef, elems[j].coef);
			swapDegs = elems[j - 1].deg;
			elems[j - 1].deg = elems[j].deg;
			elems[j].deg = swapDegs;
		}
	}

	condensePolyomial_AAU(aa);

	return aa;
}
/*
 * Given a polynomial in sorted order but with possible monomial duplicates, 
 * combine like terms and condense the polynomial.
 */

void condensePolyomial_AAU(AltArrU_t *aa)
{
	if (AA_SIZE(aa) < 1)
	{
		return;
	}

	int size = AA_SIZE(aa);
	AAElemU_t *elems = aa->elems;
	int insertIdx = 0;
	int compareIdx = 1;
	while (compareIdx < size)
	{
		if (elems[insertIdx].deg == elems[compareIdx].deg)
		{
			mpq_add(elems[insertIdx].coef, elems[insertIdx].coef, elems[compareIdx].coef);
		}
		else if (compareIdx - insertIdx > 1)
		{
			++insertIdx;
			elems[insertIdx].deg = elems[compareIdx].deg;
			mpq_swap(elems[insertIdx].coef, elems[compareIdx].coef);
		}
		else
		{
			++insertIdx;
		}
		++compareIdx;
	}

	++insertIdx;
	for (int i = insertIdx; i < size; ++i)
	{
		mpq_clear(elems[i].coef);
	}
	AA_SIZE(aa) = insertIdx;
}

void addTermSafe_AAU(AltArrU_t *aa, degree_t d, const ratNum_t coef)
{
	if (AA_SIZE(aa) >= aa->alloc)
	{
		aa->alloc *= 2;
		aa->elems = (AAElemU_t *)realloc(aa->elems, sizeof(AAElemU_t) * aa->alloc);
	}
	mpq_init(aa->elems[AA_SIZE(aa)].coef);
	mpq_set(aa->elems[AA_SIZE(aa)].coef, coef);
	aa->elems[AA_SIZE(aa)].deg = d;
	++(AA_SIZE(aa));
}

AAElemU_t *multiplyTerms_AAU(AAElemU_t *a, AAElemU_t *b)
{
	AAElemU_t *ret = (AAElemU_t *)malloc(sizeof(AAElemU_t));
	(ret->deg) = (b->deg) + (a->deg);
	mpq_init(ret->coef);
	mpq_mul(ret->coef, a->coef, b->coef);
	return ret;
}

/**
 * Negate a polynomial in place.
 */

void negatePolynomial_AAU(AltArrU_t *aa)
{
	int size = AA_SIZE(aa);
	AAElemU_t *elems = aa->elems;
	for (int i = 0; i < size; ++i)
	{
		mpq_neg(elems[i].coef, elems[i].coef);
	}
}

AltArrU_t *addPolynomials_AAU(AltArrU_t *a, AltArrU_t *b)
{
	AltArrU_t *c = makePolynomial_AAU(AA_SIZE(a) + AA_SIZE(b));

	AAElemU_t *aElems = a->elems;
	AAElemU_t *bElems = b->elems;
	AAElemU_t *cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;
	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	while (i < asize && j < bsize)
	{
		//t * b[i]
		//c* a[i]
		if (aElems[i].deg < bElems[j].deg)
		{
			//a < b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, bElems[j].coef);
			cElems[k].deg = bElems[j].deg;
			++k;
			++j;
		}
		else if (aElems[i].deg == bElems[j].deg)
		{
			// a==b
			mpq_init(cElems[k].coef);
			mpq_add(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].deg = aElems[i].deg;
			if (mpq_sgn(cElems[k].coef) == 0)
			{
				mpq_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		}
		else
		{
			//a > b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, aElems[i].coef);
			cElems[k].deg = aElems[i].deg;
			++k;
			++i;
		}
	}

	while (i < asize)
	{
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, aElems[i].coef);
		cElems[k].deg = aElems[i].deg;
		++k;
		++i;
	}
	while (j < bsize)
	{
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, bElems[j].coef);
		cElems[k].deg = bElems[j].deg;
		++k;
		++j;
	}

	AA_SIZE(c) = k;
	return c;
}

AltArrU_t *subPolynomials_AAU(AltArrU_t *a, AltArrU_t *b)
{
	AltArrU_t *c = makePolynomial_AAU(AA_SIZE(a) + AA_SIZE(b));

	AAElemU_t *aElems = a->elems;
	AAElemU_t *bElems = b->elems;
	AAElemU_t *cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;
	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	while (i < asize && j < bsize)
	{
		if (aElems[i].deg < bElems[j].deg)
		{
			//a < b
			mpq_init(cElems[k].coef);
			mpq_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].deg = bElems[j].deg;
			++k;
			++j;
		}
		else if (aElems[i].deg == bElems[j].deg)
		{
			// a==b
			mpq_init(cElems[k].coef);
			mpq_sub(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].deg = aElems[i].deg;
			if (mpq_sgn(cElems[k].coef) == 0)
			{
				mpq_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		}
		else
		{
			//a > b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, aElems[i].coef);
			cElems[k].deg = aElems[i].deg;
			++k;
			++i;
		}
	}

	while (i < asize)
	{
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, aElems[i].coef);
		cElems[k].deg = aElems[i].deg;
		++k;
		++i;
	}
	while (j < bsize)
	{
		mpq_init(cElems[k].coef);
		mpq_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].deg = bElems[j].deg;
		++k;
		++j;
	}

	AA_SIZE(c) = k;
	return c;
}

AltArrU_t *addPolynomials_AAU_inp(AltArrU_t *a, AltArrU_t *b)
{

	AltArrU_t *c = makePolynomial_AAU(AA_SIZE(a) + AA_SIZE(b));

	AAElemU_t *aElems = a->elems;
	AAElemU_t *bElems = b->elems;
	AAElemU_t *cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;
	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	while (i < asize && j < bsize)
	{
		if (aElems[i].deg < bElems[j].deg)
		{
			//a < b
			mpq_init(cElems[k].coef);
			mpq_set(cElems[k].coef, bElems[j].coef);
			cElems[k].deg = bElems[j].deg;
			++k;
			++j;
		}
		else if (aElems[i].deg == bElems[j].deg)
		{
			// a==b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			mpq_add(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpq_sgn(cElems[k].coef) == 0)
			{
				mpq_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;

			++i;
			++j;
		}
		else
		{
			//a > b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			// mpq_set(cElems[k].coef, aElems[i].coef);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if (i < asize)
	{
		memcpy(cElems + k, aElems + i, sizeof(AAElemU_t) * (asize - i));
		// mpq_init(cElems[k].coef);
		// mpq_set(cElems[k].coef, aElems[i].coef);
		// cElems[k].degs = aElems[i].degs;
		k += (asize - i);
		// ++k;
		// ++i;
	}
	while (j < bsize)
	{
		mpq_init(cElems[k].coef);
		mpq_set(cElems[k].coef, bElems[j].coef);
		cElems[k].deg = bElems[j].deg;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;
	return c;
}

AltArrU_t *subPolynomials_AAU_inp(AltArrU_t *a, AltArrU_t *b)
{

	AltArrU_t *c = makePolynomial_AAU(AA_SIZE(a) + AA_SIZE(b));

	AAElemU_t *aElems = a->elems;
	AAElemU_t *bElems = b->elems;
	AAElemU_t *cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;
	register int asize = AA_SIZE(a);
	register int bsize = AA_SIZE(b);

	while (i < asize && j < bsize)
	{
		if (aElems[i].deg < bElems[j].deg)
		{
			//a < b
			mpq_init(cElems[k].coef);
			mpq_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].deg = bElems[j].deg;
			++k;
			++j;
		}
		else if (aElems[i].deg == bElems[j].deg)
		{
			// a==b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			mpq_sub(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpq_sgn(cElems[k].coef) == 0)
			{
				mpq_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		}
		else
		{
			//a > b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			// mpq_set(cElems[k].coef, aElems[i].coef);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if (i < asize)
	{
		memcpy(cElems + k, aElems + i, sizeof(AAElemU_t) * (asize - i));
		// mpq_init(cElems[k].coef);
		// mpq_set(cElems[k].coef, aElems[i].coef);
		// cElems[k].degs = aElems[i].degs;
		k += (asize - i);
		// ++k;
		// ++i;
	}
	while (j < bsize)
	{
		mpq_init(cElems[k].coef);
		mpq_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].deg = bElems[j].deg;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;
	return c;
}

/////////////////////////////////////////////////////

void prodheapPrint_AAU(ProductHeap_AAU *h)
{
	polysize_t s = h->heapSize;
	fprintf(stderr, "size: %ld ", s);
	for (int i = 0; i < h->heapSize; ++i)
	{
		fprintf(stderr, "( %d,", h->elements[i].deg);
		// printDegs(h->elements[i]->product->degs, nvar);
		ProductHeapChain_AAU *next = h->elements[i].chain->next;
		while (next != NULL)
		{
			fprintf(stderr, "; ");
			fprintf(stderr, "%d,", h->elements[i].deg);
			// printDegs(next->product->degs, nvar);
			next = next->next;
		}
		fprintf(stderr, ")");
		fprintf(stderr, " ");
	}
	fprintf(stderr, "\n");
}

ProductHeap_AAU *prodheapInit_AAU(AltArrU_t *a, AltArrU_t *b)
{
	ProductHeap_AAU *h = prodheapCreate_AAU();
	h->elements = (ProductHeapElem_AAU *)malloc(sizeof(ProductHeapElem_AAU) * AA_SIZE(a));
	h->elements[0].chain = prodheapMakeChain_AAU(0, 0, NULL);

	(h->elements->deg) = (a->elems->deg) + (b->elems->deg);
	h->heapSize = 1;
	h->maxHeapSize = AA_SIZE(a);

	return h;
}

///////////////////////////////////// we start here
/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
void prodheapInsert_AAU(ProductHeap_AAU *h, ProductHeapChain_AAU *chain, register degree_t deg)
{

	register int s = h->heapSize;
	ProductHeapElem_AAU *elems = h->elements;

	if (s == 0)
	{
		elems[0].deg = deg;
		elems[0].chain = chain;
		h->heapSize = 1;
		return;
	}

	//first check if we can chain off the root
	if (elems[0].deg == deg)
	{
		chain->next = elems[0].chain;
		elems[0].chain = chain;
		return;
	}

	//otherwise, we must search the heap to find the new product's insertion point
	//note that since we are looking for chains we cannot use the simple swim method
	//we sort of fake the swimming, looking for a chain to be made or the eventual
	//place where the swim would stop. At this point, we insert the new elem
	//in that spot, and "push" the entire path we took down a level. Assuming
	//that we insert e and it ends up at the root, we push down the 'x' path
	//                                      //
	//     x     --->    e                  //
	//    / \           / \                 //
	//   x   o         x   o
	//                /
	//               x

	register int i = (s - 1) >> 1; //i is parent
	register int j = s;			   //j is current insertion point
	register long long unsigned int path = 1;
	while (j > 0)
	{
		if (elems[i].deg == deg)
		{
			chain->next = elems[i].chain;
			elems[i].chain = chain;
			return;
		}
		else if (elems[i].deg < deg)
		{
			path <<= 1;
			if (!(j & 1))
			{
				//set the trailing bit to 1 to distinguish left/right of path
				path += 1;
			}
			j = i;
			i = (i - 1) >> 1;
		}
		else
		{ //cmp > 0
			break;
		}
	}

	//then j is now the place we need to insert elem;
	//do so, and then push all others down the path, inserting the last
	//as the new element in elems[s];
	ProductHeapElem_AAU temp;
	ProductHeapElem_AAU elem = {deg, chain};

	//TODO use i index again to swap between elements rather than use a second temp elem.
	while (j <= s)
	{
		temp = elems[j];
		elems[j] = elem;
		elem = temp;
		j = (j << 1) + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
}

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * returns NULL if no such element in the heap.
 */
ProductHeapChain_AAU *prodheapRemoveMax_AAU(ProductHeap_AAU *h)
{
	ProductHeapElem_AAU *elems = h->elements;
	ProductHeapChain_AAU *maxElem = elems[0].chain;
	register int i = 0;
	register int j = 1;
	register int s = --(h->heapSize);

	//promote largest children
	while (j < s)
	{
		if (j + 1 < s && (elems[j].deg < elems[j + 1].deg))
		{
			++j;
		}
		elems[i] = elems[j];
		i = j;
		j = (j << 1) + 1;
	}
	//now place last element into i and swim up to make tree complete
	j = (i - 1) >> 1;
	while (i > 0)
	{
		if (elems[s].deg < elems[j].deg)
		{
			break;
		}
		elems[i] = elems[j];
		i = j;
		j = (j - 1) >> 1;
	}
	elems[i] = elems[s];

	return maxElem;
}
/**
 * Multiply two polynomials given their head Nodes, a and b.
 * This algorithm makes use of heaps as an efficient search data structure. 
 * It is assumed that both a and b have compatible exponent.
 * 
 * 
 *
 * returns a pointer to the head Node of the product polynomial.
 */
AltArrU_t *multiplyPolynomials_AAU(AltArrU_t *__restrict__ a, AltArrU_t *__restrict__ b)
{
	if (a == NULL || b == NULL)
	{
		return NULL;
	}

	// reorder to obtain smaller as a.
	if (b->size < a->size)
	{
		AltArrU_t *temp = a;
		a = b;
		b = temp;
	}

	ProductHeap_AAU *h = prodheapInit_AAU(a, b);
	ratNum_t ccoef;
	mpq_init(ccoef);

	//TODO smarter allocation here? dynamic reallocating?
	AltArrU_t *c = makePolynomial_AAU(AA_SIZE(a) * AA_SIZE(b));

	//k is c, i is a, j is b.
	register int k = 0;
	// register int i = 0;
	// register int j = 0;

	AAElemU_t *__restrict__ cElems = c->elems;
	AAElemU_t *__restrict__ bElems = b->elems;

	AAElemU_t *aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1;
	register int firstB = 0;

	ProductHeapChain_AAU *maxElem = NULL;
	ProductHeapChain_AAU *nextMaxElem = NULL;
	degree_t *nextDegs;
	while ((nextDegs = prodheapPeek_AAU(h)) != NULL)
	{
		//cache since, on RemoveMax, pointer is invalidated.
		cElems[k].deg = *nextDegs;
		mpq_init(cElems[k].coef);

		while (nextDegs != NULL && cElems[k].deg == *nextDegs)
		{
			//we will extract and accumulate the coefficents
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AAU *oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AAU(h);
			while (maxElem != NULL)
			{

				mpq_mul(ccoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);

				mpq_add(cElems[k].coef, cElems[k].coef, ccoef);

				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == firstB && maxElem->a_i != lastA)
				{
					oldMaxElem = prodheapMakeChain_AAU((maxElem->a_i) + 1, firstB, oldMaxElem);
				}

				//cache next before freeing or overwriting
				nextMaxElem = maxElem->next;

				//If the extracted term has another product in the stream,
				//update the product and push onto the oldMaxElem chain
				if (maxElem->b != lastB)
				{
					++(maxElem->b);
					maxElem->next = oldMaxElem;
					oldMaxElem = maxElem;
				}
				else
				{
					//we are done with the maxElem ProductHeapChain
					maxElem->next = NULL;
					prodheapFreeChain_AAU(maxElem);
				}

				maxElem = nextMaxElem;
			}

			//reset head of maxElem list
			maxElem = oldMaxElem;

			nextDegs = prodheapPeek_AAU(h);
		}

		//Commit new term to the product.
		if (mpq_sgn(cElems[k].coef) != 0)
		{
			++k;
		}
		else
		{
			//reset accumulator variables and do not increment k.
			//will init cElem[k] again on next loop, so clear here.
			mpq_clear(cElems[k].coef);
		}

		//Insert all successors of previously extracted products
		while (maxElem != NULL)
		{
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;

			prodheapInsert_AAU(h, maxElem, aElems[maxElem->a_i].deg + bElems[maxElem->b].deg);

			maxElem = nextMaxElem;
		}
	}

	mpq_clear(ccoef);
	prodheapFree_AAU(h);

	AA_SIZE(c) = k;

	return c;
}
///////////////////////////////////////////////
AltArrU_t *multiplyPolynomials_AAU_inp(AltArrU_t *a, AltArrU_t *b)
{
	AltArrU_t *prod = multiplyPolynomials_AAU(a, b);

	freePolynomial_AAU(a);
	return prod;
}

/*****************
 * Polynomial division
 *****************/

/**
 * Helper to determine if monomial b divides monomial a.
 *
 * Note: since we assume working in Q, do not need to check coefficients 
 * nvar: number of variables of monomials a and b
 */
// int monomialDivideTest(degrees_t adegs, degrees_t bdegs, int nvar) {
// 	return ( (adegs & FIRST_EXP) < (bdegs & SECOND_EXP)) ||
// 		   ( (adegs & SECOND_EXP) < (bdegs & SECOND_EXP)) ||
// 		   ( (adegs & THIRD_EXP) < (bdegs & THIRD_EXP));
// 	// for(int i = 0; i < nvar; ++i) {
// 	// 	if (adegs[i] < bdegs[i]) {
// 	// 		return 0;
// 	// 	}
// 	// }
// 	// return 1;
// }

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */

void divisionGetNextTerm_AAU(ProductHeap_AAU *h, const AAElemU_t *__restrict__ aElems, const AAElemU_t *__restrict__ bElems, mpq_t *retCoef)
{

	if (h->heapSize == 0)
	{
		return;
	}

	int lastB = h->lastB;

	ProductHeapChain_AAU *insertChain = NULL;
	ProductHeapChain_AAU *maxElem, *nextMaxElem;

	mpq_t prodCoef;
	mpq_init(prodCoef);
	degree_t *nextDegs = prodheapPeek_AAU(h);
	register degree_t maxDegs = *nextDegs;

	while (nextDegs != NULL && maxDegs == *nextDegs)
	{
		maxElem = prodheapRemoveMax_AAU(h);

		while (maxElem != NULL)
		{
			nextMaxElem = maxElem->next;
			mpq_mul(prodCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);

			mpq_add(*retCoef, *retCoef, prodCoef);
			if (maxElem->b != lastB)
			{
				++(maxElem->b);
				maxElem->next = insertChain;
				insertChain = maxElem;
			}
			else
			{
				maxElem->next = NULL;
				prodheapFreeChain_AAU(maxElem);
			}

			maxElem = nextMaxElem;
		}

		nextDegs = prodheapPeek_AAU(h);
	}

	while (insertChain != NULL)
	{
		maxElem = insertChain->next;
		insertChain->next = NULL;
		prodheapInsert_AAU(h, insertChain, aElems[insertChain->a_i].deg + bElems[insertChain->b].deg);

		insertChain = maxElem;
	}

	mpq_clear(prodCoef);
}

void divideBySingleTerm_AAU(AltArrU_t *c, AltArrU_t *b, AltArrU_t **res_a, AltArrU_t **res_r)
{
	if (b == NULL)
	{
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL)
	{
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	AAElemU_t *k = c->elems;
	AAElemU_t *lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0;

	AltArrU_t *a = makePolynomial_AAU(maxSize);
	AltArrU_t *r = makePolynomial_AAU(maxSize);
	AAElemU_t *bElems = b->elems;
	AAElemU_t *curA = a->elems;
	AAElemU_t *curR = r->elems;
	mpq_init(curA->coef);
	mpq_init(curR->coef);

	while (k != lenK)
	{
		if ((k->deg) >= (bElems->deg))
		{
			mpq_div(curA->coef, k->coef, bElems->coef);
			curA->deg = (k->deg) - (bElems->deg);
			++i;
			++(curA);
			mpq_init(curA->coef);
		}
		else
		{
			mpq_set(curR->coef, k->coef);
			curR->deg = k->deg;
			++j;
			++(curR);
			mpq_init(curR->coef);
		}
		++k;
	}

	mpq_clear(curA->coef);
	mpq_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;

	*res_a = a;
	*res_r = r;
}
/** 
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r 
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */

void dividePolynomials_AAU(AltArrU_t *c, AltArrU_t *b, AltArrU_t **res_a, AltArrU_t **res_r)
{

	if (b == NULL || AA_SIZE(b) == 0)
	{
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0)
	{
		//c is zero
		*res_a = makePolynomial_AAU(0);
		*res_r = NULL;
		fprintf(stderr, "c is zero!\n");
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1)
	{
		divideBySingleTerm_AAU(c, b, res_a, res_r);
		return;
	}

	AAElemU_t *__restrict__ k = c->elems;
	AAElemU_t *__restrict__ lenK = k + AA_SIZE(c);
	AAElemU_t *__restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0;

	AltArrU_t *a = makePolynomial_AAU(maxASize);
	AltArrU_t *r = makePolynomial_AAU(maxRSize);
	AAElemU_t *__restrict__ curA = a->elems;
	AAElemU_t *__restrict__ curR = r->elems;
	mpq_init(curA->coef);
	mpq_init(curR->coef);

	//init a with lt(c)/lt(b);
	register degree_t beta = b->elems->deg;
	while (k != lenK && !(k->deg >= beta))
	{
		mpq_set(curR->coef, k->coef);
		curR->deg = k->deg;
		++j;
		if (j >= maxRSize)
		{
			maxRSize <<= 1;
			r->elems = (AAElemU_t *)realloc(r->elems, maxRSize * sizeof(AAElemU_t));
			curR = r->elems + j - 1;
		}
		++(curR);
		mpq_init(curR->coef);
		++k;
	}

	if (k == lenK)
	{
		//no division to do at all!
		mpq_clear(curA->coef);
		mpq_clear(curR->coef);

		AA_SIZE(a) = i;
		AA_SIZE(r) = j;
		a->alloc = maxASize;
		r->alloc = maxRSize;

		*res_a = a;
		*res_r = r;
		return;
	}

	(curA->deg) = (k->deg) - beta;
	mpq_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	/////////////////////////////////////////////////////////
	ProductHeap_AAU *h = prodheapCreate_AAU();
	//////////////////////////////////
	prodheapResize_AAU(h, AA_SIZE(c));
	h->lastB = AA_SIZE(b) - 1;

	prodheapInsert_AAU(h, prodheapMakeChain_AAU(0, 1, NULL), curA->deg + b2Elem->deg);

	++i;
	++curA;
	mpq_init(curA->coef);

	degree_t *delta = prodheapPeek_AAU(h);
	register degree_t eps;
	register cmpExp_t cmp;
	while (k != lenK || delta != NULL)
	{

		if (k == lenK)
		{
			if (delta == NULL)
			{
				break;
			}
			cmp = 1;
		}
		else if (delta == NULL)
		{
			cmp = -1;
		}
		else
		{
			cmp = compareExponentVectors(*delta, k->deg);
		}

		if (cmp > 0)
		{
			eps = *delta;
			divisionGetNextTerm_AAU(h, a->elems, b->elems, &(curA->coef));

			if (mpq_sgn(curA->coef) == 0)
			{
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAU(h);
				continue;
			}
			else
			{
				mpq_neg(curA->coef, curA->coef);
			}
		}
		else if (cmp == 0)
		{
			eps = *delta;
			divisionGetNextTerm_AAU(h, a->elems, b->elems, &(curA->coef));

			if (mpq_sgn(curA->coef) == 0)
			{
				delta = prodheapPeek_AAU(h);
				continue; //the chains cancelled themselves out since the peek
			}
			else
			{
				mpq_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpq_sgn(curA->coef) == 0)
				{
					delta = prodheapPeek_AAU(h);
					continue;
				}
			}
		}
		else
		{
			eps = k->deg;
			mpq_set(curA->coef, k->coef);
			++k;
		}

		/*gmp_printf("curA->coef: %Qd\n", curA->coef);
   printf("eps: %d, beta: %d\n", eps, beta);*/

		if (eps >= beta)
		{
			(curA->deg) = eps - beta;
			mpq_div(curA->coef, curA->coef, b->elems->coef);
			if (i + 1 >= maxASize)
			{
				maxASize <<= 1;
				a->elems = (AAElemU_t *)realloc(a->elems, maxASize * sizeof(AAElemU_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAU(h, maxASize);
			}
			prodheapInsert_AAU(h, prodheapMakeChain_AAU(i, 1, NULL), curA->deg + b2Elem->deg);

			++i;
			++(curA);
			mpq_init(curA->coef);
		}
		else
		{

			//swap here so that curA becomes 0.
			mpq_swap(curR->coef, curA->coef);
			curR->deg = eps;
			++j;
			if (j >= maxRSize)
			{
				maxRSize <<= 1;
				r->elems = (AAElemU_t *)realloc(r->elems, maxRSize * sizeof(AAElemU_t));
				curR = r->elems + j - 1;
			}
			++(curR);
			mpq_init(curR->coef);
		}

		delta = prodheapPeek_AAU(h);
	}

	//clear since we always setup one past where we actually are.
	mpq_clear(curA->coef);
	mpq_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	if (r == NULL || !AA_SIZE(r))
	{
		//r->elems=NULL;
		/* r->size=1;
		mpq_t tempco;
		mpq_init(tempco);
		mpq_set_d(r->elems[0].coef,1);
		
		r->elems[0].deg=0; */
	}
	*res_a = a;
	*res_r = r;

	return;
}
int divideTestSingleTerm_AAU(AltArrU_t *c, AltArrU_t *b, AltArrU_t **res_a)
{
	if (b == NULL)
	{
		return 0;
	}

	if (c == NULL)
	{
		//c is zero
		*res_a = NULL;
		return 1;
	}

	AAElemU_t *k = c->elems;
	AAElemU_t *lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;

	AltArrU_t *a = makePolynomial_AAU(maxSize);
	AAElemU_t *bElems = b->elems;
	AAElemU_t *curA = a->elems;
	mpq_init(curA->coef);

	while (k != lenK)
	{
		if (monomialDivideTestU(k->deg, bElems->deg))
		{
			mpq_div(curA->coef, k->coef, bElems->coef);
			curA->deg = (k->deg) - (bElems->deg);
			++i;
			++(curA);
			mpq_init(curA->coef);
		}
		else
		{
			a->size = i;
			freePolynomial_AAU(a);
			return 0;
		}
		++k;
	}

	mpq_clear(curA->coef);

	AA_SIZE(a) = i;
	*res_a = a;

	return 1;
}
/**
 * Determine if a polynomial, c, is exactly divisble by another, b.
 * If the division is exact, returns the quoteint in res_a. 
 * Otherwise, the value of res_a is undefined.
 * returns 1 iff division is exact. 0 otherwise.
 */
int divideTest_AAU(AltArrU_t *c, AltArrU_t *b, AltArrU_t **res_a)
{

	if (b == NULL || AA_SIZE(b) == 0)
	{
		return 0;
	}

	if (c == NULL || AA_SIZE(c) == 0)
	{
		*res_a = makePolynomial_AAU(0);
		return 1;
	}
	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1)
	{
		return divideTestSingleTerm_AAU(c, b, res_a);
	}

	AAElemU_t *__restrict__ k = c->elems;
	AAElemU_t *__restrict__ lenK = k + AA_SIZE(c);
	AAElemU_t *__restrict__ b2Elem = b->elems + 1;
	register degree_t beta = b->elems->deg;

	if (!monomialDivideTestU(k->deg, beta))
	{
		return 0;
	}

	register int maxASize = AA_SIZE(c);
	register int i = 0;
	AltArrU_t *a = makePolynomial_AAU(maxASize);
	AAElemU_t *__restrict__ curA = a->elems;
	mpq_init(curA->coef);

	//init a with lt(c)/lt(b);
	curA->deg = k->deg - beta;
	mpq_div(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAU *h = prodheapCreate_AAU();
	///////////////////////////////////////////////////////////////////////////

	prodheapResize_AAU(h, AA_SIZE(c));
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAU(h, prodheapMakeChain_AAU(0, 1, NULL), curA->deg + b2Elem->deg);

	++i;
	++curA;
	mpq_init(curA->coef);

	degree_t *delta = prodheapPeek_AAU(h);
	register degree_t eps;
	register cmpExp_t cmp;
	while (k != lenK || delta != NULL)
	{

		if (k == lenK)
		{
			if (delta == NULL)
			{
				break;
			}
			cmp = 1;
		}
		else if (delta == NULL)
		{
			cmp = -1;
		}
		else
		{
			cmp = compareExponentVectors(*delta, k->deg);
		}

		if (cmp > 0)
		{
			eps = *delta;
			divisionGetNextTerm_AAU(h, a->elems, b->elems, &(curA->coef));

			if (mpq_sgn(curA->coef) == 0)
			{
				//in this case, the term with degree delta ended up
				//having its coffeicient cancelled out (i.e. 0)
				continue;
			}
			else
			{
				mpq_neg(curA->coef, curA->coef);
			}
		}
		else if (cmp == 0)
		{
			eps = *delta;
			divisionGetNextTerm_AAU(h, a->elems, b->elems, &(curA->coef));

			if (mpq_sgn(curA->coef) == 0)
			{
				continue; //the chains cancelled themselves out since the peek
			}
			else
			{
				mpq_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpq_sgn(curA->coef) == 0)
				{
					continue;
				}
			}
		}
		else
		{
			eps = k->deg;
			mpq_set(curA->coef, k->coef);
			++k;
		}

		if (monomialDivideTestU(eps, beta))
		{
			curA->deg = eps - beta;
			mpq_div(curA->coef, curA->coef, b->elems->coef);
			if (i + 1 >= maxASize)
			{
				maxASize <<= 1;
				a->elems = (AAElemU_t *)realloc(a->elems, maxASize * sizeof(AAElemU_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a
				prodheapResize_AAU(h, maxASize);
			}
			prodheapInsert_AAU(h, prodheapMakeChain_AAU(i, 1, NULL), curA->deg + b2Elem->deg);

			++i;
			++(curA);
			mpq_init(curA->coef);
		}
		else
		{

			//divide test fails
			prodheapFree_AAU(h);
			a->size = i;
			freePolynomial_AAU(a);
			return 0;
		}

		delta = prodheapPeek_AAU(h);
	}

	//clear since we always setup one past where we actually are.
	mpq_clear(curA->coef);

	AA_SIZE(a) = i;
	*res_a = a;
	return 1;
}
/**
 * Determine Integrattion of polynomial a.  
 *
*/

//  AltArrU_t * Integration(AltArrU_t *a){

//  	if (AltArrU_t->size )

// AAElemU_t* aElems = a->elems;

//}
AltArrU_t *exponentiatePoly_AAU(AltArrU_t *a, unsigned int n)
{
	if (n == 0)
	{
		AltArrU_t *ret = makePolynomial_AAU(1);
		mpq_init(ret->elems->coef);
		mpq_set_ui(ret->elems->coef, 1ul, 1ul);
		ret->elems->deg = 0;
		ret->size = 1;
		return ret;
	}
	else if (n == 1)
	{
		return deepCopyPolynomial_AAU(a);
	}

	AltArrU_t *r = NULL;
	AltArrU_t *b = deepCopyPolynomial_AAU(a);
	while (n > 1)
	{
		if (n & 1)
		{
			r = (r == NULL) ? deepCopyPolynomial_AAU(b) : multiplyPolynomials_AAU_inp(r, b);
		}
		b = multiplyPolynomials_AAU_inp(b, b);
		n >>= 1;
	}
	r = (r == NULL) ? deepCopyPolynomial_AAU(b) : multiplyPolynomials_AAU_inp(r, b);

	freePolynomial_AAU(b);
	return r;
}

/*****************
 * Content, PrimitivePart, etc.
 *****************/

void integralContent_AAU(AltArrU_t *aa, mpq_t ret)
{
	if (aa == NULL || aa->size <= 0)
	{
		mpq_set_ui(ret, 1ul, 1ul);
		return;
	}

	AAElemU_t *elems = aa->elems;
	int size = aa->size;
	mpq_set_ui(ret, 1ul, 1ul);

	mpq_t one;
	mpq_init(one);
	mpq_set_si(one, 1l, 1l);

	mpq_t cont;
	mpq_init(cont);
	mpq_abs(cont, elems->coef);

	mpq_t tempa;
	mpq_t tempb;
	mpq_init(tempa);
	mpq_init(tempb);

	for (int i = 1; i < size; ++i)
	{

		mpq_set(tempa, cont);
		mpq_set(tempb, elems[i].coef);

		mpz_mul(mpq_numref(tempa), mpq_numref(tempa), mpq_denref(tempb));
		mpz_mul(mpq_numref(tempb), mpq_numref(tempb), mpq_denref(tempa));
		mpz_mul(mpq_denref(tempa), mpq_denref(tempa), mpq_denref(tempb));

		mpz_gcd(mpq_numref(cont), mpq_numref(tempa), mpq_numref(tempb));
		mpz_set(mpq_denref(cont), mpq_denref(tempa));
		mpq_canonicalize(cont);
	}

	if (mpq_sgn(elems->coef) < 0 && mpq_sgn(cont) > 0)
	{
		mpq_neg(cont, cont);
	}

	mpq_set(ret, cont);
	mpq_clear(cont);
	mpq_clear(one);
	mpq_clear(tempa);
	mpq_clear(tempb);
}

AltArrU_t *primitivePart_AAU(AltArrU_t *aa)
{
	mpq_t content;
	mpq_init(content);

	integralContent_AAU(aa, content);

	AltArrU_t *res = deepCopyPolynomial_AAU(aa);
	AAElemU_t *elems = res->elems;
	int size = res->size;
	for (int i = 0; i < size; ++i)
	{
		mpq_div(elems[i].coef, elems[i].coef, content);
	}

	mpq_clear(content);

	return res;
}

/*****************
 * Derivative / Integral
 *****************/

AltArrU_t *derivative_AAU(AltArrU_t *aa, int k)
{
	if (aa == NULL || aa->size == 0)
	{
		return NULL;
	}

	AltArrU_t *ret = makePolynomial_AAU(aa->size);
	int insertIdx = 0;

	int size = aa->size;
	AAElemU_t *__restrict__ elems = aa->elems;
	AAElemU_t *__restrict__ retElems = ret->elems;
	degree_t deg;
	mpq_t mpqDeg;
	mpq_init(mpqDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i)
	{
		if (!(elems[i].deg))
		{
			continue;
		}
		deg = elems[i].deg;
		if (deg < k)
		{
			continue;
		}

		retElems[insertIdx].deg = (elems[i].deg);
		mpq_init(retElems[insertIdx].coef);
		mpq_set(retElems[insertIdx].coef, elems[i].coef);

		mpq_set_ui(mpqDeg, deg, 1ul);
		for (int j = 0; j < k; ++j)
		{
			mpq_mul(retElems[insertIdx].coef, retElems[insertIdx].coef, mpqDeg);
			--deg;
			mpz_sub(mpq_numref(mpqDeg), mpq_numref(mpqDeg), mpzOne);
		}
		retElems[insertIdx].deg = deg;

		++insertIdx;
	}

	mpq_clear(mpqDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	return ret;
}
/////////////////////////////////////////////// Integral

AltArrU_t *integral_AAU(AltArrU_t *aa, int k)
{
	if (aa == NULL || aa->size == 0)
	{
		return NULL;
	}
	if (k == 0)
	{
		return deepCopyPolynomial_AAU(aa);
	}
	AltArrU_t *ret = makePolynomial_AAU(aa->size);

	int insertIdx = 0;
	int size = aa->size;
	AAElemU_t *__restrict__ elems = aa->elems;
	AAElemU_t *__restrict__ retElems = ret->elems;
	degrees_t deg;
	mpq_t mpqDeg;
	mpq_init(mpqDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i)
	{

		mpq_init(retElems[insertIdx].coef);
		mpq_set(retElems[insertIdx].coef, elems[i].coef);
		deg = elems[i].deg;
		mpq_set_ui(mpqDeg, deg + 1, 1ul);
		for (int j = 0; j < k; ++j)
		{
			mpq_div(retElems[insertIdx].coef, retElems[insertIdx].coef, mpqDeg);
			++deg;
			mpz_add(mpq_numref(mpqDeg), mpq_numref(mpqDeg), mpzOne);
		}

		retElems[insertIdx].deg = deg;

		++insertIdx;
	}

	mpq_clear(mpqDeg);
	mpz_clear(mpzOne);

	ret->size = insertIdx;
	return ret;
}

/*****************
 * evaluation using horner method
 *****************/

void evaluatePoly_AAU(AltArrU_t *aa, const ratNum_t r, ratNum_t mpqRet)
{

	int size = aa->size;
	if (aa == NULL || aa->size == 0)
	{ // for Null return ?
		mpq_set_si(mpqRet, 0l, 1l);
		return;
	}

	mpq_t mpqR;

	mpq_init(mpqR);
	mpq_set(mpqR, r);

	mpq_t mpqpow;
	mpq_init(mpqpow);
	degree_t diff;

	mpq_t mpqEval;
	mpq_init(mpqEval);

	if (size == 1)
	{
		//fprintf(stderr,"i am in path 1");

		mpq_add(mpqEval, mpqEval, aa->elems[0].coef);
		mpz_pow_ui(mpq_numref(mpqpow), mpq_numref(mpqR), aa->elems[0].deg);
		mpz_pow_ui(mpq_denref(mpqpow), mpq_denref(mpqR), aa->elems[0].deg);
		mpq_canonicalize(mpqpow);
		mpq_mul(mpqEval, mpqpow, mpqEval);
		mpq_init(mpqRet);
		mpq_set(mpqRet, mpqEval);
		return;
	}

	else
	{

		//mpq_set(mpqEval, aa->elems[0].coef);
		for (int i = 0; i < size - 1; i++)
		{

			mpq_add(mpqEval, mpqEval, aa->elems[i].coef);

			diff = (aa->elems[i].deg) - (aa->elems[i + 1].deg);
			//mpq_sub(diff,aa->elems[i].deg,aa->elems[i+1].deg);

			mpz_pow_ui(mpq_numref(mpqpow), mpq_numref(mpqR), diff);
			mpz_pow_ui(mpq_denref(mpqpow), mpq_denref(mpqR), diff);
			mpq_canonicalize(mpqpow); //
			mpq_mul(mpqEval, mpqpow, mpqEval);

			// mpq_mul(mpqEval,mpqR,mpqEval);
		}

		mpq_add(mpqEval, mpqEval, aa->elems[size - 1].coef);

		mpz_pow_ui(mpq_numref(mpqpow), mpq_numref(mpqR), aa->elems[size - 1].deg);
		mpz_pow_ui(mpq_denref(mpqpow), mpq_denref(mpqR), aa->elems[size - 1].deg);
		mpq_canonicalize(mpqpow); //

		mpq_mul(mpqEval, mpqpow, mpqEval);

		mpq_init(mpqRet);
		mpq_set(mpqRet, mpqEval);
		return;
	}
}

///////////////////////////////////// pomopo /////////////////////////////
/**
 * aa = aa * c + t * b
 * TODO: re-implement by pomodo and pseudo-pomodo ideas!
 */
AltArrU_t *pomopoU(AltArrU_t *aa, AltArrU_t *c, AltArrU_t *t, AltArrU_t *b)
{

	//printf ("nvar = %d\n", nvar);                                    // TEST
	if (t == NULL || t->size == 0)
	{
		return aa;
	}
	if (b == NULL || b->size == 0)
	{
		return aa;
	}

	AltArrU_t *tb = multiplyPolynomials_AAU_inp(t, b);
	if (aa == NULL || aa->size == 0 || c == NULL || c->size == 0)
	{
		return tb;
	}
	AltArrU_t *aac = multiplyPolynomials_AAU_inp(aa, c);
	return addPolynomials_AAU_inp(aac, tb);
}
//////////////////////////////////////////////////////////////////////////////////////// lazy and pseudo division
void multiplyByRational_AAU_inp(AltArrU_t *aa, const mpq_t z)
{
	if (aa == NULL || aa->size == 0)
	{
		return;
	}

	for (int i = 0; i < aa->size; ++i)
	{
		mpq_mul(aa->elems[i].coef, aa->elems[i].coef, z);
	}
}

void univariatePseudoDividePolynomials_AAU(AltArrU_t *c, AltArrU_t *b, AltArrU_t **res_a, AltArrU_t **res_r, int *e, int lazy)
{

	if (b->elems->deg > c->elems->deg)
	{
		*res_r = deepCopyPolynomial_AAU(c);
		*res_a = NULL;
		if (e != NULL)
		{
			*e = 0;
		}
		return;
	}

	//since we're in a field, pseudo division is just division.
	dividePolynomials_AAU(c, b, res_a, res_r);

	int steps = 0;
	if (!lazy)
	{
		degree_t d = c->elems->deg - b->elems->deg + 1;
		steps = d;

		mpq_t bPow;
		mpq_init(bPow);
		mpz_pow_ui(mpq_numref(bPow), mpq_numref(b->elems->coef), d);
		mpz_pow_ui(mpq_denref(bPow), mpq_denref(b->elems->coef), d);
		mpq_canonicalize(bPow);
		multiplyByRational_AAU_inp(*res_a, bPow);
		multiplyByRational_AAU_inp(*res_r, bPow);
		mpq_clear(bPow);

		// for (int j = 0; j < d; ++j) {
		// multiplyByRational_AA_inp(*res_a, b->elems->coef);
		// multiplyByRational_AA_inp(*res_r, b->elems->coef);
		// }
	}

	if (e != NULL)
	{
		*e = steps;
	}
	return;
}
//////////////////////////////////////////////////////////////////////////////////////////   isZero()

int isZero(AltArrU_t *a)
{

	if (a == NULL || a->size == 0 ||
		((a->elems[0].deg) == 0 && mpq_sgn(a->elems[0].coef) == 0))
	{
		return 1;
	}

	return 0;
}

/////////////////////////////////////////////////////// GCD  ///////////////////////////

AltArrU_t *univariateGCD_AAU(AltArrU_t *a, AltArrU_t *b)
{

	if (a == NULL || a->size == 0)
	{
		return deepCopyPolynomial_AAU(b);
	}
	if (b == NULL || b->size == 0)
	{
		return deepCopyPolynomial_AAU(a);
	}

	AltArrU_t *r0 = NULL;
	AltArrU_t *r1 = NULL;
	AltArrU_t *r2 = NULL;

	mpq_t c0;
	mpq_t c1;
	mpq_init(c0);
	mpq_init(c1);

	if ((a->elems->deg) > (b->elems->deg))
	{

		integralContent_AAU(a, c0);

		r0 = primitivePart_AAU(a);

		//	r0 = primitivePartAndContent_AA(a, c0);

		integralContent_AAU(b, c1);

		r1 = primitivePart_AAU(b);
		//	r1 = primitivePartAndContent_AA(b, c1);
	}
	else
	{
		integralContent_AAU(b, c0);

		r0 = primitivePart_AAU(b);

		//	r0 = primitivePartAndContent_AA(b, c0);

		integralContent_AAU(a, c1);

		r1 = primitivePart_AAU(a);

		//r1 = primitivePartAndContent_AA(a, c1);
	}
    int count =0;
	AltArrU_t *quo = NULL;
	while (r1 != NULL && r1->size > 0)
	{
		dividePolynomials_AAU(r0, r1, &quo, &r2);
		freePolynomial_AAU(quo);
		quo = NULL;

		freePolynomial_AAU(r0);
		r0 = r1;
		r1 = r2;
		r1 = primitivePart_AAU(r1);
		//primitivePart_AA_inp(r1);
		r2 = NULL;
		count=count+1;
		fprintf(stderr,"\n\neterate=%d\n",count);
	}

	freePolynomial_AAU(r1);
	freePolynomial_AAU(r2);
	if (r0 != NULL && r0->size > 0 && mpq_sgn(r0->elems->coef) < 0)
	{
		negatePolynomial_AAU(r0);
	}
	mpq_clear(c0);
	mpq_clear(c1);

	return r0;
}
////////////////////////////////////////////////////////////// Square free //////////////////////////////////////////////

AltArrU_t **squareFree_AAU(AltArrU_t *A, int *factsize)
{
	mpq_t content; //=c
	mpq_init(content);
	integralContent_AAU(A, content);
	//gmp_fprintf(stderr, "\n%Qd", content);

	AltArrU_t *Sstar = makePolynomial_AAU(A->size);

	AltArrU_t *r = makePolynomial_AAU(A->size);
	AltArrU_t *Y = makePolynomial_AAU(A->size);

	AltArrU_t *S = primitivePart_AAU(A);		//  A/c->S
	AltArrU_t *dS = derivative_AAU(S, 1);		//ds/dx->S'
	AltArrU_t *Smin = univariateGCD_AAU(dS, S); //gcd(s,s')->S minus

	dividePolynomials_AAU(S, Smin, &Sstar, &r); //S/S minus->S*
	dividePolynomials_AAU(dS, Smin, &Y, &r);	//S'/Sminus->Y

	AltArrU_t *dSstar = derivative_AAU(Sstar, 1);

	AltArrU_t *Z = subPolynomials_AAU(Y, dSstar); //Y-dS*/dx->Z
	int k = 1;

	AltArrU_t **A_ = (AltArrU_t **)malloc(sizeof(AltArrU_t) * k);
	// fprintf(stderr,"\n\nZ= ")  ;printAAU(Z);

	while (!isZero(Z))

	{
		//  fprintf(stderr,"\n\nZ= ")  ;printAAU(Z);
		// fprintf(stderr,"\n\nY= ")  ;printAAU(Y);
		//fprintf(stderr,"\n\nSstar= ")  ;printAAU(Sstar);
		//  fprintf(stderr,"\n\nk=%d ",k)  ;
		if (k == 1)
		{
			A_[0] = univariateGCD_AAU(Sstar, Z);

			dividePolynomials_AAU(Sstar, A_[0], &Sstar, &r); //s*/Ak->s*

			dividePolynomials_AAU(Z, A_[0], &Y, &r); //Z/Ak->Y

			dSstar = derivative_AAU(Sstar, 1); //S*'

			Z = subPolynomials_AAU(Y, dSstar); // Y-ds*/dx->Z

			/*          for(int i = 0; i < A_[0]->size; ++i) {
	//	gmp_fprintf(stderr, "\ncontent =%Qd", content);
	//	gmp_fprintf(stderr, "\ncoef= %Qd", A_[0]->elems[i].coef );


                mpq_mul(A_[0]->elems[i].coef, A_[0]->elems[i].coef, content);
		//gmp_fprintf(stderr, "\ncoef_new= %Qd", A_[0]->elems[i].coef );
        
            } */
			k = k + 1;
		}
		else
		{
			A_[k - 1] = univariateGCD_AAU(Sstar, Z);
			dividePolynomials_AAU(Sstar, A_[k - 1], &Sstar, &r);
			dSstar = derivative_AAU(Sstar, 1); //S*'

			dividePolynomials_AAU(Z, A_[k - 1], &Y, &r);
			Z = subPolynomials_AAU(Y, dSstar);

			k = k + 1;

			A_ = (AltArrU_t **)realloc(A_, k * sizeof(AltArrU_t));
			// fprintf(stderr,"\n==================k=%d",k);
		}

		// fprintf(stderr,"\nk=%d",k);
	}
	//fprintf(stderr,"\nsstar=");

	A_[k - 1] = deepCopyPolynomial_AAU(Sstar);
	//	 A_= (AltArrU_t**) realloc(A_,(k+1)*sizeof(AltArrU_t));

	//  A_[k]->elems[0].deg=0;
	//	mpq_set( A_[k]->elems[0].coef,content);

	//fprintf(stderr,"k=%d",k);
	*factsize = k;
	A_ = (AltArrU_t **)realloc(A_, (k + 1) * sizeof(AltArrU_t));
	A_[k] = makePolynomial_AAU(10);
	A_[k]->size = 1;
	A_[k]->elems[0].deg = 0;
	mpq_set(A_[k]->elems[0].coef, content);
	//A_[k]->size=1;

	return A_;
}

////////////////////////////////////////////////////// Random Poly

AltArrU_t *buildRandomPolyU(time_t seed, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg)
{
	if (nterms <= 0)
	{
		return NULL;
	}

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand)
	{
		time_t t = seed; //time(NULL);
		srand(t);

		gmp_randinit_default(R_STATE);
		gmp_randseed_ui(R_STATE, seed);

		fprintf(stderr, "seed: %lu\n", t);
		initRand = 1;
	}

	//degree_t maxTotalDeg = sparsity * (nterms);
	//degree_t maxUniDeg = (degree_t)ceil(pow(maxTotalDeg, 1.0) - 1);
	////////////////////////////////////////////////////////
	/*Node* head;
  Node* tail;

  Node* n = (Node*) malloc(sizeof(Node));
  n->next = NULL;*/

	AltArrU_t *n = (AltArrU_t *)malloc(sizeof(AltArrU_t));

	n->size = 0;
	//int alloc = nterms;
	n->elems = (AAElemU_t *)malloc(sizeof(AAElemU_t) * nterms);
	//////////////////////////////////////////////////////////////////////////////

	degree_t degs = 0; // (degrees_t) calloc(nvar, sizeof(degree_t));
	if (sparsity == 2)
	{
		//if dense always include the constant
		degs = 0;
	}
	else
	{
		degs = rand() % 2; //50/50 chance of first term being a constant;
	}
	degree_t lastDegs = degs;

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);

	while (mpz_sgn(mpzNum) == 0)
	{
		mpz_urandomb(mpzNum, R_STATE, coefBound);
	}
	while (mpz_sgn(mpzDen) == 0)
	{
		mpz_urandomb(mpzDen, R_STATE, coefBound);
	}

	// long int coef_l = rand() % coefBound;
	// if (coef_l == 0) {
	// coef_l = 1l;
	// }

	if (includeNeg && rand() % 2)
	{
		//50/50 chance of being negative
		mpz_neg(mpzNum, mpzNum);
		// coef_l *= -1;
	}
	n->elems[(n->size)].deg = degs;
	mpq_init(n->elems[n->size].coef);

	mpz_set(mpq_numref(n->elems[n->size].coef), mpzNum);
	mpz_set(mpq_denref(n->elems[n->size].coef), mpzDen);
	mpq_canonicalize(n->elems[n->size].coef);
	// mpq_set_si(n->coef, coef_l, 1ul);
	// n->coef = coef_ul;
	n->size += 1;
	--nterms;

	degree_t step = 0;
	for (int i = 0; i < nterms; ++i)
	{
		// n = (Node*) malloc(sizeof(Node));

		// n->next = NULL;

		//n->elems;
		step = 0;
		while (step == 0)
		{
			step = rand() % sparsity;
		}
		////////////
		degs = lastDegs + step;
		//////////////////////////////////////
		mpz_urandomb(mpzNum, R_STATE, coefBound);
		while (mpz_sgn(mpzNum) == 0)
		{
			mpz_urandomb(mpzNum, R_STATE, coefBound);
		}
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while (mpz_sgn(mpzDen) == 0)
		{
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}

		// coef_l = 0l;
		// while (coef_l == 0) {
		// coef_l = rand() % coefBound;
		// }
		if (includeNeg && rand() % 2)
		{
			mpz_neg(mpzNum, mpzNum);
			// coef_l *= -1;
		}

		n->elems[n->size].deg = degs;
		lastDegs = degs;
		mpq_init(n->elems[n->size].coef);
		mpz_set(mpq_numref(n->elems[n->size].coef), mpzNum);
		mpz_set(mpq_denref(n->elems[n->size].coef), mpzDen);
		mpq_canonicalize(n->elems[n->size].coef);
		n->size += 1;
		// mpq_set_si(n->coef, coef_l, 1l);
		// n->coef = coef_ul;

		// tail->next = n;
		//tail = n;
	}
	int itera;
	//Now reverse the poly so it is in decreasing order;
	if ((n->size) % 2 == 0)
	{
		itera = ((n->size) - 1) / 2;
	}
	else
	{
		itera = ((n->size) / 2) - 1;
	}
	for (int i = 0; i <= itera; ++i)
	{

		degree_t temp_degree = n->elems[i].deg;

		ratNum_t temp_coef;
		mpq_init(temp_coef);
		mpq_set(temp_coef, n->elems[i].coef);

		n->elems[i].deg = n->elems[(n->size) - 1 - i].deg;

		mpq_init(n->elems[i].coef);
		mpq_set(n->elems[i].coef, n->elems[(n->size) - 1 - i].coef);

		n->elems[(n->size) - 1 - i].deg = temp_degree;

		mpq_init(n->elems[(n->size) - 1 - i].coef);
		mpq_set(n->elems[(n->size) - 1 - i].coef, temp_coef);
	}

	/*
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
  head = cur;*/
	n->alloc = (n->size);
	return n;
}

AltArrU_t* buildRandomPolyFromMax_seededgenU(time_t seed,  const int maxDegs, unsigned long int coefBound, float sparsity, int includeNeg) {

	static int initRand = 0;
	static gmp_randstate_t R_STATE;
	if (!initRand) {
		time_t t = seed;
		srand(t);

		gmp_randinit_default (R_STATE);
		gmp_randseed_ui(R_STATE, t);

		// fprintf(stderr, "build random poly seed: %lu\n", t);
		initRand = 1;
	}

	mpz_t mpzNum, mpzDen;
	mpz_init(mpzNum);
	mpz_init(mpzDen);

	unsigned long long int maxTerms = maxDegs+1;
	


	AltArrU_t* res = makePolynomial_AAU(maxTerms);
	AAElemU_t* elems = res->elems;
		
	float eps = 1e-6;
	if (sparsity < eps) {                   ////////////////// dens poly should comple for the case of dense
		//generate dense poly;
/* 		degree_t degs = 0;
		int curIdx = 0;
		degree_t curDegs;            //?? curDegs
		

		for (int i = nvar-1; i >= 0; --i) {
			if (i == nvar - 1) {
				for (int j = curDegs; j >= 0; --j) {
					mpz_urandomb(mpzNum, R_STATE, coefBound);
					while(mpz_sgn(mpzNum) == 0) {
						mpz_urandomb(mpzNum, R_STATE, coefBound);
					}
					mpz_urandomb(mpzDen, R_STATE, coefBound);
					while(mpz_sgn(mpzDen) == 0) {
						mpz_urandomb(mpzDen, R_STATE, coefBound);
					}
					if (includeNeg && rand() % 2) {
						mpz_neg(mpzNum, mpzNum);
					}

					mpq_init(elems[curIdx].coef);
					mpz_set(mpq_numref(elems[curIdx].coef), mpzNum);
					mpz_set(mpq_denref(elems[curIdx].coef), mpzDen);
					mpq_canonicalize(elems[curIdx].coef);

					degs = 0;
					
					elems[curIdx].deg = degs;
					
					--(curDegs);
					++curIdx;
				}
			} 
		}

		
		res->size = curIdx; */
                          //////                for dens poly

	} else if (sparsity >= 1.0f) {
		mpz_urandomb(mpzNum, R_STATE, coefBound);
		while(mpz_sgn(mpzNum) == 0) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
		}
		mpz_urandomb(mpzDen, R_STATE, coefBound);
		while(mpz_sgn(mpzDen) == 0) {
			mpz_urandomb(mpzDen, R_STATE, coefBound);
		}
		if (includeNeg && rand() % 2) {
			mpz_neg(mpzNum, mpzNum);
		}

		mpq_init(elems[0].coef);
		mpz_set(mpq_numref(elems[0].coef), mpzNum);
		mpz_set(mpq_denref(elems[0].coef), mpzDen);
		mpq_canonicalize(elems[0].coef);

		elems[0].deg = maxDegs;


		res->size = 1;
	} else {
		//we are in the general case
		unsigned long long int targetTerms = ceil(maxTerms * (1 - sparsity));
	
		//incremement maxDegs for mod operation below
		int maxDegsL= maxDegs + 1;
	

		degree_t curDeg = 0;

		for (int i = 0; i < targetTerms; ++i) {
			mpz_urandomb(mpzNum, R_STATE, coefBound);
			while(mpz_sgn(mpzNum) == 0) {
				mpz_urandomb(mpzNum, R_STATE, coefBound);
			}
			mpz_urandomb(mpzDen, R_STATE, coefBound);
			while(mpz_sgn(mpzDen) == 0) {
				mpz_urandomb(mpzDen, R_STATE, coefBound);
			}
			if (includeNeg && rand() % 2) {
				mpz_neg(mpzNum, mpzNum);
			}

			mpq_init(elems[i].coef);
			mpz_set(mpq_numref(elems[i].coef), mpzNum);
			mpz_set(mpq_denref(elems[i].coef), mpzDen);
			mpq_canonicalize(elems[i].coef);

			elems[i].deg = 0;

			if (i == 0) {
			
					elems[0].deg = ( ((degree_t)maxDegs) );
							
			} else {
				                   //?
					curDeg = rand() % (maxDegsL);
					elems[i].deg = ( curDeg );
				
			}
		}

		res->size = targetTerms;
     		sortPolynomial_AAU(res);

	}


	mpz_clear(mpzNum);
	mpz_clear(mpzDen);

	resizePolynomial_AAU(res, res->size);
	return res;

}


///////////////////////////////////////////////////////////////////// ConvertfromAlttoAr
mpz_t *ConvertfromAlttoAr(AltArrU_t *F, int *deg)
{

	mpz_t *ret;
	//printAAU(F);
	ret = (mpz_t *)malloc(sizeof(mpz_t)*(F->elems[0].deg + 1));
	for (int i = 0; i <= F->elems[0].deg; i++)
	{
		mpz_init(ret[i]);
	}
	mpz_t lcm;
	mpz_init(lcm);
	mpz_set_d(lcm,1);
   
	if (F->size < 2)
	{
		mpz_set(lcm, mpq_denref(F->elems[0].coef));
		//gmp_printf("lcm if 1= %Zd\n", lcm);
	}
	else
	{
		mpz_set(lcm, mpq_denref(F->elems[0].coef));

		mpz_lcm(lcm, lcm, mpq_denref(F->elems[1].coef));
		if (F->size > 2)
		{
			for (int i = 2; i <= F->size - 1; i++)
			{
				mpz_lcm(lcm, lcm, mpq_denref(F->elems[i].coef));
				
				//mpq_denref(F->elems[i].coef);
			}
		}
		// gmp_printf("lcm if 3= %Zd\n", lcm);
	}
	mpz_t mul;
	mpz_init(mul);
	for (int i = 0; i <= F->size - 1; i++)
	{
		mpz_mul(mul, lcm, mpq_numref(F->elems[i].coef));
		mpz_divexact(mul, mul, mpq_denref(F->elems[i].coef));

		//gmp_printf("mul iterations= %Zd\n", mul);

		mpz_set(ret[F->elems[i].deg], mul);

		//gmp_printf("the rational is: %Qd\n",mpq_numref(F->elems[i].coef));

		// gmp_printf("ret= %Zd\n", ret[F->elems[i].deg ]);
	}

	*deg = F->elems[0].deg;
	
	return ret;
}
//////////////////////////////////////////////////////////////////////////////////////////////// ConvertfromArtoAlt
AltArrU_t *ConvertfromArtoAlt(mpz_t *Ar, int deg)
{

	AltArrU_t *ret = makePolynomial_AAU(deg);

	int count = 0;
	for (int i = deg - 1; i >= 0; i--)
	{

		if (mpz_cmp_si(Ar[i], 0))
		{

			/*  fprintf(stderr, "\n count =  %d \n ", count );
				gmp_fprintf(stderr, " Ar[%d] = %Zd \n", i,  Ar[i]); */

			mpq_init(ret->elems[count].coef);
			mpz_set(mpq_numref(ret->elems[count].coef), Ar[i]);
			//	gmp_fprintf(stderr, "ret->elems[%d].coef = %Qd\n",  count, Ar[i]);

			ret->elems[count].deg = i;
			count = count + 1;
			ret->size = count;
		}
	}

	return ret;
}
////////////////////////////////////////////////////////////////Modular GCD

//////////////////////////////////////////////////////////////////////
void convertToPrimeField_SUZP_inp(SUZP_t* p, Prime_ptr* pptr, duspoly_t** pp) {
	if (p == NULL || p->alloc == 0) {
		freePolynomial_spX (&*pp);
		*pp = NULL;
		return;
	}

	if ((*pp)->alloc < p->lt+1) {
		freePolynomial_spX (&*pp);
		*pp = makePolynomial_spX(p->lt + 1);
	}
	polysize_t lt = (*pp)->lt = p->lt;
	mpz_t* coefs = p->coefs;
	elem_t* primeElems = (*pp)->elems;

	unsigned long uPrime = (unsigned long) pptr->prime;

	for (polysize_t i = 0; i <= lt; ++i) {
		primeElems[i] = smallprimefield_convert_in(mpz_fdiv_ui(coefs[i], uPrime), pptr);
	}

}
SUZP_t* deepCopyPolynomial_SUZP(SUZP_t* p) {
	SUZP_t* ret = makePolynomial_SUZP(p->lt + 1);
	polysize_t lt = ret->lt = p->lt;
	mpz_t* retcoefs = ret->coefs;
	mpz_t* coefs = p->coefs;
	for (polysize_t i = 0; i <= lt; ++i) {
		mpz_init_set(retcoefs[i], coefs[i]);
	}

	return ret;
}
void content_SUZP(SUZP_t* p, mpz_t c) {
	if (p == NULL || p->alloc == 0){
		mpz_set_ui(c, 0ul);
		return;
	}

	mpz_t* coefs = p->coefs;
	polysize_t lt = p->lt;

	polysize_t i = 0;
	for (; i <= lt; ++i) {
		if (mpz_sgn(coefs[i]) != 0) {
			mpz_set(c, coefs[i]);
			break;
		}
	}

	for (; i <= lt; ++i) {
		if (mpz_sgn(coefs[i]) != 0) {
			mpz_gcd(c, c, coefs[i]);
			if (mpz_cmp_si(c, 1l) == 0) {
				return;
			}
		}
	}

	if (mpz_sgn(coefs[lt]) < 0 && mpz_sgn(c) > 0) {
		mpz_neg(c, c);
	}

	return;
}
SUZP_t* primitivePart_SUZP(SUZP_t* p) {

	SUZP_t* pp = deepCopyPolynomial_SUZP(p);
	mpz_t cont;
	mpz_init(cont);
	content_SUZP(pp, cont);
	if (mpz_cmp_si(cont, 1l) == 0) {
		mpz_clear(cont);
		return pp;
	}

	divideByIntegerExact_SUZP_inp(pp, cont);
	mpz_clear(cont);
	return pp;
}
///////////////////////////// divide test
int divideTest_SUZP(SUZP_t* a, SUZP_t* b, SUZP_t** q) {
	if (a == NULL) {
		return 1;
	}
	if (b == NULL) {
		fprintf(stderr, "Tried to divide by 0 in divideTest_DUZP!\n");
		return 0;
	}

	polysize_t lt_a = a->lt;
	polysize_t lt_b = b->lt;

	if (lt_b > lt_a) {
		return 0;
	}

	polysize_t diff_deg = lt_a - lt_b;

	mpz_t* lcb = &(b->coefs[lt_b]);

	SUZP_t* rr = deepCopyPolynomial_SUZP(a);
	mpz_t* rr_coefs = rr->coefs;

	SUZP_t* quo = makePolynomial_SUZP(diff_deg + 1);

	polysize_t i = 0, j = 0;
    for (i = diff_deg; i >= 0; --i) {
    	if (mpz_sgn(rr_coefs[lt_b + i]) != 0) {
    		if (!mpz_divisible_p(rr_coefs[lt_b + i], *lcb)) {
			    rr->lt = lt_a;
			    freePolynomial_SUZP(rr);
			    for (j = i+1; j <= diff_deg; ++j) {
			    	mpz_clear(quo->coefs[j]);
			    }
			    quo->lt = -1;
			    freePolynomial_SUZP(quo);
    			return 0;
    		}

    		mpz_init(quo->coefs[i]);
    		mpz_divexact(quo->coefs[i], rr_coefs[lt_b + i], *lcb);
    	    
  			mpz_set_si(rr_coefs[lt_b + i], 0l);
    	    --(rr->lt);
    	    for (j = 0; j < lt_b; ++j) {
    	    	mpz_submul(rr_coefs[i+j], b->coefs[j], quo->coefs[i]);
    	    }
    	} else {
    		mpz_init(quo->coefs[i]); //quo[i] = 0;
    	}
    }

    for (i = rr->lt; i >= 0; --i) {
    	if (mpz_sgn(rr_coefs[i]) != 0) {
		    rr->lt = lt_a;
		    freePolynomial_SUZP(rr);
		   	freePolynomial_SUZP(quo);
    		return 0;
    	}
    }

    rr->lt = lt_a;
    freePolynomial_SUZP(rr);
	quo->lt = diff_deg;

	if (q != NULL) {
		*q = quo;
	} else {
		freePolynomial_SUZP(quo);
	}
    return 1;
}
////////////////////////////////////////////////////modular gcd
SUZP_t* primitiveGCD_SUZP(SUZP_t* a, SUZP_t* b) {

	mpz_t lcg;
	mpz_init(lcg);
	mpz_gcd(lcg, a->coefs[a->lt], b->coefs[b->lt]);

	polysize_t d = MIN(a->lt, b->lt);

	SUZP_t* g_work = makePolynomial_SUZP(d+1); 
	mpz_t* g_work_coefs = g_work->coefs;
	for (polysize_t i = 0; i <= d; ++i) {
		mpz_init(g_work_coefs[i]);
	}
	g_work->lt = d;

	mpz_t m; 
	mpz_init_set_ui(m, 1l);
	mpz_t halfm;
	mpz_init(halfm);
	mpz_t newm;
	mpz_init(newm);
	mpz_t lcgm; 
	mpz_init_set(lcgm, lcg);

	mpz_t s, t, mpz_prime; //for CRT
	mpz_inits(s, t, mpz_prime, NULL);

	int primeIdx = 0;
	Prime_ptr* Pptr = (Prime_ptr*) malloc(sizeof(Pptr));

	duspoly_t* moda = makePolynomial_spX(a->lt + 1);
	duspoly_t* modb = makePolynomial_spX(b->lt + 1);
	duspoly_t* modg = NULL;
	
	int fibs[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584};
	int fibIdx = 0, nFibs = 17;
  
	for(; primeIdx < n_prime64_ptr1; ++primeIdx) {
		 *Pptr = prime64_ptr1[primeIdx];
	
		if (mpz_divisible_ui_p(lcgm, (unsigned long long) Pptr->prime)) {

			continue;
		}

		convertToPrimeField_SUZP_inp(a, Pptr, &moda);
		convertToPrimeField_SUZP_inp(b, Pptr, &modb);
		plainGCDInForm_spX(moda, modb, &modg, Pptr);

		if (modg->lt == 0) {
			SUZP_t* ret = makePolynomial_SUZP(1);
			mpz_init(ret->coefs[0]);
			mpz_set_si(ret->coefs[0], 1l);
			ret->lt = 0;
			return ret;
		} else if (modg->lt < g_work->lt) {
			polysize_t lt = modg->lt;
			polysize_t i = 0;
			elem_t* primeElems = modg->elems;
			long long int halfP = (Pptr->prime - 1) >> 1;
			long long int tmp;
			long long int lcg_pf = smallprimefield_convert_in(mpz_fdiv_ui(lcg, Pptr->prime), Pptr);
			for (; i <= lt; ++i) {
				tmp = smallprimefield_mul(primeElems[i], lcg_pf, Pptr);
				tmp = smallprimefield_convert_out(tmp, Pptr);
				if (tmp > halfP) {
					tmp -= Pptr->prime;
				}
				mpz_set_si(g_work_coefs[i], tmp);
			}
			lt = g_work->lt;
			for (; i <= lt; ++i) {
				mpz_clear(g_work_coefs[i]);
			}
			g_work->lt = modg->lt;

			mpz_set_si(m, Pptr->prime);
			mpz_mul(lcgm, m, lcg);

			freePolynomial_spX (&modg);
			modg = NULL;
			continue;
		} else if (modg->lt > g_work->lt) {

			freePolynomial_spX (&modg);
			modg = NULL;
			continue;
		}	

		long long int lcg_pf = smallprimefield_convert_in(mpz_fdiv_ui(lcg, Pptr->prime), Pptr);
		scalarMulPolynomialInForm_spX_inp (&modg, lcg_pf, Pptr);

		mpz_set_si(mpz_prime, Pptr->prime);
		mpz_gcdext(lcgm, s, t, mpz_prime, m);

		mpz_mul(newm, m, mpz_prime);
		mpz_sub_ui(halfm, newm, 1ul);
		mpz_fdiv_q_2exp(halfm, halfm, 1ul);

		//CRT
		polysize_t lt = modg->lt;
		unsigned long long modg_out;
		for (int i = 0; i <= lt; ++i) {
 
			modg_out = (unsigned long long) smallprimefield_convert_out(modg->elems[i], Pptr);
			mpz_sub_ui(g_work_coefs[i], g_work_coefs[i], modg_out);
			mpz_mul(g_work_coefs[i], g_work_coefs[i], s);
			mpz_mul(g_work_coefs[i], g_work_coefs[i], mpz_prime);
			mpz_add_ui(g_work_coefs[i], g_work_coefs[i], modg_out);
			mpz_mod(g_work_coefs[i], g_work_coefs[i], newm);

			
			if(mpz_cmp(g_work_coefs[i], halfm) > 0) {
				mpz_sub(g_work_coefs[i], g_work_coefs[i], newm);
			}
		}

		freePolynomial_spX (&modg);

		int doDivideTest = fibIdx < nFibs ? (primeIdx + 1 == fibs[fibIdx++]) : 1;
		if (doDivideTest) {
			SUZP_t* primG = primitivePart_SUZP(g_work);
			if (divideTest_SUZP(a, primG, NULL) && divideTest_SUZP(b, primG, NULL)) {
				mpz_clears(m, halfm, newm, lcgm, s, t, mpz_prime, NULL);
				free(Pptr);
				freePolynomial_spX (&moda);
				freePolynomial_spX (&modb);
				freePolynomial_SUZP(g_work);

				return primG;
			}
		}

		mpz_set(m, newm);
		mpz_mul(lcgm, m, lcg);
	}
     

	//all primes failed.... 
	SUZP_t* ret = makePolynomial_SUZP(1);
	mpz_init(ret->coefs[0]);
	mpz_set_si(ret->coefs[0], 1l);
	ret->lt = 0;

	mpz_clears(m, halfm, newm, lcgm, s, t, mpz_prime, NULL);
	free(Pptr);
	freePolynomial_spX (&moda);
	freePolynomial_spX (&modb);
	freePolynomial_SUZP(g_work);    

	return ret;
}

AltArrU_t *modularGCD(AltArrU_t *F, AltArrU_t *G){
	int curdf;
	int curdg;
    
	AltArrU_t* res;
	mpz_t *f = ConvertfromAlttoAr(F, &curdf);
	mpz_t *g = ConvertfromAlttoAr(G, &curdg);
	SUZP_t* a=makePolynomial_SUZP(curdf+1);
	SUZP_t* b=makePolynomial_SUZP(curdg+1);


 


	for(int i =0;i<=curdf;i++){
		a->lt=i;
		mpz_init(a->coefs[i]);
		mpz_set(a->coefs[i],f[i]);

	} 
	for(int i =0;i<=curdg;i++){
		b->lt=i;
		mpz_init(b->coefs[i]);
		mpz_set(b->coefs[i],g[i]);
		
	}
    SUZP_t*gcd= primitiveGCD_SUZP(a,b);
	mpz_t* ret;
	ret= (mpz_t*)malloc( sizeof(mpz_t)*(curdf+curdg));
	for(int i=0;i<gcd->alloc;i++){
		mpz_init(ret[i]);
		}

	for(int i=0;i<gcd->lt+1;i++){
		gmp_fprintf(stderr,"\ni=%d  coef=%Zd",i,gcd->coefs[i]);
		mpz_set(ret[i],gcd->coefs[i]);
	}	
	res=ConvertfromArtoAlt(ret,gcd->lt+1);


    return res;
}