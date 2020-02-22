
#include "IntegerPolynomial/SMZP_Support_Recursive_Unpacked.h"
#include "IntegerPolynomial/SMZP_Support_Recursive.h"
#include "IntegerPolynomial/SMZP_Support_Unpacked.h"

/**
 * Convert a polynomial to a recursively viewed polynomial using the variable
 * of highest precedence (index 0) as the main variable.
 *
 * Note this conversion is done INPLACE. And therefore the input Node* is 
 * invalidated by this operation.
 *
 * returns teh recursively viewed polynomial.
 */
RecArrZ_t* convertToRecursiveArrayZ_unpk(AltArrZ_t* aa) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (!aa->unpacked) {
		return convertToRecursiveArrayZ(aa);
	}

	AAZElem_t* elems = aa->elems;
	degree_t* degs = (degree_t*) aa->elems->degs;
	degree_t mvarDeg =  degs[0];
	RecArrZ_t* poly = (RecArrZ_t*) malloc(sizeof(RecArrZ_t));
	poly->alloc = mvarDeg;
	RecArrElemZ_t* recElems = (RecArrElemZ_t*) malloc(sizeof(RecArrElemZ_t)*(mvarDeg+1)); 
	poly->elems = recElems;

	int nvar = aa->nvar;
	int curIdx = 0, lastSize = 0;
	degree_t curDeg;
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		curDeg = degs[i*nvar];
		
		if (curDeg != mvarDeg) {
			recElems[curIdx].exp = mvarDeg;
			recElems[curIdx].coefSize = i - lastSize;
			recElems[curIdx].coef = elems + lastSize;
			++curIdx;
			lastSize = i;
			mvarDeg = curDeg;
		}

		//throw away exp of mvar for coef.
		degs[i*nvar] = 0;
	}

	//one more to commit after loop
	recElems[curIdx].exp = mvarDeg;
	recElems[curIdx].coefSize = AA_SIZE(aa) - lastSize;
	recElems[curIdx].coef = elems + lastSize;
	++curIdx;

	poly->size = curIdx;
	poly->origAA = aa;
	poly->unpacked = 1;

	return poly;

}

/**
 * Convert a recursively viewed polynomial to it's distributed form. 
 * This conversion is done INPLACE and therefore the original recursively
 * viewed polynomial is then invalidated.
 * 
 * nvar is the number of variables the returning AltArr_t should have. That is, 
 * the number of variables in the RecArr's coefficients + 1.
 *
 * returns the distributed polynomial.
 */
AltArrZ_t* convertFromRecursiveArrayZ_unpk(RecArrZ_t* poly, int nvar) {
	if (poly == NULL || poly->size == 0) {
		return NULL;
	}

	if (!poly->unpacked) {
		return convertFromRecursiveArrayZ(poly, nvar);
	}

	AltArrZ_t* aa = NULL;
	if (poly->origAA != NULL) {
		aa = poly->origAA;
	} else {
		aa = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	}


	RecArrElemZ_t* recElems = poly->elems;
	AAZElem_t* elems = poly->elems->coef;
	degree_t* degs = (degree_t*) elems->degs;
	aa->elems = elems;
	register int curSize = 0;
	for (int i = 0; i < poly->size; ++i) {
		degrees_t curDeg = recElems[i].exp;
		for (int j = 0; j < recElems[i].coefSize; ++j) {
			degs[(curSize + j)*nvar] = curDeg;
		}

		curSize += recElems[i].coefSize;
	}

	freeRecArrayZ(poly);

	aa->size = curSize;
	aa->alloc = curSize;
	aa->nvar = nvar;
	aa->unpacked = 1;
	return aa;
}

