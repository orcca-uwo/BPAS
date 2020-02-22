
#ifndef _EXP_PACKING_DEFS_H_
#define _EXP_PACKING_DEFS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int degree_t;

typedef long long unsigned int degrees_t; 

typedef int cmpExp_t;


#include "ExponentPackingDefinesValues.h"
#include "PackedMonomialDivisionDefines.h"


/**
 * Get the nth exponent from a packed exponent vector, given the mask
 * and offset for that exponent.
 */
#define GET_NTH_EXP(A, M, O) ( ((A) & (M)) >> (O))

/**
 * Compare two exponent vectors for lexicographical term order.
 * a: int array representing exponent vector
 * b: int array representing exponent vector
 * returns negative int if a < b, 0 if equal, a positive int if a > b
 */
// #define compareExponentVectors(A, B) ((A) - (B))
#define compareExponentVectors(A, B) ((A) < (B) ? -1 : ((A) > (B) ? 1 : 0) )

static inline cmpExp_t compareExponentVectors_unpk(degrees_t degs_a, degrees_t degs_b, int nvar) {
	if (nvar == 0) {
		return 0;
	}
	if (degs_a == 0 || degs_b == 0) {
		fprintf(stderr, "compareExponentVectors_unpk: exponent pointer is NULL!\n");
		exit(1);
	}

	degree_t* a = (degree_t*) degs_a;
	degree_t* b = (degree_t*) degs_b;
	degree_t diff = a[0] - b[0];
	int i = 1;
	while (diff == 0 && i < nvar) {
		diff = a[i] - b[i];
		++i;
	}
	return diff;
}

#define isLessExponentVectors(A, B) ((A) < (B))

static inline cmpExp_t isLessExponentVectors_unpk(degrees_t degs_a, degrees_t degs_b, int nvar) {
	if (nvar == 0) {
		return 0;
	}
	if (degs_a == 0 || degs_b == 0) {
		fprintf(stderr, "isLessExponentVectors_unpk: exponent pointer is NULL!\n");
		exit(1);
	}

	degree_t* a = (degree_t*) degs_a;
	degree_t* b = (degree_t*) degs_b;
	for (int k = 0; k < nvar; ++k) {
		if (a[k] < b[k]) {
			return 1;
		}
		if (a[k] > b[k]) {
			return 0;
		}
	}

	return 0;
}

#define isEqualExponentVectors(A, B) ((A) == (B))

static inline cmpExp_t isEqualExponentVectors_unpk(degrees_t degs_a, degrees_t degs_b, int nvar) {
	if (nvar == 0) {
		return 1;
	}
	if (degs_a == 0 || degs_b == 0) {
		fprintf(stderr, "isLessExponentVectors_unpk: exponent pointer is NULL!\n");
		exit(1);
	}

	degree_t* a = (degree_t*) degs_a;
	degree_t* b = (degree_t*) degs_b;
	for (int k = 0; k < nvar; ++k) {
		if (a[k] != b[k]) {
			return 0;
		}
	}

	return 1;
}

#define isGreaterExponentVectors(A, B) ((A) > (B))

static inline cmpExp_t isGreaterExponentVectors_unpk(degrees_t degs_a, degrees_t degs_b, int nvar) {
	if (nvar == 0) {
		return 0;
	}
	if (degs_a == 0 || degs_b == 0) {
		fprintf(stderr, "isGreaterExponentVectors_unpk: exponent pointer is NULL!\n");
		exit(1);
	}

	degree_t* a = (degree_t*) degs_a;
	degree_t* b = (degree_t*) degs_b;
	for (int k = 0; k < nvar; ++k) {
		if (a[k] > b[k]) {
			return 1;
		}
		if (a[k] < b[k]) {
			return 0;
		}
	}

	return 0;
}

/**
 * Check if an exponent vector has all 0 exponents.
 * a: the exponent vector to check
 */ 
#define isZeroExponentVector(A) ((A) == 0)

static inline cmpExp_t isZeroExponentVector_unpk(degrees_t degs, int nvar) {
	if (nvar == 0) {
		return 1;
	}
	if (degs == 0) {
		fprintf(stderr, "isZeroExponentVector_unpk: exponent pointer is NULL!\n");
		free((int*)-1);
		exit(1);
	}
	degree_t* degs_unpk = (degree_t*) degs;
	for (int k = 0; k < nvar; ++k) {
		if (degs_unpk[k] > 0) {
			return 0;
		}
	}
	return 1;
}

static inline void setExponentVector_unpk(degrees_t dest, degrees_t src, int nvar) {
	if (nvar == 0) {
		return;
	}
	degree_t* dest_p = (degree_t*) dest;
	degree_t* src_p = (degree_t*) src;
	memcpy(dest_p, src_p, nvar*sizeof(degree_t));
}

/**
 * Add two exponent vectors, a and b, and store the result in c.
 * a,b,c: degree_t  
 */
#define addExponentVectors(A, B, C) ((C) = (A) + (B)) 

static inline void addExponentVectors_unpk(degrees_t degs_a, degrees_t degs_b, degrees_t degs_c, int nvar) {
	if (nvar == 0) {
		return;
	}
	if (degs_a == 0 || degs_b  == 0 || degs_c == 0) {
		fprintf(stderr, "addExponentVectors_unpk: exponent pointer is NULL!\n");
		exit(1);
	}

	degree_t* a = (degree_t*) degs_a;
	degree_t* b = (degree_t*) degs_b;
	degree_t* c = (degree_t*) degs_c;
	for (int k = 0; k < nvar; ++k) {
		c[k] = a[k] + b[k];
	}
}


/**
 * Subtract two exponent vectors, b from a, and store the result in c.
 * It is assumed that this operation is valid.
 * see also, monomialDivideTest.
 *
 * a,b,c: degree_tarrays
 */
#define subtractExponentVectors(A, B, C) ((C) = (A) - (B))

static inline void subtractExponentVectors_unpk(degrees_t degs_a, degrees_t degs_b, degrees_t degs_c, int nvar) {
	if (nvar == 0) {
		return;
	}
	if (degs_a == 0 || degs_b  == 0 || degs_c == 0) {
		fprintf(stderr, "subtractExponentVectors_unpk: exponent pointer is NULL!\n");
		exit(1);
	}

	degree_t* a = (degree_t*) degs_a;
	degree_t* b = (degree_t*) degs_b;
	degree_t* c = (degree_t*) degs_c;
	for (int k = 0; k < nvar; ++k) {
		c[k] = a[k] - b[k];
	}
}

//forward declare
static const degrees_t* getMaxExpArray (int nvar);
static const degrees_t* getExpMaskArray(int nvar);
static const int* getExpOffsetArray(int nvar);

static inline int isPackableExponentVector(degree_t* degs, int nvar) {
	const degrees_t* maxExps = getMaxExpArray(nvar);
	int ret = 1;
	for (int k = 0; k < nvar; ++k) {
		if ((degrees_t) degs[k] > maxExps[k]) {
			ret = 0;
			break;
		}
	}
	return ret;
}

static inline void unpackExponentVector(degrees_t packedExp, degree_t* unpackedExp, int nvar) {
	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);

	for (int k = 0; k < nvar; ++k) {
		unpackedExp[k] = GET_NTH_EXP(packedExp, masks[k], sizes[k]);
	}

}


static int checkValidMonomialMult(degrees_t maxA, degrees_t maxB, int nvar) {
	if (nvar <= 0) {
		return 1;
	}

	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* maxExps = getMaxExpArray(nvar);

	degree_t a, b;
	for (int j = 0; j < nvar; ++j) {
		a = GET_NTH_EXP(maxA, masks[j], sizes[j]);
		b = GET_NTH_EXP(maxB, masks[j], sizes[j]);
		if ((degrees_t) (a + b) > maxExps[j]) {
			fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for mult at index %d; %d + %d.\n", j, a, b);
			free((void*) -1);

			exit(1); //TODO
		}
	}

	return 1;
}

static int monomialMultFitsPacked(degrees_t maxA, degrees_t maxB, int nvar) {
	if (nvar <= 0) {
		return 1;
	}

	const degrees_t* masks = getExpMaskArray(nvar);
	const int* sizes = getExpOffsetArray(nvar);
	const degrees_t* maxExps = getMaxExpArray(nvar);

	int fits = 1;
	degree_t a, b;
	for (int j = 0; j < nvar; ++j) {
		a = GET_NTH_EXP(maxA, masks[j], sizes[j]);
		b = GET_NTH_EXP(maxB, masks[j], sizes[j]);
		if ((degrees_t) (a + b) > maxExps[j]) {
			fits = 0;
			break;
		}
	}

	return fits;
}



static int getMVarExpOffset(int nvar) {
	switch(nvar) {
    	case 1: {
			return EXP_OFFSET_1_V1;    		
    	}
    	case 2: {
			return EXP_OFFSET_1_V2;
    	}
    	case 3: {
			return EXP_OFFSET_1_V3;
    	}
    	case 4: {
			return EXP_OFFSET_1_V4;
    	}
    	case 5: {
			return EXP_OFFSET_1_V5;
    	}
    	case 6: {
			return EXP_OFFSET_1_V6;
    	}
    	case 7: {
			return EXP_OFFSET_1_V7;
		}
    	case 8: {
    		return EXP_OFFSET_1_V8;
    	}
    	case 9: {
    		return EXP_OFFSET_1_V9;
    	}
    	case 10: {
			return EXP_OFFSET_1_V10;
    	}
    	case 11: {
			return EXP_OFFSET_1_V11;
    	}
    	case 12: {
			return EXP_OFFSET_1_V12;
    	}
    	case 13: {
			return EXP_OFFSET_1_V13;
    	}
    	case 14: {
			return EXP_OFFSET_1_V14;
    	}
    	case 15: {
			return EXP_OFFSET_1_V15;
    	}
    	case 16: {
			return EXP_OFFSET_1_V16;
    	}
    	case 17: {
			return EXP_OFFSET_1_V17;
    	}
    	case 18: {
			return EXP_OFFSET_1_V18;
		}
    	case 19: {
			return EXP_OFFSET_1_V19;
    	}
    	case 20: {
			return EXP_OFFSET_1_V20;
    	}
    	case 21: {
			return EXP_OFFSET_1_V21;
    	}
    	case 22: {
			return EXP_OFFSET_1_V22;
    	}
    	case 23: {
			return EXP_OFFSET_1_V23;
    	}
    	case 24: {
			return EXP_OFFSET_1_V24;
    	}
    	case 25: {
			return EXP_OFFSET_1_V25;
    	}
    	case 26: {
			return EXP_OFFSET_1_V26;
    	}
    	case 27: {
			return EXP_OFFSET_1_V27;
		}
    	case 28: {
			return EXP_OFFSET_1_V28;
    	}
    	case 29: {
			return EXP_OFFSET_1_V29;
    	}
    	case 30: {
			return EXP_OFFSET_1_V30;
    	}
    	case 31: {
			return EXP_OFFSET_1_V31;
    	}
    	case 32: {
			return EXP_OFFSET_1_V32;
    	}
    	default : {
    		fprintf(stderr, "EXP SIZES NOT DEFINED FOR NVAR %d\n",nvar);
    		return -1;
    	}
    }
}



static const degrees_t* getExpMaskArray (int nvar) {
	switch(nvar) {
		case 1: {
			return (EXP_MASK_LIST);
		}
		case 2: {
			return (EXP_MASK_LIST + 1);
		}
		case 3: {
			return (EXP_MASK_LIST + 3);
		}
		case 4: {
			return (EXP_MASK_LIST + 6);
		}
		case 5: {
			return (EXP_MASK_LIST + 10);
		}
		case 6: {
			return (EXP_MASK_LIST + 15);
		}
		case 7: {
			return (EXP_MASK_LIST + 21);
		}
		case 8: {
			return (EXP_MASK_LIST + 28);
		}
		case 9: {
			return (EXP_MASK_LIST + 36);
		}
		case 10: {
			return (EXP_MASK_LIST + 45);
		}
		case 11: {
			return (EXP_MASK_LIST + 55);
		}
		case 12: {
			return (EXP_MASK_LIST + 66);
		}
		case 13: {
			return (EXP_MASK_LIST + 78);
		}
		case 14: {
			return (EXP_MASK_LIST + 91);
		}
		case 15: {
			return (EXP_MASK_LIST + 105);
		}
		case 16: {
			return (EXP_MASK_LIST + 120);
		}
		case 17: {
			return (EXP_MASK_LIST + 136);
		}
		case 18: {
			return (EXP_MASK_LIST + 153);
		}
		case 19: {
			return (EXP_MASK_LIST + 171);
		}
		case 20: {
			return (EXP_MASK_LIST + 190);
		}
		case 21: {
			return (EXP_MASK_LIST + 210);
		}
		case 22: {
			return (EXP_MASK_LIST + 231);
		}
		case 23: {
			return (EXP_MASK_LIST + 253);
		}
		case 24: {
			return (EXP_MASK_LIST + 276);
		}
		case 25: {
			return (EXP_MASK_LIST + 300);
		}
		case 26: {
			return (EXP_MASK_LIST + 325);
		}
		case 27: {
			return (EXP_MASK_LIST + 351);
		}
		case 28: {
			return (EXP_MASK_LIST + 378);
		}
		case 29: {
			return (EXP_MASK_LIST + 406);
		}
		case 30: {
			return (EXP_MASK_LIST + 435);
		}
		case 31: {
			return (EXP_MASK_LIST + 465);
		}
		case 32: {
			return (EXP_MASK_LIST + 496);
		}
    	default : {
    		fprintf(stderr, "EXP MASK ARRAY NOT DEFINED FOR NVAR %d\n",nvar);
    		return NULL;
    	}
    }
}


static const int* getExpOffsetArray (int nvar) {
	switch(nvar) {
		case 1: {
			return (EXP_OFFSET_LIST);
		}
		case 2: {
			return (EXP_OFFSET_LIST + 1);
		}
		case 3: {
			return (EXP_OFFSET_LIST + 3);
		}
		case 4: {
			return (EXP_OFFSET_LIST + 6);
		}
		case 5: {
			return (EXP_OFFSET_LIST + 10);
		}
		case 6: {
			return (EXP_OFFSET_LIST + 15);
		}
		case 7: {
			return (EXP_OFFSET_LIST + 21);
		}
		case 8: {
			return (EXP_OFFSET_LIST + 28);
		}
		case 9: {
			return (EXP_OFFSET_LIST + 36);
		}
		case 10: {
			return (EXP_OFFSET_LIST + 45);
		}
		case 11: {
			return (EXP_OFFSET_LIST + 55);
		}
		case 12: {
			return (EXP_OFFSET_LIST + 66);
		}
		case 13: {
			return (EXP_OFFSET_LIST + 78);
		}
		case 14: {
			return (EXP_OFFSET_LIST + 91);
		}
		case 15: {
			return (EXP_OFFSET_LIST + 105);
		}
		case 16: {
			return (EXP_OFFSET_LIST + 120);
		}
		case 17: {
			return (EXP_OFFSET_LIST + 136);
		}
		case 18: {
			return (EXP_OFFSET_LIST + 153);
		}
		case 19: {
			return (EXP_OFFSET_LIST + 171);
		}
		case 20: {
			return (EXP_OFFSET_LIST + 190);
		}
		case 21: {
			return (EXP_OFFSET_LIST + 210);
		}
		case 22: {
			return (EXP_OFFSET_LIST + 231);
		}
		case 23: {
			return (EXP_OFFSET_LIST + 253);
		}
		case 24: {
			return (EXP_OFFSET_LIST + 276);
		}
		case 25: {
			return (EXP_OFFSET_LIST + 300);
		}
		case 26: {
			return (EXP_OFFSET_LIST + 325);
		}
		case 27: {
			return (EXP_OFFSET_LIST + 351);
		}
		case 28: {
			return (EXP_OFFSET_LIST + 378);
		}
		case 29: {
			return (EXP_OFFSET_LIST + 406);
		}
		case 30: {
			return (EXP_OFFSET_LIST + 435);
		}
		case 31: {
			return (EXP_OFFSET_LIST + 465);
		}
		case 32: {
			return (EXP_OFFSET_LIST + 496);
		}
		default : {
	    		fprintf(stderr, "EXP SIZES NOT DEFINED FOR NVAR %d\n",nvar);
	    		return NULL;
	    }
	}
}

static degrees_t getMVarExpMask(int nvar) {
	switch(nvar) {
    	case 1: {
			return EXP_1_V1;    		
    	}
    	case 2: {
			return EXP_1_V2;
    	}
    	case 3: {
			return EXP_1_V3;
    	}
    	case 4: {
			return EXP_1_V4;
    	}
    	case 5: {
			return EXP_1_V5;
    	}
    	case 6: {
			return EXP_1_V6;
    	}
    	case 7: {
			return EXP_1_V7;
		}
    	case 8: {
    		return EXP_1_V8;
    	}
    	case 9: {
    		return EXP_1_V9;
    	}
    	case 10: {
			return EXP_1_V10;
    	}
    	case 11: {
			return EXP_1_V11;
    	}
    	case 12: {
			return EXP_1_V12;
    	}
    	case 13: {
			return EXP_1_V13;
    	}
    	case 14: {
			return EXP_1_V14;
    	}
    	case 15: {
			return EXP_1_V15;
    	}
    	case 16: {
			return EXP_1_V16;
    	}
    	case 17: {
			return EXP_1_V17;
    	}
    	case 18: {
			return EXP_1_V18;
		}
    	case 19: {
			return EXP_1_V19;
    	}
    	case 20: {
			return EXP_1_V20;
    	}
    	case 21: {
			return EXP_1_V21;
    	}
    	case 22: {
			return EXP_1_V22;
    	}
    	case 23: {
			return EXP_1_V23;
    	}
    	case 24: {
			return EXP_1_V24;
    	}
    	case 25: {
			return EXP_1_V25;
    	}
    	case 26: {
			return EXP_1_V26;
    	}
    	case 27: {
			return EXP_1_V27;
		}
    	case 28: {
			return EXP_1_V28;
    	}
    	case 29: {
			return EXP_1_V29;
    	}
    	case 30: {
			return EXP_1_V30;
    	}
    	case 31: {
			return EXP_1_V31;
    	}
    	case 32: {
			return EXP_1_V32;
    	}
    	default : {
    		fprintf(stderr, "EXP MASKS NOT DEFINED FOR NVAR %d\n",nvar);
    		return 0;
    	}
    }
}


static const degrees_t* getMaxExpArray (int nvar) {

	switch(nvar) {
		case 1: {
			return (MAX_EXP_LIST);
		}
		case 2: {
			return (MAX_EXP_LIST + 1);
		}
		case 3: {
			return (MAX_EXP_LIST + 3);
		}
		case 4: {
			return (MAX_EXP_LIST + 6);
		}
		case 5: {
			return (MAX_EXP_LIST + 10);
		}
		case 6: {
			return (MAX_EXP_LIST + 15);
		}
		case 7: {
			return (MAX_EXP_LIST + 21);
		}
		case 8: {
			return (MAX_EXP_LIST + 28);
		}
		case 9: {
			return (MAX_EXP_LIST + 36);
		}
		case 10: {
			return (MAX_EXP_LIST + 45);
		}
		case 11: {
			return (MAX_EXP_LIST + 55);
		}
		case 12: {
			return (MAX_EXP_LIST + 66);
		}
		case 13: {
			return (MAX_EXP_LIST + 78);
		}
		case 14: {
			return (MAX_EXP_LIST + 91);
		}
		case 15: {
			return (MAX_EXP_LIST + 105);
		}
		case 16: {
			return (MAX_EXP_LIST + 120);
		}
		case 17: {
			return (MAX_EXP_LIST + 136);
		}
		case 18: {
			return (MAX_EXP_LIST + 153);
		}
		case 19: {
			return (MAX_EXP_LIST + 171);
		}
		case 20: {
			return (MAX_EXP_LIST + 190);
		}
		case 21: {
			return (MAX_EXP_LIST + 210);
		}
		case 22: {
			return (MAX_EXP_LIST + 231);
		}
		case 23: {
			return (MAX_EXP_LIST + 253);
		}
		case 24: {
			return (MAX_EXP_LIST + 276);
		}
		case 25: {
			return (MAX_EXP_LIST + 300);
		}
		case 26: {
			return (MAX_EXP_LIST + 325);
		}
		case 27: {
			return (MAX_EXP_LIST + 351);
		}
		case 28: {
			return (MAX_EXP_LIST + 378);
		}
		case 29: {
			return (MAX_EXP_LIST + 406);
		}
		case 30: {
			return (MAX_EXP_LIST + 435);
		}
		case 31: {
			return (MAX_EXP_LIST + 465);
		}
		case 32: {
			return (MAX_EXP_LIST + 496);
		}
		default : {
			fprintf(stderr, "MAX EXPS NOT DEFINED FOR NVAR %d\n",nvar);
			return NULL;
		}
	}
}


#endif
