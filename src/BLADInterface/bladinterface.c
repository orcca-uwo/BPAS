

#include "BLADInterface/bladinterface.h"
#include <string.h>
#include "Utils/Unix_Timer.h"

/********************************************
 ** This file is littered with commented out
 ** debug statements and the likes...
 ** Might come in handy in future.
 ********************************************/

// static void printDegs_AAZ(FILE* fp, degrees_t degs, const char** syms, int nvar, degrees_t* masks, int* sizes) {
//         if (degs == 0) {
//                 return;
//         }
//         degree_t deg = GET_NTH_EXP(degs, masks[0], sizes[0]);
//         fprintf(fp, "%s^%d", syms[0], deg);
//         for (int k = 1; k < nvar; ++k) {
//                 deg = GET_NTH_EXP(degs, masks[k], sizes[k]);
//                 fprintf(fp, "*%s^%d", syms[k], deg);
//         }
// }

// static void printPoly_AAZ(FILE* fp, const AltArrZ_t* aa, const char** syms, int nvar) {
//         if (aa == NULL || aa->size == 0) {
//                 fprintf(fp,"0");
//                 return;
//         }

//         degrees_t* masks = getExpMaskArray(nvar);
//         int* sizes = getExpOffsetArray(nvar);

//         gmp_fprintf(fp, "%Zd", aa->elems[0].coef);
//         if (!isZeroExponentVector(aa->elems[0].degs)) {
//                 fprintf(fp, "*");
//                 printDegs_AAZ(fp, aa->elems[0].degs, syms, nvar, masks, sizes);
//         }
//         for (int i = 1; i < AA_SIZE(aa)-1; ++i) {
//                 if (mpz_sgn(aa->elems[i].coef) > 0) {
//                         fprintf(fp, " + ");
//                 } else {
//                         fprintf(fp, " ");
//                 }
//                 gmp_fprintf(fp, "%Zd", aa->elems[i].coef);
//                 if (!isZeroExponentVector(aa->elems[i].degs)) {
//                         fprintf(fp, "*");
//                         printDegs_AAZ(fp, aa->elems[i].degs, syms, nvar, masks, sizes);
//                 }
//         }
//         if (aa->size > 1) {
//                 if (mpz_sgn(aa->elems[aa->size-1].coef) > 0) {
//                         fprintf(fp, " + ");
//                 } else {
//                         fprintf(fp, " ");
//                 }
//                 gmp_fprintf(fp, "%Zd", aa->elems[aa->size-1].coef);
//                 if (!isZeroExponentVector(aa->elems[aa->size-1].degs)) {
//                         fprintf(fp, "*");
//                         printDegs_AAZ(fp, aa->elems[aa->size-1].degs, syms, nvar, masks, sizes);
//                 }
//         }

//         //fprintf(fp, "\n");

//         free(masks);
//         free(sizes);

// }

void pushOrderingBLAD(const char** syms, int nvar) {
	bav_Iordering r;

	char blocksCstr[8192];
	blocksCstr[0] = '\0';
	char* blocksPtr = &blocksCstr[0];
	char varCstr[8192];
	for (int i = 0; i < nvar-1; ++i) {
		// sprintf(varCstr, "x_%d,", i);
		sprintf(varCstr, "%s,", syms[i]);
		strcat(blocksPtr, varCstr);
		blocksPtr += strlen(varCstr);
	}
	// sprintf(varCstr, "x_%d", nvar-1);
	sprintf(varCstr, "%s", syms[nvar-1]);
	strcat(blocksPtr, varCstr);

	char orderCstr[8192];
	sprintf(orderCstr, "ordering (derivations= [], blocks=[%s])", blocksCstr);
	// fprintf(stderr, "\n\n\nPUSHING ORDERING\n%s\n", orderCstr);
	ba0_sscanf2(orderCstr, "%ordering", &r);
	bav_R_push_ordering(r);

}

blad_poly_z convertToBLAD_AAZ(const AltArrZ_t* aa, const char** syms) {

	int nvar = aa->nvar;

	blad_var* blad_vars = (blad_var*) malloc(sizeof(blad_var)*nvar);
	for (int i = 0; i < nvar; ++i) {
		ba0_sscanf2((char*) syms[i], "%v", &blad_vars[i]);
		// printf("var[%d]:", i);	
		// ba0_printf("%v\n", blad_vars[i]);bloblocksPtrcksPtr
	}

	struct bav_term T;
	struct bav_term resultTerm;
	bav_init_term(&T);
	bav_init_term(&resultTerm);
	for (int i = 0; i < nvar; ++i) {
		bav_mul_term_variable(&resultTerm, &T, blad_vars[i], 1);
		bav_set_term(&T, &resultTerm);
	}
	// ba0_printf("Term: %term\n\n", &T);

	enum bap_typeof_total_rank creatorRankType = bap_exact_total_rank;

	degree_t maxDegs[nvar];
	partialDegrees_AAZ(aa, maxDegs);
	for (int i = 0; i < nvar; ++i) {
		// fprintf(stderr, "maxDegs[%d]: %ld\n",i, GET_NTH_EXP(maxDegs, masks[i], sizes[i]));
		T.rg[i].deg = maxDegs[i];
		if (T.rg[i].deg == 0) {
			T.rg[i].deg = 1;
			creatorRankType = bap_approx_total_rank;
		}
	}
	
	blad_poly_z A = bap_new_polynom_mpz();
	struct bap_creator_mpz crea;
	// ba0_printf("Term before create: %term\n\n", &T);

	bap_begin_creator_mpz(&crea, A, &T, creatorRankType, aa->size);
	// ba0_printf("creators clot rang totall: %term\n\n", &crea.crea.iter.clot->tgest.rang_total);
	// fprintf(stderr, "got bap greator!\n");
	// ba0_printf("A rang total after creator: %term\n\n", &A->rang_total);

	AAZElem_t* elems = aa->elems;
	degree_t curDegs[aa->nvar];
	for (int i = 0; i < aa->size; ++i) {
		partialDegreesTerm_AAZ(aa, i, curDegs);
		for (int k = 0; k < nvar; ++k) {
			// fprintf(stderr, "curDegs[%d]: %ld\n",k, GET_NTH_EXP(elems[i].degs, masks[k], sizes[k]));
			T.rg[k].deg = curDegs[k];
		}
		// ba0_printf("Term before write: %term\n\n", &T);
		// ba0_printf("A rang total before write: %term\n\n", &A->rang_total);
		// ba0_printf("creators clot before write: %term\n\n", &crea.crea.iter.clot->tgest.rang_total);

		// int isFactor = bav_is_factor_term (&crea.crea.iter.clot->tgest.rang_total, &T, NULL);
		// fprintf(stderr, "is factor: %d\n", isFactor);
		
		// gmp_fprintf(stderr, "writing coef[%d]: %Zd\n", i, elems[i].coef);
		bap_write_creator_mpz(&crea, &T, elems[i].coef);

		// fprintf(stderr, "wrote factor!\n");
		// ba0_printf("current A: %Az\n\n", A);
	}
	bap_close_creator_mpz(&crea);
	// fprintf(stderr, "closing factor cretor\n");
	// ba0_printf("current A: %Az\n\n", A);

	return A;
}


// typedef struct bav_block {
// bav_subranking subr;
// struct ba0_tableof_int_p indices;
// struct ba0_tableof_string strs;
// } *bav_block;

// typedef struct bav_ordering {
// struct bav_tableof_symbol ders;
// struct bav_tableof_block blocks;
// struct bav_block operator_block;
// struct bav_tableof_variable varmax;
// } *bav_ordering;


// enum bav_typeof_symbol {
// bav_independent_symbol, bav_dependent_symbol,
// bav_operator_symbol, bav_temporary_symbol };

AltArrZ_t* convertFromBLAD_AAZ(const blad_poly_z A) {
	struct bap_itermon_mpz iter;
	bap_begin_itermon_mpz(&iter, A);

	long long Asize = bap_nbmon_polynom_mpz(A);

	mpz_t* coefPtr;
	struct bav_term T;
	bav_init_term(&T);
	blad_term degsTerm = &(A->rang_total);
	int nvar = degsTerm->size;

	AltArrZ_t* aa = makePolynomial_AAZ(Asize, T.size);
	AAZElem_t* elems = aa->elems;
	int i = 0;
	degrees_t degs;
	int* sizes = getExpOffsetArray(nvar);

	degree_t degsList[nvar];
	while (!bap_outof_itermon_mpz(&iter)) {
		bap_term_itermon_mpz(&T, &iter);
		coefPtr = bap_coeff_itermon_mpz(&iter);
		mpz_init(elems[i].coef);
		mpz_set(elems[i].coef, *coefPtr);

		for (int k = 0; k < nvar; ++k) {
			degsList[k] = 0;
		}	
		for (int l = 0; l < T.size; ++l) {
			for (int k = 0; k < nvar; ++k) {
				if (T.rg[l].var == degsTerm->rg[k].var) {
					degsList[k] = T.rg[l].deg;
					// break;
				}
			}
		}
		setDegrees_AAZ_inp(aa, i, degsList, nvar);

		bap_next_itermon_mpz(&iter);
		++i;
	}

	bap_close_itermon_mpz(&iter);

	aa->size = i;
	return aa;
}

AltArrZ_t* convertFromBLADWithVars_AAZ(const blad_poly_z A, const char** syms, int nvar) {
	struct bap_itermon_mpz iter;
	bap_begin_itermon_mpz(&iter, A);

	long long Asize = bap_nbmon_polynom_mpz(A);

	mpz_t* coefPtr;
	struct bav_term T;
	bav_init_term(&T);
	blad_term degsTerm = &(A->rang_total);
	int bladNvar = degsTerm->size;

	AltArrZ_t* aa = makePolynomial_AAZ(Asize, T.size);
	AAZElem_t* elems = aa->elems;
	aa->nvar = nvar;
	int i = 0;
	degrees_t degs;
	int* sizes = getExpOffsetArray(nvar);

	degree_t degsList[nvar];
	while (!bap_outof_itermon_mpz(&iter)) {
		bap_term_itermon_mpz(&T, &iter);
		coefPtr = bap_coeff_itermon_mpz(&iter);
		mpz_init(elems[i].coef);
		mpz_set(elems[i].coef, *coefPtr);

		for (int k = 0; k < nvar; ++k) {
			degsList[k] = 0;
		}
		for (int l = 0; l < T.size; ++l) {
			for (int k = 0; k < nvar; ++k) {
				if (strcmp(T.rg[l].var->root->ident, syms[k]) == 0) {
				// if (T.rg[l].var == degsTerm->rg[k].var) {
					degsList[k] = T.rg[l].deg;
					// break;
				}
			}
		}
		setDegrees_AAZ_inp(aa, i, degsList, nvar);

		bap_next_itermon_mpz(&iter);
		++i;
	}

	bap_close_itermon_mpz(&iter);

	aa->size = i;
	return aa;
}

void factorPolynomialBLAD_AAZ(const AltArrZ_t* aa, const char** syms, int* numFacts, AltArrZ_t*** retFacts, int** retExps, mpz_t numericFact) {

	initBLAD();
	pushOrderingBLAD(syms, aa->nvar);

	// fprintf(stderr, "FACTORING POLY:\n" );
	// printPoly_AAZ(stderr, aa, syms, aa->nvar);
	// fprintf(stderr, "\n");
	blad_poly_z A = convertToBLAD_AAZ(aa, syms);
	blad_product_z facts = bap_new_product_mpz();
	// ba0_fprintf(stderr, "got blad poly: %Az\n\n", A);
	// fprintf(stderr, "calling factor polynom\n");
	bap_factor_polynom_mpz(facts, A);
	// fprintf(stderr, "factored polynom\n");

//	mpz_t nnn;
//	mpz_init(nnn);
//	bap_signed_numeric_content_polynom_mpz (nnn, A);
//	gmp_fprintf(stderr,"%Zd\n",nnn);
	
	AltArrZ_t** facts_AAZ = (AltArrZ_t**) malloc(sizeof(AltArrZ_t)*facts->size);
	int* exps_AAZ = (int*) malloc(sizeof(int)*facts->size);
    for (int i = 0; i < facts->size; ++i) {
        facts_AAZ[i] = convertFromBLADWithVars_AAZ(&facts->tab[i].factor, syms, aa->nvar);
        exps_AAZ[i] = facts->tab[i].exponent;
    }

    if (numFacts != NULL) {
    	*numFacts = facts->size;
    }
    if (retExps != NULL) {
    	*retExps = exps_AAZ;
    } else {
    	free(exps_AAZ);
    }
    if (retFacts != NULL) {
    	*retFacts = facts_AAZ;
    } else {
    	for (int i = 0; i < facts->size; ++i) {
    		freePolynomial_AAZ(facts_AAZ[i]);
    	}
    	free(facts_AAZ);
    }
    mpz_set(numericFact, facts->num_factor);



    popOrderingBLAD();
    freeBLAD();

}

void gcdBLAD_AAZ(const AltArrZ_t* a, const AltArrZ_t* b, const char** syms, AltArrZ_t** retG) {
	if (retG == NULL) {
		return;
	}

	// timer_id id = start_timer();

	initBLAD();
	pushOrderingBLAD(syms, a->nvar);
	// timer_time elapsed1 = elapsed_time(&id);
 //    double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Init blad time: %f\n", time );



	// id = start_timer();
	blad_poly_z A = convertToBLAD_AAZ(a, syms);
	blad_poly_z B = convertToBLAD_AAZ(b, syms);
	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Convert to blad time: %f\n", time );


	// id = start_timer();
	blad_poly_z G = bap_new_polynom_mpz();
	bap_gcd_polynom_mpz(G, NULL, NULL, A, B);
	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "gcd time: %f\n", time );



	// id = start_timer();
	*retG = convertFromBLADWithVars_AAZ(G, syms, a->nvar);

    popOrderingBLAD();
    freeBLAD();
	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Convert from blad time: %f\n", time );

}