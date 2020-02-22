


#include "../../include/RationalNumberPolynomial/SMQP_Support_Unpacked.h"
#include "../../include/RationalNumberPolynomial/SMQP_Support_Recursive_Unpacked.h"
#include "../../include/RationalNumberPolynomial/SMQP_Support_Recursive-AA.h"
#include "../../include/RationalNumberPolynomial/SMQP_Support_Test-AA.h"

long nvar = 3;
long numTerms = 10;
int coefBound = 5ul;
degree_t sparsity = 10;
int includeNeg = 1;

char* syms[] = {"x", "y", "z", "s", "t", "u", "v", "w", "p", "q", "r", "f", "g", "h", "j", "k", "l", "m", "n", "o"};

int recArrayMonomialsEqual(RecArr_t* a, RecArr_t* b) {
    if (a->size != b->size) {
        return 0;
    }

    int size = a->size;
    for (int i = 0 ; i < size; ++i) {
        if(a->elems[i].coefSize != b->elems[i].coefSize ||
            a->elems[i].exp != b->elems[i].exp) {
            return 0;
        } 
    }

    return 1;
}


int main(int argc, char *argv[]) {
    if (argc > 1 && atol(argv[1]) >= 0) {
        nvar = atol(argv[1]);
    }
    if (argc > 2 && atol(argv[2]) > 0) { 
        numTerms = atol(argv[2]);
    }
    if (argc > 3 && atoi(argv[3]) > 0) {
        coefBound = atoi(argv[3]);
    }
    if (argc > 4 && atol(argv[4]) > 1) {
        sparsity = atol(argv[4]);
    }
    if (argc > 5 && atoi(argv[5]) >= 0) {
        includeNeg = atoi(argv[5]);
    }

    if (nvar == 0) {
        numTerms = 1;
    }

    Node* node = buildRandomPoly(nvar, numTerms, coefBound, sparsity, includeNeg);

    fprintf(stderr, "Testing deepCopy, equal, etc. ================================\n");

    AltArr_t* aa = deepCopyPolynomial_AAFromNode(node, nvar);
    // printPoly_AA(stderr, aa, syms, aa->nvar);
    unpackExponentVectors_AA_inp(aa);
    // fprintf(stderr, "\n");
    // printPoly_AA_unpk(stderr, aa, syms, aa->nvar);
    // fprintf(stderr, "\n");

    AltArr_t* ba = deepCopyPolynomial_AAFromNode(node, nvar);
    int equal1 = isExactlyEqual_AA(aa, ba);
    // unpackExponentVectors_AA_inp(ba);
    int equal2 = isExactlyEqual_AA(aa, ba);

    fprintf(stderr, "equal1: %d, equal2: %d\n", equal1, equal2);

    AltArr_t* packed = deepCopyPolynomial_AAFromNode(node, nvar);
    AltArr_t* ca = deepCopyPolynomial_AA(ba);
    AltArr_t* da = deepCopyPolynomial_AA_unpk(packed);

    int equal3 = isExactlyEqual_AA(packed, ca);
    int equal4 = isExactlyEqual_AA(ca, da);
    fprintf(stderr, "equal3: %d, equal4: %d\n", equal3, equal4);



    fprintf(stderr, "\n\nTesting unpack, merge sort, etc. ================================\n" );

	//reverse the node to test sort
    Node* prev = NULL, *cur = node, *next;
    while(cur != NULL) {
    	next = cur->next;
    	cur->next = prev;
    	prev = cur;
    	cur = next;
    }
    node = prev;

    AltArr_t* unsorted = deepCopyPolynomial_AAFromNode(node, nvar);
    unpackExponentVectors_AA_inp(unsorted);
    // printPoly_AA(stderr, unsorted, syms, unsorted->nvar);
    // fprintf(stderr, "\n" );
    mergeSortPolynomial_AA(unsorted);
    // printPoly_AA(stderr, unsorted, syms, unsorted->nvar);
    // fprintf(stderr, "\n" );
    
		
    fprintf(stderr, "\n\nTesting eval poly etc. ================================\n");

	mpq_t vals[nvar];
	for (int k = 0; k < nvar; ++k) {
		mpq_init(vals[k]);
		mpz_set_si(mpq_numref(vals[k]), rand() % (1 << coefBound));
		while (mpz_cmp_si(mpq_denref(vals[k]), 0l) == 0) {
			mpz_set_si(mpq_denref(vals[k]), rand() % (1 << coefBound));
		}
	}

	mpq_t res;
	mpq_init(res);
	mpq_t res_unpk;
	mpq_init(res_unpk);

	// fprintf(stderr, "packed->unpacked: %d\n", packed->unpacked);
	evalPolyToVal_AA(packed, vals, nvar, res);
	evalPolyToVal_AA_unpk(da, vals, nvar, res_unpk);

	gmp_fprintf(stderr, "%Qd == %Qd ?? %d\n", res, res_unpk, mpq_cmp(res, res_unpk) == 0);

	if (nvar > 2) {
		int active[nvar]; 
		for (int i = 0; i < nvar; ++i) {
			active[i] = 1;
		}
		active[2] = 0;
		AltArr_t* evalPoly1 = evaluatePoly_AA(packed, active, vals, nvar);
		AltArr_t* evalPoly2 = evaluatePoly_AA_unpk(da, active, vals, nvar);
		int equal5 = isExactlyEqual_AA(evalPoly1, evalPoly2);
		// if (evalPoly2->unpacked) { 
		// 	printPoly_AA_unpk(stderr, evalPoly2, syms, evalPoly2->nvar);
		// } else {
		// 	printPoly_AA(stderr, evalPoly2, syms, evalPoly2->nvar);
		// }
		fprintf(stderr, "\neval poly test: %d\n", equal5);

	}

    fprintf(stderr, "\n\nTesting add poly etc. ================================\n");

    // fprintf(stderr, "da->unpacked: %d, ca->unpacked: %d\n", da->unpacked, ca->unpacked);
    AltArr_t* sum_unpk = addPolynomials_AA_unpk(da, da, nvar);
    AltArr_t* sum_pk = addPolynomials_AA(ca, ca, nvar);
    // fprintf(stderr, "sum:\n");
    // printPoly_AA(stderr, sum_pk, syms, nvar);
    int equal6 = isExactlyEqual_AA(sum_unpk, sum_pk);
    fprintf(stderr, "\nsummations equal: %d\n", equal6);

    fprintf(stderr, "\n\nTesting sub poly etc. ================================\n");

    // fprintf(stderr, "da->unpacked: %d, ca->unpacked: %d\n", da->unpacked, ca->unpacked);
    sum_unpk = subPolynomials_AA_unpk(da, da, nvar);
    sum_pk = subPolynomials_AA(ca, ca, nvar);
    // fprintf(stderr, "difference:\n");
    // printPoly_AA(stderr, sum_pk, syms, nvar);





// TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//equal 7 not equal!!!!






    equal6 = isExactlyEqual_AA(sum_unpk, sum_pk);
    fprintf(stderr, "\ndifferences equal: %d\n", equal6);

    Node* node2 = buildRandomPoly(nvar, numTerms, coefBound, sparsity, includeNeg);
    AltArr_t* packed2 = deepCopyPolynomial_AAFromNode(node2, nvar);
    AltArr_t* diff1 = subPolynomials_AA(ca, packed2, nvar);
    AltArr_t* diff2 = subPolynomials_AA(da, packed2, nvar);

    fprintf(stderr, "ca:\n");
    printPoly_AA(stderr, ca, syms, nvar);
    fprintf(stderr, "\n");

    fprintf(stderr, "da:\n");
    printPoly_AA(stderr, da, syms, nvar);
    fprintf(stderr, "\n");

    fprintf(stderr, "packed2:\n");
    printPoly_AA(stderr, packed2, syms, nvar);
    fprintf(stderr, "\n");

    fprintf(stderr, "diff1:\n");
    printPoly_AA(stderr, diff1, syms, nvar);
    fprintf(stderr, "\n");

    fprintf(stderr, "diff2:\n");
    printPoly_AA(stderr, diff2, syms, nvar);
    fprintf(stderr, "\n");

	int equal7 = isExactlyEqual_AA(diff1, diff2);
    fprintf(stderr, "\ndifferences equal: %d\n", equal7);
    if (!equal7) {
        exit(1);
    }


    fprintf(stderr, "\n\nTesting add poly inplace etc. ================================\n");

    AltArr_t* da2 = deepCopyPolynomial_AA_unpk(da);
	da2 = addPolynomials_AA_inp_unpk(da2, da, nvar);
	AltArr_t* ca2 = deepCopyPolynomial_AA(ca);
	ca2 = addPolynomials_AA_inp(ca2, ca, nvar);
	int equal8 = isExactlyEqual_AA(da2, ca2);
    fprintf(stderr, "\n summations in place equal: %d\n", equal8);


    fprintf(stderr, "\n\nTesting sub poly inplace etc. ================================\n");
	// printPoly_AA_unpk(stderr, da2, syms, nvar);
	// fprintf(stderr, "\n" );

	da2 = subPolynomials_AA_inp_unpk(da2, da, nvar);
	int equal9 = isExactlyEqual_AA(da2, da);
	// printPoly_AA_unpk(stderr, da, syms, nvar);
    fprintf(stderr, "\n differences in place equal: %d\n", equal9);

    unpackExponentVectors_AA_inp(da2);
    da2 = subPolynomials_AA_inp_unpk(da2, packed2, nvar);
    AltArr_t* properDiff = subPolynomials_AA(ca, packed2, nvar);
    int equal10 = isExactlyEqual_AA(da2, properDiff);
    fprintf(stderr, "\n differences in place equal: %d\n", equal10);


    fprintf(stderr, "\n\nTesting multiply poly etc. ================================\n");


    // fprintf(stderr, "packed:\n" );
    // printPoly_AA(stderr, packed2, syms, nvar);
    // fprintf(stderr, "\nda:\n");
    // printPoly_AA(stderr, da, syms, nvar);
    AltArr_t* prod_unpk = multiplyPolynomials_AA_unpk(da, packed2, nvar);

    // fprintf(stderr, "\npacked:\n" );
    // printPoly_AA(stderr, packed2, syms, nvar);
    // fprintf(stderr, "\nca:\n");
    // printPoly_AA(stderr, ca, syms, nvar);
    AltArr_t* prod_pk = multiplyPolynomials_AA(ca, packed2, nvar);

    // fprintf(stderr, "\n\nprod_unpk:" );
    // printPoly_AA_unpk(stderr, prod_unpk, syms, nvar);
    // fprintf(stderr, "\n\nprod_pk:");
    // printPoly_AA(stderr, prod_pk, syms, nvar);
    int equal11 = isExactlyEqual_AA(prod_unpk, prod_pk);
    fprintf(stderr, "\n products equal: %d\n", equal11);


    fprintf(stderr, "\n\nTesting exponentiate poly etc. ================================\n");

    // int testNvar = 10;
    // AltArr_t* testExp = makePolynomial_AA(1, testNvar);
    // int* sizes = getExpOffsetArray(testNvar);
    // testExp->elems[0].degs = 0;
   	// for (int i = 0; i < testNvar; ++i) {
   	// 	testExp->elems[0].degs |= 1ll << sizes[i];
   	// }
   	// mpq_init(testExp->elems->coef);
   	// mpq_set_si(testExp->elems->coef, 1ll, 1ll);
   	// testExp->size = 1;

   	// testExp = exponentiatePoly_AA(testExp, 12200, testNvar);

   	// fprintf(stderr, "exp poly unpacked?: %d\n", testExp->unpacked);
   	// printPoly_AA(stderr, testExp, syms, testNvar);
   	// fprintf(stderr, "\n");





    fprintf(stderr, "\n\nTesting divide poly etc. ================================\n");

    AltArr_t* singleTerm = makePolynomial_AA(1, nvar);
    mpq_init(singleTerm->elems[0].coef);
    mpq_set_si(singleTerm->elems[0].coef, 5l, 1l); 
    singleTerm->elems[0].degs = numTerms > 5 ? prod_pk->elems[4].degs : prod_pk->elems[0].degs;
    singleTerm->size = 1;
    unpackExponentVectors_AA_inp(prod_pk);

    // fprintf(stderr, "prod_pk dividend:\n");
    // printPoly_AA(stderr, prod_pk, syms, nvar);
    // fprintf(stderr, "\n");

    // fprintf(stderr, "single term divisor:\n");
    // printPoly_AA(stderr, singleTerm, syms, nvar);
    // fprintf(stderr, "\n");

    AltArr_t* quo, *rem;
    // divideBySingleTerm_AA_unpk(prod_pk, singleTerm, &quo, &rem, nvar);
	// fprintf(stderr, "single term quo:\n");
 //    printPoly_AA(stderr, quo, syms, nvar);
 //    fprintf(stderr, "\n");

 //    fprintf(stderr, "single term rem:\n");
 //    printPoly_AA(stderr, rem, syms, nvar);
 //    fprintf(stderr, "\n");


    //prod_pk = ca * packed2
    // fprintf(stderr, "prod_pk:\n");
    // printPoly_AA(stderr, prod_pk, syms, nvar);
    // fprintf(stderr, "\n");

    // fprintf(stderr, "ca:\n");
    // printPoly_AA(stderr, ca, syms, nvar);
    // fprintf(stderr, "\n");

    // dividePolynomials_AA_unpk(prod_pk, ca, &quo, &rem, nvar);

	// fprintf(stderr, "quo:\n");
	// printPoly_AA(stderr, quo, syms, nvar);
	// fprintf(stderr, "\n");

	// fprintf(stderr, "packed2:\n");
	// printPoly_AA(stderr, packed2, syms, nvar);
	// fprintf(stderr, "\n");

    // int equal12 = isExactlyEqual_AA(quo, packed2);
    // int equal13 = (rem == NULL || rem->size == 0);

    // fprintf(stderr, "divide polynomails constructive quo correct: %d\n", equal12);
    // fprintf(stderr, "divide polynomails constructive rem correct: %d\n", equal13);

    Node* node4 = buildRandomPoly(nvar, numTerms, coefBound, sparsity, includeNeg);
    AltArr_t* randDividend = deepCopyPolynomial_AAFromNode(node4, nvar);
    unpackExponentVectors_AA_inp(randDividend);

    Node* node3 = buildRandomPoly(nvar, (numTerms < 4 ? numTerms : numTerms - 3), coefBound, sparsity, includeNeg);
    AltArr_t* randDivisor = deepCopyPolynomial_AAFromNode(node3, nvar);

    dividePolynomials_AA_unpk(prod_pk, randDivisor, &quo, &rem, nvar);

    packExponentVectors_AA_inp(prod_pk);
	AltArr_t* quo_pk, * rem_pk;
    dividePolynomials_AA(prod_pk, randDivisor, &quo_pk, &rem_pk, nvar);

    int equal14 = isExactlyEqual_AA(quo, quo_pk);
    int equal15 = isExactlyEqual_AA(rem, rem_pk);

	// fprintf(stderr, "rem_unpk:=\n");
	// printPoly_AA(stderr, rem, syms, nvar); `
	// fprintf(stderr, ";\n");
	
	// fprintf(stderr, "rem_pk:=\n");
	// printPoly_AA(stderr, rem_pk, syms, nvar);
	// fprintf(stderr, ";\n");

    

	fprintf(stderr, "divide polynomails random quo correct: %d\n", equal14);
    fprintf(stderr, "divide polynomails random rem correct: %d\n", equal15);
    if (!equal14 || !equal15) {
    	exit(1);
    }

    //force unpack of prod_unpk;
    unpackExponentVectors_AA_inp(prod_unpk);

    exactDividePolynomials_AA_unpk(prod_unpk, da, &quo, nvar);
    int equal16 = isExactlyEqual_AA(quo, packed2);

    fprintf(stderr, "exact divide polynomails quo correct: %d\n", equal16);

    AltArr_t* divTestQuo = NULL;
    int divisible = divideTest_AA_unpk(prod_unpk, da, &divTestQuo, nvar);
    int equal17 = isExactlyEqual_AA(quo, divTestQuo);
    fprintf(stderr, "\ndivide test polynomials correct: %d\n", divisible && equal17 );


    fprintf(stderr, "\n\nTesting derivative etc. ================================\n");

    AltArr_t* div0_unpk = derivative_AA_unpk(prod_unpk, 0, 1);
    AltArr_t* div1_unpk = derivative_AA_unpk(prod_unpk, nvar < 2 ? 0 : 1, 1);

    packExponentVectors_AA_inp(prod_unpk);
    AltArr_t* div0_pk = derivative_AA(prod_unpk, 0, 1);
    AltArr_t* div1_pk = derivative_AA(prod_unpk, nvar < 2 ? 0 : 1, 1);

    int equal18 = isExactlyEqual_AA(div0_pk, div0_unpk);
    int equal19 = isExactlyEqual_AA(div1_pk, div1_unpk);

    fprintf(stderr, "derivative0 correct: %d\n", equal18);
    fprintf(stderr, "derivative1 correct: %d\n", equal19);


    fprintf(stderr, "\n\nTesting integral etc. ================================\n");

    unpackExponentVectors_AA_inp(prod_unpk);
    AltArr_t* int0_unpk = integral_AA_unpk(prod_unpk, -1, 1);
    AltArr_t* int1_unpk = integral_AA_unpk(prod_unpk, 0, 1);

    packExponentVectors_AA_inp(prod_unpk);
    AltArr_t* int0_pk = integral_AA(prod_unpk, -1, 1);
    AltArr_t* int1_pk = integral_AA(prod_unpk, 0, 1);

    int equal20 = isExactlyEqual_AA(int0_pk, int0_unpk);
    int equal21 = isExactlyEqual_AA(int1_pk, int1_unpk);

    fprintf(stderr, "integral0 correct: %d\n", equal20);
    fprintf(stderr, "integral1 correct: %d\n", equal21);


    fprintf(stderr, "\n\nTesting common factor etc. ================================\n");

    unpackExponentVectors_AA_inp(prod_unpk);
    AltArr_t* fact_unpk = NULL;
    AltArr_t* comFact_unpk = commonFactor_AA_unpk(prod_unpk, &fact_unpk);
    packExponentVectors_AA_inp(prod_unpk);
    AltArr_t* fact_pk = NULL;
    AltArr_t* comFact_pk = commonFactor_AA(prod_unpk, &fact_pk);

    int equal22 = isExactlyEqual_AA(fact_unpk, fact_pk);
    int equal23 = isExactlyEqual_AA(comFact_pk, comFact_unpk);
    fprintf(stderr, "common factors equal: %d\n", equal23);
    fprintf(stderr, "factored polys equal: %d\n", equal22);
    fprintf(stderr, "non-trivial common factor: %d\n", !isZeroExponentVector(comFact_pk->elems->degs) );


    fprintf(stderr, "\n\nTesting variable reorderings etc. ================================\n");

    AltArr_t* v1 = deepCopyPolynomial_AA_unpk(prod_unpk);
    AltArr_t* v2 = deepCopyPolynomial_AA(prod_unpk);
    unpackExponentVectors_AA_inp(v1);

    expandNumVars_AA_unpk(v1, nvar + 2);
    expandNumVars_AA(v2, nvar + 2);
    int equal24 = isExactlyEqual_AA(v1, v2);

    expandNumVarsLeft_AA_unpk(v1, nvar + 4);
    expandNumVarsLeft_AA(v2, nvar + 4);
    int equal25 = isExactlyEqual_AA(v1, v2);

    fprintf(stderr, "expand variables equal: %d\n", equal24);
    fprintf(stderr, "expand variables left equal: %d\n", equal25);
    
    // fprintf(stderr, "before shrink:=\n");
    // printPoly_AA(stderr, v1, syms, nvar+4);
    // fprintf(stderr, ";\n");
    shrinkNumVarsAtIdx_AA_unpk(v1, 2);
    // fprintf(stderr, "after shrink:=\n");
    // printPoly_AA(stderr, v1, syms+1, nvar+4-1);
    // fprintf(stderr, ";\n");

    shrinkNumVarsAtIdx_AA(v2, 2);
    int equal26 = isExactlyEqual_AA(v1, v2);
    fprintf(stderr, "shrink vars at idx equal: %d\n", equal26);

    int newNvar = nvar+3;
    int varMap[newNvar];
    for (int j = 0; j < newNvar; ++j) {
        varMap[j] = 0;
    }
    varMap[1] = 3;
    varMap[3] = 2;
    varMap[2] = 1;

    AltArr_t* v11 = deepCopyPolynomial_AA(v1);
    AltArr_t* v22 = deepCopyPolynomial_AA(v2);

    unpackExponentVectors_AA_inp(v1);
    reorderVars_AA_unpk(v1, varMap, newNvar);
    reorderVars_AA_unpk(v2, varMap, newNvar);
    int equal27 = isExactlyEqual_AA(v1, v2);
    fprintf(stderr, "reoerder vars equal: %d\n", equal27);

    varMap[0] = -1;
    varMap[4] = -1;

    unpackExponentVectors_AA_inp(v11);
    shrinkAndReorderVars_AA_unpk(v11, varMap, newNvar);
    shrinkAndReorderVars_AA(v22, varMap, newNvar);
    int equal28 = isExactlyEqual_AA(v11, v22);

    fprintf(stderr, "shrink and reoerder vars equal: %d\n", equal28);




    fprintf(stderr, "\n\nTesting main shift etc. ================================\n");

    AltArr_t* toShift = deepCopyPolynomial_AA(prod_unpk);
    unpackExponentVectors_AA_inp(toShift);
    AltArr_t* shifted_unpk = mainLShiftPolynomial_AA_unpk(toShift, 5);
    packExponentVectors_AA_inp(toShift);
    AltArr_t* shifted_pk = mainLShiftPolynomial_AA(toShift,5);
    unpackExponentVectors_AA_inp(toShift);
    mainLShiftPolynomial_AA_inp_unpk(toShift, 5);

    int equal29 = isExactlyEqual_AA(shifted_unpk, toShift);
    int equal30 = isExactlyEqual_AA(shifted_unpk, shifted_pk);

    fprintf(stderr, "shifted and shifted_inp equal: %d\n", equal29 && equal30);


    fprintf(stderr, "\n\nTesting leading variable etc. ================================\n");

    unpackExponentVectors_AA_inp(prod_unpk);
    AltArr_t* lt_unpk = leadingTerm_AA_unpk(prod_unpk, prod_unpk->nvar);
    packExponentVectors_AA_inp(prod_unpk);
    AltArr_t* lt_pk = leadingTerm_AA(prod_unpk, prod_unpk->nvar);
    int equal31 = isExactlyEqual_AA(lt_unpk, lt_pk);

    fprintf(stderr, "leading terms equal: %d\n", equal31);

    unpackExponentVectors_AA_inp(prod_unpk);
    int lv_unpk = leadingVariable_AA_unpk(prod_unpk);
    degree_t ld_unpk = mainLeadingDegree_AA_unpk(prod_unpk);
    AltArr_t* mlc_unpk = mainLeadingCoefficient_AA_unpk(prod_unpk);
    AltArr_t* coefAt1_unpk = mainCoefficientAtIdx_AA_unpk(prod_unpk, 2);

    packExponentVectors_AA_inp(prod_unpk);
    int lv_pk = leadingVariable_AA(prod_unpk);
    degree_t ld_pk = mainLeadingDegree_AA(prod_unpk);
    AltArr_t* mlc_pk = mainLeadingCoefficient_AA(prod_unpk);
    AltArr_t* coefAt1_pk = mainCoefficientAtIdx_AA(prod_unpk, 2);

    int equal32 = isExactlyEqual_AA(mlc_unpk, mlc_pk);
    int equal33 = isExactlyEqual_AA(coefAt1_unpk, coefAt1_pk);
    fprintf(stderr, "leading variables equal: %d\n", lv_unpk == lv_pk);
    fprintf(stderr, "main leading degree equal: %d\n", ld_unpk == ld_pk);
    fprintf(stderr, "main leading coefficient equal: %d\n", equal32);
    fprintf(stderr, "main leading coefficient at idx equal: %d\n", equal33);


    fprintf(stderr, "\n\nTesting max polynomials etc. ================================\n");
    AltArr_t* max1 = maxPolynomials_AA(packed, packed2);
    unpackExponentVectors_AA_inp(packed);
    AltArr_t* max2 = maxPolynomials_AA(packed, packed2);
    unpackExponentVectors_AA_inp(packed2);
    AltArr_t* max3 = maxPolynomials_AA(packed, packed2);
    AltArr_t* max4 = maxPolynomials_AA_inp_unpk(packed, packed2);
    packExponentVectors_AA_inp(packed);
    packExponentVectors_AA_inp(packed2);
    AltArr_t* max5 = maxPolynomials_AA_inp(packed, packed2);

    int equal34 = isExactlyEqual_AA(max1, max2) && isExactlyEqual_AA(max2, max3);
    int equal35 = isExactlyEqual_AA(max4, max5);
    fprintf(stderr, "maxPolynomials equal: %d\n", equal34);
    fprintf(stderr, "maxPolynomials inp equal: %d\n", equal35);





    fprintf(stderr, "\n\nTesting recursive conversion etc. ================================\n");
    prod_pk = deepCopyPolynomial_AA(prod_unpk);
    unpackExponentVectors_AA_inp(prod_unpk);
    RecArr_t* recPoly_unpk = convertToRecursiveArray(prod_unpk);
    packExponentVectors_AA_inp(prod_pk);
    RecArr_t* recPoly_pk = convertToRecursiveArray(prod_pk);

    RecArr_t* recCopy_unpk = deepCopyRecArrayPolynomial(recPoly_unpk, prod_pk->nvar);
    AltArr_t* prodCopy_unpk = convertFromRecursiveArray(recCopy_unpk, prod_unpk->nvar);

    RecArr_t* recCopy_pk = deepCopyRecArrayPolynomial(recPoly_pk, prod_pk->nvar);
    unpackExponentVectors_RecArray_inp(recCopy_pk, prod_pk->nvar);
    AltArr_t* prodCopy_pk = convertFromRecursiveArray(recCopy_pk, prod_pk->nvar);

    int equal36 = recArrayMonomialsEqual(recPoly_unpk, recPoly_pk);

    prod_unpk = convertFromRecursiveArray(recPoly_unpk, prod_unpk->nvar);
    prod_pk = convertFromRecursiveArray(recPoly_pk, prod_pk->nvar);
    int equal37 = isExactlyEqual_AA(prod_unpk, prod_pk);

    fprintf(stderr, "convert in rec equal: %d\n", equal36);
    fprintf(stderr, "convert out rec equal: %d\n", equal37);


    int equal38 = isExactlyEqual_AA(prod_unpk, prodCopy_unpk);

    fprintf(stderr, "deepCopyRecArrayPolynomial equal: %d\n", equal38);

    int equal39 = isExactlyEqual_AA(prodCopy_unpk, prodCopy_pk);

    fprintf(stderr, "unpackExponentVectors_RecArray_inp equal: %d\n", equal39);


    fprintf(stderr, "\n\nTesting psuedodivide etc. ================================\n");

    packExponentVectors_AA_inp(prod_pk);
    // if you do this with PDIVIDE_DIVISBLE_CHECK you're gonna be waiting forever and a day
    // prod_pk->elems->degs |= (0xfffffll) << 28;

    unpackExponentVectors_AA_inp(prod_pk);
    unpackExponentVectors_AA_inp(packed2);

    recPoly_unpk = convertToRecursiveArray(prod_pk);
    RecArr_t* recPoly2_unpk = convertToRecursiveArray(packed2);

    AltArr_t* pquo, * prem;
    int e = 0;
    AltArr_t* hPow;
    int lazy = 0;
    pesudoDivide_RecArray(recPoly_unpk, recPoly2_unpk, &pquo, &prem, &e, &hPow, nvar, lazy);


    prod_pk = convertFromRecursiveArray(recPoly_unpk, nvar);
    packed2 = convertFromRecursiveArray(recPoly2_unpk, nvar);

    // fprintf(stderr, "c:= ");
    // printPoly_AA(stderr, prod_pk, syms, nvar);
    // fprintf(stderr, ";\n\nb:= ");
    // printPoly_AA(stderr, packed2, syms, nvar);
    // fprintf(stderr, ";\n\nq:= ");
    // printPoly_AA(stderr, pquo, syms, nvar);
    // fprintf(stderr, ";\n\nr:= ");
    // printPoly_AA(stderr, prem, syms, nvar);
    // fprintf(stderr, "\n");

    recPoly_unpk = convertToRecursiveArray(prod_pk);
    RecArr_t* recRandDiv = convertToRecursiveArray(randDivisor);
    pesudoDivide_RecArray(recPoly_unpk, recRandDiv, &pquo, &prem, &e, &hPow, nvar, lazy);

    prod_pk = convertFromRecursiveArray(recPoly_unpk, nvar);
    randDivisor = convertFromRecursiveArray(recRandDiv, nvar);

    AltArr_t* quo1 = pquo;
    AltArr_t* rem1 = prem;

    // fprintf(stderr, "c:= ");
    // printPoly_AA(stderr, prod_pk, syms, nvar);
    // fprintf(stderr, ";\n\nb:= ");
    // printPoly_AA(stderr, randDivisor, syms, nvar);
    // fprintf(stderr, ";\nq:= ");
    // printPoly_AA(stderr, pquo, syms, nvar);
    // fprintf(stderr, ";\n\nr:= ");
    // printPoly_AA(stderr, prem, syms, nvar);
    // fprintf(stderr, "\n\n");


    packExponentVectors_AA_inp(prod_pk);
    packExponentVectors_AA_inp(randDivisor);

    recPoly_unpk = convertToRecursiveArray(prod_pk);
    recRandDiv = convertToRecursiveArray(randDivisor);
    pesudoDivide_RecArray(recPoly_unpk, recRandDiv, &pquo, &prem, &e, &hPow, nvar, lazy);

    prod_pk = convertFromRecursiveArray(recPoly_unpk, nvar);
    randDivisor = convertFromRecursiveArray(recRandDiv, nvar);

    AltArr_t* quo2 = pquo;
    AltArr_t* rem2 = prem;

    int equal40 = isExactlyEqual_AA(quo1, quo2);
    int equal41 = isExactlyEqual_AA(rem1, rem2);    

    fprintf(stderr, "pseudodiv quo equal: %d\n", equal40 );
    fprintf(stderr, "pseudodiv rem equal: %d\n", equal41 );

    randDivisor->size = 1;
    recPoly_unpk = convertToRecursiveArray(prod_pk);
    recRandDiv = convertToRecursiveArray(randDivisor);
    pesudoDivide_RecArray(recPoly_unpk, recRandDiv, &pquo, &prem, &e, &hPow, nvar, lazy);

    prod_pk = convertFromRecursiveArray(recPoly_unpk, nvar);
    randDivisor = convertFromRecursiveArray(recRandDiv, nvar);

    // fprintf(stderr, "quo1:= ");
    // printPoly_AA(stderr, quo1, syms, nvar);
    // fprintf(stderr, ";\n\nquo2:= ");
    // printPoly_AA(stderr, quo2, syms, nvar);
    // fprintf(stderr, "c:= ");
    // printPoly_AA(stderr, prod_pk, syms, nvar);
    // fprintf(stderr, ";\n\nb:= ");
    // printPoly_AA(stderr, randDivisor, syms, nvar);
    // fprintf(stderr, ";\n\nq:= ");
    // printPoly_AA(stderr, pquo, syms, nvar);
    // fprintf(stderr, ";\n\nr:= ");
    // printPoly_AA(stderr, prem, syms, nvar);
    // fprintf(stderr, ";\n\nmapler:= prem(c,b,x,\'m\',\'mapleq\'):\nexpand(mapler - r);\nexpand(mapleq - q);\n");



    fprintf(stderr, "\n\nTesting variables, degrees, etc. etc. ================================\n");

    unpackExponentVectors_AA_inp(randDivisor);
    int vars_unpk[nvar];
    nonZeroVariables_AA_unpk(randDivisor, vars_unpk);
    int vars_pk[nvar];
    packExponentVectors_AA_inp(randDivisor);
    nonZeroVariables_AA_unpk(randDivisor, vars_pk);
    int equal42 = 1;
    for (int i = 0; i < nvar; ++i) {
        equal42 = equal42 && vars_pk[i] == vars_unpk[i];
    }

    fprintf(stderr, "nonZeroVariables_AA equal: %d\n", equal42);

    degree_t totalDeg_pk = totalDegree_AA(randDivisor);
    unpackExponentVectors_AA_inp(randDivisor);
    degree_t totalDeg_unpk = totalDegree_AA_unpk(randDivisor);
    int equal43 = totalDeg_unpk == totalDeg_pk;
    fprintf(stderr, "totalDegree_AA equal: %d\n", equal43);

    int equal44 = 1;
    for (int k = 0; k < nvar; ++k) {
        unpackExponentVectors_AA_inp(randDivisor);
        degree_t pDeg_unpk = partialDegree_AA_unpk(randDivisor, k);
        packExponentVectors_AA_inp(randDivisor);
        degree_t pDeg_pk = partialDegree_AA(randDivisor, k);
        equal44 = equal44 && (pDeg_pk == pDeg_unpk);
    }
    fprintf(stderr, "partialDegree_AA equal: %d\n", equal44);

    degree_t mainDeg_pk = mainDegree_AA(randDivisor);
    unpackExponentVectors_AA_inp(randDivisor);
    degree_t mainDeg_unpk = mainDegree_AA_unpk(randDivisor);
    int equal45 = mainDeg_pk == mainDeg_unpk;
    fprintf(stderr, "mainDegree_AA equal: %d\n", equal45);

    packExponentVectors_AA_inp(randDivisor);
    int mainVar_pk = mainVariable_AA(randDivisor);
    unpackExponentVectors_AA_inp(randDivisor);
    degree_t mainVar_unpk = mainVariable_AA_unpk(randDivisor);
    int equal46 = mainVar_pk == mainVar_unpk;
    fprintf(stderr, "mainVariable_AA equal: %d\n", equal46);


    mpq_t coef1;
    mpq_init(coef1);
    mpq_t coef2;
    mpq_init(coef2);
    degree_t degs[nvar];
    for (int i = 0; i < nvar; ++i) {
        degs[i] = rand() % 4;
    }

    packExponentVectors_AA_inp(randDivisor);
    coefficient_AA(randDivisor, degs, nvar, coef1);
    unpackExponentVectors_AA_inp(randDivisor);
    coefficient_AA_unpk(randDivisor, degs, nvar, coef2);
    int equal47 = (mpq_cmp(coef1, coef2) == 0);
    fprintf(stderr, "coefficient_AA equal: %d\n", equal47);

    AltArr_t* randDivisor2 = deepCopyPolynomial_AA(randDivisor);
    mpq_set_si(coef1, 456789876l, 7l);

    unpackExponentVectors_AA_inp(randDivisor);
    setCoefficient_AA_unpk(randDivisor, degs, nvar, coef1);
    packExponentVectors_AA_inp(randDivisor2);
    setCoefficient_AA_unpk(randDivisor2, degs, nvar, coef1);
    int equal48 = isExactlyEqual_AA(randDivisor, randDivisor2);
    fprintf(stderr, "setCoefficient_AA equal: %d\n", equal48);

    AltArr_t* test1 = makePolynomial_AA(10, 2);
    AltArr_t* test2 = makePolynomial_AA(10, 2);

    //test1 = x*y^2*z^0 + x*y*z^0 + x^0*y*z^0;
    //test2 = a^0*x*y^2 + a^0*x*y + a^0*x^0*y;
    degs[0] = 1;
    degs[1] = 2;
    degs[2] = 0;
    setCoefficient_AA(test1, degs, nvar, coef1);
    degs[1] = 1;
    setCoefficient_AA(test1, degs, nvar, coef1);
    degs[0] = 0;
    setCoefficient_AA(test1, degs, nvar, coef1);

    degs[0] = 0;
    degs[1] = 1;
    degs[2] = 2;
    setCoefficient_AA(test2, degs, nvar, coef1);
    degs[2] = 1;
    setCoefficient_AA(test2, degs, nvar, coef1);
    degs[1] = 0;
    setCoefficient_AA(test2, degs, nvar, coef1);

    int xs[2*nvar];
    xs[0] = 0; //a not in 1
    xs[1] = 1; //a in 1
    xs[2] = 1; //x in index 0 of 1
    xs[3] = 2; //x in index 1 of 2
    xs[4] = 2; //y in index 1 of 1
    xs[5] = 3; //y in index 2 of 2 
    xs[6] = 3; //z in index 2 of 1
    xs[7] = 0; //z not in 2.

    int equal49 = isEqualWithVariableOrdering_AA(test1, test2, xs, 8);
    unpackExponentVectors_AA_inp(test1);
    unpackExponentVectors_AA_inp(test2);
    int equal50 = isEqualWithVariableOrdering_AA_unpk(test1, test2, xs, 8);
    int equal51 = (equal49 == equal50) && equal49;

    mpq_set_si(coef2, 1l, 1l);
    setCoefficient_AA(test2, degs, nvar, coef2);

    int equal52 = isEqualWithVariableOrdering_AA_unpk(test1, test2, xs, 8);
    packExponentVectors_AA_inp(test1);
    packExponentVectors_AA_inp(test2);
    int equal53 = isEqualWithVariableOrdering_AA(test1, test2, xs, 8);
    int equal54 = (equal52 == equal53) && equal52;

    fprintf(stderr, "isEqualWithVariableOrdering_AA equal: %d\n", equal51 && equal54);




int isEqualWithVariableOrdering_AA_unpk(AltArr_t* a, AltArr_t* b, const int* xs, int xsSize);


    

    
    return 0;
}


