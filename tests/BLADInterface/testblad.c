
#include <iostream>
#include <gmpxx.h>

#include "../../include/BLADInterface/bladinterface.h"
#include "../../include/IntegerPolynomial/SMZP_Support_Test.h"
#include "../../include/IntegerPolynomial/SMZP_CppSupport.hpp"
#include "../../include/Parser/bpas_parser.h"

// struct bav_term {

//     ba0_int_p alloc;
//     ba0_int_p size;
//     struct bav_rank* rg;
// };

// struct bav_rank {
//     struct bav_variable* var;
//     bav_Idegree deg;
// };


// struct bav_variable {
//     struct bav_symbol* root;
//     ba0_int_p hack;
//     ba0_int_p index;
//     struct bav_tableof_Inumber number;
//     struct bav_tableof_Iorder order;
//     struct bav_tableof_variable derivative;
// };


int nvar = 3;
int numTerms = 10;
int coefBound = 5;
degree_t sparsity = 10;
int includeNeg = 1;

void testBLADConversion() {

	char* syms[nvar];
	char var[16];
	for (int i = 0; i < nvar; ++i) 
	{
		sprintf(var, "x_%d", i);
		syms[i] = (char*) malloc(sizeof(char)*16);
		strcpy(syms[i], var);
	}

	const char** names = (const char**) &(syms[0]);

	initBLAD();

	//Construct random poly
	Node* node = buildRandomZPoly(nvar, numTerms, coefBound, sparsity, includeNeg);
	AltArrZ_t* aaz = deepCopyPolynomial_AAZFromNode(node, nvar);
    freePolynomial(node);

    blad_poly_z A = convertToBLAD_AAZ(aaz, names);

	std::string aazStr = polyToString_AAZ(aaz, names);
	char* bladCStr = (char*) malloc(sizeof(char)*aazStr.length()*2);
	char formatString[16] = "%Az";
	ba0_sprintf(bladCStr, formatString, A);

	// parse blad string
	AltArr_t* aaq_blad = generate_altarr_var_defined(bladCStr, (char**)names, nvar);
	AltArrZ_t* aaz_blad = deepCopyPolynomial_AAZFromAA(aaq_blad);
	freePolynomial_AA(aaq_blad);
	// for (int i = 0; i < pack->numVars; ++i) {
	// 	fprintf(stderr, "vars[%d]: %s\n", i, pack->vars[i]);
	// }
	// AltArrZ_t* aaz_blad = deepCopyPolynomial_AAZFromAA(pack->altarr_t_data);
	// freePolynomial_AA(pack->altarr_t_data);
	// free(pack->vars);
	// free(pack);


	//As of November 2018 the parsing is broken if some variables have degree 0 
	//So check strings directly as well
	std::string parsedStr = polyToString_AAZ(aaz, names);
	if(!isExactlyEqual_AAZ(aaz, aaz_blad) && strcmp(bladCStr, parsedStr.c_str()) != 0) {
		std::cerr << "Convert to BLAD: ERROR! Polynomials did not match!" << std::endl;
		std::cerr << "BLAD Print: " << std::endl;
		char formatString[16] = "%Az\n";
		ba0_fprintf(stderr, formatString, A);
		std::cerr << "BLAD String: " << std::endl;
		fprintf(stderr, "%s\n", bladCStr);
		std::cerr << "BPAS Parse From BLAD: " << std::endl << polyToString_AAZ(aaz_blad, names) << std::endl;
		std::cerr << "BPAS Original: " << std::endl << aazStr << std::endl;
		exit(1);
	}

	std::cerr << "Convert to BLAD: PASSED" << std::endl;

	degree_t maxDegs[nvar];
	partialDegrees_AAZ(aaz, maxDegs);
	int hasAllVars = 1;
	for (int i = 0; i < nvar; ++i) {
		if (maxDegs[i] == 0) {
			hasAllVars = 0;
			break;
		}
	}

	AltArrZ_t* aaz2;
	if (hasAllVars) {
		aaz2 = convertFromBLAD_AAZ(A);
	} else {
		aaz2 = convertFromBLADWithVars_AAZ(A, names, nvar);
	}

	if (!isExactlyEqual_AAZ(aaz, aaz2)) {
		std::cerr << "Convert from BLAD: ERROR! Polynomials did not match!" << std::endl;
		std::cerr << "BLAD: " << std::endl;
		char formatString[16] = "%Az\n";
		ba0_fprintf(stderr, formatString, A);
		std::cerr << "Expected: " << std::endl << aazStr << std::endl;
		std::string aaz2Str = polyToString_AAZ(aaz2, names);
		std::cerr << "Got: " << std::endl << aaz2Str << std::endl;
		exit(1);		
	}

	std::cerr << "Convert from BLAD: PASSED" << std::endl;

	freeBLAD();
}

void testBLADFactor() {
	Node* node = buildRandomZPoly(nvar, numTerms, coefBound, sparsity, includeNeg);
	AltArrZ_t* aaz = deepCopyPolynomial_AAZFromNode(node, nvar);
    freePolynomial(node);
    unpackExponentVectors_AAZ_inp(aaz);

    char* syms[] = {"x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7"};
    mpz_t numericFact;
    mpz_init(numericFact);
    AltArrZ_t** facts = NULL;
    int* exps;
    int numFacts;
    factorPolynomialBLAD_AAZ(aaz, (const char**) syms, &numFacts, &facts, &exps, numericFact);

    fprintf(stderr, "aaz: \n");
    printPoly_AAZ(stderr, aaz, syms, nvar);
    fprintf(stderr, "\n");
    gmp_fprintf(stderr, "numfact: %Zd\n", numericFact);
    for (int i = 0;i  < numFacts; ++i) {
    	fprintf(stderr, "factor %d, exp: %d\n", i, exps[i]);
	    printPoly_AAZ(stderr, facts[i], syms, nvar);
    	fprintf(stderr, "\n");
    }
}


int main(int argc, char** argv) {

	if (argc > 1 && atol(argv[1]) >= 0) {
        nvar = atol(argv[1]);
    }
    if (argc > 2 && atol(argv[2]) > 0) { 
        numTerms = atol(argv[2]);
    }
    if (argc > 3 && atol(argv[3]) > 0) {
        coefBound = atol(argv[3]);
    }
    if (argc > 4 && atol(argv[4]) > 1) {
        sparsity = atol(argv[4]);
    }
    if (argc > 5 && atoi(argv[5]) >= 0) {
        includeNeg = atoi(argv[5]);
    }

    testBLADFactor();
    // testBLADConversion();

    // initBLAD();
    // freeBLAD();
    // initBLAD();
 //    //Construct random poly
	// Node* node = buildRandomZPoly(nvar, numTerms, coefBound, sparsity, includeNeg);
	// AltArrZ_t* aaz = deepCopyPolynomial_AAZFromNode(node, nvar);
 //    freePolynomial(node);

	// node = buildRandomZPoly(nvar, numTerms, coefBound, sparsity, includeNeg);
 //    AltArrZ_t* aaz2 = deepCopyPolynomial_AAZFromNode(node, nvar);
 //    freePolynomial(node);

 //    AltArrZ_t* prod = multiplyPolynomials_AAZ(aaz, aaz2, nvar);
 //    aaz = prod;

	// char* syms[nvar];
	// char var[16];
	// for (int i = 0; i < nvar; ++i) 
	// {
	// 	sprintf(var, "x_%d", i);
	// 	syms[i] = (char*) malloc(sizeof(char)*16);
	// 	strcpy(syms[i], var);
	// }

	// const char** names = (const char**) &(syms[0]);

 //    blad_product_z facts = factorPolynomialBLAD_AAZ(aaz, names);
	// std::cerr << "Original: \n" << polyToString_AAZ(aaz, names) << std::endl;

 //    fprintf(stderr, "facts size: %d\n", facts->size );
 //    gmp_fprintf(stderr, "facts coef: %Zd\n", facts->num_factor);

 //    AltArrZ_t** aaFacts = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*)*facts->size);
 //    long long* exps = (long long*) malloc(sizeof(long long)*facts->size);

 //    for (int i = 0; i < facts->size; ++i) {
 //    	exps[i] = facts->tab[i].exponent;
 //    	aaFacts[i] = convertFromBLADWithVars_AAZ(&facts->tab[i].factor, names, nvar);
 //    	fprintf(stderr, "exp: %lld\n", facts->tab[i].exponent);
 //    	ba0_fprintf(stderr, "fact:\n%Az\n", &facts->tab[i].factor);
 //    	// ba0_fprintf(stderr, "(%Az)^%d\n", facts->tab[i].factor, facts->tab[i].exponent);
 //    }

 //    for (int i = 0; i < facts->size; ++i) {
 //    	std::cerr << "AA nvar: " << aaFacts[i]->nvar << std::endl;
 //    	std::cerr << "AAZ factor: " << std::endl << polyToString_AAZ(aaFacts[i], names) << std::endl;
 //    }

    // freeBLAD();

    return 0;


}
