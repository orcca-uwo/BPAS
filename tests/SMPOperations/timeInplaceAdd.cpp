

#include <gmpxx.h>
#include "RationalNumberPolynomial/SMQP_Support_Test-AA.h"
#include "RationalNumberPolynomial/SMQP_SparseInterpolation-AA.h"
#include <iostream>

#include "Utils/Unix_Timer.h"

int nvar = 3;
unsigned long int coefBound = 4;
long sparsity = 5;
int includeNeg = 1;
double eps = 1e-20;
char* syms[] = {"x", "y", "z", "s", "t", "u", "v", "w"};
time_t seed = 1001l;
long numTerms = 10;
int ntrials = 50;

//USAGE: timeInterpolation.bin nvar sparsity eps maxD T seed
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
	if (argc > 6 && atoi(argv[6]) >= 0) {
		// seed = atoi(argv[6]);
		ntrials = atoi(argv[6]);
	}
	if (argc > 7 && atoi(argv[7]) >= 0) {
	}

    if (nvar == 0) {
        numTerms = 1;
    }

	std::cerr << "nvar " << nvar << " sparsity " << sparsity << " coefbound " << coefBound << " numterms " << numTerms << " ntrials: " << ntrials << std::endl;

	double totalOut = 0;
	double totalIn = 0;


	for (int i = 0; i < ntrials; ++i) {


	Node* a = buildRandomPoly(nvar,  numTerms,  coefBound, sparsity, includeNeg);
	Node* b = buildRandomPoly(nvar,  numTerms,  coefBound, sparsity, includeNeg);

	AltArr_t* aa = deepCopyPolynomial_AAFromNode(a, nvar);
	AltArr_t* ba = deepCopyPolynomial_AAFromNode(b, nvar);
	freePolynomial(a);
	freePolynomial(b);



	timer_id id2 = start_timer();
	aa = addPolynomials_AA_inp(aa,ba,nvar);
	timer_time elapsed2 = elapsed_time(&id2);
	timer_id id = start_timer();
	AltArr_t* ca = addPolynomials_AA(aa,ba,nvar);
	timer_time elapsed1 = elapsed_time(&id);



	double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	double time2 = (elapsed2.tv_sec + ((double)elapsed2.tv_usec / 1000000));


	totalOut += time;
	totalIn += time2;

	freePolynomial_AA(aa);
	freePolynomial_AA(ba);
	freePolynomial_AA(ca);

	}

	std::cerr << "Out Time: " << std::fixed << (totalOut / ntrials) << std::endl; 
	std::cerr << "In  Time: " << std::fixed << (totalIn / ntrials) << std::endl; 

	return 0;

	// //univariate timings
	// AltArr_t* aa = buildRandomPolyFromMax(1, degreeBounds, coefBound, sparsity, includeNeg);
	// std::cerr << "degree: " << degreeBounds[0] << std::endl;

	// mpq_t points[degreeBounds[0]+1];
	// mpq_t vals[degreeBounds[0]+1];

	// time_t seed = time(NULL);

	// MPQRandomGenerator_t* rgen = initMPQRandomGenerator(seed, coefBound*2, coefBound*2);
	// for (int i = 0; i < degreeBounds[0] + 1; ++i) {
	// 	mpq_init(points[i]);
	// 	mpq_init(vals[i]);
	// 	while(1) {
	// 		MPQRandomGeneratorGetQ(rgen, points[i]);
	// 		int unique = 1;
	// 		for (int j = 0; j < i; ++j) {
	// 			if (mpq_cmp(points[i], points[j]) == 0) {
	// 				unique = 0;
	// 				break;
	// 			}
	// 		}
	// 		if (unique) {
	// 			break;
	// 		}
	// 	}
	// 	evalPolyToVal_AA(aa, &(points[i]), 1, vals[i]);
	// }

	// timer_id id = start_timer();
	// AltArr_t* ba = univarInterpolate_AA(points, vals, degreeBounds[0]);
	// timer_time elapsed1 = elapsed_time(&id);

	// id = start_timer();
	// AltArr_t* ca = univariateInterpolation_AA(degreeBounds[0], points, vals);
	// timer_time elapsed2 = elapsed_time(&id);

	// double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// std::cerr << "Lagrange Time: " << time << std::endl;
	// time = (elapsed2.tv_sec + ((double)elapsed2.tv_usec / 1000000));
	// std::cerr << "LinearEquations Time: " << time << std::endl;

	// for (int i = 0; i < degreeBounds[0] + 1; ++i) {
	// 	mpq_clear(points[i]);
	// 	mpq_clear(vals[i]);
	// }




}









