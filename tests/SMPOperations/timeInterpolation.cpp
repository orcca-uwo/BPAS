

#include <gmpxx.h>
#include "RationalNumberPolynomial/SMQP_Support_Test-AA.h"
#include "RationalNumberPolynomial/SMQP_SparseInterpolation-AA.h"
#include <iostream>

#include "Utils/Unix_Timer.h"

int nvar = 3;
unsigned long int coefBound = 4;
float sparsity = 0.8f;
int includeNeg = 1;
double eps = 1e-20;
char* syms[] = {"x", "y", "z", "s", "t", "u", "v", "w"};
int maxD = 3;
time_t seed = 1001l;

//USAGE: timeInterpolation.bin nvar sparsity eps maxD T seed
int main(int argc, char** argv) {

	int T = 3;

	if (argc > 1 && atoi(argv[1]) > 0) {
		nvar = atoi(argv[1]);
	}
	if (argc > 2 && atof(argv[2]) >= 0.0) {
		sparsity = atof(argv[2]);
	}
	if (argc > 3 && atof(argv[3]) > 0) {
		eps = atof(argv[3]);
	}
	if (argc > 4 && atoi(argv[4]) > 0) {
		maxD = atoi(argv[4]);
	}
	if (argc > 5 && atoi(argv[5])) {
		T = atoi(argv[5]);
	}
	if (argc > 6 && atoi(argv[6]) >= 0) {
		seed = atoi(argv[6]);
	}

	int degreeBounds[nvar];
	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] = maxD;
	}

	std::cerr << "nvar " << nvar << " sparsity " << sparsity << " eps " << eps << " maxD " << maxD << " seed " << seed << std::endl;

	//test tight degree bounds


	AltArr_t* aa = buildRandomPolyFromMax_seeded(seed, nvar, degreeBounds, coefBound, sparsity, includeNeg);
	// printPoly_AA(stderr, aa, syms, nvar);
	SMQPBlackBox_t* bb = initSMQPBlackBox(aa);
	int numTerms = T < 0 ? T : aa->size + T;
	// std::cerr << "numTerms: " << numTerms << std::endl;

	timer_id id = start_timer();
	AltArr_t* probInterp = multivariateInterpolate_AA(bb, degreeBounds, nvar, numTerms, eps);
	timer_time elapsed1 = elapsed_time(&id);

	timer_id id2 = start_timer();
	AltArr_t* detInterp = multivariateInterpolateDeterministic_AA(bb, degreeBounds, nvar, numTerms);
	timer_time elapsed2 = elapsed_time(&id2);

	double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	double time2 = (elapsed2.tv_sec + ((double)elapsed2.tv_usec / 1000000));

	std::cerr << "Probabilistic Time: " << std::fixed << time << std::endl; 
	std::cerr << "Deterministic Time: " << std::fixed << time2 << std::endl; 

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(probInterp);
	freePolynomial_AA(detInterp);

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









