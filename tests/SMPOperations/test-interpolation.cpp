

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


void testUnivariateInterpolation (int degreeBound, int coefBound, float sparsity, int includeNeg) {
	AltArr_t* aa = buildRandomPolyFromMax(1, &degreeBound, coefBound, sparsity, includeNeg);
	std::cerr << "aa size: " << aa->size << std::endl;

	mpq_t points[degreeBound+1];
	mpq_t vals[degreeBound+1];

	time_t seed = time(NULL);

	MPQRandomGenerator_t* rgen = initMPQRandomGenerator(seed, coefBound*2, coefBound*2);
	for (int i = 0; i < degreeBound + 1; ++i) {
		mpq_init(points[i]);
		mpq_init(vals[i]);
		while(1) {
			MPQRandomGeneratorGetQ(rgen, points[i]);
			int unique = 1;
			for (int j = 0; j < i; ++j) {
				if (mpq_cmp(points[i], points[j]) == 0) {
					unique = 0;
					break;
				}
			}
			if (unique) {
				break;
			}
		}
		evalPolyToVal_AA(aa, &(points[i]), 1, vals[i]);
	}

	timer_id id = start_timer();
	AltArr_t* ba = univarInterpolate_AA(points, vals, degreeBound+1);
	timer_time elapsed1 = elapsed_time(&id);

	id = start_timer();
	AltArr_t* ca = univariateInterpolation_AA(degreeBound, points, vals);
	timer_time elapsed2 = elapsed_time(&id);

	if (!isExactlyEqual_AA(aa, ba) || !isExactlyEqual_AA(aa, ca)) {
		std::cerr << "Univariate Interpolation Test: FAILED";
		std::cerr << "aa: ";
		printPoly_AA(stderr, aa, syms, nvar);
		std::cerr << std::endl;
		std::cerr << "ba: ";
		printPoly_AA(stderr, ba, syms, nvar);
		std::cerr << std::endl;
		std::cerr << "ca: ";
		printPoly_AA(stderr, ca, syms, nvar);
		std::cerr << std::endl;

	}

	double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	std::cerr << "Lagrange Time: " << time << std::endl;
	time = (elapsed2.tv_sec + ((double)elapsed2.tv_usec / 1000000));
	std::cerr << "LinearEquations Time: " << time << std::endl;

	std::cerr << "Univariate Interpolation Test: PASSED" << std::endl;


	for (int i = 0; i < degreeBound + 1; ++i) {
		mpq_clear(points[i]);
		mpq_clear(vals[i]);
	}


}



int main(int argc, char** argv) {

	std::cerr << "argc: " << argc << std::endl;
	if (argc > 1 && atoi(argv[1]) > 0) {
		nvar = atoi(argv[1]);
	}
	if (argc > 2 && atof(argv[2]) >= 0.0) {
		sparsity = atof(argv[2]);
	}
	if (argc > 3 && atof(argv[3]) > 0) {
		eps = atof(argv[3]);
	}

	int degreeBounds[nvar];
	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] = 3;
	}


	//univar
	testUnivariateInterpolation(10, 32, sparsity, includeNeg);

	std::cerr << "nvar: " << nvar << std::endl;
	std::cerr << "sparsity: " << sparsity << std::endl;
	std::cerr << "eps: " << eps << std::endl << std::endl;

	//test tight degree bounds
	AltArr_t* aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);
	// printPoly_AA(stderr, aa, syms, nvar);
	// fprintf(stderr, "\n");	

	SMQPBlackBox_t* bb = initSMQPBlackBox(aa);
	AltArr_t* interp;
	interp = multivariateInterpolate_AA(bb, degreeBounds, nvar, aa->size + 3, eps);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Interpolation, tight bounds, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Interpolation FAILED!\nInterpolated: ");
		printPoly_AA(stderr, interp, syms, nvar);
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);


	//test loose degreebounds
	aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);
	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] += 2;
	}

	bb = initSMQPBlackBox(aa);
	interp = multivariateInterpolate_AA(bb, degreeBounds, nvar, aa->size + 3, eps);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Interpolation, loose bounds, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Interpolation FAILED!\nInterpolated: ");
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);

	//test uneven degree bounds
	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] -= 2;
	}
	aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);
	if (nvar > 2) {
		degreeBounds[2] += 3;
	} else if (nvar == 2) {
		degreeBounds[1] += 3;
	} else {
		degreeBounds[0] += 3;
	}

	bb = initSMQPBlackBox(aa);
	interp = multivariateInterpolate_AA(bb, degreeBounds, nvar, aa->size + 3, eps);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Interpolation, uneven bounds, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Interpolation FAILED!\nInterpolated: ");
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);

	//test probabilistic with T < 0
	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] = 3;
	}
	aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);

	bb = initSMQPBlackBox(aa);
	interp = multivariateInterpolate_AA(bb, degreeBounds, nvar, -1, eps);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Interpolation, T = -1, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Interpolation FAILED!\nInterpolated: ");
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);

	//////////////////////////////////  Test deterministic one ////////////////

	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] = 3;
	}

	//test tight degree bounds
	aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);
	// printPoly_AA(stderr, aa, syms, nvar);
	// fprintf(stderr, "\n");	

	bb = initSMQPBlackBox(aa);
	interp = multivariateInterpolateDeterministic_AA(bb, degreeBounds, nvar, aa->size + 2);
	// interp = multivariateInterpolateDeterministic_AA(bb, degreeBounds, nvar, -1);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Deterministic Interpolation, tight bounds, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Deterministic Interpolation FAILED!\nInterpolated: ");
		printPoly_AA(stderr, interp, syms, nvar);
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);

	//test loose degreebounds
	aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);
	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] += 2;
	}

	bb = initSMQPBlackBox(aa);
	interp = multivariateInterpolateDeterministic_AA(bb, degreeBounds, nvar, aa->size+5);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Deterministic Interpolation, loose bounds, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Deterministic Interpolation FAILED!\nInterpolated: ");
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);

	//test uneven degree bounds
	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] -= 2;
	}
	aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);
	if (nvar > 2) {
		degreeBounds[2] += 3;
	} else if (nvar == 2) {
		degreeBounds[1] += 3;
	} else {
		degreeBounds[0] += 3;
	}

	bb = initSMQPBlackBox(aa);
	interp = multivariateInterpolateDeterministic_AA(bb, degreeBounds, nvar, aa->size + 5);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Deterministic Interpolation, uneven bounds, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Deterministic Interpolation FAILED!\nInterpolated: ");
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);

	//test -1 T size;

	for (int i = 0; i < nvar; ++i) {
		degreeBounds[i] = 3;
	}

	//test tight degree bounds
	aa = buildRandomPolyFromMax(nvar, degreeBounds, coefBound, sparsity, includeNeg);
	
	bb = initSMQPBlackBox(aa);
	interp = multivariateInterpolateDeterministic_AA(bb, degreeBounds, nvar, -1 );
	// interp = multivariateInterpolateDeterministic_AA(bb, degreeBounds, nvar, -1);

	if (isExactlyEqual_AA(aa, interp)) {
		fprintf(stderr, "Sparse Deterministic Interpolation, tight bounds & -1 T, PASSED!\n\n");
	} else {
		fprintf(stderr, "Sparse Deterministic Interpolation FAILED!\nInterpolated: ");
		printPoly_AA(stderr, interp, syms, nvar);
		fprintf(stderr, "\n");
		exit(1);
	}

	freeSMQPBlackBox(bb);
	freePolynomial_AA(aa);
	freePolynomial_AA(interp);


}









