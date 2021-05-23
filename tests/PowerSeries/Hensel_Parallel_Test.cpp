
#include <stdlib.h>

#include <gmpxx.h>
#include "../../include/PowerSeries/PowerSeries.h"
#include "../../include/PowerSeries/UnivariatePolynomialOverPowerSeries.h"
#include "../../include/PowerSeries/UPOPS_Weierstrass.h"
#include "../../include/PowerSeries/UPOPS_Hensel.hpp"
#include "../../include/PowerSeries/UPOPS_Parallel.hpp"
#include "../../include/Utils/Unix_Timer.h"
#include "../../include/Parser/bpas_parser.h"
#include <time.h>

#include <vector>

int degree =50;
int nthreads = 4;
int family = 0;
int upopDeg = 8;

std::string family0[] = {
	"0", //0
	"0", //1
	"0", //2
	"0", //3
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)*(z-11)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)*(z-11)*(z-12)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)*(z-11)*(z-12)*(z-13)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)*(z-11)*(z-12)*(z-13)*(z-14)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)*(z-11)*(z-12)*(z-13)*(z-14)*(z-15)",
	"y*(z^3+z)+(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)*(z-11)*(z-12)*(z-13)*(z-14)*(z-15)*(z-16)"
};

std::string family1[] = {
	"0", //0
	"0", //1
	"0", //2
	"0", //3
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9*(z-10)^10",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9*(z-10)^10*(z-11)^11",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9*(z-10)^10*(z-11)^11*(z-12)^12",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9*(z-10)^10*(z-11)^11*(z-12)^12*(z-13)^13",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9*(z-10)^10*(z-11)^11*(z-12)^12*(z-13)^13*(z-14)^14",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9*(z-10)^10*(z-11)^11*(z-12)^12*(z-13)^13*(z-14)^14*(z-15)^15",
	"y*(z^3+z)+(z-1)*(z-2)^2*(z-3)^3*(z-4)^4*(z-5)^5*(z-6)^6*(z-7)^7*(z-8)^8*(z-9)^9*(z-10)^10*(z-11)^11*(z-12)^12*(z-13)^13*(z-14)^14*(z-15)^15*(z-16)^16"
};


std::string family2[] = {
	"0", //0
	"0", //1
	"0", //2
	"0", //3
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)*(z+x+y-10)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)*(z+x+y-10)*(z+x+y-11)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)*(z+x+y-10)*(z+x+y-11)*(z+x+y-12)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)*(z+x+y-10)*(z+x+y-11)*(z+x+y-12)*(z+x+y-13)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)*(z+x+y-10)*(z+x+y-11)*(z+x+y-12)*(z+x+y-13)*(z+x+y-14)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)*(z+x+y-10)*(z+x+y-11)*(z+x+y-12)*(z+x+y-13)*(z+x+y-14)*(z+x+y-15)",
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)*(z+x+y-5)*(z+x+y-6)*(z+x+y-7)*(z+x+y-8)*(z+x+y-9)*(z+x+y-10)*(z+x+y-11)*(z+x+y-12)*(z+x+y-13)*(z+x+y-14)*(z+x+y-15)*(z+x+y-16)"
};

std::string family3[] = {
	"0", //0
	"0", //1
	"0", //2
	"0", //3
	"x*y*(z^3+z)+(z+x-1)*(z+x-2)*(z+x-3)*(z+x-4)"
};

std::string family4[] = {
	"0", //0
	"0", //1
	"0", //2
	"0", //3
	"x*y*(z^3+z)+(z+x+y-1)*(z+x+y-2)*(z+x+y-3)*(z+x+y-4)"
};



std::vector<int> famThreads0[] = {
	{0}, //0
	{0}, //1
	{0}, //2
	{0}, //3
	{7, 3, 0, 2},
	{0}, //5
	{5, 3, 2, 0, 0, 2},
	{0}, //7
	{3, 2, 2, 2, 2, 0, 0, 1},
	{0}, //9
	{3, 2, 2, 2, 1, 1, 0, 0, 0, 1}
};

std::vector<int> famThreads1[] = {
	{0}, //0
	{0}, //1
	{0}, //2
	{0}, //3
	{4, 4, 0, 4},
	{0}, //5
	{1, 3, 3, 3, 0, 2},
	{0}, //7
	{0, 1, 2, 2, 3, 2, 0, 2},
	{0}, //9
	{0, 0, 2, 2, 2, 2, 2, 1, 0, 1}
};


std::vector<int> famThreads2[] = {
	{0}, //0
	{0}, //1
	{0}, //2
	{0}, //3
	{7, 3, 0, 2},
	{0}, //5
	{5, 3, 2, 0, 0, 2},
	{0}, //7
	{4, 2, 2, 2, 1, 0, 0, 1},
	{0}, //9
	{3, 2, 2, 2, 1, 1, 0, 0, 0, 2}
};

#if PS_PARALLEL_TIME
extern float g_shiftTime;
float g_totalShift;
#endif


/**
 * f : the input upop
 * C : the given roots
 * bound : the upper bound for the accuracy
 */
void TestingFactorizationViaHensel(Upops_t* f, int bound, const char** sym, int nvar, int nthreads) {

#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
	fprintf(stderr, "F:\n");
	print_UPOPS(stderr, f , sym);
	fprintf(stderr, "\n");
#endif

	int n;
	Upops_t** facts_out = NULL;
	HenselFactorization_UPOPS(f, &facts_out, &n);

	// Upops_t* prod = facts_out[0];
	// reserve_UPOPS(prod);
	// for (int i = 1; i <= n-1; ++i) {
	// 	Upops_t* tmp = multiplyUnivariatePolyOverPowerSeries_UPOPS(prod, facts_out[i]);
	// 	destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);
	// 	prod = tmp;
	// }
	// Upops_t* sub =  subUnivariatePolyOverPowerSeries_UPOPS(prod,f);

	float time[n];
	for (int j = 0; j < n; ++j) {
		time[j] = 0.f;
	}


	unsigned long long start;
	// for (int d = 1; d <= bound; ++d) {
#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
		fprintf(stderr, "\n");
		fprintf(stderr, "In degree: %d\n ", d);
		fprintf(stderr, "\n");
#endif


if (nthreads > 1) {

	updateHenselFactsParallel_UPOPS(bound, facts_out, n, nthreads);//n-1);
	// updateToDeg_UPOPS(bound, facts_out[n-1]);

	for (int j = 0; j < n; ++j) {
		gmp_fprintf(stderr, "upops[%d]: %Qd, %d\n",j, facts_out[j]->data[0]->polys[0]->elems->coef, facts_out[j]->deg);
	}
}
else if (nthreads <= 0) {

	std::vector<int> threadsPer;
	if (family == 0) {
		threadsPer = famThreads0[upopDeg];
	}
	if (family == 1) {
		threadsPer = famThreads1[upopDeg];
	}
	if (family == 2) {
		threadsPer = famThreads2[upopDeg];
	}

	for (size_t i = 0; i < threadsPer.size(); ++i) {
		fprintf(stderr, "facts[%lu]: %d\n", i, threadsPer[i]);
	}
	updateHenselFactsParallel_UPOPS(bound, facts_out, n, threadsPer);
	// updateHenselFactsParallel_UPOPS(bound, facts_out, n, nthreads);//n-1);
	// updateToDeg_UPOPS(bound, facts_out[n-1]);

	for (int j = 0; j < n; ++j) {
		gmp_fprintf(stderr, "upops[%d]: %Qd, %d\n",j, facts_out[j]->data[0]->polys[0]->elems->coef, facts_out[j]->deg);
	}
}
else {
/////////////////////////////////////////////////////////////////////////////////
		// for (int j = 0; j < 1; ++j) {
		for (int j = 0; j < n; ++j) {
			fprintf(stderr, "starting facts[%d]\n", j );
#if PS_PARALLEL_TIME
			_startTimer(&start);
			updateToDeg_UPOPS(bound, facts_out[j]);
			_stopTimerAddElapsed(&start, time + j);
			fprintf(stderr, "shiftTime[%d] = %10f\n", j, g_shiftTime );
			g_totalShift += g_shiftTime;
			// time[j] -= g_shiftTime;
			g_shiftTime = 0.f;
#else
			updateToDeg_UPOPS(bound, facts_out[j]);
#endif


#if defined(BPAS_POWERSERIES_DEBUG) && BPAS_POWERSERIES_DEBUG
			fprintf(stderr, "\n");
			fprintf(stderr, "The upops at index: %d\n",j);
			fprintf(stderr, "\n");
			fprintf(stderr, "\n");
			print_UPOPS(stderr, facts_out[j] , sym);
			fprintf(stderr, "\n");
#endif
		}
/////////////////////////////////////////////////////////////////////////////////
	for (int j = 0; j < n; ++j) {
		gmp_fprintf(stderr, "upops[%d]: %Qd, %d\n",j, facts_out[j]->data[0]->polys[0]->elems->coef, facts_out[j]->deg);
	}

#if PS_PARALLEL_TIME
	for (int j = 0; j < n; ++j) {
		fprintf(stderr, "time[%d]: %10f\n",j, time[j]);
	}
	fprintf(stderr, "total shift: %10f\n", g_totalShift);
#endif
}

	// Poly_ptr polyPart = polynomialPart_UPOPS(bound, sub);
	// if (!isZero_AA(polyPart)) {
	// 	fprintf(stderr, "Testing Factorization Via Hensel: ERROR\n");
	// 	exit(1);
	// } else {
	// 	fprintf(stderr, "Testing Factorization Via Hensel: PASSED\n");
	// }
	// freePolynomial_AA(polyPart);
	// destroyUnivariatePolynomialOverPowerSeries_UPOPS(sub);
	// destroyUnivariatePolynomialOverPowerSeries_UPOPS(prod);

	for (int i = 0; i < n; ++i) {
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(facts_out[i]);
	}
	free(facts_out);
}


int main(int argc, char** argv) {
	if (argc > 1) {
        degree = atoi(argv[1]);
    }

    if (argc > 2) {
        nthreads = atoi(argv[2]);
    }

    if (argc > 3 && atoi(argv[3]) >= 0) {
    	family = atoi(argv[3]);
    }

    if (argc > 4 && atoi(argv[4]) > 0) {
    	upopDeg = atoi(argv[4]);
    }

    fprintf(stderr, "degree: %d, nthreads: %d, family: %d, upopDeg: %d\n", degree, nthreads, family, upopDeg);

	const char* sym[] = {"z", "y", "x"};

	AltArr_t* a;
	int nvar;
	switch (family) {
		case 0 : {
			a = generate_altarr_var_defined(family0[upopDeg].c_str(), sym, 2);
			nvar = 2;
			break;
		}
		case 1 : {
			a = generate_altarr_var_defined(family1[upopDeg].c_str(), sym, 2);
			nvar = 2;
			break;
		}
		case 2 : {
			a = generate_altarr_var_defined(family2[upopDeg].c_str(), sym, 3);
			nvar = 3;
			break;
		}
	}
	// if (family == 3) {
	// 	a = generate_altarr_var_defined(family3[upopDeg].c_str(), sym, 3);
	// } else if (family == 4) {
	// 	a = generate_altarr_var_defined(family4[upopDeg].c_str(), sym, 3);
	// }

	// fprintf(stderr, "a: \n" );
	// printPoly_AA(stderr, a, sym, a->nvar);

	// AltArr_t* a = generate_altarr_var_defined("y*z^2+z^3+y*z-6*z^2+11*z-6", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("y*z^3+z^4-10*z^3+y*z+35*z^2-50*z+24", sym, 2); //(y-1)(y-2)(y-3)(y-4) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("y*z^3+z^4-31*z^3+y*z+305*z^2-1025*z+750", sym, 2); //(z-1)(z-5)(z-10)(z-15) + y(z^3+z)
	// AltArr_t* a = generate_altarr_var_defined("y*z^3+z^4-1610*z^3+y*z+666000*z^2-56500000*z+500000000", sym, 2); //(y-10)(y-100)(y-500)(y-1000) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("y*z^3+z^4-1111*z^3+y*z+112110*z^2-1111000*z+1000000", sym, 2); //(y-1)(y-10)(y-100)(y-1000) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("y*z^3+z^4-46*z^3+y*z+791*z^2-6026*z+17160", sym, 2); //(y-10)(y-11)(y-12)(y-13) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("y*z^3+z^4-406*z^3+y*z+61811*z^2-4182206*z+106110600", sym, 2); //(y-100)(y-101)(y-102)(y-103) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("z^5+y*z^3-510*z^4+104035*z^3+y*z-10610550*z^2+541060024*z-11035502400", sym, 2); //(y-100)(y-101)(y-102)(y-103)*(y-104) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("z^5+y*z^3-20*z^4+155*z^3+y*z-580*z^2+1044*z-720", sym, 2); //(y-2)(y-3)(y-4)(y-5)*(y-6) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("z^6-615*z^5+y*z^3+157585*z^4-21534225*z^3+y*z+1655167774*z^2-67846804920*z+1158727752000", sym, 2); //(y-100)(y-101)(y-102)(y-103)*(y-104)(y-105) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("z^7-721*z^6+222775*z^5+y*z^3-38238235*z^4+3937795624*z^3+y*z-243294588964*z^2+8350489073520*z-122825141712000", sym, 2); //(y-100)(y-101)(y-102)(y-103)*(y-104)(y-105)(y-106) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-16)*(z-128)*(z-1024)*(z-8192)*(z-98304)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-16)*(z-128)*(z-1024)*(z-8192)*(z-98304)*(z-2)*(z-16)*(z-128)*(z-1024)*(z-8192)*(z-98304)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-16)*(z-16)*(z-16)*(z-128)*(z-128)*(z-8192)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-16)*(z-16)*(z-16)*(z-128)*(z-128)*(z-65536)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-16)*(z-16)*(z-16)*(z-128)*(z-128)*(z-65536)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-16)*(z-16)*(z-16)*(z-128)*(z-128)*(z-65536)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-128)*(z-128)*(z-128)*(z-8192)*(z-8192)*(z-524288)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-16)*(z-16)*(z-16)*(z-16)*(z-128)*(z-128)*(z-128)*(z-1024)*(z-1024)*(z-8192)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-8)*(z-8)*(z-8)*(z-32)*(z-32)*(z-128)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-8)*(z-8)*(z-8)*(z-8)*(z-32)*(z-32)*(z-32)*(z-128)*(z-128)*(z-512)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-5)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-6)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-8192)*(z-8192)*(z-8192)*(z-8192)*(z-128)*(z-128)*(z-128)*(z-16)*(z-16)*(z-2)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-100)*(z-101)*(z-102)*(z-103)*(z-104)*(z-105)*(z-106)*(z-107)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-100)*(z-101)*(z-102)*(z-103)*(z-104)*(z-105)*(z-106)*(z-107)*(z-108)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)*(z-11)*(z-12)*(z-13)*(z-14)*(z-15)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-100)*(z-200)*(z-300)*(z-400)*(z-500)*(z-600)*(z-700)*(z-800)*(z-900)*(z-1000)*(z-1100)*(z-1200)*(z-1300)*(z-1400)*(z-1500)*(z-1600)*(z-1700)*(z-1800)*(z-1900)*(z-2000)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-1)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-1)*(z-1)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-4)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-5)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-6)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-6)*(z-7)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-6)*(z-7)*(z-8)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)*(z-5)*(z-6)*(z-6)*(z-6)*(z-7)*(z-7)*(z-8)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-400)*(z-400)*(z-400)*(z-400)*(z-400)*(z-500)*(z-500)*(z-500)*(z-500)*(z-600)*(z-600)*(z-600)*(z-700)*(z-700)*(z-800)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-10004)*(z-10004)*(z-10004)*(z-10004)*(z-10004)*(z-10003)*(z-10003)*(z-10003)*(z-10003)*(z-10002)*(z-10002)*(z-10002)*(z-10001)*(z-10001)*(z-10000)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-1000000004)*(z-1000000004)*(z-1000000004)*(z-1000000004)*(z-1000000003)*(z-1000000003)*(z-1000000003)*(z-1000000003)*(z-1000000002)*(z-1000000002)*(z-1000000002)*(z-1000000002)*(z-1000000001)*(z-1000000001)*(z-1000000001)*(z-1000000001)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-1000000001)*(z-1000000001)*(z-1000000001)*(z-1000000001)*(z-1000000002)*(z-1000000002)*(z-1000000002)*(z-1000000003)*(z-1000000003)*(z-1000000004)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-3)*(z-11)*(z-11)*(z-11)*(z-11)*(z-12)*(z-12)*(z-12)*(z-13)*(z-13)*(z-14)+y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2 + y*z)*(z-3 + 2*y*z)*(z-4 - y*z)*(z-5 - 2*y*z) + y*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-2 + y^5*z + 3*y^4 + 2*y^3*z - 1*y^2 + y)*(z-3 + 2*y*z + y^5*z + 3*y^4 + 2*y^3*z - 1*y^2 + y)*(z-4 - y*z + y^5*z + 3*y^4 + 2*y^3*z - 1*y^2 + y)*(z-5 - 2*y*z + y^5*z + 3*y^4 + 2*y^3*z - 1*y^2 + y) + y*(z^3+z)", sym, 2);


	// AltArr_t* a = generate_altarr_var_defined("z^6-27*z^5+y*z^3+295*z^4-1665*z^3+y*z+5104*z^2-8028*z+5040", sym, 2); //(y-2)(y-3)(y-4)(y-5)*(y-6)(y-7) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("z^7-35*z^6+511*z^5+y*z^3-4025*z^4+18424*z^3+y*z-48860*z^2+69264*z-40320", sym, 2); //(y-2)(y-3)(y-4)(y-5)*(y-6)(y-7)(y-8) + x(y^3+y)
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)+(x^2 + y)*(z^3+z)", sym, 3);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)+(y)*(z^3+z)", sym, 2); //3 threads for first Weierstrass
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-2)*(z-2)*(z-3)*(z-3)*(z-4)*(z-4)*(z-5)*(z-5)+(y)*(z^3+z)", sym, 2); //3 threads for first Weierstrass
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)*(z-5)+(y)*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)*(z-5)+(y)*(z^3+z)", sym, 2);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-1)*(z-2)*(z-2)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)*(z-5)*(z-5)+(y+x^2)*(z^3+z)", sym, 3);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)*(z-5)*(z-5)+(y)*(z^3+z)", sym, 3);
	// AltArr_t* a = generate_altarr_var_defined("(z-1)*(z-2)*(z-2)*(z-3)*(z-3)*(z-3)*(z-4)*(z-4)*(z-4)*(z-4)*(z-5)*(z-5)*(z-5)*(z-5)*(z-5)+(y+x^2)*(z^3+z)", sym, 3);

	// AltArr_t* a = generate_altarr_var_defined("(z-1+x)*(z-2+x)*(z-3+x)*(z-4+x) + x*y*(z^3+z)", sym, 3);

	// fprintf(stderr, "a:= \n" );
	// printPoly_AA(stderr, a, sym, 2);
	// fprintf(stderr, "\n");


	Upops_t* f = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(a);
	unsigned long long start;
    _startTimer(&start);
	TestingFactorizationViaHensel(f, degree, sym, nvar, nthreads);
	float elapsedTime;
	_stopTimer(&start, &elapsedTime);
	fprintf(stderr, "Time total: %f\n", elapsedTime);

	destroyUnivariatePolynomialOverPowerSeries_UPOPS(f);
	freePolynomial_AA(a);

	return 0;


}