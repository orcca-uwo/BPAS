

#include <gmpxx.h>
#include "../../include/Parser/bpas_parser.h"
#include "../../include/Utils/Unix_Timer.h"

#include <sstream>

#include "../../include/PowerSeries/PowerSeries_Parallel.hpp"

int main(int argc, char** argv) {

	int d = 40;
	if (argc > 1 && argv[1] > 0) {
		d = atoi(argv[1]);
	}

	const char* vars[] = {"x", "y", "z"};

    Poly_ptr  polynom = generate_altarr_var_defined("1+x+y+z", vars, 3);
    PowerSeries_t*  a1 = convertPolyToPowerSeries_PS(polynom);
    PowerSeries_t*  a = inversePowerSeries_PS(a1);
    PowerSeries_t*  b = geometricSeries_PS(3);
    PowerSeries_t* c = dividePowerSeries_PS(a, b);
    // PowerSeries_t* c = multiplyPowerSeries_PS(a, b);
    // PowerSeries_t* c = addPowerSeries_PS(b, b);
    // PowerSeries_t* c = subPowerSeries_PS(b, b);
    // PowerSeries_t* c = negatePowerSeries_PS(b);
    updateToDeg_PS(d, b);
    updateToDeg_PS(d, a);

    fprintf(stderr, "starting...\n" );
    unsigned long long start;
    _startTimer(&start);

    // updateToDeg_PS(d, c);
	updateToDegParallel_PS(d,c,4);
	// updateToDeg_PS(d,b);
	// homogPartParallel_PS(d,b,4);


	float elapsed;
    _stopTimer(&start, &elapsed);

    // printPoly_AA(stderr, b->polys[100], vars, 3);
    // fprintf(stderr, "\n");
    // printPoly_AA(stderr, c->polys[100], vars, 3);
    // fprintf(stderr, "\np->size: %d\n", p->size);
    fprintf(stderr, "\ntime: %9.4f\n", elapsed);
    // int longCount = 0;
    // for (int i = 0; i < p->size; ++i) {
    // 	if (mpz_sizeinbase(mpq_numref(p->elems[i].coef), 2) > 63) {
    // 		++longCount;
    // 	}
    // }
    // fprintf(stderr, "\nlongCounf: %d\n", longCount);


    freePolynomial_AA(polynom);

    // destroyPowerSeries_PS (c);
    destroyPowerSeries_PS (a);
    destroyPowerSeries_PS (b);

	return 0;
}

