
#include <gmpxx.h>
#include <stdlib.h>

#include "../../include/Parser/bpas_parser.h"
#include "../../include/Utils/Unix_Timer.h"
#include "../../include/PowerSeries/UPOPS_Weierstrass.h"
#include "../../include/PowerSeries/UPOPS_Parallel.hpp"


#include "../../include/Utils/RandomHelpers.h"
#include "../../include/Utils/Parallel/ExecutorThreadPool.hpp"

void weierstrassUpdateParallel_UPOPS(int d, Upops_t* p, int nthreads);

std::string familyA[] = {
	"0", //0
	"0", //1
	"0", //2
	"0", //3
	"(4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	"0", //5
	"(6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	"0", //7
	"(8+y+x^2)*z^8 + (7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	"0", //9
	"(10+y+x^2)*z^10 + (9+y+x^2)*z^9 + (8+y+x^2)*z^8 + (7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	"0", //11
	"(12+y+x^2)*z^12 + (11+y+x^2)*z^11 + (10+y+x^2)*z^10 + (9+y+x^2)*z^9 + (8+y+x^2)*z^8 + (7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2"
};


std::string familyB[] = {
	"0", //0
	"0", //1
	"0", //2
	"0", //3
	"(4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	"0", //5
	// "(5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	// "(6)*z^6 + (5)*z^5 + (4)*z^4 + (3)*z^3 + (y)*z^2 + (y)*z + x^2",
	"(6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (y+x^2)*z^2 + (y+x^2)*z + x^2 + y^2",
	"0", //7
	// "(7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (y+x^2)*z^3 + (y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	// "(8)*z^8 + (7)*z^7 + (6)*z^6 + (5)*z^5 + (4)*z^4 + (y)*z^3 + (y)*z^2 + (y)*z + x^2",
	"(8+y+x^2)*z^8 + (7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (y+x^2)*z^3 + (y+x^2)*z^2 + (y+x^2)*z + x^2 + y^2",
	"0", //9
	"(10+y+x^2)*z^10 + (9+y+x^2)*z^9 + (8+y+x^2)*z^8 + (7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (y+x^2)*z^4 + (y+x^2)*z^3 + (y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2",
	"0", //11
	"(12+y+x^2)*z^12 + (11+y+x^2)*z^11 + (10+y+x^2)*z^10 + (9+y+x^2)*z^9 + (8+y+x^2)*z^8 + (7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (y+x^2)*z^5 + (y+x^2)*z^4 + (y+x^2)*z^3 + (y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2"
};


int main(int argc, char** argv) {
	int testDeg = 50;
	if (argc > 1 && atoi(argv[1])) {
		testDeg = atoi(argv[1]);
	}

	int nthreads = 4;
	if (argc > 2 && atoi(argv[2])) {
		nthreads = atoi(argv[2]);
	}

	int family = 0;
	if (argc > 3 && atoi(argv[3]) >= 0) {
		family = atoi(argv[3]);
	}

	int deg = 4;
	if (argc > 4 && atoi(argv[4]) > 3) {
		deg = atoi(argv[4]);
	}

	const char* sym[] = {"z", "y", "x"};
    // Poly_ptr  poly_lc = generate_altarr_var_defined("1+x+y", sym, 3);
    // Poly_ptr  poly_y = generate_altarr_var_defined("y", sym, 3);
    // Poly_ptr  poly_4 = generate_altarr_var_defined("y+1", sym, 3);
    // Poly_ptr  poly_x = generate_altarr_var_defined("x", sym, 3);
    // PowerSeries_t* ps_lc = convertPolyToPowerSeries_PS(poly_lc);

    // PowerSeries_t* ps_1 = inversePowerSeries_PS(ps_lc);
    // PowerSeries_t* ps_y = convertPolyToPowerSeries_PS(poly_y);
    // PowerSeries_t* ps_x = convertPolyToPowerSeries_PS(poly_x);
    // PowerSeries_t* ps_one = onePowerSeries_PS(3);
    // PowerSeries_t* ps_array[100];

    // ps_array[0] = ps_y;
    // ps_array[1] = ps_x;
    // ps_array[2] = ps_x;
    // ps_array[3] = ps_x;
    // ps_array[4] = ps_one;
    // ps_array[5] = ps_one;
    // ps_array[6] = ps_one;
    // ps_array[7] = ps_one;
    // ps_array[8] = ps_1;

    // ps_array[0] = ps_y;
    // ps_array[1] = ps_y;
    // ps_array[2] = ps_y;
    // ps_array[3] = ps_y;
    // ps_array[4] = ps_one;
    // ps_array[5] = ps_one;
    // ps_array[6] = ps_one;
    // ps_array[7] = ps_one;
    // ps_array[8] = ps_one;

    // ps_array[0] = ps_y;
    // ps_array[1] = ps_x;
    // ps_array[2] = ps_x;
    // ps_array[3] = ps_x;
    // ps_array[4] = ps_x;
    // ps_array[5] = ps_x;
    // ps_array[6] = ps_one;
    // ps_array[7] = ps_one;
    // ps_array[8] = ps_one;
    // ps_array[9] = ps_one;
    // ps_array[10] = ps_one;
    // ps_array[11] = ps_one;
    // ps_array[12] = ps_one;//ps_1;
    // Upops_t* upop = convertArrayOfPSToUPOPS_UPOPS(ps_array,9);

    Poly_ptr wpCase;
	// Poly_ptr wpCase = generate_altarr_var_defined("x*z^8 + z^7 + z^6 + z^5 + z^4 + z^3 + z^2 + (x+y)*z + x+y", sym, 3);
	// Poly_ptr wpCase = generate_altarr_var_defined("x*z^4 + z^3 + z^2 + (x+y)*z + x+y", sym, 3);
	// Poly_ptr wpCase = generate_altarr_var_defined("z^3 + z^2 + y*z + y", sym, 3);
	// Poly_ptr wpCase = generate_altarr_var_defined("(6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (y+x^2)*z^2 + (y+x^2)*z + y^2 + x^2 + x*y", sym, 3);
	// Poly_ptr wpCase = generate_altarr_var_defined("(6+y+y^2)*z^6 + (5+y+y^2)*z^5 + (4+y+y^2)*z^4 + (3+y+y^2)*z^3 + (y+y^2)*z^2 + (y+y^2)*z + y^2 + y", sym, 2);
	// Poly_ptr wpCase = generate_altarr_var_defined("(6+y+y^2)*z^6 + (5+y+y^2)*z^5 + (4+y+y^2)*z^4 + (3+y+y^2)*z^3 + (2+y+y^2)*z^2 + (y+y^2)*z + y^2 + y", sym, 2);
	// Poly_ptr wpCase = generate_altarr_var_defined("(6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2", sym, 3);
	// Poly_ptr wpCase = generate_altarr_var_defined("(8+y+x^2)*z^8 + (7+y+x^2)*z^7 + (6+y+x^2)*z^6 + (5+y+x^2)*z^5 + (4+y+x^2)*z^4 + (3+y+x^2)*z^3 + (2+y+x^2)*z^2 + (y+x^2)*z + x^2 + x*y + y^2", sym, 3);
	// Poly_ptr wpCase = generate_altarr_var_defined("z^5+z^4+(4*y-2)*z^3+(6*y^2-6*y+1)*z^2+(4*y^3-6*y^2+2*y)*z+y^4-2*y^3+y^2", sym, 3);

    if (family == 0) {
    	wpCase = generate_altarr_var_defined(familyA[deg].c_str(), sym, 3);
    } else {
    	wpCase = generate_altarr_var_defined(familyB[deg].c_str(), sym, 3);
    }

	Upops_t* upop = convertPolyToUnivariatePolyOverPowerSeries_UPOPS(wpCase);

	print_UPOPS(stderr, upop, sym);

    //update input fully first
    updateToDeg_UPOPS(testDeg, upop);

    Upops_t* p_out = NULL; Upops_t* alpha_out = NULL;
    weierstrassPreparation_UPOPS(upop, &p_out, &alpha_out);
    fprintf(stderr,"p->deg: %d, a->deg: %d\n", p_out->deg, alpha_out->deg);

    unsigned long long start;
    _startTimer(&start);
	weierstrassUpdateParallel_UPOPS(testDeg, p_out, nthreads);
	// updateToDeg_UPOPS(testDeg, p_out);
	float elapsed;
	_stopTimer(&start, &elapsed);

    // print_UPOPS(stderr,  p_out, sym);
    // fprintf(stderr, "\n");

    // print_UPOPS(stderr,  alpha_out, sym);
    // fprintf(stderr, "\n");


	fprintf(stderr, "Time: %10f\n", elapsed);

    // TestingWeierstrassPreparation(upop, testDeg, sym);

    return 0;
}
