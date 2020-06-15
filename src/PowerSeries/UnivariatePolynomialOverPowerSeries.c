#include "PowerSeries/PowerSeries.h"
#include "PowerSeries/UnivariatePolynomialOverPowerSeries.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))

/**
 * alloc : the size of the array of power series associated with the Upops
 * allocate a Upops
 */
Upops_t* allocateUnivariatePolynomialOverPowerSeries_UPOPS(int alloc) {
    Upops_t* upops = (Upops_t*) malloc(sizeof(Upops_t));
    upops->refCount = 1;
    if (alloc > 0) {
    	upops->data = (PowerSeries_t**) malloc(sizeof(PowerSeries_t*)*alloc);
    	upops->alloc = alloc;
    } else {
    	upops->data = NULL;
    	upops->alloc = 0;
    }
    upops->weierstrassFData = NULL;
    upops->fDataSize = 0;
    upops->deg = -1;
    return upops;
}


/**
 * upops : a given univariate polynomial over power series
 * destroy a Upops
 */
void destroyUnivariatePolynomialOverPowerSeries_UPOPS(Upops_t* upops) {

    --(upops->refCount);

    if (upops->refCount <= 0 ) {

    	//If this is the p from weierstrass
        if (upops->weierstrassFData != NULL) {
        	for (int i = 0; i < upops->fDataSize; ++i) {
        		destroyPowerSeries_PS(upops->weierstrassFData[i]);
        	}
        	free(upops->weierstrassFData);

	        for (int i = 0; i < upops->deg; ++i) {
	        	if (upops->data[i]->refCount > 1) {
	        		upops->data[i]->genOrder = -1;
	        		upops->data[i]->genParam1 = NULL;
	        		upops->data[i]->genParam2 = NULL;
	        		destroyUnivariatePolynomialOverPowerSeries_UPOPS((Upops_t*) (upops->data[i]->genParam3));
	        		upops->data[i]->genParam3 = NULL;
	        		upops->data[i]->paramType1 = PLAIN_DATA;
	        		upops->data[i]->paramType2 = PLAIN_DATA;
	        		upops->data[i]->paramType3 = PLAIN_DATA;
	        		upops->data[i]->gen.tertiaryGen = NULL;
	        	}
	        	destroyPowerSeries_PS(upops->data[i]);
	        }
        	destroyPowerSeries_PS(upops->data[upops->deg]);

        } else {
	        for (int i = 0; i <= upops->deg; ++i) {
	        	destroyPowerSeries_PS(upops->data[i]);
	        }
        }


        free(upops->data);
        free(upops);
    }
}


/**
 * upops : a given univariate polynomial over power series
 * reserve a Upops
 */
void reserve_UPOPS(Upops_t* upops) {
    ++(upops->refCount);
}


/**
 * Construct zero Upops
 */
Upops_t* zero_UPOPS() {

	Upops_t* zero = allocateUnivariatePolynomialOverPowerSeries_UPOPS(1);
	zero->data[0] = zeroPowerSeries_PS();
	zero->deg = 0;
	return zero;
}


/**
 * Construct one Upops
 */
Upops_t* one_UPOPS(int nvar) {

	Upops_t* one = allocateUnivariatePolynomialOverPowerSeries_UPOPS(1);
	one->data[0] = onePowerSeries_PS(nvar);
	one->deg = 0;
	return one;
}


void  updateToDeg_UPOPS(int d, Upops_t* f) {

	int deg = f->deg;
	for (int i = 0; i <= deg; ++i) {
	     if (d + 1 > f->data[i]->alloc) {
		     int newAlloc = (2*(f->data[i])->alloc < d+1) ? d + 1 : 2*f->data[i]->alloc;
			 f->data[i]->polys = (Poly_ptr*) realloc(f->data[i]->polys, sizeof(Poly_ptr)*newAlloc);
			 f->data[i]->alloc = newAlloc;
		}
	}

	for (int j = 0; j <= deg; ++j) {
		updateToDeg_PS(d, f->data[j]);
	}

}


/**
 * upops : a Upops
 * sym : a list of variables
 * print a Upops
 */
void print_UPOPS(FILE* fp, Upops_t* upops, const char** sym) {

    int First = 1;
    for (int i = 0; i <= upops->deg; ++i) {

        if (!isZero_PS(upops->data[i])) {
            if (!First) {
                fprintf(fp, " + (");
                print_PS(fp, upops->data[i], sym);
                fprintf(fp, ")*%s^%d", sym[0], i);
              } else {
                 print_PS(fp, upops->data[i], sym);
                 First = 0;
             }
        }
    }
}


/********************
 * Conversion Helpers
 ********************/


/**
 * Convert a polynomial to a univariate polynomial with power series coefficients.
 * The variable of index 0 in p is used as the univariate polynomial variable.
 * Then, the polynomials in the power series coefficients remain
 * to be defined in their original polynomial ring, but have degree 0 in their main variable.
 * @param p : a polynomial to convert.
 * @return p converted to a upops.
 */
Upops_t* convertPolyToUnivariatePolyOverPowerSeries_UPOPS(const Poly_ptr p) {

	//const char* syms[] = {"y", "z", "x", "t"};
	//AltArr_t*  b = generate_altarr_var_defined("0", syms, 4);	//ToDo : how to convert NULL to a zero polynomial without knowing the variables
	degree_t coefflistsize;
	if (isConstant_AA(p)) {
		coefflistsize = 0;
	} else {
		coefflistsize = partialDegree_AA(p, 0);
	}

	Upops_t* arrayofps = allocateUnivariatePolynomialOverPowerSeries_UPOPS(coefflistsize+1);
	int i;
	AltArr_t** cList;
	int sz = 0;

	if (isZero_AA(p))  {                                /*it is not always equal to 1 but 5, 7,...we need a function that creates a constant and converts it to a power series*/
		arrayofps->data[0] = zeroPowerSeries_PS();
	} else if (isConstant_AA(p)){
		arrayofps->data[0] = constPowerSeries_PS(p->elems->coef, p->nvar);
	} else {
	    mainCoefficientListAtIdx_AA(p, 0, &cList, &sz);

	    for (i = 0; i <=coefflistsize; ++i) {

		    if (cList[i]== NULL){

			//arrayofps->data[i] = convertPolyToPowerSeries_PS(b);
		   	    arrayofps->data[i] = zeroPowerSeries_PS();

		    } else {

                arrayofps->data[i] = convertPolyToPowerSeries_PS(cList[i]);
                freePolynomial_AA(cList[i]);
		    }

	    }
	    free(cList);
	}

	//arrayofps->deg = coefflistsize - 1;//updated on March 3
    arrayofps->deg = coefflistsize;
	return arrayofps;
}



/**
 * Convert an array of power series to a univariate polynomial over power series (upops).
 * Where the index of the power series implies its associated monomial's degree.
 * @note the power series should be defined to exist in Q[[Y,X_1,...,X_n]]
 * but have degree 0 in Y, the eventual polynomial variable of the upops.
 *
 * @see convertPolyToUnivariatePolyOverPowerSeries_UPOPS
 *
 * @param ps : an array of power series
 * @param size : the size of the array
 * @return a upops whose coefficients are that of the array ps.
 */
Upops_t* convertArrayOfPSToUPOPS_UPOPS(PowerSeries_t** ps, int size) {

	Upops_t* upops = allocateUnivariatePolynomialOverPowerSeries_UPOPS(size);
	upops->deg = size - 1;
	for (int i = 0; i <= size - 1; ++i) {
		reserve_PS(ps[i]);
		upops->data[i] = ps[i];
	}
	return upops;
}


Poly_ptr polynomialPart_UPOPS(int d, Upops_t* f) {
	if ( f == NULL) {
		return NULL;
	}

	Poly_ptr finalPolyPart = NULL;
	Poly_ptr powerSeriesPolyPart;

	for (int i = 0; i <= f->deg; ++i) {
		powerSeriesPolyPart = polynomialPart_PS(d, f->data[i]);
		for (int k = 0; powerSeriesPolyPart != NULL && k < powerSeriesPolyPart->size; ++k) {
			//set term k's 0th variable to degree i
			setPartialDegreeTerm_AA(powerSeriesPolyPart, k, 0, i);
		}

		finalPolyPart = addPolynomials_AA_inp(finalPolyPart, powerSeriesPolyPart, powerSeriesPolyPart == NULL ? 0 : powerSeriesPolyPart->nvar);
	}

	return finalPolyPart;
}



/***********************
 * UPOPS Arithmetic
 ***********************/

/**
 * f : a Upops
 * g : a Upops
 * add f and g
 */
Upops_t*  addUnivariatePolyOverPowerSeries_UPOPS(Upops_t* f, Upops_t* g) {

	int fDeg = f->deg;
	int gDeg = g->deg;
	int d;
	int e;
	if (gDeg == -1 && fDeg == -1) {
	   d = 0;
	   e = 0;
	} else if (gDeg == -1) {
	    d = 0;
	    e = fDeg;
	} else if (fDeg == -1){
		d = 0;
		e = gDeg;
	} else {
	    d = MIN(fDeg, gDeg);
	    e = MAX(fDeg, gDeg);
	}

	Upops_t* sum = allocateUnivariatePolynomialOverPowerSeries_UPOPS(e+1);
	sum->deg = e;

	for (int i = 0; i <= d; ++i) {
		sum->data[i] = addPowerSeries_PS(f->data[i], g->data[i]);

	}

	if (fDeg > gDeg){
		for (int j = d + 1; j <= fDeg; ++j) {
			reserve_PS(f->data[j]);
			sum->data[j] = f->data[j];
		}

	} else if (fDeg < gDeg) {
		for (int k = d + 1; k <= gDeg; ++k) {
			reserve_PS(g->data[k]);
			sum->data[k] = g->data[k];
		}
	}
	return sum;
}


/**
 * f : a Upops
 * g : a Upops
 * subtract g from  f
 */
Upops_t*  subUnivariatePolyOverPowerSeries_UPOPS(Upops_t* f, Upops_t* g) {


	int fDeg = f->deg;
	int gDeg = g->deg;
	int d;
	int e;
	if (gDeg == -1 && fDeg == -1) {
        d = 0;
		e = 0;
	} else if (gDeg == -1) {
		d = 0;
		e = fDeg;
	} else if (fDeg == -1){
		d = 0;
		e = gDeg;
	} else {
	    d = MIN(fDeg, gDeg);
		e = MAX(fDeg, gDeg);
	}


	Upops_t* sub = allocateUnivariatePolynomialOverPowerSeries_UPOPS(e+1);

	sub->deg = e;

	for (int i = 0; i <= d; ++i) {
		sub->data[i] = subPowerSeries_PS(f->data[i], g->data[i]);
	}

	if (fDeg > gDeg){
		for (int j = d + 1; j <= fDeg; ++j) {
			reserve_PS(f->data[j]);
			sub->data[j] = f->data[j];
		}

	} else if (fDeg < gDeg) {
		for (int k = d + 1; k <= gDeg; ++k) {
			sub->data[k] = negatePowerSeries_PS(g->data[k]);
		}
	}
	return sub;

}


/**
 * f : a Upops
 * g : a Upops
 * multiply f by g
 */
Upops_t* multiplyUnivariatePolyOverPowerSeries_UPOPS(Upops_t* f, Upops_t* g) {

	int fDeg = f->deg;
	int gDeg = g->deg;
	if (fDeg == -1 && gDeg == -1) {
		fDeg = 0;
		gDeg = 0;
	} else if (fDeg == -1) {
		fDeg = 0;
	} else if (gDeg == -1) {
		gDeg = 0;
	}

    int prodDeg = fDeg + gDeg;

    Upops_t* prod = allocateUnivariatePolynomialOverPowerSeries_UPOPS(prodDeg+1);

    prod->deg = prodDeg;

	PowerSeries_t* middleProd;
	for (int m = 0; m <= prodDeg; ++m) {
		prod->data[m] = NULL;
	}
	for (int j = 0; j <= fDeg; ++j) {

		for (int k=0; k <= gDeg; ++k) {

			middleProd = multiplyPowerSeries_PS(f->data[j], g->data[k]);

			if (prod->data[j+k] == NULL) {
				prod->data[j+k] = middleProd;
			} else {
				PowerSeries_t* temp = addPowerSeries_PS(prod->data[j+k], middleProd);
				destroyPowerSeries_PS(middleProd);
				destroyPowerSeries_PS(prod->data[j+k]);
				prod->data[j+k] = temp;
			}

		}
	}
	return prod;
}



/** ToDo: it is not an efficient code
 *
 */

//Poly_ptr sumMonoms(int d, const char** sym) {

//	int sizeVars = sizeof(sym)/sizeof(sym[0]);
//	Poly_ptr finalPoly;
//	Poly_ptr newMonom;
//	Poly_ptr multiplyResults = generate_altarr_var_defined("1",sym, sizeVars);

//	for (int i = 0; i <= sizeVars; ++i) {
//		finalPoly = generate_altarr_var_defined("1",sym, sizeVars);
//		for (int j = 1; j <= d; ++j) {
//			newMonom = generate_altarr_var_defined("sym[i]^j",sym, sizeVars);
//		    finalPoly = addPolynomials_AA_inp(finalPoly,newMonom,sizeVars);
//		}
//		multiplyResults = multiplyPolynomials_AA(multiplyResults, finalPoly,sizeVars);
//	}
//	return multiplyResults;
//}


/** ToDo: it is not an efficient code
 *
 */
//Poly_ptr homogeneousPart_sum_of_monoms_PS(int d,  const char** sym) {

//	Poly_ptr sumWithExtraMonoms =  sumMonoms(d, sym);
//	Poly_ptr sumOfRequestedDeg =  truncatePolynomial_PS(sumWithExtraMonoms, d);
//    return sumOfRequestedDeg;
//}


/** ToDo: it is not an efficient code
 *
 */
//Poly_ptr homogeneousPartVoid_sum_of_monoms_PS(int deg,  const char** sym) {
//        return homogeneousPart_sum_of_monoms_PS(deg,  sym);
//}



/** ToDo : it is not an efficient code
 * Given a list of variables it generates a series of monomials up to d
 */

//PowerSeries_t*  sumOfAllMonomials(const char** sym) {

//	PowerSeries_t* sum_of_monoms = allocatePowerSeries_PS(1);
//	sum_of_monoms->deg = 0;
//	sum_of_monoms->polys[0] = makeConstIntPolynomial_AA(1,sizeof(sym)/sizeof(sym[0]), 1, 1);
//	sum_of_monoms->genOrder = 1;
//	sum_of_monoms->genParam1 = sym;
//	sum_of_monoms->gen.unaryGen = &(homogeneousPartVoid_sum_of_monoms_PS);
//	return sum_of_monoms;
//}








