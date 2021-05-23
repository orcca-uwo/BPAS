#include "PowerSeries/PowerSeriesZ.h"
#include "PowerSeries/UPOPS_Z.h"



#define MAX(x, y) ((x) > (y) ? (x) : (y))

/**
 * alloc : the size of the array of power series associated with the Upops
 * allocate a Upops
 */
UpopsZ_t* allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(int alloc) {
    UpopsZ_t* upops = (UpopsZ_t*) malloc(sizeof(UpopsZ_t));
    upops->refCount = 1;
    if (alloc > 0) {
    	upops->data = (PowerSeriesZ_t**) malloc(sizeof(PowerSeriesZ_t*)*alloc);
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
void destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(UpopsZ_t* upops) {

    --(upops->refCount);

    if (upops->refCount <= 0 ) {

    	//If this is the p from weierstrass
        if (upops->weierstrassFData != NULL) {
        	for (int i = 0; i < upops->fDataSize; ++i) {
        		destroyPowerSeries_PSZ(upops->weierstrassFData[i]);
        	}
        	free(upops->weierstrassFData);

	        for (int i = 0; i < upops->deg; ++i) {
	        	if (upops->data[i]->refCount > 1) {
	        		upops->data[i]->genOrder = -1;
	        		upops->data[i]->genParam1 = NULL;
	        		upops->data[i]->genParam2 = NULL;
	        		destroyUnivariatePolynomialOverPowerSeries_UPOPSZ((UpopsZ_t*) (upops->data[i]->genParam3));
	        		upops->data[i]->genParam3 = NULL;
	        		upops->data[i]->paramType1 = PLAIN_DATA;
	        		upops->data[i]->paramType2 = PLAIN_DATA;
	        		upops->data[i]->paramType3 = PLAIN_DATA;
	        		upops->data[i]->gen.tertiaryGen = NULL;
	        	}
	        	destroyPowerSeries_PSZ(upops->data[i]);
	        }
        	destroyPowerSeries_PSZ(upops->data[upops->deg]);

        } else {
	        for (int i = 0; i <= upops->deg; ++i) {
	        	destroyPowerSeries_PSZ(upops->data[i]);
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
void reserve_UPOPSZ(UpopsZ_t* upops) {
    ++(upops->refCount);
}


/**
 * Construct zero Upops
 */
UpopsZ_t* zero_UPOPSZ() {

	UpopsZ_t* zero = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(1);
	zero->data[0] = zeroPowerSeries_PSZ();
	zero->deg = 0;
	return zero;
}


/**
 * Construct one Upops
 */
UpopsZ_t* one_UPOPSZ(int nvar) {

	UpopsZ_t* one = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(1);
	one->data[0] = onePowerSeries_PSZ(nvar);
	one->deg = 0;
	return one;
}


void  updateToDeg_UPOPSZ(int d, UpopsZ_t* f) {

	int deg = f->deg;
	for (int i = 0; i <= deg; ++i) {
	     if (d + 1 > f->data[i]->alloc) {
		     int newAlloc = (2*(f->data[i])->alloc < d+1) ? d + 1 : 2*f->data[i]->alloc;
			 f->data[i]->polys = (PolyZ_ptr*) realloc(f->data[i]->polys, sizeof(PolyZ_ptr)*newAlloc);
			 f->data[i]->alloc = newAlloc;
		}
	}

	for (int j = 0; j <= deg; ++j) {
		updateToDeg_PSZ(d, f->data[j]);
	}

}


/**
 * upops : a Upops
 * sym : a list of variables
 * print a Upops
 */
void print_UPOPSZ(FILE* fp, UpopsZ_t* upops, const char** sym) {

    int First = 1;
    for (int i = 0; i <= upops->deg; ++i) {

        if (!isZero_PSZ(upops->data[i])) {
            if (!First) {
                fprintf(fp, " + (");
                print_PSZ(fp, upops->data[i], sym);
                fprintf(fp, ")*%s^%d", sym[0], i);
              } else {
                 print_PSZ(fp, upops->data[i], sym);
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
UpopsZ_t* convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ(const PolyZ_ptr p) {

	//const char* syms[] = {"y", "z", "x", "t"};
	//AltArrZ_t*  b = generate_altarr_var_defined("0", syms, 4);	//ToDo : how to convert NULL to a zero polynomial without knowing the variables
	degree_t coefflistsize;
	if (isConstant_AAZ(p)) {
		coefflistsize = 0;
	} else {
		coefflistsize = partialDegree_AAZ(p, 0);
	}

	UpopsZ_t* arrayofps = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(coefflistsize+1);
	int i;
	AltArrZ_t** cList;
	int sz = 0;

	if (isZero_AAZ(p))  {                                /*it is not always equal to 1 but 5, 7,...we need a function that creates a constant and converts it to a power series*/
		arrayofps->data[0] = zeroPowerSeries_PSZ();
	} else if (isConstant_AAZ(p)){
		arrayofps->data[0] = constPowerSeries_PSZ(p->elems->coef, p->nvar);
	} else {
	    mainCoefficientListAtIdx_AAZ(p, 0, &cList, &sz);

	    for (i = 0; i <=coefflistsize; ++i) {

		    if (cList[i]== NULL){

			//arrayofps->data[i] = convertPolyToPowerSeries_PSZ(b);
		   	    arrayofps->data[i] = zeroPowerSeries_PSZ();

		    } else {

                arrayofps->data[i] = convertPolyToPowerSeries_PSZ(cList[i]);
                freePolynomial_AAZ(cList[i]);
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
 * @see convertPolyToUnivariatePolyOverPowerSeries_UPOPSZ
 *
 * @param ps : an array of power series
 * @param size : the size of the array
 * @return a upops whose coefficients are that of the array ps.
 */
UpopsZ_t* convertArrayOfPSToUPOPS_UPOPSZ(PowerSeriesZ_t** ps, int size) {

	UpopsZ_t* upops = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(size);
	upops->deg = size - 1;
	for (int i = 0; i <= size - 1; ++i) {
		reserve_PSZ(ps[i]);
		upops->data[i] = ps[i];
	}
	return upops;
}


PolyZ_ptr polynomialPart_UPOPSZ(int d, UpopsZ_t* f) {
	if ( f == NULL) {
		return NULL;
	}

	PolyZ_ptr finalPolyPart = NULL;
	PolyZ_ptr powerSeriesPolyPart;

	for (int i = 0; i <= f->deg; ++i) {
		powerSeriesPolyPart = polynomialPart_PSZ(d, f->data[i]);
		for (int k = 0; powerSeriesPolyPart != NULL && k < powerSeriesPolyPart->size; ++k) {
			//set term k's 0th variable to degree i
			setPartialDegreeTerm_AAZ(powerSeriesPolyPart, k, 0, i);
		}

		finalPolyPart = addPolynomials_AAZ_inp(finalPolyPart, powerSeriesPolyPart, powerSeriesPolyPart == NULL ? 0 : powerSeriesPolyPart->nvar);
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
UpopsZ_t*  addUnivariatePolyOverPowerSeries_UPOPSZ(UpopsZ_t* f, UpopsZ_t* g) {

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

	UpopsZ_t* sum = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(e+1);
	sum->deg = e;

	for (int i = 0; i <= d; ++i) {
		sum->data[i] = addPowerSeries_PSZ(f->data[i], g->data[i]);

	}

	if (fDeg > gDeg){
		for (int j = d + 1; j <= fDeg; ++j) {
			reserve_PSZ(f->data[j]);
			sum->data[j] = f->data[j];
		}

	} else if (fDeg < gDeg) {
		for (int k = d + 1; k <= gDeg; ++k) {
			reserve_PSZ(g->data[k]);
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
UpopsZ_t*  subUnivariatePolyOverPowerSeries_UPOPSZ(UpopsZ_t* f, UpopsZ_t* g) {


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


	UpopsZ_t* sub = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(e+1);

	sub->deg = e;

	for (int i = 0; i <= d; ++i) {
		sub->data[i] = subPowerSeries_PSZ(f->data[i], g->data[i]);
	}

	if (fDeg > gDeg){
		for (int j = d + 1; j <= fDeg; ++j) {
			reserve_PSZ(f->data[j]);
			sub->data[j] = f->data[j];
		}

	} else if (fDeg < gDeg) {
		for (int k = d + 1; k <= gDeg; ++k) {
			sub->data[k] = negatePowerSeries_PSZ(g->data[k]);
		}
	}
	return sub;

}


/**
 * f : a Upops
 * g : a Upops
 * multiply f by g
 */
UpopsZ_t* multiplyUnivariatePolyOverPowerSeries_UPOPSZ(UpopsZ_t* f, UpopsZ_t* g) {

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

    UpopsZ_t* prod = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(prodDeg+1);

    prod->deg = prodDeg;

	PowerSeriesZ_t* middleProd;
	for (int m = 0; m <= prodDeg; ++m) {
		prod->data[m] = NULL;
	}
	for (int j = 0; j <= fDeg; ++j) {

		for (int k=0; k <= gDeg; ++k) {

			middleProd = multiplyPowerSeries_PSZ(f->data[j], g->data[k]);

			if (prod->data[j+k] == NULL) {
				prod->data[j+k] = middleProd;
			} else {
				PowerSeriesZ_t* temp = addPowerSeries_PSZ(prod->data[j+k], middleProd);
				destroyPowerSeries_PSZ(middleProd);
				destroyPowerSeries_PSZ(prod->data[j+k]);
				prod->data[j+k] = temp;
			}

		}
	}
	return prod;
}

