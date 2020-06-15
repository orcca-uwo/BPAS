#include "PowerSeries/PowerSeries.h"
#include "PowerSeries/UnivariatePolynomialOverPowerSeries.h"
#include "PowerSeries/UPOPS_Weierstrass.h"

/**
 * requestedDeg : the requested degree
 * coefDeg : the index of the power series in data array
 * p : a upops
 * update all coefficients (power series) in p and alpha up to degree requestedDeg
 * and returns the updated power series at index coefDeg in p
 */
Poly_ptr weierstrassCoefPUpdate_UPOPS(int requestedDeg, long long coefDeg, Upops_t* p, Upops_t* alpha) {

	if (p == NULL) {
		return NULL;
	}
	int needUpdate = 0;
	int curdP = -1;
	for (int i = 0; i < p->deg; ++i) {
		if (p->data[i]->deg < requestedDeg) {
			needUpdate = 1;
			curdP = curdP == -1 ? p->data[i]->deg : MIN(curdP, p->data[i]->deg);
		}
	}
	if (needUpdate) {
		for (int j = curdP + 1; j <= requestedDeg; ++j) {
			weierstrassUpdate_UPOPS(p, alpha);
		}
	}
	return p->data[coefDeg]->polys[requestedDeg];
}


/**
 *
 */
Poly_ptr weierstrassCoefVoidPUpdate_UPOPS(int i, void* coefDeg, void* p, void* alpha) {
	return weierstrassCoefPUpdate_UPOPS(i, (long long) coefDeg, (Upops_t*) p, (Upops_t*) alpha);
}


// This was the original for updateing the upops alpha...
// but it's not really needed since alpha is designed to updated
// lazily based on PS arithmetic. This, eventually,
// we call one of p's coefficient's generators to trigger
// weierstrassUpdate through weierstrassCoefVoidPUpdate_UPOPS
#if 0
/**
 * requestedDeg : the requested degree
 * coefDeg : the index of the power series in data array
 * alpha : a upops
 * update all coefficients (power series) in p and alpha up to degree requestedDeg
 * and returns the updated power series at index coefDeg in alpha
 */
Poly_ptr weierstrassCoefAUpdate_UPOPS(int requestedDeg, long long coefDeg, Upops_t* alpha) {

	if (alpha == NULL) {
		return NULL;
	}
	int needUpdate = 0;
	int curdA = -1;
	for (int i = 0; i < alpha->deg; ++i) {
		if (alpha->data[i]->deg < requestedDeg) {
			needUpdate = 1;
			curdA = curdA == -1 ? alpha->data[i]->deg : MIN(curdA, alpha->data[i]->deg);
		}
	}
	if (needUpdate) {
		for (int j = curdA + 1; j <= requestedDeg; ++j) {
			weierstrassUpdate_UPOPS((Upops_t*) alpha->genParam1, (Upops_t*) alpha->genParam2, alpha);
		}
	}
	return alpha->data[coefDeg]->polys[requestedDeg];
}


/**
 *
 */
Poly_ptr weierstrassCoefVoidAUpdate_UPOPS(int i, void* coefDeg, void* alpha) {
	return weierstrassCoefAUpdate_UPOPS(i, (long long) coefDeg, (Upops_t*) alpha);
}
#endif

/**
 * F, G, H : power series
 * r : an integer
 * nvar : the number of variables
 *
 */
Poly_ptr lemmaForWeiestrass_UPOPS(PowerSeries_t* F, PowerSeries_t* G, PowerSeries_t* H, int r) {
	if (F == NULL || G == NULL || H == NULL) {
		return NULL;
	}

	int nvar = F->polys[0] == NULL ? H->polys[0]->nvar : F->polys[0]->nvar;

	Poly_ptr s = NULL;
	// fprintf(stderr, "lemma for r=%d\n", r);
	for (int i = 1; i < r; ++i) {

		//Here we need homogPart, because sometimes alpha or p may be truncated.
	    Poly_ptr multiply = multiplyPolynomials_AA(homogPart_PS(r-i, G), homogPart_PS(i, H), nvar);
	    // const char* syms[] = {"z", "x", "y"};
	    // fprintf(stderr, "multiply: ");
	    // printPoly_AA(stderr, multiply, syms, multiply == NULL ? 0 : multiply->nvar);
	    // fprintf(stderr, "\n ");
	    s = addPolynomials_AA_inp(s,  multiply, nvar);
	    freePolynomial_AA(multiply);
	}

	if(r > F->deg) {
		if (!isZero_AA(s)) {
			negatePolynomial_AA(s);
		}
	} else {
		s = subPolynomials_AA_inp(s, F->polys[r], nvar);
		negatePolynomial_AA(s);
	}

	divideByRational_AA_inp(s, H->polys[0]->elems[0].coef);

	return s;
}



/**
 * upops : the input upops
 * p : the upops p
 * alpha : the upops alpha
 * update all power series of p and alpha
 */
void weierstrassUpdate_UPOPS(Upops_t* p, Upops_t* alpha) {

	int curdP = -1;
	int curdA = -1;
	for (int i = 0; i < p->deg; ++i) {
		curdP = curdP == -1 ? p->data[i]->deg : MIN(curdP, p->data[i]->deg);
	}
	for (int i = 0; i < alpha->deg; ++i) {
		curdA = curdA == -1 ? alpha->data[i]->deg : MIN(curdA, alpha->data[i]->deg);
	}

	//we look to increment each by one herein.
	curdP += 1;
	curdA += 1;

	int d = p->deg;
	int m = alpha->deg;
	Poly_ptr polyP;


	PowerSeries_t** F = p->weierstrassFData;

    for (int lP = 0; lP <= d-1; ++lP) {
		// fprintf(stderr, "curdP %d, alloc: %d\n\n", curdP, p->data[lP]->alloc);
     	if (curdP + 1 > p->data[lP]->alloc) {
		    int newAlloc = (2*(p->data[lP])->alloc < curdP+1) ? curdP + 1 : 2*p->data[lP]->alloc;
		    p->data[lP]->polys = (Poly_ptr*) realloc(p->data[lP]->polys, sizeof(Poly_ptr)*newAlloc);
			p->data[lP]->alloc = newAlloc;
		}
	}
	//alpha's allocation automatically managed by its power series structure.


	for (int lP = 0; lP <= d - 1; ++lP) {
		// fprintf(stderr, "before degree of p->data[%d]: %d\n", lP, p->data[lP]->deg);
		// fprintf(stderr, "before degree of alpha->Data[0]: %d\n", alpha->data[0]->deg);
		// fprintf(stderr, "before degree of F[%d]: %d\n", lP, F[lP]->deg);

		if (F[lP]->deg == -1){
			 polyP = NULL;
		}else{
		   updateToDeg_PS(curdP,F[lP]); //update F based on new terms of p->data[lP-1]
			// fprintf(stderr, "after degree of F[%d]: %d\n", lP, F[lP]->deg);
		   polyP = lemmaForWeiestrass_UPOPS(F[lP], p->data[lP], alpha->data[0], curdP);
		}

		// const char* syms[] = {"z", "x", "y"};
		// fprintf(stderr, "\n\ncomputed lemma poly for degree %d:\n", curdP);
		// printPoly_AA(stderr, polyP, syms, polyP == NULL ? 0 : polyP->nvar);
		// fprintf(stderr, "\n" );
		p->data[lP]->polys[curdP] = polyP;
		//manually set new degree since this updated PS is needed for next iteration of the loop.
		p->data[lP]->deg = curdP;
		// fprintf(stderr, "after degree of p->data[%d]: %d\n", lP, p->data[lP]->deg);
	}

	for (int lA = m; lA >= 0; --lA) {
		// fprintf(stderr, "before degree of alpha->data[%d]: %d\n", lA, alpha->data[lA]->deg);
		updateToDeg_PS(curdA,  alpha->data[lA]);
		// fprintf(stderr, "after degree of alpha->data[%d]: %d\n", lA, alpha->data[lA]->deg);
	}
}



/**
 * upops : the input upops
 * p_out : the output p
 * alpha_out : the output alpha
 * It takes a upops as input and computes p and alpha
 */
void weierstrassPreparation_UPOPS(Upops_t* upops, Upops_t** p_out, Upops_t** alpha_out) {

	/* Check if the input is valid or not */
	int n = upops->deg;
	PowerSeries_t** a = upops->data;
	int d = -1;
	for (int i = 0; i <= n; i++) {
		if (isUnit_PS(a[i])) {
			d = i;
			break;
		}
	}

	if (d == -1) {
		fprintf(stderr, "ERROR in weierstrassPreparation_UPOPS: Invalid input!\n");
		exit(1);
	}

	int nvar = upops->data[d]->polys[0] == NULL ? 0 : upops->data[d]->polys[0]->nvar;
	if (d == 0) {
		if (p_out != NULL) {
			*p_out = one_UPOPS(nvar);
		}
		if (alpha_out != NULL) {
			reserve_UPOPS(upops);
			*alpha_out = upops;
		}
		return;
	}

	long long m = n - d;

	Upops_t* p = allocateUnivariatePolynomialOverPowerSeries_UPOPS(d+1);
	Upops_t* alpha = allocateUnivariatePolynomialOverPowerSeries_UPOPS(m+1);

	PowerSeries_t** b = p->data;
	PowerSeries_t** c = alpha->data;
	PowerSeries_t** F = (PowerSeries_t**) malloc(sizeof(PowerSeries_t*) * (d));

	/*  Initializing b  */

    for (int j = 0; j <= d - 1; j++) {
		b[j]= allocatePowerSeries_PS(0);
		b[j]->polys= copyUpTo_PS(upops->data[j], 0);
		b[j]->alloc = 1;
		b[j]->deg = 0;

	}

    b[d] = onePowerSeries_PS(nvar);


	/*  Initializing c */
    reserve_PS(a[d+m]);
    c[m] = a[d+m];
	for (int lA = m-1; lA >= 0; --lA) {

		c[lA] = upops->data[d + lA];
		reserve_PS(c[lA]);

		for (int jA = d-1; jA >= 0; --jA) {
		    int k = d + lA - jA;

			if (k <= m) {

		        PowerSeries_t* mulpsA = multiplyPowerSeries_PS(b[jA], c[k]);
				PowerSeries_t* tmp = subPowerSeries_PS(c[lA], mulpsA);
				destroyPowerSeries_PS(c[lA]); //c[lA] is reserved inside subPowerSeries_PS
				destroyPowerSeries_PS(mulpsA); //mulpsA is reserved inside subPowerSeries_PS
				c[lA] = tmp;
			}
	    }
	}

	/*  Initializing F */
	reserve_PS(upops->data[0]);
	F[0] = upops->data[0];
	for (int lP = 1; lP <= d - 1; ++lP) {
	    F[lP] = upops->data[lP];
	    reserve_PS(F[lP]);
	    if (m > 0) {
	        for (int jP = 1; jP <= lP; ++jP) {
	    	    PowerSeries_t* mulpsP = multiplyPowerSeries_PS(b[lP - jP], c[jP]);
	    	    PowerSeries_t* tmp = subPowerSeries_PS(F[lP], mulpsP);
	    	    destroyPowerSeries_PS(F[lP]);
	    	    destroyPowerSeries_PS(mulpsP); //mulpsP is reserved inside sub
	    		F[lP] = tmp;
	    	}
	    }
	}
	p->weierstrassFData = F;
	p->fDataSize = d;

	/*  p */
	p->deg = d;
	/*  generators of each b[i] in p */
	for (long long int k = 0; k < d; k++) {
	    b[k]->genOrder = 3;
	    b[k]->gen.tertiaryGen = &(weierstrassCoefVoidPUpdate_UPOPS);
	    b[k]->genParam1 = (void*) k; //the index of the coefficient
	    b[k]->genParam2 = (void*) p; //yes, a bit recrusive in nature. But it should be fine with the correct number of reserves and destroys.
	   	b[k]->genParam3 = (void*) alpha;
	    // reserve_UPOPS(p); // use a weak reference to p
	    reserve_UPOPS(alpha);
		b[k]->paramType1 = PLAIN_DATA;
		b[k]->paramType2 = WEAK_UPOPS;
		b[k]->paramType3 = UPOPS;
	}

	// p->BinaryGen = &(weierstrassCoefVoidPUpdate_UPOPS);
	// p->genParam1 = upops;
 //    p->genParam2 = alpha;
 //    p->paramType1 = UPOPS;
 //    p->paramType2 = UPOPS;
 //    reserve_UPOPS(upops);
 //    reserve_UPOPS(alpha);

	/*  alpha */
	alpha->deg = m;

	// alpha->BinaryGen = &(weierstrassCoefVoidAUpdate_UPOPS);
	// alpha->genParam1 = upops;
	// alpha->genParam2 = p;
    // alpha->paramType1 = UPOPS;
    // alpha->paramType2 = UPOPS;
    // reserve_UPOPS(upops);
    // reserve_UPOPS(p);


	if (p_out != NULL) {
		*p_out = p;
	} else {
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(p);
	}

	if (alpha_out != NULL) {
		*alpha_out = alpha;
	} else {
		destroyUnivariatePolynomialOverPowerSeries_UPOPS(alpha);
	}

}


