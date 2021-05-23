#include "PowerSeries/PowerSeriesZ.h"
#include "PowerSeries/UPOPS_Z.h"
#include "PowerSeries/UPOPS_Z_Weierstrass.h"

/**
 * requestedDeg : the requested degree
 * coefDeg : the index of the power series in data array
 * p : a upops
 * update all coefficients (power series) in p and alpha up to degree requestedDeg
 * and returns the updated power series at index coefDeg in p
 */
PolyZ_ptr weierstrassCoefPUpdate_UPOPSZ(int requestedDeg, long long coefDeg, UpopsZ_t* p, UpopsZ_t* alpha) {

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
			weierstrassUpdate_UPOPSZ(p, alpha);
		}
	}
	return p->data[coefDeg]->polys[requestedDeg];
}


/**
 *
 */
PolyZ_ptr weierstrassCoefVoidPUpdate_UPOPSZ(int i, void* coefDeg, void* p, void* alpha) {
	return weierstrassCoefPUpdate_UPOPSZ(i, (long long) coefDeg, (UpopsZ_t*) p, (UpopsZ_t*) alpha);
}

/**
 * F, G, H : power series
 * r : an integer
 * nvar : the number of variables
 *
 */
PolyZ_ptr lemmaForWeierstrass_UPOPSZ(PowerSeriesZ_t* F, PowerSeriesZ_t* G, PowerSeriesZ_t* H, int r) {
	if (F == NULL || G == NULL || H == NULL) {
		return NULL;
	}

	int nvar = F->polys[0] == NULL ? H->polys[0]->nvar : F->polys[0]->nvar;

	PolyZ_ptr s = NULL;
	for (int i = 1; i < r; ++i) {

		//Here we need homogPart, because sometimes alpha or p may be truncated.
		PolyZ_ptr multiply = multiplyPolynomials_AAZ(homogPart_PSZ(r-i, G), homogPart_PSZ(i, H), nvar);
		s = addPolynomials_AAZ_inp(s,  multiply, nvar);
		freePolynomial_AAZ(multiply);
	}

	if(r > F->deg) {
		if (!isZero_AAZ(s)) {
			negatePolynomial_AAZ(s);
		}
	} else {
		s = subPolynomials_AAZ_inp(s, F->polys[r], nvar);
		negatePolynomial_AAZ(s);
	}

	divideByIntegerExact_AAZ_inp(s, H->polys[0]->elems[0].coef);

	return s;
}



/**
 * upops : the input upops
 * p : the upops p
 * alpha : the upops alpha
 * update all power series of p and alpha
 */
void weierstrassUpdate_UPOPSZ(UpopsZ_t* p, UpopsZ_t* alpha) {


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
	PolyZ_ptr polyP;


	PowerSeriesZ_t** F = p->weierstrassFData;

	for (int lP = 0; lP <= d-1; ++lP) {
		// fprintf(stderr, "curdP %d, alloc: %d\n\n", curdP, p->data[lP]->alloc);
		if (curdP + 1 > p->data[lP]->alloc) {
			int newAlloc = (2*(p->data[lP])->alloc < curdP+1) ? curdP + 1 : 2*p->data[lP]->alloc;
			p->data[lP]->polys = (PolyZ_ptr*) realloc(p->data[lP]->polys, sizeof(PolyZ_ptr)*newAlloc);
			p->data[lP]->alloc = newAlloc;
		}
	}
	//alpha's allocation automatically managed by its power series structure.


	for (int lP = 0; lP <= d - 1; ++lP) {
		if (F[lP]->deg == -1){
			polyP = NULL;
		}else{
		   updateToDeg_PSZ(curdP,F[lP]); //update F based on new terms of p->data[lP-1]
		   polyP = lemmaForWeierstrass_UPOPSZ(F[lP], p->data[lP], alpha->data[0], curdP);
		}

		p->data[lP]->polys[curdP] = polyP;
		//manually set new degree since this updated PS is needed for next iteration of the loop.
		p->data[lP]->deg = curdP;
	}

	for (int lA = m; lA >= 0; --lA) {
		updateToDeg_PSZ(curdA,  alpha->data[lA]);
	}
}



/**
 * upops : the input upops
 * p_out : the output p
 * alpha_out : the output alpha
 * It takes a upops as input and computes p and alpha
 */
void weierstrassPreparation_UPOPSZ(UpopsZ_t* upops, UpopsZ_t** p_out, UpopsZ_t** alpha_out) {

	/* Check if the input is valid or not */
	int n = upops->deg;
	PowerSeriesZ_t** a = upops->data;
	int d = -1;
	for (int i = 0; i <= n; i++) {
		if (isUnit_PSZ(a[i])) {
			d = i;
			break;
		}
	}

	if (d == -1) {
		fprintf(stderr, "ERROR in weierstrassPreparation_UPOPSZ: Invalid input!\n");
		exit(1);
	}

	int nvar = upops->data[d]->polys[0] == NULL ? 0 : upops->data[d]->polys[0]->nvar;
	if (d == 0) {
		if (p_out != NULL) {
			*p_out = one_UPOPSZ(nvar);
		}
		if (alpha_out != NULL) {
			reserve_UPOPSZ(upops);
			*alpha_out = upops;
		}
		return;
	}

	long long m = n - d;

	UpopsZ_t* p = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(d+1);
	UpopsZ_t* alpha = allocateUnivariatePolynomialOverPowerSeries_UPOPSZ(m+1);

	PowerSeriesZ_t** b = p->data;
	PowerSeriesZ_t** c = alpha->data;
	PowerSeriesZ_t** F = (PowerSeriesZ_t**) malloc(sizeof(PowerSeriesZ_t*) * (d));

	/*  Initializing b  */

	for (int j = 0; j <= d - 1; j++) {
		b[j]= allocatePowerSeries_PSZ(0);
		b[j]->polys= copyUpTo_PSZ(upops->data[j], 0);
		b[j]->alloc = 1;
		b[j]->deg = 0;

	}

	b[d] = onePowerSeries_PSZ(nvar);

	/*  Initializing c */
	reserve_PSZ(a[d+m]);
	c[m] = a[d+m];
	for (int lA = m-1; lA >= 0; --lA) {

		c[lA] = upops->data[d + lA];
		reserve_PSZ(c[lA]);

		for (int jA = d-1; jA >= 0; --jA) {
			int k = d + lA - jA;

			if (k <= m) {

				PowerSeriesZ_t* mulpsA = multiplyPowerSeries_PSZ(b[jA], c[k]);
				PowerSeriesZ_t* tmp = subPowerSeries_PSZ(c[lA], mulpsA);
				destroyPowerSeries_PSZ(c[lA]); //c[lA] is reserved inside subPowerSeries_PSZ
				destroyPowerSeries_PSZ(mulpsA); //mulpsA is reserved inside subPowerSeries_PSZ
				c[lA] = tmp;
			}
		}
	}

	/*  Initializing F */
	reserve_PSZ(upops->data[0]);
	F[0] = upops->data[0];
	for (int lP = 1; lP <= d - 1; ++lP) {
		F[lP] = upops->data[lP];
		reserve_PSZ(F[lP]);
		if (m > 0) {
			for (int jP = 1; jP <= lP && jP <= m; ++jP) {
				PowerSeriesZ_t* mulpsP = multiplyPowerSeries_PSZ(b[lP - jP], c[jP]);
				PowerSeriesZ_t* tmp = subPowerSeries_PSZ(F[lP], mulpsP);
				destroyPowerSeries_PSZ(F[lP]);
				destroyPowerSeries_PSZ(mulpsP); //mulpsP is reserved inside sub
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
		b[k]->gen.tertiaryGen = &(weierstrassCoefVoidPUpdate_UPOPSZ);
		b[k]->genParam1 = (void*) k; //the index of the coefficient
		b[k]->genParam2 = (void*) p; //yes, a bit recrusive in nature. But it should be fine with the correct number of reserves and destroys.
		b[k]->genParam3 = (void*) alpha;
		// reserve_UPOPSZ(p); // use a weak reference to p
		reserve_UPOPSZ(alpha);
		b[k]->paramType1 = PLAIN_DATA;
		b[k]->paramType2 = WEAK_UPOPS;
		b[k]->paramType3 = UPOPS;
	}

	/*  alpha */
	alpha->deg = m;


	if (p_out != NULL) {
		*p_out = p;
	} else {
		destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(p);
	}

	if (alpha_out != NULL) {
		*alpha_out = alpha;
	} else {
		destroyUnivariatePolynomialOverPowerSeries_UPOPSZ(alpha);
	}

}
