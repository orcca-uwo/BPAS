

#include "LinearAlgebra/Vandermonde.h"
#include "IntegerPolynomial/DUZP_Support.h"


int _solveFFVandermondeSystem_1_T(const elem_t* t, const elem_t* b, int size, elem_t** x, const Prime_ptr* Pptr) {
	if (x == NULL) {
		return 0;
	}

	//t values must be pairwise distinct.
	for (int i = 0; i < size; ++i) {
		for (int j = i+1; j < size; ++j) {
			if (t[i] == t[j] || t[i] == 0) {
				return 0;
			}
		}
	}

	duspoly_t* P = makePolynomial_spX(size + 1);
	P->elems[0] = smallprimefield_convert_in(1, Pptr);

	for (int i = 0; i < size; ++i) {
		multiplyByBinomialInForm_spX_inp(&P, t[i], Pptr);
	}

	if (*x == NULL) {
		*x = (elem_t*) malloc(sizeof(elem_t)*size);
	}

	elem_t* xVals = *x;

	duspoly_t* q = NULL;
	elem_t eval = 0;

	for (int i = 0; i < size; ++i) {
		divideByMonicLinearInForm_spX(P, t[i], &q, NULL, Pptr);
		divideByMonicLinearInForm_spX(q, t[i], NULL, &eval, Pptr);
		eval = smallprimefield_mul(eval, t[i], Pptr);
	
		scalarMulPolynomialInForm_spX_inp(&q, smallprimefield_inv(eval, Pptr), Pptr);
		xVals[i] = smallprimefield_convert_in(0, Pptr);
		for (int j = 0; j < size; ++j) {
			xVals[i] = smallprimefield_add(xVals[i], smallprimefield_mul(b[j], q->elems[j], Pptr), Pptr);
		}
	}

	return 1;
}

int _solveFFVandermondeSystem_0_T(const elem_t* t, const elem_t* b, int size, elem_t** x, const Prime_ptr* Pptr) {
	if (x == NULL) {
		return 0;
	}

	//t values must be pairwise distinct.
	for (int i = 0; i < size; ++i) {
		for (int j = i+1; j < size; ++j) {
			if (t[i] == t[j] || t[i] == 0) {
				return 0;
			}
		}
	}

	duspoly_t* P = makePolynomial_spX(size + 1);
	P->elems[0] = smallprimefield_convert_in(1, Pptr);

	for (int i = 0; i < size; ++i) {
		multiplyByBinomialInForm_spX_inp(&P, t[i], Pptr);

	}

	if (*x == NULL) {
		*x = (elem_t*) malloc(sizeof(elem_t)*size);
	}

	elem_t* xVals = *x;

	duspoly_t* q = NULL;
	elem_t eval = 0;

	for (int i = 0; i < size; ++i) {
		divideByMonicLinearInForm_spX(P, t[i], &q, NULL, Pptr);
		divideByMonicLinearInForm_spX(q, t[i], NULL, &eval, Pptr);
	
		scalarMulPolynomialInForm_spX_inp(&q, smallprimefield_inv(eval, Pptr), Pptr);
		xVals[i] = smallprimefield_convert_in(0, Pptr);
		for (int j = 0; j < size; ++j) {
			xVals[i] = smallprimefield_add(xVals[i], smallprimefield_mul(b[j], q->elems[j], Pptr), Pptr);
		}
	}

	return 1;
}


int _solveIntegerVandermondeSystem_1_T(const mpz_t* t, const mpz_t* b, int size, mpz_t** x) {
	if (x == NULL) {
		return 0;
	}

	//begin by creating the master polynomial P.
	DUZP_t* P = makePolynomial_DUZP(size + 1);
	mpz_set_ui(P->coefs[0],1ul);

	for (int i = 0; i < size; ++i) {
		multiplyByBinomial_DUZP_inp(&P, t[i]);
	}


	if (*x == NULL) {
		*x = (mpz_t*) malloc(sizeof(mpz_t)*size);
		for (int j = 0; j < size; ++j) {
			mpz_init((*x)[j]);
		}
	}
	mpz_t* xVals = *x;

	DUZP_t* q = NULL;
	mpz_t eval;
	mpz_init(eval);

	for (int i = 0; i < size; ++i) {
		//TODO inplace!!
		fprintf(stderr, " i = %d\n", i);
		divideByMonicLinear_DUZP(P, t[i], &q, NULL);
		printPoly_DUZP(q, "x");
		evaluate_DUZP(q, t[i], eval);
		mpz_mul(eval, eval, t[i]);
		divideByIntegerExact_DUZP_inp(q, eval);
		mpz_set_ui(xVals[i], 0ul);
		for (int j = 0; j < size; ++j) {
			mpz_addmul(xVals[i], b[j], q->coefs[j]);
		}
	}


	return 1;
}


int solveFFVandermondeSystemInForm(const elem_t* t, const elem_t* b, int size, int startExp, int transposed, elem_t** x, const Prime_ptr* Pptr) {
	if (x == NULL) {
		return 0;
	}

	if (transposed) {
		if (startExp == 0) {
			return _solveFFVandermondeSystem_0_T(t, b, size, x, Pptr);
		} else {
			return _solveFFVandermondeSystem_1_T(t, b, size, x, Pptr);
		}
	}

	//TODO more configs
	return 0;
}

int solveFFVandermondeSystem(const elem_t* t, const elem_t* b, int size, int startExp, int transposed, elem_t** x, const Prime_ptr* Pptr) {
	if (x == NULL) {
		return 0;
	}

	elem_t t_in[size];
	memcpy(t_in, t, sizeof(elem_t)*size);
	for (int i = 0; i < size; ++i) {
		t_in[i] = smallprimefield_convert_in(t_in[i], Pptr);
	}
	elem_t b_in[size];
	memcpy(b_in, b, sizeof(elem_t)*size);
	for (int i = 0; i < size; ++i) {
		t_in[i] = smallprimefield_convert_in(t_in[i], Pptr);
	}

	int ret;
	if (transposed) {
		if (startExp == 0) {
			ret = _solveFFVandermondeSystem_0_T(t, b, size, x, Pptr);
		} else {
			ret = _solveFFVandermondeSystem_1_T(t, b, size, x, Pptr);
		}
	}

	for (int i = 0; i < size; ++i) {
		(*x)[i] = smallprimefield_convert_out((*x)[i], Pptr);
	}

	//TODO more configs
	return 0;
}