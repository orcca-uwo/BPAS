#include "RationalNumberPolynomial/urpolynomial.h"

void DenseUnivariateRationalPolynomial::taylorShift(int ts) {
	int B = 16, m = curd + 1;
	if (m % B) { B = 20; }
	if (!ts)
		taylorShiftIncrementalCilkFor(coef, m, B);
	else if (ts == 1)
		taylorShiftTableau(B);
	else if (ts == 2)
		taylorShiftDnC(URPTSDNCBASE);
	else if (ts == 3)
		taylorShiftCVL();
	else
		taylorShiftDnC(URPTSDNCBASE, B);
}

void DenseUnivariateRationalPolynomial::taylorShiftDnC(int base, int B) {
	int m = curd + 1;
	if (m < base) {
		taylorShiftIncrementalCilkFor(coef, m, B);
	}
	else {
		int pos = m, k = 0;
		while (pos) {
			pos >>= 1;
			k++;
		}

		k = 1 << (k - 1);
		if (m > k)
			k <<= 1;
		lfixz* bi = new lfixz[k];
		binomials(bi, k);
		if (m != k) {
			DenseUnivariateRationalPolynomial A(k);
			A.curd = k - 1;
			for (int i = 0; i < k; ++i) {
				if (i < m)
					A.coef[i] = coef[i];
				else
					A.coef[i] = 0;
			}

			A.taylorShiftBasePower2(bi, URPTSDNCBASE);
			for (int i = 0; i < m; ++i)
				coef[i] = A.coef[i];
		}
		else
			taylorShiftBasePower2(bi, URPTSDNCBASE);
		delete [] bi;
	}
}

void ts_modulo (lfixz* r, lfixz a, lfixz p, lfixz p2, int k) {
	//*r = a % p;
	lfixz highest = a >> k;
	lfixz least = (highest << k) ^ a;
	if (a >= p2) {
		lfixz tmp = highest >> k;
		least += tmp;
		highest = (tmp << k) ^ highest;
	}
	*r = least - highest;
	if (*r < 0)
		*r += p;
	else if (*r > p)
		*r -= p;
}


void DenseUnivariateRationalPolynomial::taylorShiftCVL() {
	int m = curd, k = 0;
	while (m) {
		m >>= 1;
		k++;
	}
	k = 1 << (k + 1);
	lfixz p = 1;
	p <<= k;
	p++;
	lfixz p2 = p * p;

	m = curd + 1;
	lfixz* h = new lfixz[m];
	for (int i = 0; i < m; ++i)
		h[i] = coef[i];

	lfixz* f = new lfixz[m];
	f[1] = f[0] = 1;
	for (int i = 2; i <= curd; ++i) {
		ts_modulo(&f[i], f[i-1] * i, p, p2, k);
		ts_modulo(&h[i], h[i] * f[i], p, p2, k);
	}

	lfixz* g = new lfixz[m];
	g[0] = 1;
	for (int i = 1; i < m; ++i)
		ts_modulo(&g[i], g[i-1] * (curd - i + 1), p, p2, k);

	lfixz* r = new lfixz[m<<1];
	univariateMultiplication(r, h, m, g, m);

	lfixz inv;
	mpz_invert(inv.get_mpz_t(), f[curd].get_mpz_t(), p.get_mpz_t());
	ts_modulo(&inv, inv * inv, p, p2, k);
	for (int i = 0; i <= curd; ++i) {
		lfixz tmp;
		ts_modulo(&tmp, r[i+curd], p, p2, k);
		ts_modulo(&tmp, r[i+curd] * g[curd-i], p, p2, k);
		ts_modulo(&tmp, tmp * inv, p, p2, k);
		coef[i] = tmp;
	}

	delete [] f;
	delete [] h;
	delete [] g;
	delete [] r;
}

/***************************** Divide & Conquer *******************************/

void DenseUnivariateRationalPolynomial::binomials(lfixz* bi, int m) {
        for (int i = 1; i < m; i <<= 1) {
                bi[i-1] = 1;
                for (int j = 0; j < i/2; ++j) {
                        bi[i+j] = bi[i+j-1] * (i-j) / (j+1);
                        if (j+1 < i/2)
                                bi[2*i-j-2] = bi[i+j];
                }
        }
        bi[m-1] = 0;
}

void DenseUnivariateRationalPolynomial::taylorShiftBasePower2(lfixz* bi, int base) {
	int m = curd + 1;
	int d = (base > m)? m : base;

	// Base case
	cilk_for (int i = 0; i < m/d; ++i)
		taylorShiftIncrementalCilkFor(&coef[i*base], d, 16);

	lfixq* mul = new lfixq[m];
        for (int i = base; i < m; i <<= 1) {
		cilk_for (int k = 0; k < m/(2*i); ++k)
			univariateMultiplication(&mul[2*k*i], &coef[2*k*i+i], i, &bi[i-1], i);
                for (int k = 0; k < m; ++k)
                        coef[k] += mul[k];
        }
	delete [] mul;
}

void DenseUnivariateRationalPolynomial::taylorShiftDnC(int base) {
	std::vector<int> degbin;

	int m = curd + 1;
	int pos = m, k = 0;
	while (pos) {
		if (pos % 2) { degbin.push_back(k); }
		pos >>= 1;
		k++;
	}

	k = 1 << (k - 1);
	lfixz* bi = new lfixz[k];
	binomials(bi, k);

	pos = degbin.size();
	DenseUnivariateRationalPolynomial* A = new DenseUnivariateRationalPolynomial[pos];
#if __GNUC__ == 4
	cilk_for (int i = 0; i < pos; i++) {
#else
	for (int i = 0; i < pos; i++) {
#endif
		int d = (int) pow(2, degbin[i]);
		int offset = 0;
		for (int j = 0; j < i; ++j)
			offset += (int) pow(2, degbin[j]);

		A[i].coef = new lfixq[d];
		A[i].n = d;
		A[i].curd = d - 1;
		for (int j = 0; j < d; ++j) {
			A[i].coef[j] = coef[offset+j];
			coef[offset+j] = 0;
		}
		A[i].taylorShiftBasePower2(bi, base);
	}

	k = 0;
	lfixz* b = new lfixz[m+1];
	for (int i = 0; i < pos; i++) {
		int d = A[i].curd + 1;
		if (i) {
			lfixq* mul = new lfixq[d+k+1];
			// Multiply (1+x)^(k+d), k can be any integer
			univariateMultiplication(mul, A[i].coef, d, b, k+1);
			for (int j = d+k-1; j > -1; --j)
				coef[j] += mul[j];
			delete [] mul;
			if (i < pos - 1) {
				lfixz* prevb = new lfixz[k+1];
				for (int j = 0; j <= k; ++j)
					prevb[j] = b[j];
				univariateMultiplication(b, prevb, k+1, &bi[d-1], d);
				for (int j = k+d; j >= d; --j)
					b[j] += prevb[j-d];
				k += d;
				delete [] prevb;
			}
		}
		else {
			for (int j = 0; j < d; j++)
				coef[j] += A[i].coef[j];

			if (i < pos - 1) {
				for (int j = 0; j < d; ++j)
					b[j] = bi[d-1+j];
				b[d] = 1;
				k += d;
			}
		}
	}

	delete [] b;
	delete [] bi;
	delete [] A;
}


/***************************** IncrementalCilkFor *****************************/

void DenseUnivariateRationalPolynomial::tableauBase(lfixq* a, int B) {
	int m = B << 1;

	// First row
	for (int i = B; i < m; i++)
		a[i] += a[i-1];

	// Middle rows
	for (int i = B-2; i > 0; i--)
		for (int j = 0; j < B; j++)
			a[i+j+1] = a[i+j] + a[i+j+2];

	// First element of last row
	a[0] += a[2];

	//The rest elements of last row
	for (int i = 1; i < B; i++)
		a[i] = a[i-1] + a[i+2];

	// First element of first row
	a[B] = a[B-1];
}

void DenseUnivariateRationalPolynomial::polygonBase(lfixq* a, int B, int r, int k) {
	int m = B << 1;

	// First row
	for (int i = B; i < m; i++)
		a[i] += a[i-1];

	// The rectangle rows, same as tableau
	for (int i = B-2; i >= k; i--)
		for (int j = 0; j < B; j++)
			a[i+j+1] = a[i+j] + a[i+j+2];

	// The trapezoid starts
	for (int i = k-1; i > 0; i--) {
		a[k+1] += a[i];
		for (int j = k+2; j <= B+i; j++)
			a[j] += a[j-1];
	}

	// The first element of last row
	a[0] += a[k+1];

	// The rest elements of last row
	for (int i = 1; i <= r; i++)
		a[i] = a[i-1] + a[i+k+1];
}

void DenseUnivariateRationalPolynomial::taylorShiftBase(lfixq* a, int B) {
	for (int i = 0; i < B; i++) {
		a[B-1] += a[B+i];
		for (int j = B-2; j >= i; j--)
			a[j] += a[j+1];
	}
}

void DenseUnivariateRationalPolynomial::taylorShift3RegularCilkFor(lfixq* p, int d, int B) {
	int m = d << 1;

	lfixq* b = new lfixq[m];
	for (int i = 0; i < d; i++)
		b[i] = p[i];
	for (int i = d; i < m; i++)
		b[i] = 0;

	int k = d / B;
	lfixq* bn = b + d;

	if (k > 1) {
		tableauBase(bn-B, B);
		for (int i = 2; i < k; i++) {
			cilk_for (int j = 0; j < i; j++)
				tableauBase(bn+(2*j-i)*B, B);
		}
	}
	cilk_for (int i = 0; i < k; i++)
		taylorShiftBase(b+2*i*B, B);

	// Copy elements
	for (int i = 0; i < k; i++) 
		for (int j = 0; j < B; j++)
			p[j+i*B] = b[j+2*i*B];

	delete [] b;
}

void DenseUnivariateRationalPolynomial::taylorShiftIncrementalCilkFor(lfixq* p, int d, int B) {
	int r = d % B;
	if (!r) {
		taylorShift3RegularCilkFor(p, d, B);
		return;
	}

	int q = d / B;
	int m = d << 1;
	int k = B - r - 1;
	k = (k > 0)? k : 1;

	lfixq* b = new lfixq[m];
	cilk_for (int i = 0; i < d; i++)
		b[i] = p[i];
	cilk_for (int i = d; i < m; i++)
		b[i] = 0;

	// First tableau
	lfixq* bn = b + d;
	if (q > 1) {
		tableauBase(bn-B, B);
		for (int i = 2; i < q; i++) {
			cilk_for (int j = 0; j < i; j++)
				tableauBase(bn+(2*j-i)*B, B);
		}
	}

	// Second polygon
	cilk_for (int i = 0; i < q; i++)
		polygonBase(bn+(2*i-q)*B, B, r, k);

	// Final triangle
	cilk_for (int i = 0; i <= q; i++)
		taylorShiftBase(b+2*i*B, r);


	// Copy the first triangle
	cilk_for (int j = 0; j < r; j++)
		p[j] = b[j];

	cilk_for (int i = 1; i <= q; i++) {
		// Copy the triangle
		for (int j = 0; j < r; j++)
			p[j+i*B] = b[j+2*i*B];
		// Copy the polygon
		for (int j = 1; j < B-r; j++)
			p[i*B-j] = b[2*i*B-j];
		p[(i-1)*B+r] = b[2*((i-1)*B+r)];
	}

	delete [] b;
}

/******************************** Tableau *************************************/

void DenseUnivariateRationalPolynomial::taylorShiftBase(lfixq* p, lfixq* q, int B) {
	for (int i = 0; i < B; i++) {
		p[0] += q[i];
		for (int j = 1; j < B-i; j++)
			p[j] += p[j-1];
	}
}

void DenseUnivariateRationalPolynomial::tableauBaseInplace(lfixq* p, lfixq* q, int l, int r) {
	for (int i = 0; i < r; i++) {
		p[0] += q[i];
		for (int j = 1; j < l; j++)
			p[j] += p[j-1];
		q[i] = p[l-1];
	}
}

void DenseUnivariateRationalPolynomial::tableauConstruction(lfixq* p, lfixq* q, int B, int l, int r) {
	int k = (l > r)? l : r;

	if (k <= B) 
		tableauBaseInplace(p, q, l, r);
	else {
		int i = (l+1) >> 1;
		int j = (r+1) >> 1;

		tableauConstruction(p, q, B, i, j);
		cilk_spawn tableauConstruction(p, q+j, B, i, r-j);
		tableauConstruction(p+i, q, B, l-i, j);
		cilk_sync;
		tableauConstruction(p+i, q+j, B, l-i, r-j);    
	}
}

void DenseUnivariateRationalPolynomial::taylorShiftGeneral(lfixq* p, lfixq* q, int d, int B) {
	if (d <= B)
		taylorShiftBase(p, q, d);
	else {
		int m = (d+1) >> 1;
		tableauConstruction(p, q, B, m, m);
		cilk_spawn taylorShiftGeneral(p, q+m, d-m, B);
		taylorShiftGeneral(p+m, q, d-m, B);
		cilk_sync;
	}
}

void DenseUnivariateRationalPolynomial::taylorShiftTableau(int B) {
	int m = curd + 1;
	lfixq* q = new lfixq[m];
	for (int i = 0; i < m; i++)
		q[i] = 0;

	reciprocal();
	taylorShiftGeneral(coef, q, m, B);
	reciprocal();

	delete [] q;
}
