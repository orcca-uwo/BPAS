#include "DyadicRationalNumber/Multiplication/multiplication.h"

unsigned long long getMaxBits(DyadicRationalNumber* a, int n) {
        unsigned long long maxbits = 0;
        for (int j = 0; j < n; j++) {
                if (a[j].getDen() > maxbits)
                        maxbits = a[j].getDen();
        }
	return maxbits;
}

unsigned long long getMaxBits(mpz_class* a, int n) {
	lfixz maxelem = 0;
	for (int j = 0; j < n; j++) {
		if (a[j] > maxelem)
			maxelem = a[j];
	}

	unsigned long long maxbits = 0;
	while (maxelem > 0) {
		maxbits++;
		maxelem >>= 1;
	}
	return maxbits;
}

unsigned long long getMaxBits(mpq_class* a, int n) {
	lfixz maxelem = 0;
	for (int j = 0; j < n; j++) {
		if (a[j].get_den() > maxelem)
			maxelem = a[j].get_den();
	}

	if (maxelem == 1) { return 0; }

	unsigned long long maxbits = 0;
	while (maxelem > 0) {
		maxbits++;
		maxelem >>= 1;
	}
	return maxbits;
}

void univariateMultiplication(lfixz* mul, lfixz* a, int n, lfixz* b, int m) {
	//try {
		int size = n + m - 1;
		//unsigned long long ak = getMaxBits(a, n);
		//unsigned long long bk = getMaxBits(b, m);
		UnivariateIntegerPolynomial aPoly(n, a);
		UnivariateIntegerPolynomial bPoly(m, b);
		UnivariateIntegerPolynomial rPoly(size, mul);

		//Multiplication f;
		//f.multiply(&rPoly, &aPoly, &bPoly);

		if (size < 64) {
			MulNaive naive;
			naive.multiply(&aPoly, &bPoly, &rPoly);
		}
		else if (size < 2048) {
//			std::cout << "Sch\"onaghe-Strassen via Kronecker's substitution";
			MulKS ks;
			ks.multiply(&aPoly, &bPoly, &rPoly);
		}
		else if (size < 4096) {
			int worker = __cilkrts_get_nworkers();
			if (worker <= 6) {
//				std::cout << "4-way Toom-Cook";
				MulToom4 toom4;
				toom4.multiply(&aPoly, &bPoly, &rPoly, 4);
			}
			else {
//				std::cout << "8-way Toom-Cook";
				MulToom8 toom8;
				toom8.multiply(&aPoly, &bPoly, &rPoly, 8);
			}
		}
		else {
//			std::cout << "two convolution";
			MulSSA ssa;
			ssa.multiply(&aPoly, &bPoly, &rPoly);
		}

	//}
	//catch (std::exception& e) {
	//	std::cout << "Multiplicaion Error: using naive multiplication instead!" << std::endl;
	//	naiveUnivariateMultiplication(mul, a, n, b, m);
	//}
}

void univariateMultiplication(lfixq* mul, DyadicRationalNumber* a, int n, lfixz* b, int m) {
        unsigned long long k = getMaxBits(a, n);

        lfixz* c = new lfixz[n];
        for (int i = 0; i < n; i++) {
		DyadicRationalNumber elem = a[i] << k;
		c[i] = elem.getNum();
	}

        lfixz* d = new lfixz[n+m];
        univariateMultiplication(d, c, n, b, m);

        for (int i = 0; i < n+m-1; i++) {
                mul[i] = d[i];
		mul[i] >>= k;	
	}
        mul[n+m-1] = 0;

        delete [] c;
        delete [] d;
}

void univariateMultiplication(lfixq* mul, mpq_class* a, int n, lfixz* b, int m) {
        unsigned long long k = getMaxBits(a, n);

        lfixz* c = new lfixz[n];
        for (int i = 0; i < n; i++)
                c[i] = a[i] << k;

        lfixz* d = new lfixz[n+m];
        univariateMultiplication(d, c, n, b, m);

        for (int i = 0; i < n+m-1; i++) {
                mul[i] = d[i];
                mul[i] >>= k;
        }
        mul[n+m-1] = 0;

        delete [] c;
        delete [] d;
}


void univariateMultiplication(lfixq* mul, DyadicRationalNumber* a, int n, DyadicRationalNumber* b, int m) {
	unsigned long long ak = getMaxBits(a, n);
	unsigned long long bk = getMaxBits(b, m);

	lfixz* c = new lfixz[n];
	for (int i = 0; i < n; i++) {
		DyadicRationalNumber elem = a[i] << ak;
		c[i] = elem.getNum();
	}
	lfixz* d = new lfixz[m];
	for (int i = 0; i < m; i++) {
		DyadicRationalNumber elem = b[i] << bk;
		d[i] = elem.getNum();
	}

	lfixz* f = new lfixz[n+m];
	univariateMultiplication(f, c, n, d, m);

	for (int i = 0; i < n+m-1; i++) {
		mul[i] = f[i];
		mul[i] >>= ak+bk;
	}
	mul[n+m-1] = 0;

	delete [] c;
	delete [] d;
	delete [] f;
}

void univariateMultiplication(lfixq* mul, mpq_class* a, int n, mpq_class* b, int m) {
        unsigned long long ak = getMaxBits(a, n);
        unsigned long long bk = getMaxBits(b, m);

        lfixz* c = new lfixz[n];
        for (int i = 0; i < n; i++)
                c[i] = a[i] << ak;
        lfixz* d = new lfixz[m];
        for (int i = 0; i < m; i++)
                d[i] = b[i] << bk;

        lfixz* f = new lfixz[n+m];
        univariateMultiplication(f, c, n, d, m);

        for (int i = 0; i < n+m-1; i++) {
                mul[i] = f[i];
                mul[i] >>= ak+bk;
        }
        mul[n+m-1] = 0;

        delete [] c;
        delete [] d;
        delete [] f;
}

void naiveUnivariateMultiplication(lfixz* mul, lfixz* a, int n, lfixz* b, int m) {
        for (int i = 0; i < n+m-1; i++)
                mul[i] = 0;

        for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++)
                        mul[i+j] += a[i] * b[j];
        }
}
