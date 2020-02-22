#include "../../include/IntegerPolynomial/primes.h"
#include "../../include/IntegerPolynomial/modpoly.h"
//#include <cilktools/cilkview.h>

sfixn nextprime(int* s, mpz_class m) {
	if (*s >= BPASPRIMESNUM) {
		std::cout << "BPAS: error, not enough primes for modular GCD of DUZP.." << std::endl;
		exit(1);
	}
	sfixn p = BPAS_PRIMES.primes[*s];
	if (m > p) {
		while (m % p == 0) {
			if (*s >= BPASPRIMESNUM) {
				std::cout << "BPAS: error, not enough primes for modular GCD of DUZP.." << std::endl;
				exit(1);
			}
			p = BPAS_PRIMES.primes[*s];
			++*s;
		}
	}
	else { ++*s; }
	return p;
}

// forces the leading coefficient of the gcd to be lc
void monicGCD(sfixn* gp, int* gd, sfixn* fp, int fd, sfixn p, sfixn lc) {
	sfixn *g = gp, *f = fp;
	while (fd > 0) {
		int k = 0;
		sfixn invlc = inverseMod(fp[fd], p);

		// monic division
		while (*gd >= fd) {
			bool isFirst = 0;
			k = *gd - 1;
			sfixn e = MulMod(gp[*gd], invlc, p);
			for (int i = *gd-1, j = fd-1 ; j > -1; --i, --j) {
				gp[i] = SubMod(gp[i], MulMod(fp[j], e, p), p);
				// check deg(gp)
				if (!isFirst && gp[i]) {
					k = i;
					isFirst = 1;
				}
			}
			// continue to check deg(gp)
			for (int i = k; i > 0; --i) {
				if (gp[i] != 0) { break; }
				else { k = i - 1; }
			}
			*gd = k;
		}

		// swap gp and fp
		sfixn* t = gp;
		gp = fp;
		fp = t;
		*gd = fd;
		fd = k;
	}

	if (fp[0]) {
		*gd = 0;
		g[0] = 1;
	}
	else {
		sfixn e = MulMod(lc, inverseMod(gp[*gd], p), p);
		for (int i = 0; i <= *gd; ++i)
			g[i] = MulMod(gp[i], e, p);
	}
	fp = f;
}

// Pseudo remainder g := prem(g, f)
// deg(g) >= deg(f)
// prem (g, -f) = (-1)^(gd-fd)*prem(g, f)
void prem (sfixn* gp, int* gd, sfixn* fp, int fd, sfixn p) {
	bool delta = (*gd - fd + 1) & 1;
	while (*gd >= fd) {
		for (int i = *gd-1, j = fd-1; i > -1; --i, --j) {
			if (j > -1)
				gp[i] = SubMod(MulMod(gp[i], fp[fd], p), MulMod(fp[j], gp[*gd], p), p);
			else { gp[i] = MulMod(gp[i], fp[fd], p); }
		}
		--*gd;
	}

	if (*gd < 0) { *gd = 0; }
	if (delta) {
		for (int i = 0; i < fd; ++i) {
			if (gp[i]) {
				gp[i] = NegMod(gp[i], p);
				*gd = i; 
			}
		}
	}
	else {
		for (int i = fd-1; i > -1; --i)
			if (gp[i]) { *gd = i; break; }
	}
}

// Brown's subresultant algorithm
// deg(gp) >= deg(fp)
void brownSubresultant(UniModularSubresultantChain* S, sfixn* gp, int gd, sfixn* fp, int fd, sfixn p) {
	sfixn *g = gp, *f = fp;

	sfixn s = PowerMod(fp[fd], gd-fd, p);
	prem (gp, &gd, fp, fd, p);
	sfixn *t = gp;
	gp = fp;
	fp = t;
	int k = gd;
	gd = fd;
	fd = k;

	while (fd > 0 || fp[0]) {
		int delta = gd - fd;
		if (gd > 0) {
			S->deg[gd-1] = fd;
			std::copy(fp, fp+fd+1, S->coef[gd-1]);
		}
		if (delta > 1) {
			S->deg[fd] = fd;
			sfixn e = PowerMod(DivMod(fp[fd], s, p), delta-1, p);
			for (int i = 0; i <= fd; ++i)
				S->coef[fd][i] = MulMod(fp[i], e, p);
		}

		if (!fd) { break; }
		sfixn a = gp[gd];
		prem (gp, &gd, fp, fd, p);
		t = gp;
		gp = fp;
		fp = t;
		k = gd;
		gd = fd;
		fd = k;
		a = inverseMod(MulMod(a, PowerMod(s, delta, p), p), p);
		for (int i = 0; i <= fd; ++i)
			fp[i] = MulMod(fp[i], a, p);
		if (delta > 1) {
			std::copy(S->coef[gd], S->coef[gd]+S->deg[gd]+1, gp);
			gd = S->deg[gd];
		}
		s = gp[gd];
	}
	gp = g, fp = f;
}

// g := \prod (x-t[i]) mod p
// deg(g) = n
// t stores the evaluation points
void lagrangeBasis (sfixn* g, int n, sfixn p, sfixn* t) {
	g[0] = NegMod(t[0], p), g[1] = 1;
	for (int i = 2; i <= n; ++i) {
		g[i] = 0;
		for (int j = i; j > -1; --j) {
			g[j] = MulMod(g[j], NegMod(t[i-1], p), p);
			if (j) { g[j] = AddMod(g[j], g[j-1], p); }
		}
	}
}

// q = g / (x - u) mod p
// n is the size of q
// q/a is the result
sfixn longDivision(sfixn* q, sfixn u, sfixn* g, int n, sfixn p) {
	sfixn invu = inverseMod(u, p);
	sfixn a = q[0] = NegMod(MulMod(g[0], invu, p), p), x = u;
	for (int i = 1; i < n; ++i) {
		q[i] = MulMod(SubMod(q[i-1], g[i], p), invu, p);
		// also evaluated at point x = u
		a = AddMod(a, MulMod(q[i], x, p), p);
		x = MulMod(x, u, p);
	}
	return a;
}


// S is a sequence of subresultant chains: each of those was obtained
// atfer (1) specializing one variable in two bivariate polynomials, and
// (2) computing the subresultant chain of the specialized polynomials
// n is the length of the sequence, that is, the degree+1 of the targeted resultant
// t is the sequence of the evaluation points
// s is the output sequence of bivariate polynomials (still modulo p)
// Evaluation ppints stored in t
void interpolateSubresultant (BiModularSubresultantChain<sfixn>* s, UniModularSubresultantChain* S, int n, sfixn p, int alpha, sfixn* t) {
	sfixn* g = new sfixn[s->n*(n+1)], *q = new sfixn[s->n*n];
	// for each row of the subresultant chain
	#pragma cilk_grainsize = 128;
	cilk_for (int k = 0; k < s->n; ++k) {
		int delta = n - k * alpha, offset = n * k;
		lagrangeBasis (&g[offset+k], delta, p, t);

		// for each specialized regular chain at row k
		for (int i = 0; i < delta; ++i) {
			sfixn e = longDivision(&q[offset], t[i], &g[offset+k], delta, p);
			e = inverseMod(e, p);
                        // compute the corresponding Lagrange basis element
			for (int j = 0; j < delta; ++j)
				q[offset+j] = MulMod(q[offset+j], e, p);

			// \sum v[l] * q[j]

                        // for the l-th coefficient of Row k, add the contribution 
                        // v[l] * q[j] to it
			#pragma cilk_grainsize = 1024;
			cilk_for (int l = 0; l <= S[i].deg[k]; ++l) {
				for (int j = 0; j < delta; ++j)
					s->coef[k][l*delta+j] = AddMod(s->coef[k][l*delta+j], MulMod(S[i].coef[k][l], q[offset+j], p), p);
			}
		}
	}

	delete [] q;
	delete [] g;
}

// Subresultant chain modular a prime
// Bivariate polynomials f and g with their partial degrees
// Assuming deg(g, x) > deg(f, x)
BiModularSubresultantChain<sfixn> modularSubresultantChain(mpz_class* g, int* gd, mpz_class* f, int* fd, sfixn p) {
        // n2 is the total degree of g
        // n3 is the total degree of f
        int n, n2 = 0, n3 = 0;
	n = gd[0] + 1;
	for (int i = 0; i <= gd[1]; ++i) {
		for (int j = 0; j <= gd[0]; ++j) {
			if (g[n*i+j] != 0 && (i+j) > n2)
				n2 = i + j;
		}
	}
	n = fd[0] + 1;
	for (int i = 0; i <= fd[1]; ++i) {
		for (int j = 0; j <= fd[0]; ++j) {
			if (f[n*i+j] != 0 && (i+j) > n3)
				n3 = i + j;
		}
	}
	n2 = n2 * n3 + 1;
	// n is the size of the resultant based on the Sylvester matrix
	n = fd[1] * gd[0] + fd[0] * gd[1] + 1;
	if (n2 < n) { n = n2; }
	if (n < 2) { n = 2; }

        // n2 is the number of subresultants to be computed
	n2 = (fd[1] > 2)? fd[1] : 2;
        // degree in x of f(t,x) and g(t,x) mod <t -t_0>
	int as = gd[1] + 1, bs = fd[1] + 1;
	sfixn *a = new sfixn[as*n];
	sfixn *b = new sfixn[bs*n];
	// Evaluation points
	sfixn *t = new sfixn[n];

	// allocate copies of subresultant chain over a machine prime
	UniModularSubresultantChain *S = new UniModularSubresultantChain[n];
	for (int i = 0; i < n; ++i) {
		S[i].n = n2;
		S[i].deg = new int[S[i].n];
		S[i].size = new int[S[i].n];
		for (int j = 0; j < S[i].n; ++j) {
			S[i].deg[j] = 0;
			S[i].size[j] = j+2;
		}
		S[i].coef = new sfixn*[S[i].n];
		for (int j = 0; j < S[i].n; ++j) {
			S[i].coef[j] = new sfixn[j+2];
			for (int k = 0; k <= j; ++k)
				S[i].coef[j][k] = 0;
		}
	}

	//unsigned long long start, end;
	//float elapsed;

	//start = __cilkview_getticks();
	// evaluate variable t and compute subresultant chain of variable x
	mpz_class m (p);
	#pragma cilk_grainsize = 1024;
	cilk_for (int k = 0; k < n; ++k) {
		// evaluation point l * n + k + 1
		t[k] = k + 1;
		int ad = 0, bd = 0, at = gd[0] + 1, bt = fd[0] + 1, l = 0;
		while (ad < gd[1] || bd < fd[1]) {
			t[k] = AddMod(t[k], MulMod(l, n, p), p);
			for (int i = 0; i <= gd[1]; ++i) {
				mpz_class e = g[at*i+gd[0]];
				for (int j = gd[0]-1; j > -1; --j)
					e = e * t[k] + g[at*i+j];
				mpz_mod (e.get_mpz_t(), e.get_mpz_t(), m.get_mpz_t());
				a[k*as+i] = e.get_mpz_t()->_mp_d[0];
				if (a[k*as+i]) { ad = i; }
			}
			for (int i = 0; i <= fd[1]; ++i) {
				mpz_class e = f[bt*i+fd[0]];
				for (int j = fd[0]-1; j > -1; --j)
					e = e * t[k] + f[bt*i+j];
				mpz_mod (e.get_mpz_t(), e.get_mpz_t(), m.get_mpz_t());
				b[k*bs+i] = e.get_mpz_t()->_mp_d[0];
				if (b[k*bs+i]) { bd = i; }
			}
			l++;
		}

		// Brown's subresultant algorithm
		brownSubresultant(&S[k], &a[k*as], ad, &b[k*bs], bd, p);
	}
	//end = __cilkview_getticks();
	//elapsed = (end - start) / 1000.f;
	//std::cout << "Brown's subresultant with " << n << " images:\t" << elapsed << std::endl;

	// Lagrange interpolate
	int k = gd[0] + fd[0];
	BiModularSubresultantChain<sfixn> s;
	s.n = n2;
	s.deg = new int*[n2];
	s.size = new int[n2];
	s.coef = new sfixn*[n2];
	for (int i = 0; i < n2; ++i) {
		s.deg[i] = new int[2];
		s.deg[i][1] = i, s.deg[i][0] = n - i * k - 1;
		s.size[i] = (s.deg[i][1] + 1) * (s.deg[i][0] + 1);
		s.coef[i] = new sfixn[s.size[i]];
		for (int j = 0; j < s.size[i]; ++j)
			s.coef[i][j] = 0;
	}
	//start = __cilkview_getticks();
	interpolateSubresultant(&s, S, n, p, k, t);
	//end = __cilkview_getticks();
	//elapsed = (end - start) / 1000.f;
	//std::cout << "Lagrange interpolation:\t\t\t" << elapsed << std::endl;

	delete [] a;
	delete [] b;
	delete [] t;
	delete [] S;

	return s;
}

// Subresultant chain modular a set of primes
// Bivariate polynomials f and g with their partial degrees
// Assuming deg(g, x) > deg(f, x)
BiModularSubresultantChain<mpz_class> modularSetSubresultantChain(mpz_class* g, int* gd, mpz_class* f, int* fd, const sfixn* p, int s, mpz_class* m) {
	// The first modular subresultatn chain
	BiModularSubresultantChain<sfixn> S = modularSubresultantChain (g, gd, f, fd, p[0]);

	if (!S.n) { return BiModularSubresultantChain<mpz_class>(); }
	BiModularSubresultantChain<mpz_class> a (S.n);
	*m = p[0];
	mpz_class halfp = (p[0] - 1) >> 1;
	for (int j = 0; j < S.n; ++j) {
		a.deg[j][1] = S.deg[j][1], a.deg[j][0] = S.deg[j][0];
		a.size[j] = S.size[j];
		a.coef[j] = new mpz_class[a.size[j]];
		for (int l = 0; l < a.size[j]; ++l) {
			a.coef[j][l] = S.coef[j][l];
			if (a.coef[j][l] > halfp) { a.coef[j][l]-= p[0]; }
		}
	}

	for (int i = 1; i < s; ++i) {
		// next modular subresultant chain
		S = modularSubresultantChain (g, gd, f, fd, p[i]);
		mpz_class t, x(p[i]), e, o;
		mpz_gcdext (t.get_mpz_t(), e.get_mpz_t(), o.get_mpz_t(), x.get_mpz_t(), m->get_mpz_t());
		t = *m * x;
		halfp = (t - 1) >> 1;
		// Merge two modular subresultant chain
		for (int j = 0; j < a.n; ++j) {
			for (int l = 0; l < a.size[j]; ++l) {
				//a.coef[j][l] = S.coef[j][l] * e * p[i] + a.coef[j][l] * o * *m;
				mpz_class tmp = (S.coef[j][l] - a.coef[j][l]) * o;
				mpz_mod (tmp.get_mpz_t(), tmp.get_mpz_t(), x.get_mpz_t());
				a.coef[j][l] += tmp * *m;
				if (a.coef[j][l] > halfp) { a.coef[j][l] -= t; }
			}
		}
		*m = t;
	}
	return a;
}

BiModularSubresultantChain<mpz_class> modularStableSubresultantChain (mpz_class* g, int* gd, mpz_class* f, int* fd, int v) {
	if (v != 2) { return BiModularSubresultantChain<mpz_class>(); }

	int offset = 0, k = 1;
	mpz_class m;
	// The fist set of modular subresultant chains
	BiModularSubresultantChain<mpz_class> s = modularSetSubresultantChain (g, gd, f, fd, &BPAS_PRIMES.primes[offset], k, &m);
	BiModularSubresultantChain<mpz_class> z (s.n);
	for (int i = 0; i < z.n; ++i) {
		z.deg[i][1] = s.deg[i][1], z.deg[i][0] = s.deg[i][0];
		z.size[i] = s.size[i];
		z.coef[i] = new mpz_class[z.size[i]];
		for (int j = 0; j < z.size[i]; ++j)
			z.coef[i][j] = 0;
	}

	// Until the chain is stable
	while (1) {
		offset += k;
		mpz_class x, e, o, t, halft;
		// Next set of modular subresultant chains
		BiModularSubresultantChain<mpz_class> h = modularSetSubresultantChain (g, gd, f, fd, &BPAS_PRIMES.primes[offset], k, &x);
		mpz_gcdext (t.get_mpz_t(), e.get_mpz_t(), o.get_mpz_t(), x.get_mpz_t(), m.get_mpz_t());

		t = x * m;
		halft = (t - 1) >> 1;
		// Merge two set of modular subresultant chains
		for (int i = 0; i < z.n; ++i) {
			for (int j = 0; j < z.size[i]; ++j) {
				//z.coef[i][j] = s.coef[i][j] * e * m + h.coef[i][j] * o * x;
				z.coef[i][j] = (h.coef[i][j] - s.coef[i][j]) * o;
				mpz_mod (z.coef[i][j].get_mpz_t(), z.coef[i][j].get_mpz_t(), x.get_mpz_t());
				z.coef[i][j] = s.coef[i][j] + z.coef[i][j] * m;
				if (z.coef[i][j] > halft) { z.coef[i][j] -= t; }
			}
		}
		if (s == z) { break; }
		s = z;
		m = t;
	}

	return s;
}
