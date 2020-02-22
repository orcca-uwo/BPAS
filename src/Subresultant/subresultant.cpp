#include "RationalNumberPolynomial/urpolynomial.h"
#include "RingPolynomial/upolynomial.h"
#include "IntegerPolynomial/modpoly.h"

mpz_class SUPDUQP_to_mpq (SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial>& a, std::vector<mpq_class>& q, int* n) {
	mpz_class x = 1;
	n[1] = a.degree().get_si();
	n[0] = 0;
	for (int i = 0; i <= n[1]; ++i) {
		DenseUnivariateRationalPolynomial e = a.coefficient(i);
		int d = e.degree().get_si();
		if (d > n[0]) { n[0] = d; }
		q.push_back(mpq_class(d));
		for (int j = 0; j <= d; ++j) {
			mpq_class t = e.coefficient(j).get_mpq();
			q.push_back(t);
			x *= t.get_den();
		}
	}
	return x;
}

void SUPDUQP_to_mpz (mpz_class* g, int* n, std::vector<mpq_class> q, mpq_class x) {
	int k = 0, s = n[0] + 1;
	for (int i = 0; i <= n[1]; ++i) {
		mpz_class t = mpz_class(q[k]);
		for (int j = 0; j <= n[0]; ++j) {
			if (t >= j)
				g[s*i+j] = mpz_class(x * q[k+j+1]);
			else { g[s*i+j] = 0; }
		}
		k += t.get_si() + 2;
	}
}

std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > bivariateRationalSubresultantChain (SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial>& a, SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial>& b) {
	Symbol av[2], bv[2];
	av[1] = a.variable(), bv[1] = b.variable();
	av[0] = a.coefficient(a.degree().get_si()).variable(), bv[0] = b.coefficient(b.degree().get_si()).variable();
	if (av[0] != bv[0] || av[1] != bv[1]) {
		std::cout << "BPAS: error, trying to compute subresultant chain between SUP<DUQP>(" << av[1] << ", " << av[0] << ") and SUP<DUQP>(" << bv[1] << ", " << bv[0] << ")." << std::endl;
		exit(1);
	}

	int fd[2], gd[2], fs, gs;
	mpz_class *f, *g, x, y;
	std::vector< mpq_class > q;
	if (a.degree() >= b.degree()) {
		x =  SUPDUQP_to_mpq (a, q, gd);
		gs = (gd[0] + 1) * (gd[1] + 1);
		g = new mpz_class[gs];
		SUPDUQP_to_mpz (g, gd, q, x);
		q.clear();
		y =  SUPDUQP_to_mpq (b, q, fd);
		fs = (fd[0] + 1) * (fd[1] + 1);
		f = new mpz_class[fs];
		SUPDUQP_to_mpz (f, fd, q, y);
		q.clear();
	}
	else {
		x =  SUPDUQP_to_mpq (b, q, gd);
		gs = (gd[0] + 1) * (gd[1] + 1);
		g = new mpz_class[gs];
		SUPDUQP_to_mpz (g, gd, q, x);
		q.clear();
		y =  SUPDUQP_to_mpq (a, q, fd);
		fs = (fd[0] + 1) * (fd[1] + 1);
		f = new mpz_class[fs];
		SUPDUQP_to_mpz (f, fd, q, y);
		q.clear();
	}

	BiModularSubresultantChain<mpz_class> src = modularStableSubresultantChain (g, gd, f, fd);

	delete [] f;
	delete [] g;

	std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > r;
	for (int i = 0; i < src.n; ++i) {
		SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> s;
		s.setVariableName(av[1]);
		for (int j = 0; j <= src.deg[i][1]; ++j) {
			int k = src.deg[i][0] + 1;
			DenseUnivariateRationalPolynomial e (k);
			e.setVariableName(av[0]);
			for (int l = 0; l <= src.deg[i][0]; ++l)
				e.setCoefficient(l, src.coef[i][j*k+l]);
			if (!e.isZero()) { s.setCoefficient(j, e); }
		}
		r.push_back(s);
	}
	int k = r.size() - 1;
	while (r[k].isZero()) {
		r.pop_back();
		k = r.size() - 1;
	}
	
	return r;
}
