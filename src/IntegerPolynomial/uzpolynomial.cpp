#include "../../include/IntegerPolynomial/uzpolynomial.h"
#include "../../include/DyadicRationalNumber/Multiplication/multiplication.h"

bool DenseUnivariateIntegerPolynomial::isEqual(const DenseUnivariateIntegerPolynomial& q) const {
    if (curd && q.curd && (name != q.name))
        return 0;
    if (curd != q.curd)
        return 0;
    for (int i = 0; i <= curd; ++i) {
        if (coef[i] != q.coef[i])
            return 0;
    }
    return 1;
}

void DenseUnivariateIntegerPolynomial::pomopo(const mpz_class c, const mpz_class t, const DenseUnivariateIntegerPolynomial& b) {
    if (c == 1) {
        for (int i = curd, j = b.curd; j > -1; --i, --j) {
            coef[i] += t * b.coef[j];
	}
    }
    else {
        for (int i = curd, j = b.curd; i > -1; --i, --j) {
            mpz_class elem = coef[i] * c;
            if (j > -1)
                elem += t * b.coef[j];
            coef[i] = elem;
        }
    }
    resetDegree();
}

void DenseUnivariateIntegerPolynomial::resetDegree() {
    for (int i = curd; i > 0; --i) {
        if (coef[i] != 0) { break;}
        else { curd = i - 1; }
    }
}

void DenseUnivariateIntegerPolynomial::print(std::ostream &out) const {
    bool isFirst = 0;
    for (int i = 0; i <= this->curd; ++i) {
        if (this->coef[i] != 0) {
            if (isFirst && this->coef[i] > 0)
                out << "+";
            if (i) {
                if (this->coef[i] != 1 && this->coef[i] != -1)
                    out << this->coef[i] << "*";
                else if (this->coef[i] < 0)
                    out << "-";
                out << this->name;
                if (i > 1)
                    out << "^" << i;
            }
            else { out << this->coef[i]; }
            isFirst = 1;
        }
    }
    if (!isFirst) { out << "0"; }
}

ExpressionTree DenseUnivariateIntegerPolynomial::convertToExpressionTree() const {
    //TODO 
    std::cerr << "BPAS ERROR: DUZP::convertToExpressionTree NOT YET IMPLEMENTED" << std::endl;
    return ExpressionTree(); 
}


DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::operator^ (long long int e) const {
    DenseUnivariateIntegerPolynomial res;
    res.name = name;
    res.one();
    unsigned long int q = e / 2, r = e % 2;
    DenseUnivariateIntegerPolynomial power2 = *this * *this;
    for (int i = 0; i < q; ++i)
        res *= power2;
    if (r) { res *= *this; }
    return res;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::operator<< (int k) const {
    int s = curd + k + 1;
    DenseUnivariateIntegerPolynomial r;
    r.n = s;
    r.curd = s - 1;
    r.name = name;
    r.coef = new mpz_class[s];
    for (int i = 0; i < k; ++i)
        r.coef[i] = 0;
    for (int i = k; i < s; ++i)
        r.coef[i] = coef[i-k];
    return r;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::operator>> (int k) const {
    DenseUnivariateIntegerPolynomial r;
    r.name = name;
    int s = curd - k + 1;
    if (s > 0) {
        r.n = s;
        r.curd = s - 1;
        delete [] r.coef;
        r.coef = new mpz_class[s];
        for (int i = 0; i < s; ++i)
            r.coef[i] = coef[i+k];
        return r;
    }
    return r;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::operator+ (const DenseUnivariateIntegerPolynomial& b) const {
    if (!curd) { return (b + coef[0]); }
    if (!b.curd) { return (*this + b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to add between Z[" << name << "] and Z[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int size = (curd > b.curd)? curd+1 : b.curd+1;
    DenseUnivariateIntegerPolynomial res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new mpz_class[size];
    for (int i = 0; i < size; ++i) {
        mpz_class elem = 0;
        if (i <= curd)
            elem += coef[i];
        if (i <= b.curd)
            elem += b.coef[i];
        res.coef[i] = elem;
    }
    res.resetDegree();
    return res;
}

void DenseUnivariateIntegerPolynomial::add(const DenseUnivariateIntegerPolynomial& b) {
    for (int i = curd; i >= 0; --i) {
        if (i <= b.curd)
            coef[i] += b.coef[i];
    }
    resetDegree();
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::operator- (const DenseUnivariateIntegerPolynomial& b) const {
    if (!curd) { return (coef[0] - b); }
    if (!b.curd) { return (*this - b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to subtract between Z[" << name << "] and Z[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int size = (curd > b.curd)? curd+1 : b.curd+1;
    DenseUnivariateIntegerPolynomial res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new mpz_class[size];
    for (int i = 0; i < size; ++i) {
        mpz_class elem = 0;
        if (i <= curd)
            elem = coef[i];
        if (i <= b.curd)
            elem -= b.coef[i];
        res.coef[i] = elem;
    }
    res.resetDegree();
    return res;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::operator- () const {
    DenseUnivariateIntegerPolynomial res(curd+1);
    res.name = name;
    res.curd = curd;
    for (int i = 0; i <= curd; ++i)
        res.coef[i] = -coef[i];
    return res;
}

void DenseUnivariateIntegerPolynomial::subtract(const DenseUnivariateIntegerPolynomial& b) {
    for (int i = curd; i >= 0; --i) {
        if (i <= b.curd)
            coef[i] -= b.coef[i];
    }
    resetDegree();
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::operator* (const DenseUnivariateIntegerPolynomial& b) const {
    if (!curd) { return (b * coef[0]); }
    if (!b.curd) { return (*this * b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to multiply between Z[" << name << "] and Z[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int d = curd + 1;
    int m = b.curd + 1;
    int size = curd + m;
    DenseUnivariateIntegerPolynomial res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new mpz_class[size];
    
	univariateMultiplication(res.coef, coef, d, b.coef, m);

    res.resetDegree();
    return res;
}

bool DenseUnivariateIntegerPolynomial::isDivide(DenseUnivariateIntegerPolynomial f, DenseUnivariateIntegerPolynomial g, const DenseUnivariateIntegerPolynomial& h) const {
	while (f.curd >= h.curd) {
		mpz_class lc = f.coef[f.curd] / h.coef[h.curd];
		f.pomopo(1, -lc, h);
	}
	if (!f.isZero()) { return 0; }
	while (g.curd >= h.curd) {
		mpz_class lc = g.coef[g.curd] / h.coef[h.curd];
		g.pomopo(1, -lc, h);
	}
	return g.isZero();
}

DenseUnivariateIntegerPolynomial& DenseUnivariateIntegerPolynomial::operator/= (const DenseUnivariateIntegerPolynomial& b) {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUZP." << std::endl;
        exit(1);
    }
    if (!b.curd)
        return (*this /= b.coef[0]);
    if (!curd) {
        coef[0] = 0;
        return *this;
    }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to exact divide between Z[" << name << "] and Z[" << b.name << "]." << std::endl;
        exit(1);
    }

    DenseUnivariateIntegerPolynomial rem (*this);
    zeros();
    curd -= b.curd;
    while (rem.curd >= b.curd) {
        mpz_class lc = rem.coef[rem.curd] / b.coef[b.curd];
        int diff = rem.curd - b.curd;
        rem.pomopo(1, -lc, b);
        coef[diff] = lc;
    }
    resetDegree();
    if (!rem.isZero()) {
        std::cout << "BPAS: error, not exact division from DUZP." << std::endl;
        exit(1);
    }
    return *this;    
}

DenseUnivariateIntegerPolynomial& DenseUnivariateIntegerPolynomial::operator/= (const mpz_class& e) {
    if (e == 0) {
        std::cout << "BPAS: error, dividend is zero from DUZP." << std::endl;
        exit(1);
    }
    else if (e != 1) {
    	for (int i = 0; i <= curd; ++i)
    	    coef[i] /= e;
    	resetDegree();
    }
    return *this;
}

DenseUnivariateIntegerPolynomial operator/ (const mpz_class& e, const DenseUnivariateIntegerPolynomial& p) {
    if (p.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUZP." << std::endl;
        exit(1);
    }
    
    DenseUnivariateIntegerPolynomial q;
    q.name = p.name;
    q.curd = 0;
    q.n = 1;
    q.coef = new mpz_class[1];
    if (p.isConstant())
        q.coef[0] = e / q.coef[0];
    else
        q.coef[0] = 0;
    return q;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::monicDivide(const DenseUnivariateIntegerPolynomial& b) {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUZP." << std::endl;
        exit(1);
    }
    else if (b.coef[b.curd] != 1) {
        std::cout << "BPAS: error, leading coefficient is not one in monicDivide() from DUZP." << std::endl;
        exit(1);
    }
    if (!b.curd) {
        DenseUnivariateIntegerPolynomial r (*this);
        zero();
        return r;
    }
    if (!curd) {
        DenseUnivariateIntegerPolynomial r;
        r.name = name;
        return r;
    }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to monic divide between Z[" << name << "] and Z[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int size = curd - b.curd + 1;
    DenseUnivariateIntegerPolynomial quo(size);
    quo.curd = size - 1;
    quo.name = name;
    while (curd >= b.curd) {
        mpz_class lc = coef[curd];
        int diff = curd - b.curd;
        pomopo(1, -lc, b);
        quo.coef[diff] = lc;
    }
    quo.resetDegree();
    return quo;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::lazyPseudoDivide (const DenseUnivariateIntegerPolynomial& b, Integer* c, Integer* d) {
    if (d == NULL)
        d = new Integer;
    int da = curd, db = b.curd;
    if (b.isZero() || !db) {
        std::cout << "BPAS: error, dividend is zero or constant." << std::endl;
        exit(1);
    }
    *c = Integer(1), *d = Integer(1);
    if (!curd) {
        DenseUnivariateIntegerPolynomial r;
        r.name = name;
        return r;
    }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to pseudo divide between Z[" << name << "] and Z[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    if (da < db) {
        DenseUnivariateIntegerPolynomial r;
        r.name = name;
        return r;
    }
    
    int size = curd - b.curd + 1;
    DenseUnivariateIntegerPolynomial quo(size);
    quo.curd = size - 1;
    quo.name = name;
    int e = 0, diff = da - db;
    mpz_class blc = b.coef[b.curd];
    while (curd >= b.curd) {
        mpz_class lc = coef[curd];
        int k = curd - b.curd;
        *c *= Integer(blc);
        e++;
        pomopo(blc, -coef[curd], b);
        quo.coef[k] = lc;
    }
    quo.resetDegree();
    for (int i = e; i <= diff; ++i)
        *d *= Integer(blc);
    return quo;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::pseudoDivide (const DenseUnivariateIntegerPolynomial& b, Integer* d) {
    Integer c;
    if (d == NULL)
        d = new Integer;
    DenseUnivariateIntegerPolynomial quo = lazyPseudoDivide(b, &c, d);
    quo *= *d;
    *this *= *d;
    *d *= c;
    return quo;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::pseudoDivide (const DenseUnivariateIntegerPolynomial& b, DenseUnivariateIntegerPolynomial* rem, Integer* d) const {
    Integer c;
    DenseUnivariateIntegerPolynomial quo = lazyPseudoDivide(b, rem, &c, d);
    quo *= *d;
    *rem *= *d;
    *d *= c;
    return quo;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::euclideanGCD (const DenseUnivariateIntegerPolynomial& q) const {
    int da, db, maxs = (curd < q.curd)? q.curd+1 : curd+1;
    mpq_class* a = new mpq_class[maxs];
    mpq_class* b = new mpq_class[maxs];
    if (curd < q.curd) {
        for (int i = 0; i <= q.curd; ++i)
            a[i] = mpq_class(q.coef[i]);
        for (int i = 0; i <= curd; ++i)
            b[i] = mpq_class(coef[i]);
        da = q.curd, db = curd;
    }
    else {
        for (int i = 0; i <= curd; ++i)
            a[i] = mpq_class(coef[i]);
        for (int i = 0; i <= q.curd; ++i)
            b[i] = mpq_class(q.coef[i]);
        da = curd, db = q.curd;
    }

    while (db > 0) {
	int k = 0;
        while (da >= db) {
	    k = da - 1;
	    bool isFirst = 0;
	    mpq_class e = a[da] / b[db];
            for (int i = da-1, j = db-1; j > -1; --i, --j) {
		a[i] -= e * b[j];
		if (!isFirst && a[i] != 0) {
			k = i;
			isFirst = 1;
		}
	    }
	    for (int i = k; i > 0; --i) {
		if (a[i] != 0) { break; }
		else { k = i - 1; }
	    }
	    da = k;
        }

        // swap a and b
	mpq_class* t = a;
	a = b;
	b = t;
        da = db;
	db = k;
    }

    if (b[0] != 0) {
	delete [] a;
	delete [] b;
	DenseUnivariateIntegerPolynomial r (1);
	r.name = name;
	r.coef[0] = 1;
	return r;
    }

    mpz_class lc = 1;
    for (int i = 0; i <= da; ++i)
    	lc *= a[i].get_den();

    DenseUnivariateIntegerPolynomial r (da + 1);
    r.name = name;
    r.curd = da;
    for (int i = 0; i <= da; ++i)
        r.coef[i] = a[i] * lc;
    r /= r.content();
    
    delete [] a;
    delete [] b;
    return r;
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::modularGCD (const DenseUnivariateIntegerPolynomial& g) const {
	int d, s = 0;
	sfixn *fp, *gp;
	if (curd > g.curd) {
		d = g.curd;
		fp = new sfixn[g.curd+1];
		gp = new sfixn[curd+1];
	}
	else {
		d = curd;
		fp = new sfixn[curd+1];
		gp = new sfixn[g.curd+1];
	}

	mpz_class b, m = 1;
	mpz_gcd (b.get_mpz_t(), coef[curd].get_mpz_t(), g.coef[g.curd].get_mpz_t());

	DenseUnivariateIntegerPolynomial hp (1);
	hp.name = name;

	while(1) {
		int k;
		// next prime
		sfixn p = nextprime(&s, m * b);

		mpz_class t (p), e, a;
		mpz_mod (a.get_mpz_t(), b.get_mpz_t(),  t.get_mpz_t());
		// to guarantee deg(gp) >= deg(fp)
		if (curd > g.curd) {
			// f' = g mod p
			for (int i = 0; i <= g.curd; ++i) {
				if (g.coef[i] != 0) {
					mpz_mod (e.get_mpz_t(), g.coef[i].get_mpz_t(), t.get_mpz_t());
					fp[i] = e.get_mpz_t()->_mp_d[0];
				}
				else { fp[i] = 0; }
			}
			// g' = f mod p
			for (int i = 0; i <= curd; ++i) {
				if (coef[i] != 0) {
					mpz_mod (e.get_mpz_t(), coef[i].get_mpz_t(), t.get_mpz_t());
					gp[i] = e.get_mpz_t()->_mp_d[0];
				}
				else { gp[i] = 0; }
			}
			k = curd;
			monicGCD(gp, &k, fp, g.curd, p, a.get_mpz_t()->_mp_d[0]);
		}
		else {
			// f' = f mod p
			for (int i = 0; i <= curd; ++i) { 
				if (coef[i] != 0) {
					mpz_mod (e.get_mpz_t(), coef[i].get_mpz_t(), t.get_mpz_t());
					fp[i] = e.get_mpz_t()->_mp_d[0];
				}
				else { fp[i] = 0; }
			}
			// g' = g mod p
			for (int i = 0; i <= g.curd; ++i) {
				if (g.coef[i] != 0) {
					mpz_mod (e.get_mpz_t(), g.coef[i].get_mpz_t(), t.get_mpz_t());
					gp[i] = e.get_mpz_t()->_mp_d[0];
				}
				else { gp[i] = 0; }
			}
			k = g.curd;
			monicGCD(gp, &k, fp, curd, p, a.get_mpz_t()->_mp_d[0]);
		}

		if (!k) {
			// deg(gp) == 0, then return 1
			DenseUnivariateIntegerPolynomial r (1);
			r.coef[0] = 1;
			r.name = name;
            delete [] fp;
            delete [] gp;
			return r;
		}
		else if (k < d) {
			// deg(gp) < d, then re-initialize
			m = 1;
			d = k;
			hp.zero();
		}
		if (k == d) {
			// deg(gp) == d, then do the work!
			DenseUnivariateIntegerPolynomial h (k+1);
			h.curd = k;

			bool isEqual = 1;
			if (hp.isZero()) {
				// hp is zero, h = gp
				// coefficients of h in [-(p-1)/2, (p-1)/2]
				sfixn halfp = (p-1) >> 1;
				for (int i = 0; i <= d; ++i) {
					h.coef[i] = gp[i];
					if (h.coef[i] > halfp) { h.coef[i] -= p; }
				}
				m *= p;
			}
			else {
				// t = e * p + a * m = 1
				mpz_gcdext (t.get_mpz_t(), e.get_mpz_t(), a.get_mpz_t(), t.get_mpz_t(), m.get_mpz_t());
				t = m;
				m *= p;
				// h = combine(p,m)(gp, hp)
				// coefficients of h in [-(mp-1)/2, (mp-1)/2]
				k = 0;
				mpz_class halfm = (m-1) >> 1;
				for (int i = 0; i <= d; ++i) {
					//h.coef[i] = gp[i] * e * p + hp.coef[i] * a * t;
					h.coef[i] = (hp.coef[i] - gp[i]) * a;
					mpz_mod (h.coef[i].get_mpz_t(), h.coef[i].get_mpz_t(), t.get_mpz_t());
					h.coef[i] = gp[i] + h.coef[i] * p;
					if (h.coef[i] > halfm)
						h.coef[i] -= m;
					if (h.coef[i] != 0) { k = i; }

					if (h.coef[i] != hp.coef[i]) { isEqual = 0; }
				}
				h.curd = k;
			}

			if (isEqual) {
				DenseUnivariateIntegerPolynomial r = h / h.content();
				r.name = name;
				if (isDivide(*this, g, r)) { 
                    delete [] fp;
                    delete [] gp;
                    return r;
                }
			}
			hp = h;
		}
	}

	delete [] fp;
	delete [] gp;
    DenseUnivariateIntegerPolynomial r (1);
    r.coef[0] = 1;
    r.name = name;
    return r;    
}

DenseUnivariateIntegerPolynomial DenseUnivariateIntegerPolynomial::gcd (const DenseUnivariateIntegerPolynomial& b, int type) const {
	if (isZero()) { return b; }
	if (b.isZero()) { return *this; }
	if (!curd || !b.curd) {
		DenseUnivariateIntegerPolynomial h (1);
		h.coef[0] = 1;
		h.name = name;
		return h;
	}
	if (name != b.name) {
		std::cout << "BPAS: error, trying to compute GCD between Z[" << name << "] and Z[" << b.name << "]." << std::endl;
		exit(1);
	}

    DenseUnivariateIntegerPolynomial p(*this);
    DenseUnivariateIntegerPolynomial q(b);
	mpz_class e = p.content().get_mpz();
	p /= e;
	mpz_class o = q.content().get_mpz();
	q /= o;

	DenseUnivariateIntegerPolynomial r;
	if (!type)
		r = euclideanGCD(q);
	else 
		r = modularGCD(q);

	r *= e * o;

	return r;
}

void DenseUnivariateIntegerPolynomial::differentiate(int k) {
    if (k <= 0) { return; }
    for (int i = k; i <= curd; ++i) {
        coef[i-k] = coef[i];
        for (int j = 0; j < k; ++j)
            coef[i-k] *= (i - j);
    }
    curd -= k;
    resetDegree();
}

Integer DenseUnivariateIntegerPolynomial::evaluate(const Integer& x) const {
    if (curd) {
        mpz_class px = coef[curd];
        for (int i = curd-1; i > -1; --i)
            px = px * x.get_mpz() + coef[i];
        return px;
    }
    return coef[0];
}

Factors<DenseUnivariateIntegerPolynomial> DenseUnivariateIntegerPolynomial::squareFree() const {
    std::vector<DenseUnivariateIntegerPolynomial> sf;
    if (!curd)
        sf.push_back(*this);
    else if (curd == 1) {
        DenseUnivariateIntegerPolynomial t;
        t.name = name;
        t.coef[0] = coef[curd];
        sf.push_back(t);
        t = *this / t.coef[0];
        sf.push_back(t);
    }
    else {
        DenseUnivariateIntegerPolynomial a (*this), b(*this);
        b.differentiate(1);
        DenseUnivariateIntegerPolynomial g = a.gcd(b);
        g /= g.content();
        DenseUnivariateIntegerPolynomial x = a / g;
        DenseUnivariateIntegerPolynomial y = b / g;
        DenseUnivariateIntegerPolynomial z = -x;
        z.differentiate(1);
        z += y;

        while (!z.isZero()) {
            g = x.gcd(z);
            g /= g.content();
            sf.push_back(g);
            x /= g;
            y = z / g;
            z = -x;
            z.differentiate(1);
            z += y;
        }
        sf.push_back(x);
        
        mpz_class e = 1;
        for (int i = 0; i < sf.size(); ++i) {
            e *= sf[i].coef[sf[i].curd];
            sf[i] /= sf[i].coef[sf[i].curd];
        }
        DenseUnivariateIntegerPolynomial t;
        t.name = name;
        t.coef[0] = e;
        sf.insert(sf.begin(), t);
    }
    Factors<DenseUnivariateIntegerPolynomial> f;
    f.setRingElement(sf[0]);
    for (int i = 1; i < sf.size(); ++i) {
        f.addFactor(sf[i], i);
    }
    return f;
}
