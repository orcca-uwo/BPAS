#include "RationalNumberPolynomial/urpolynomial.h"
#include "../../include/IntegerPolynomial/uzpolynomial.h"

typedef mpq_class lfixq;

RationalNumber DenseUnivariateRationalPolynomial::evaluate(const RationalNumber& x) const {
	if (curd) {
		RationalNumber px = coef[curd];
		for (int i = curd-1; i > -1; --i)
			px = px * x + coef[i];
		return px;
	}
	return coef[0];
}

void DenseUnivariateRationalPolynomial::reciprocal() {
	for (int i = 0; i < (curd+1)/2; ++i) {
		lfixq elem = coef[i];
		coef[i] = coef[curd-i];
		coef[curd-i] = elem;
	}
	resetDegree();
}

void DenseUnivariateRationalPolynomial::homothetic(int k) {
        for (int i = 0; i <= curd; ++i)
                coef[i] <<= (curd - i) * k;
}

void DenseUnivariateRationalPolynomial::scaleTransform(int k) {
        for (int i = 0; i <= curd; ++i)
                coef[i] <<= k * i;
}

void DenseUnivariateRationalPolynomial::negativeVariable() {
        for (int i = 0; i <= curd; ++i) {
                if (i%2)
                        coef[i] = -coef[i];
        }
}

void DenseUnivariateRationalPolynomial::negate() {
	for (int i = 0; i <= curd; ++i)
		coef[i] = -coef[i];
}

bool DenseUnivariateRationalPolynomial::isConstantTermZero() const {
    return (coef[0] == 0);
}

bool DenseUnivariateRationalPolynomial::divideByVariableIfCan() {
        if (coef[0] != 0)
                return 0;
        else {
                curd--;
                for (int i = 0; i <= curd; ++i)
                        coef[i] = coef[i+1];
                return 1;
        }
}

bool DenseUnivariateRationalPolynomial::isEqual(const DenseUnivariateRationalPolynomial& q) const {
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

void DenseUnivariateRationalPolynomial::resetDegree() {
    for (int i = curd; i > 0; --i) {
        if (coef[i] == 0)
            curd = i - 1;
        else { break; }
    }
}

void DenseUnivariateRationalPolynomial::pomopo(const lfixq c, const lfixq t, const DenseUnivariateRationalPolynomial& b) {
    if (c == 1) {
        for (int i = curd, j = b.curd; j > -1; --i, --j) {
            lfixq elem = coef[i] + t * b.coef[j];
            coef[i]  = elem;
        }
    }
    else {
        for (int i = curd, j = b.curd; i > -1; --i, --j) {
            lfixq elem = coef[i] * c;
            if (j > -1)
                elem += t * b.coef[j];
            coef[i] = elem;
        }
    }
    resetDegree();
}

void DenseUnivariateRationalPolynomial::print(std::ostream &out) const {
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

ExpressionTree DenseUnivariateRationalPolynomial::convertToExpressionTree() const {
    //TODO
    std::cerr << "BPAS ERROR: DUQP::convertToExpressionTree NOT YET IMPLEMENTED" << std::endl;
    return ExpressionTree();
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::operator^ (long long int e) const {
    DenseUnivariateRationalPolynomial res;
    res.name = name;
    res.one();
    unsigned long int q = e / 2, r = e % 2;
    DenseUnivariateRationalPolynomial power2 = *this * *this;
    for (int i = 0; i < q; ++i)
        res *= power2;
    if (r) { res *= *this; }
    return res;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::operator<< (int k) const {
    int s = curd + k + 1;
    DenseUnivariateRationalPolynomial r;
    r.n = s;
    r.curd = s - 1;
    r.name = name;
    r.coef = new lfixq[s];
    for (int i = 0; i < k; ++i)
        r.coef[i] = 0;
    for (int i = k; i < s; ++i)
        r.coef[i] = coef[i-k];
    return r;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::operator>> (int k) const {
    DenseUnivariateRationalPolynomial r;
    r.name = name;
    int s = curd - k + 1;
    if (s > 0) {
        r.n = s;
        r.curd = s - 1;
        delete [] r.coef;
        r.coef = new lfixq[s];
        for (int i = 0; i < s; ++i)
            r.coef[i] = coef[i+k];
    }
    return r;
}
    

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::operator+ (const DenseUnivariateRationalPolynomial& b) const {
    if (!curd) { return (b + coef[0]); }
    if (!b.curd) { return (*this + b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to add between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int size = (curd > b.curd)? curd+1 : b.curd+1;
    DenseUnivariateRationalPolynomial res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new lfixq[size];
    for (int i = 0; i < size; ++i) {
        lfixq elem = 0;
        if (i <= curd)
            elem += coef[i];
        if (i <= b.curd)
            elem += b.coef[i];
        res.coef[i] = elem;
    }
    res.resetDegree();
    return res;
}

void DenseUnivariateRationalPolynomial::add(const DenseUnivariateRationalPolynomial& b) {
    for (int i = curd; i >= 0; --i) {
        if (i <= b.curd)
            coef[i] += b.coef[i];
    }
    resetDegree();
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::operator- (const DenseUnivariateRationalPolynomial& b) const {
    if (!curd) { return (coef[0] - b); }
    if (!b.curd) { return (*this - b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to subtract between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    
    int size = (curd > b.curd)? curd+1 : b.curd+1;
    DenseUnivariateRationalPolynomial res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new lfixq[size];
    for (int i = 0; i < size; ++i) {
        lfixq elem = 0;
        if (i <= curd)
            elem = coef[i];
        if (i <= b.curd)
            elem -= b.coef[i];
        res.coef[i] = elem;
    }
    res.resetDegree();
    return res;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::operator- () const {
    DenseUnivariateRationalPolynomial res(curd+1);
    res.name = name;
    res.curd = curd;
    for (int i = 0; i <= curd; ++i)
        res.coef[i] = -coef[i];
    return res;
}

void DenseUnivariateRationalPolynomial::subtract(const DenseUnivariateRationalPolynomial& b) {
    for (int i = curd; i >= 0; --i) {
        if (i <= b.curd)
            coef[i] -= b.coef[i];
    }
    resetDegree();
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::operator* (const DenseUnivariateRationalPolynomial& b) const {
    if (!curd) { return (b * coef[0]); }
    if (!b.curd) { return (*this * b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to multiply between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int d = curd + 1;
    int m = b.curd + 1;
    int size = curd + m;
    DenseUnivariateRationalPolynomial res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new lfixq[size];
    
    lfixz aden = 1, bden = 1;
    for (int i = 0; i < d; ++i)
        aden *= coef[i].get_den();
    for (int i = 0; i < m; ++i)
        bden *= b.coef[i].get_den();
    
    lfixz* acoef = new lfixz[d];
    for (int i = 0; i < d; ++i)
        acoef[i] = aden * coef[i];
    lfixz* bcoef = new lfixz[m];
    for (int i = 0; i < m; ++i)
        bcoef[i] = bden * b.coef[i];
    
    lfixz* mul = new lfixz[size];
    univariateMultiplication(mul, acoef, d, bcoef, m);
    lfixz den = aden * bden;
    for (int i = 0; i < size; ++i) {
        lfixq elem = mul[i];
        elem /= den;
        res.coef[i] = elem;
    }
    res.resetDegree();
    
    delete [] acoef;
    delete [] bcoef;
    delete [] mul;
    return res;
}

DenseUnivariateRationalPolynomial& DenseUnivariateRationalPolynomial::operator*= (const RationalNumber& e) {
    lfixq c (e.get_mpq());
    *this *= c;
    return *this;
}

DenseUnivariateRationalPolynomial& DenseUnivariateRationalPolynomial::operator*= (const mpq_class& e) {
    if (e != 0 && e != 1) {
        for (int i = 0; i <= curd; ++i)
            coef[i] *= e;
    }
    else if (e == 0) { zero(); }
    return *this;
}

DenseUnivariateRationalPolynomial& DenseUnivariateRationalPolynomial::operator*= (const sfixn& e) {
    if (e != 0 && e != 1) {
        for (int i = 0; i <= curd; ++i)
            coef[i] *= e;
    }
    else if (e == 0) { zero(); }
    return *this;
}

DenseUnivariateRationalPolynomial& DenseUnivariateRationalPolynomial::operator/= (const DenseUnivariateRationalPolynomial& b) {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUQP." << std::endl;
        exit(1);
    }
    if (!b.curd)
        return (*this /= b.coef[0]);
    if (!curd) {
        coef[0] = 0;
        return *this;
    }
    
    if (name != b.name) {
        std::cout << "BPAS: error, trying to exact divide between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    DenseUnivariateRationalPolynomial rem(*this);
    zeros();
    while (rem.curd >= b.curd) {
        lfixq lc = rem.coef[rem.curd] / b.coef[b.curd];
        int diff = rem.curd - b.curd;
        rem.pomopo(1, -lc, b);
        coef[diff] = lc;
    }
    resetDegree();
    if (!rem.isZero()) {
        std::cout << "BPAS: error, not exact division from DUQP." << std::endl;
        exit(1);
    }
    return *this;
}

DenseUnivariateRationalPolynomial& DenseUnivariateRationalPolynomial::operator/= (const RationalNumber& e) {
    lfixq c (e.get_mpq());
    *this /= c;
    return *this;
}

DenseUnivariateRationalPolynomial& DenseUnivariateRationalPolynomial::operator/= (const mpq_class& e) {
    if (e == 0) {
        std::cout << "BPAS: error, dividend is zero from DUQP." << std::endl;
        exit(1);
    }
    else if (e != 1) {
    	for (int i = 0; i <= curd; ++i)
    	    coef[i] /= e;
    }
    return *this;
}

DenseUnivariateRationalPolynomial operator/ (const mpq_class& c, const DenseUnivariateRationalPolynomial& p) {
    if (p.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUQP." << std::endl;
        exit(1);
    }
    
    DenseUnivariateRationalPolynomial q;
    q.name = p.name;
    q.curd = 0;
    q.n = 1;
    q.coef = new lfixq[1];
    
    if (p.isConstant())
        q.coef[0] = c / p.coef[0];
    else
        q.coef[0] = 0;
    return q;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::monicDivide(const DenseUnivariateRationalPolynomial& b) {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUQP." << std::endl;
        exit(1);
    }
    else if (b.coef[b.curd] != 1) {
        std::cout << "BPAS: error, leading coefficient is not one in monicDivide() from DUQP." << std::endl;
        exit(1);
    }
    if (!b.curd) {
        DenseUnivariateRationalPolynomial r (*this);
        zero();
        return r;
    }
    if (!curd) {
        DenseUnivariateRationalPolynomial r;
        r.name = name;
        return r;
    }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to monic divide between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int size = curd - b.curd + 1;
    DenseUnivariateRationalPolynomial quo(size);
    quo.curd = size - 1;
    quo.name = name;
    while (curd >= b.curd) {
        lfixq lc = coef[curd];
        int diff = curd - b.curd;
        pomopo(1, -lc, b);
        quo.coef[diff] = lc;
    }
    quo.resetDegree();
    return quo;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::monicDivide(const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* rem) const {
    *rem = *this;
    return rem->monicDivide(b);
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::lazyPseudoDivide (const DenseUnivariateRationalPolynomial& b, RationalNumber* c, RationalNumber* d) {
    if (d == NULL)
        d = new RationalNumber;
    int da = curd, db = b.curd;
    if (b.isZero() || !db) {
        std::cout << "BPAS: error, dividend is zero or constant." << std::endl;
        exit(1);
    }
    *c = 1, *d = 1;
    if (!curd) {
        DenseUnivariateRationalPolynomial r;
        r.name = name;
        return r;
    }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to pseudo divide between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    if (da < db) {
        DenseUnivariateRationalPolynomial r;
        r.name = name;
        return r;
    }
    
    int size = curd - b.curd + 1;
    DenseUnivariateRationalPolynomial quo(size);
    quo.curd = size - 1;
    quo.name = name;
    int e = 0, diff = da - db;
    mpq_class blc = b.coef[b.curd];
    while (curd >= b.curd) {
        mpq_class lc = coef[curd];
        int k = curd - b.curd;
        *c *= RationalNumber(blc);
        e++;
        pomopo(blc, -coef[curd], b);
        quo.coef[k] = lc;
    }
    quo.resetDegree();
    for (int i = e; i <= diff; ++i)
        *d *= RationalNumber(blc);
    return quo;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::lazyPseudoDivide (const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* rem, RationalNumber* c, RationalNumber* d) const {
    *rem = *this;
    return rem->lazyPseudoDivide(b, c, d);
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::pseudoDivide (const DenseUnivariateRationalPolynomial& b, RationalNumber* d) {
    RationalNumber c;
    if (d == NULL)
        d = new RationalNumber;
    DenseUnivariateRationalPolynomial quo = lazyPseudoDivide(b, &c, d);
    quo *= *d;
    *this *= *d;
    *d *= c;
    return quo;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::pseudoDivide (const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* rem, RationalNumber* d) const {
    RationalNumber c;
    DenseUnivariateRationalPolynomial quo = lazyPseudoDivide(b, rem, &c, d);
    quo *= *d;
    *rem *= *d;
    *d *= c;
    return quo;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::halfExtendedEuclidean (const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* g) const {
    if (g == NULL)
        g = new DenseUnivariateRationalPolynomial;
    *g = *this;
    
    DenseUnivariateRationalPolynomial a1, b1, b2(b);
    a1.name = name;
    b1.name = b.name;
    a1.coef[0] = 1;
    while (!b2.isZero()) {
        DenseUnivariateRationalPolynomial q, r;
        q.name = r.name = name;
        lfixq e = b2.coef[b2.curd];
        if (e != 1) {
            b2 /= e;
            *g /= e;
        }
        q = g->monicDivide(b2, &r);
        if (e != 1) {
            *g = b2 * e;
            b2 = r * e;
        }
        else {
            *g = b2;
            b2 = r;
        }
        
        r = a1;
        r -= q * b1;
        a1 = b1;
        b1 = r;
    }
    
    a1 /= g->coef[g->curd];
    *g /= g->coef[g->curd];
    
    return a1;
}

void DenseUnivariateRationalPolynomial::diophantinEquationSolve(const DenseUnivariateRationalPolynomial& a, const DenseUnivariateRationalPolynomial& b, DenseUnivariateRationalPolynomial* s, DenseUnivariateRationalPolynomial* t) const {
    DenseUnivariateRationalPolynomial f(*this), g, q, r;
    *s = a.halfExtendedEuclidean(b, &g);
    if (g.coef[g.curd] != 1) {
        f /= g.coef[g.curd];
        g /= g.coef[g.curd];
    }
    q = f.monicDivide(g, &r);
    if (!r.isZero()) {
        std::cout << "BPAS: error, " << *this << " is not in the ideal (" << a << ", " << b << ") from DUQP." << std::endl;
        exit(1);
    }
    *s *= q;
    
    DenseUnivariateRationalPolynomial b2;
    lfixq e = b.coef[b.curd];
    if (e != 1) { b2 = b / e; }
    if (s->curd >= b.curd) {
        *s /= e;
        s->monicDivide(b2, &r);
        *s = r * e;
    }
    
    g = *this;
    g -= *s * a;
    if (e != 1) { g /= e; }
    *t = g.monicDivide(b2);
}

void DenseUnivariateRationalPolynomial::differentiate(int k) {
    if (k <= 0) { return; }
    for (int i = k; i <= curd; ++i) {
        coef[i-k] = coef[i];
        for (int j = 0; j < k; ++j)
            coef[i-k] *= (i - j);
    }
    curd -= k;
    resetDegree();
}

void DenseUnivariateRationalPolynomial::integrate() {
    DenseUnivariateRationalPolynomial b;
    b.name = name;
    b.n = curd+2;
    b.coef = new lfixq[b.n];
    b.coef[0] = 0;
    for (int i = 0; i <= curd; ++i)
        b.coef[i+1] = coef[i] / (i + 1);
    b.resetDegree();
    *this = b;
}

Factors<DenseUnivariateRationalPolynomial> DenseUnivariateRationalPolynomial::squareFree() const {
    std::vector<DenseUnivariateRationalPolynomial> sf;
    if (!curd)
        sf.push_back(*this);
    else if (curd == 1) {
        DenseUnivariateRationalPolynomial t;
        t.name = name;
        t.coef[0] = coef[curd];
        sf.push_back(t);
        t = *this / t.coef[0];
        sf.push_back(t);
    }
    else {
        DenseUnivariateRationalPolynomial a (*this), b(*this);
        b.differentiate(1);
        DenseUnivariateRationalPolynomial g = a.gcd(b);
        DenseUnivariateRationalPolynomial x = a / g;
        DenseUnivariateRationalPolynomial y = b / g;
        DenseUnivariateRationalPolynomial z = -x;
        z.differentiate(1);
        z += y;
        
        while (!z.isZero()) {
            g = x.gcd(z);
            sf.push_back(g);
            x /= g;
            y = z / g;
            z = -x;
            z.differentiate(1);
            z += y;
        }
        sf.push_back(x);
        
        lfixq e = 1;
        for (int i = 0; i < sf.size(); ++i) {
            e *= sf[i].coef[sf[i].curd];
            sf[i] /= sf[i].coef[sf[i].curd];
        }
        DenseUnivariateRationalPolynomial t;
        t.name = name;
        t.coef[0] = e;
        sf.insert(sf.begin(), t);
    }

    Factors<DenseUnivariateRationalPolynomial> f;
    f.setRingElement(sf[0]);
    for (int i = 1; i < sf.size(); ++i) {
        f.addFactor(sf[i], i);
    }
    return f;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::euclideanGCD (const DenseUnivariateRationalPolynomial& q) const {
    DenseUnivariateRationalPolynomial a, b;
    if (curd < q.curd) {
        a = q;
        b = *this;
    }
    else {
        a = *this;
        b = q;
    }
    
    while (!b.isZero()) {
        lfixq lc = b.coef[b.curd];
        b /= lc;
        DenseUnivariateRationalPolynomial r;
        r.name = name;
        a.monicDivide(b, &r);
        a = b * lc;
        b = r;
    }
    return a / a.coef[a.curd];
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::modularGCD (const DenseUnivariateRationalPolynomial& q) const {
	DenseUnivariateIntegerPolynomial a (curd+1), b (q.curd+1);
	mpz_class e = 1, o = 1;
	for (int i = 0; i <= curd; ++i) {
		if (coef[i] != 0)
			e *= coef[i].get_den();
	}
	for (int i = 0; i <= curd; ++i)
		a.setCoefficient(i, mpz_class(coef[i]) * e);
	for (int i = 0; i <= q.curd; ++i) {
		if (q.coef[i] != 0)
			o *= q.coef[i].get_den();
	}
	for (int i = 0; i <= q.curd; ++i)
		b.setCoefficient(i, mpz_class(q.coef[i]) * o);

	DenseUnivariateIntegerPolynomial c = a.gcd(b,1);

	int d = c.degree().get_si();
	DenseUnivariateRationalPolynomial r (d+1);
	r.curd = d;
	r.name = name;
	for (int i = 0; i < d; ++i)
		r.coef[i] = mpq_class(c.coefficient(i).get_mpz()) / c.coefficient(d).get_mpz();
	r.coef[d] = 1;

	return r;
}

DenseUnivariateRationalPolynomial DenseUnivariateRationalPolynomial::gcd (const DenseUnivariateRationalPolynomial& q, int type) const {
    if (isZero()) { return q; }
    if (q.isZero()) { return *this; }
    if (!curd || !q.curd) {
        DenseUnivariateRationalPolynomial h (1);
        h.coef[0] = 1;
        h.name = name;
        return h;
    }

    if (name != q.name) {
        std::cout << "BPAS: error, trying to compute GCD between Q[" << name << "] and Q[" << q.name << "]." << std::endl;
        exit(1);
    }

	DenseUnivariateRationalPolynomial r;
	if (!type)
		r = euclideanGCD(q);
	else 
		r = modularGCD(q);

	return r;
}
