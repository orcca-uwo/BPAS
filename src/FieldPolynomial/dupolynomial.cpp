#include "../../include/RingPolynomial/dupolynomial.h"

template <class Field>
Field DenseUnivariatePolynomial<Field>::evaluate(Field x) {
	if (curd) {
		Field px = coef[curd];
		for (int i = curd-1; i > -1; --i)
			px = px * x + coef[i];
		return px;
	}
	return coef[0];
}

template <class Field>
void DenseUnivariatePolynomial<Field>::resetDegree() {
    for (int i = curd; i > 0; --i) {
        if (coef[i] == 0)
            curd = i - 1;
        else { break; }
    }
}

template <class Field>
void DenseUnivariatePolynomial<Field>::reciprocal() {
	for (int i = 0; i < (curd+1)/2; ++i) {
		Field elem = coef[i];
		coef[i] = coef[curd-i];
		coef[curd-i] = elem;
	}
	resetDegree();
}

template <class Field>
void DenseUnivariatePolynomial<Field>::homothetic(int k) {
        for (int i = 0; i <= curd; ++i)
                coef[i] <<= (curd - i) * k;
}

template <class Field>
void DenseUnivariatePolynomial<Field>::scaleTransform(int k) {
        for (int i = 0; i <= curd; ++i)
                coef[i] <<= k * i;
}

template <class Field>
void DenseUnivariatePolynomial<Field>::negativeVariable() {
        for (int i = 0; i <= curd; ++i) {
                if (i%2)
                        coef[i] = -coef[i];
        }
}

template <class Field>
void DenseUnivariatePolynomial<Field>::negate() {
	for (int i = 0; i <= curd; ++i)
		coef[i] = -coef[i];
}

template <class Field>
bool DenseUnivariatePolynomial<Field>::isTrailingCoefficientZero() {
    return (coef[0] == 0);
}

template <class Field>
bool DenseUnivariatePolynomial<Field>::divideByVariableIfCan() {
        if (coef[0] != 0)
                return 0;
        else {
                curd--;
                for (int i = 0; i <= curd; ++i)
                        coef[i] = coef[i+1];
                return 1;
        }
}

template <class Field>
bool DenseUnivariatePolynomial<Field>::isEqual(DenseUnivariatePolynomial<Field>& q) {
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

template <class Field>
void DenseUnivariatePolynomial<Field>::pomopo(Field c, Field t, DenseUnivariatePolynomial<Field>& b) {
    if (c == 1) {
        for (int i = curd, j = b.curd; j > -1; --i, --j) {
            Field elem = coef[i] + t * b.coef[j];
            coef[i]  = elem;
        }
    }
    else {
        for (int i = curd, j = b.curd; i > -1; --i, --j) {
            Field elem = coef[i] * c;
            if (j > -1)
                elem += t * b.coef[j];
            coef[i] = elem;
        }
    }
    resetDegree();
}

template <class Field>
std::ostream& operator<< (std::ostream &out, DenseUnivariatePolynomial<Field> b) {
	bool isFirst = 0;
	for (int i = 0; i <= b.curd; ++i) {
		if (b.coef[i] != 0) {
			if (isFirst && b.coef[i] > 0)
				out << "+";
			if (i) {
				if (b.coef[i] != 1 && b.coef[i] != -1)
					out << b.coef[i] << "*";
				else if (b.coef[i] < 0)
					out << "-";
				out << b.name;
				if (i > 1)
					out << "^" << i;
			}
			else { out << b.coef[i]; }
			isFirst = 1;
		}
	}
	if (!isFirst) { out << "0"; }
	return out;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::operator^ (int e) {
    DenseUnivariatePolynomial<Field> res;
    res.name = name;
    res.one();
    unsigned long int q = e / 2, r = e % 2;
    DenseUnivariatePolynomial<Field> power2 = *this * *this;
    for (int i = 0; i < q; ++i)
        res *= power2;
    if (r) { res *= *this; }
    return res;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::operator<< (int k) {
    int s = curd + k + 1;
    DenseUnivariatePolynomial<Field> r;
    r.n = s;
    r.curd = s - 1;
    r.name = name;
    r.coef = new Field[s];
    for (int i = 0; i < k; ++i)
        r.coef[i] = 0;
    for (int i = k; i < s; ++i)
        r.coef[i] = coef[i-k];
    return r;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::operator>> (int k) {
    DenseUnivariatePolynomial<Field> r;
    r.name = name;
    int s = curd - k + 1;
    if (s > 0) {
        r.n = s;
        r.curd = s - 1;
        delete [] r.coef;
        r.coef = new Field[s];
        for (int i = 0; i < s; ++i)
            r.coef[i] = coef[i+k];
    }
}
    
template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::operator+ (DenseUnivariatePolynomial<Field> b) {
    if (!curd) { return (b + coef[0]); }
    if (!b.curd) { return (*this + b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to add between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int size = (curd > b.curd)? curd+1 : b.curd+1;
    DenseUnivariatePolynomial<Field> res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new Field[size];
    for (int i = 0; i < size; ++i) {
        Field elem = 0;
        if (i <= curd)
            elem += coef[i];
        if (i <= b.curd)
            elem += b.coef[i];
        res.coef[i] = elem;
    }
    res.resetDegree();
    return res;
}

template <class Field>
void DenseUnivariatePolynomial<Field>::add(DenseUnivariatePolynomial<Field> b) {
    for (int i = curd; i >= 0; --i) {
        if (i <= b.curd)
            coef[i] += b.coef[i];
    }
    resetDegree();
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::operator- (DenseUnivariatePolynomial<Field> b) {
    if (!curd) { return (coef[0] - b); }
    if (!b.curd) { return (*this - b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to subtract between Q[" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    
    int size = (curd > b.curd)? curd+1 : b.curd+1;
    DenseUnivariatePolynomial<Field> res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new Field[size];
    for (int i = 0; i < size; ++i) {
        Field elem = 0;
        if (i <= curd)
            elem = coef[i];
        if (i <= b.curd)
            elem -= b.coef[i];
        res.coef[i] = elem;
    }
    res.resetDegree();
    return res;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::operator- () {
    DenseUnivariatePolynomial<Field> res(curd+1);
    res.name = name;
    res.curd = curd;
    for (int i = 0; i <= curd; ++i)
        res.coef[i] = -coef[i];
    return res;
}

template <class Field>
void DenseUnivariatePolynomial<Field>::subtract(DenseUnivariatePolynomial<Field> b) {
    for (int i = curd; i >= 0; --i) {
        if (i <= b.curd)
            coef[i] -= b.coef[i];
    }
    resetDegree();
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::operator* (DenseUnivariatePolynomial<Field> b) {
    if (!curd) { return (b * coef[0]); }
    if (!b.curd) { return (*this * b.coef[0]); }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to multiply between Field[" << name << "] and Field[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int d = curd + 1;
    int m = b.curd + 1;
    int size = curd + m;
    DenseUnivariatePolynomial<Field> res;
    res.n = size;
    res.curd = size - 1;
    res.name = name;
    res.coef = new Field[size];
    
    /*
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
    univariateMultiplication(mul, acoef, d, bcoef, m);  //need to change , comes from Multiplication.h class
    lfixz den = aden * bden;
    for (int i = 0; i < size; ++i) {
        Field elem = mul[i];
        elem /= den;
        res.coef[i] = elem;
    }
    */
    for(int i=0; i<d; ++i){
    	for(int j=0; i<m; ++j){
    		res[i+j] += coef[i] * b.coef[j];
    	}
    }
    res.resetDegree();
    
    //delete [] acoef;
    //delete [] bcoef;
    //delete [] mul;
    return res;
}

template <class Field>
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator*= (Field e) {
    if (e != 0 && e != 1) {
        for (int i = 0; i <= curd; ++i)
            coef[i] *= e;
    }
    else if (e == 0) { zero(); }
    return *this;
}
/*
template <class Field>
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator*= (RationalNumber e) {
    Field c (e.get_mpq_t());
    *this *= c;
    return *this;
}

template <class Field>
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator*= (mpq_class e) {
    if (e != 0 && e != 1) {
        for (int i = 0; i <= curd; ++i)
            coef[i] *= e;
    }
    else if (e == 0) { zero(); }
    return *this;
}

template <class Field>
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator*= (sfixn e) {
    if (e != 0 && e != 1) {
        for (int i = 0; i <= curd; ++i)
            coef[i] *= e;
    }
    else if (e == 0) { zero(); }
    return *this;
}
*/
template <class Field>
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator/= (DenseUnivariatePolynomial<Field> b) {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
        exit(1);
    }
    if (!b.curd)
        return (*this /= b.coef[0]);
    if (!curd) {
        coef[0].zero();
        return *this;
    }
    
    if (name != b.name) {
        std::cout << "BPAS: error, trying to exact divide between Field[" << name << "] and Field[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    DenseUnivariatePolynomial<Field> rem(*this);
    zeros();
    while (rem.curd >= b.curd) {
        Field lc = rem.coef[rem.curd] / b.coef[b.curd];
        int diff = rem.curd - b.curd;
        rem.pomopo(1, -lc, b);
        coef[diff] = lc;
    }
    resetDegree();
    if (!rem.isZero()) {
        std::cout << "BPAS: error, not exact division from DUFP." << std::endl;
        exit(1);
    }
    return *this;
}

template <class Field>
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator/= (Field e) {
    if (e == 0) {
        std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
        exit(1);
    }
    else if (e != 1) {
    	for (int i = 0; i <= curd; ++i)
    	    coef[i] /= e;
    }
    return *this;
}
/*
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator/= (RationalNumber e) {
    Field c (e.get_mpq_t());
    *this /= c;
    return *this;
}

template <class Field>
DenseUnivariatePolynomial<Field>& DenseUnivariatePolynomial<Field>::operator/= (mpq_class e) {
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
*/
template <class Field>
DenseUnivariatePolynomial<Field> operator/ (Field c, DenseUnivariatePolynomial<Field> p) {
    if (p.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
        exit(1);
    }
    
    DenseUnivariatePolynomial<Field> q;
    q.name = p.name;
    q.curd = 0;
    q.n = 1;
    q.coef = new Field[1];
    
    if (p.isConstant())
        q.coef[0] = c / p.coef[0];
    else
        q.coef[0].zero();   
    return q;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::monicDivide(DenseUnivariatePolynomial<Field>& b) {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from DUFP." << std::endl;
        exit(1);
    }
    else if (b.coef[b.curd] != 1) {
        std::cout << "BPAS: error, leading coefficient is not one in monicDivide() from DUFP." << std::endl;
        exit(1);
    }
    if (!b.curd) {
        DenseUnivariatePolynomial<Field> r (*this);
        zero();
        return r;
    }
    if (!curd) {
        DenseUnivariatePolynomial<Field> r;
        r.name = name;
        return r;
    }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to monic divide between [" << name << "] and Q[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    int size = curd - b.curd + 1;
    DenseUnivariatePolynomial<Field> quo(size);
    quo.curd = size - 1;
    quo.name = name;
    while (curd >= b.curd) {
        Field lc = coef[curd];
        int diff = curd - b.curd;
        pomopo(1, -lc, b);
        quo.coef[diff] = lc;
    }
    quo.resetDegree();
    return quo;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::monicDivide(DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* rem) {
    *rem = *this;
    return rem->monicDivide(b);
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::lazyPseudoDivide (DenseUnivariatePolynomial<Field>& b, Field* c, Field* d) {
    if (d == NULL)
        d = new Field;
    int da = curd, db = b.curd;
    if (b.isZero() || !db) {
        std::cout << "BPAS: error, dividend is zero or constant." << std::endl;
        exit(1);
    }
    c->one(), d->one();
    if (!curd) {
        DenseUnivariatePolynomial<Field> r;
        r.name = name;
        return r;
    }
    if (name != b.name) {
        std::cout << "BPAS: error, trying to pseudo divide between Field[" << name << "] and Field[" << b.name << "]." << std::endl;
        exit(1);
    }
    
    if (da < db) {
        DenseUnivariatePolynomial<Field> r;
        r.name = name;
        return r;
    }
    
    int size = curd - b.curd + 1;
    DenseUnivariatePolynomial<Field> quo(size);
    quo.curd = size - 1;
    quo.name = name;
    int e = 0, diff = da - db;
    Field blc = b.coef[b.curd];
    while (curd >= b.curd) {
        Field lc = coef[curd];
        int k = curd - b.curd;
        *c *= blc;
        e++;
        pomopo(blc, -coef[curd], b);
        quo.coef[k] = lc;
    }
    quo.resetDegree();
    for (int i = e; i <= diff; ++i)
        *d *= blc;
    return quo;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::lazyPseudoDivide (DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* rem, Field* c, Field* d) {
    *rem = *this;
    return rem->lazyPseudoDivide(b, c, d);
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::pseudoDivide (DenseUnivariatePolynomial<Field>& b, Field* d) {
    Field c;
    if (d == NULL)
        d = new Field;
    DenseUnivariatePolynomial<Field> quo = lazyPseudoDivide(b, &c, d);
    quo *= *d;
    *this *= *d;
    *d *= c;
    return quo;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::pseudoDivide (DenseUnivariatePolynomial<Field>& b, DenseUnivariatePolynomial<Field>* rem, Field* d) {
    Field c;
    DenseUnivariatePolynomial<Field> quo = lazyPseudoDivide(b, rem, &c, d);
    quo *= *d;
    *rem *= *d;
    *d *= c;
    return quo;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::halfExtendedEuclidean (DenseUnivariatePolynomial<Field> b, DenseUnivariatePolynomial<Field>* g) {
    if (g == NULL)
        g = new DenseUnivariatePolynomial<Field>;
    *g = *this;
    
    DenseUnivariatePolynomial<Field> a1, b1;
    a1.name = name;
    b1.name = b.name;
    a1.coef[0].one();
    while (!b.isZero()) {
        DenseUnivariatePolynomial<Field> q, r;
        q.name = r.name = name;
        Field e = b.coef[b.curd];
        if (e != 1) {
            b /= e;
            *g /= e;
        }
        q = g->monicDivide(b, &r);
        if (e != 1) {
            *g = b * e;
            b = r * e;
        }
        else {
            *g = b;
            b = r;
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

template <class Field>
void DenseUnivariatePolynomial<Field>::diophantinEquationSolve(DenseUnivariatePolynomial<Field> a, DenseUnivariatePolynomial<Field> b, DenseUnivariatePolynomial<Field>* s, DenseUnivariatePolynomial<Field>* t) {
    DenseUnivariatePolynomial<Field> f(*this), g, q, r;
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
    
    Field e = b.coef[b.curd];
    if (e != 1) { b /= e; }
    if (s->curd >= b.curd) {
        *s /= e;
        s->monicDivide(b, &r);
        *s = r * e;
    }
    
    g = *this;
    g -= *s * a;
    if (e != 1) { g /= e; }
    *t = g.monicDivide(b);
}

template <class Field>
void DenseUnivariatePolynomial<Field>::differentiate(int k) {
    if (k <= 0) { return; }
    for (int i = k; i <= curd; ++i) {
        coef[i-k] = coef[i];
        for (int j = 0; j < k; ++j)
            coef[i-k] *= (i - j);
    }
    curd -= k;
    resetDegree();
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::integrate() {
    DenseUnivariatePolynomial<Field> b;
    b.name = name;
    b.n = curd+2;
    b.coef = new Field[b.n];
    b.coef[0].zero();
    for (int i = 0; i <= curd; ++i)
        b.coef[i+1] = coef[i] / (i + 1);
    b.resetDegree();
    return b;
}

template <class Field>
std::vector<DenseUnivariatePolynomial<Field> > DenseUnivariatePolynomial<Field>::squareFree() {
    std::vector<DenseUnivariatePolynomial<Field> > sf;
    if (!curd)
        sf.push_back(*this);
    else if (curd == 1) {
        DenseUnivariatePolynomial<Field> t;
        t.name = name;
        t.coef[0] = coef[curd];
        sf.push_back(t);
        t = *this / t.coef[0];
        sf.push_back(t);
    }
    else {
        DenseUnivariatePolynomial<Field> a (*this), b(*this);
        b.differentiate(1);
        DenseUnivariatePolynomial<Field> g = a.gcd(b);
        DenseUnivariatePolynomial<Field> x = a / g;
        DenseUnivariatePolynomial<Field> y = b / g;
        DenseUnivariatePolynomial<Field> z = -x;
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
        
        Field e;
        e.one();
        for (int i = 0; i < sf.size(); ++i) {
            e *= sf[i].coef[sf[i].curd];
            sf[i] /= sf[i].coef[sf[i].curd];
        }
        DenseUnivariatePolynomial<Field> t;
        t.name = name;
        t.coef[0] = e;
        sf.insert(sf.begin(), t);
    }
    return sf;
}

template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::euclideanGCD (DenseUnivariatePolynomial<Field> q) {
    DenseUnivariatePolynomial<Field> a, b;
    if (curd < q.curd) {
        a = q;
        b = *this;
    }
    else {
        a = *this;
        b = q;
    }
    
    while (!b.isZero()) {
        Field lc = b.coef[b.curd];
        b /= lc;
        DenseUnivariatePolynomial<Field> r;
        r.name = name;
        a.monicDivide(b, &r);
        a = b * lc;
        b = r;
    }
    return a / a.coef[a.curd];
}
/*
template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::modularGCD (DenseUnivariatePolynomial<Field> q) {
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

	DenseUnivariateIntegerPolynomial c = a.gcd(b);

	int d = c.degree();
	DenseUnivariatePolynomial<Field> r (d+1);
	r.curd = d;
	r.name = name;
	for (int i = 0; i < d; ++i)
		r.coef[i] = c.coefficient(i) / c.coefficient(d);
	r.coef[d].one();

	return r;
}
*/
template <class Field>
DenseUnivariatePolynomial<Field> DenseUnivariatePolynomial<Field>::gcd (DenseUnivariatePolynomial<Field> q, int type) {
    if (isZero()) { return q; }
    if (q.isZero()) { return *this; }
    if (!curd || !q.curd) {
        DenseUnivariatePolynomial<Field> h (1);
        h.coef[0].one();
        h.name = name;
        return h;
    }

    if (name != q.name) {
        std::cout << "BPAS: error, trying to compute GCD between Q[" << name << "] and Q[" << q.name << "]." << std::endl;
        exit(1);
    }

	DenseUnivariatePolynomial<Field> r;
	if (!type)
		r = euclideanGCD(q);
	else 
		r = modularGCD(q);

	return r;
}
//to avoid linking errors.
template void DenseUnivariatePolynomial<SmallPrimeField>::resetDegree();
template void DenseUnivariatePolynomial<SmallPrimeField>::differentiate(int);
template bool DenseUnivariatePolynomial<SmallPrimeField>::isTrailingCoefficientZero();
//template std::ostream& operator<< (std::ostream &out, DenseUnivariatePolynomial<SmallPrimeField> b);


