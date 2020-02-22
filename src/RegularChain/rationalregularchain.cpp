#include "RegularChain/rationalregularchain.h"


RationalRegularChain& RationalRegularChain::operator+= (DenseUnivariateRationalPolynomial up) {
    if (up.isConstant()) {
        std::cout << "BPAS: error, cannot adding a constant to RationalRegularChain" << std::endl;
        exit(1);
    }
    if (names[0] != up.variable()) {
        std::cout << "BPAS: error, adding a univariate polynomial of " << up.variable() << " to RationalRegularChain[";
        for (int i = var-1; i > -1; --i) {
            std::cout << names[i];
            if (i)
                std::cout << ", ";
        }
        std::cout << "]." << std::endl;
        exit(1);
    }
    
    p += up;
    return *this;
}

RationalRegularChain& RationalRegularChain::operator+= (SparseMultivariateRationalPolynomial mp) {
    if (mp.isConstant()) {
        std::cout << "BPAS: error, cannot adding a constant to RationalRegularChain" << std::endl;
        exit(1);
    }
    int v = mp.numberOfVariables();
    if (v > var || v < 2) {
        std::cout << "BPAS: error, RationalRegularChain(" << var << "), but adding a multivariate polynomial of " << v << " variables." << std::endl;
        exit(1);
    }
    std::vector<Symbol> xs = mp.variables();
    for (int i = 0, j = v-1; j > -1; ++i, --j) {
        if (names[i] != xs[j]) {
            std::cout << "BPAS: error, adding Q[";
            for (int k = 0; k < v; ++k) {
                std::cout << xs[k];
                if (k < v-1) { std::cout << ", "; }
            }
            std::cout << "] in RationalRegularChain[";
            for (int k = var-1; k > -1; --k) {
                std::cout << names[k];
                if (k) { std::cout << ", "; }
            }
            std::cout << "], where is obtained by Triangularize." << std::endl;
            exit(1);
        }
    }
    
    mps[v-2] += mp;
    
    return *this;
}

SparseMultivariateRationalPolynomial RationalRegularChain::select(const Symbol& x) const {
    int k = -1;
    for (int i = 0; i < var; ++i) {
        if (x == names[i]) {
            k = i;
            break;
        }
    }
    
    if (k < 0) {
        SparseMultivariateRationalPolynomial r;
        return r;
    }
    else if (!k) {
        SparseMultivariateRationalPolynomial r(p);
        return r;
    }
    else
        return mps[k-1];
    
}

void RationalRegularChain::lower(const Symbol& x, RationalRegularChain& ts) const {
    int k = -1;
    for (int i = 0; i < var; ++i) {
        if (x == names[i]) {
            k = i;
            break;
        }
    }
    
    if (k <= 0)
        ts = RationalRegularChain();
    else {
        RationalRegularChain r;
        r.var = k;
        r.names = new Symbol[r.var];
        std::copy(names, names+k, r.names);
        
        r += p;
        for (int i = 0; i < k-1; ++i)
            r += mps[i];
        ts = r;
    }
}

std::ostream& operator<< (std::ostream &out, RationalRegularChain& rc) {
    if (rc.var) {
        out << rc.p << "\n";
        for (int i = 0; i < rc.var-1; ++i) {
            out << rc.mps[i];
            if (i < rc.var-2)
                out << "\n";
        }
    }
    else { out << "rational regular chain doesn't exist."; }
    return out;
}
