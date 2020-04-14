#include "ModularPolynomial/sdmpolynomial.h"

mpz_class SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::characteristic(101);

void SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::zeros() {
	for (int i = 0; i < n; ++i)
		coefs[i] = 0;
}

bool SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::isSameRing(const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const {
	if (var != b.var || names[0] != b.names[0])
		return 0;
	if (names[0] == "9") {
		for (int i = 1; i <= var; ++i) {
			if (names[i] != b.names[i])
				return 0;
		}
	}
	return 1;
}

bool SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::shrink() {
	int* ds = new int[var];
	for (int i = 0; i < var; ++i) { ds[i] = 0; }
	for (int i = 0; i < n; ++i) {
		if (coefs[i] != 0) {
		 	for (int j = 0; j < var; ++j) {
				int e = (i / acc[j]) % (degs[j] + 1);
				if (e > ds[j])
					ds[j] = e;
			}
		}
	}
	bool isShrink = 0;
	for (int i = 0; i < var; ++i) {
		if (ds[i] < degs[i]) {
			isShrink = 1;
			break;
		}
	}
	if (isShrink) {
		int* acs = new int[var];
		acs[0] = 1;
		for (int i = 0; i < var-1; ++i)
			acs[i+1] = acs[i] * ds[i];
		int s = acs[var-1] * ds[var-1];
		sfixn* cfs = new sfixn[s];

		for (int i = 0; i < s; ++i) {
			int k = 0;
			for (int j = 0; j < var; ++j) {
				int e = (i / acs[j]) % (ds[j] + 1);
				k += acc[j] * e;
			}
			cfs[i] = coefs[k];
		}

		delete [] coefs;
		delete [] degs;
		delete [] acc;
		coefs = cfs;
		degs = ds;
		acc = acs;
		return 1;
	}
	else {
		delete [] ds;
		return 0;
	}
}

bool SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator== (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const{
	if (var != b.var || p != b.p)
		return 0;

	int prev = 0;
	for (int i = 0; i < n; ++i) {
		int k = 0;
		for (int j = 0; j < var; ++j) {
			int e = (i / acc[j]) % (degs[j] + 1);
			if (e <= b.degs[j])
				k += e * b.acc[j];
			else if (coefs[i] != 0)
				return 0;
		}
		for (int j = prev+1; j < k; ++j) {
			if (b.coefs[j] != 0)
				return 0;
		}
		if (coefs[i] != b.coefs[k])
			return 0;
		prev = k;
	}
	for (int i = n; i < b.n; ++i) {
		if (b.coefs[i] != 0)
			return 0;
	}
	return 1;
}

SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator+= (const sfixn& e) {
	coefs[0] = AddMod(coefs[0], e, p);
	return *this;
}

SmallPrimeFieldDistributedDenseMultivariateModularPolynomial SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator+ (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const {
	if (p != b.p) {
		std::cout << "BPAS: error, trying to add between Z/" << p << "Z and Z/" << b.p << "Z from SFDDMMP." << std::endl;
		exit(1);
	}

	if (isConstant()) { return b + coefs[0]; }
	if (b.isConstant()) { return *this + b.coefs[0]; }

	bool isSame = isSameRing(b);
	if (!isSame) {
		std::cout << "BPAS: error, trying to add between Z/" << p << "Z[";
		for (int i = 1; i <= var; ++i) {
			std::cout << names[i];
			if (i < var)
				std::cout << ", ";
		}
		std::cout << "] and Z/" << b.p << "Z[";
		for (int i = 1; i <= b.var; ++i) {
			std::cout << b.names[i];
			if (i < b.var)
				std::cout << ", ";
		}
		std::cout << "] from SFDDMMP." << std::endl;
		exit(1);
	}

	int* ds = new int[var];
	for (int i = 0; i < var; ++i)
		ds[i] = (degs[i] >= b.degs[i])? degs[i] : b.degs[i];

	SmallPrimeFieldDistributedDenseMultivariateModularPolynomial res (var, ds, p);
	std::copy(names, names+var+1, res.names);
	//#pragma cilk_grainsize = 1024;
	for (int i = 0; i < res.n; ++i) {
		sfixn elem = 0;
		int offseta = 0, offsetb = 0;
		for (int j = 0; j < var; ++j) {
			int k = (i / res.acc[j]) % (res.degs[j] + 1);
			if (offseta >= 0 && k <= degs[j])
				offseta += k * acc[j];
			else
				offseta = -1;
			if (offsetb >= 0 && k <= b.degs[j])
				offsetb += k * b.acc[j];
			else
				offsetb = -1;
		}
		if (offseta >= 0 && offsetb >= 0)
			elem = AddMod(coefs[offseta], b.coefs[offsetb], p);
		else if (offseta >= 0)
			elem = coefs[offseta];
		else if (offsetb >= 0)
			elem = b.coefs[offsetb];
		res.coefs[i] = elem;
	}
	delete [] ds;
	return res;
}

SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator-= (const sfixn& e) {
	coefs[0] = SubMod(coefs[0], e, p);
	return *this;
}

SmallPrimeFieldDistributedDenseMultivariateModularPolynomial SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator- (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const {
	if (p != b.p) {
		std::cout << "BPAS: error, trying to subtract between Z/" << p << "Z and Z/" << b.p << "Z from SFDDMMP." << std::endl;
		exit(1);
	}
	if (isConstant()) { return -b + coefs[0]; }
	if (b.isConstant()) { return *this - b.coefs[0]; }

	bool isSame = isSameRing(b);
	if (!isSame) {
		std::cout << "BPAS: error, trying to subtract between Z/" << p << "Z[";
		for (int i = 1; i <= var; ++i) {
			std::cout << names[i];
			if (i < var)
				std::cout << ", ";
		}
		std::cout << "] and Z/" << b.p << "Z[";
		for (int i = 1; i <= b.var; ++i) {
			std::cout << b.names[i];
			if (i < b.var)
				std::cout << ", ";
		}
		std::cout << "] from SFDDMMP." << std::endl;
		exit(1);
	}

	int* ds = new int[var];
	for (int i = 0; i < var; ++i)
		ds[i] = (degs[i] >= b.degs[i])? degs[i] : b.degs[i];
	SmallPrimeFieldDistributedDenseMultivariateModularPolynomial res (var, ds, p);
	std::copy(names, names+var+1, res.names);
	for (int i = 0; i < res.n; ++i) {
		sfixn elem = 0;
		int offseta = 0, offsetb = 0;
		for (int j = 0; j < var; ++j) {
			int k = (i / res.acc[j]) % (res.degs[j] + 1);
			if (offseta >= 0 && k <= degs[j])
				offseta += k * acc[j];
			else { offseta = -1; }
			if (offsetb >= 0 && k <= b.degs[j])
				offsetb += k * b.acc[j];
			else { offsetb = -1; }
		}
		if (offseta >= 0 && offsetb >= 0)
			elem = SubMod(coefs[offseta], b.coefs[offsetb], p);
		else if (offseta >= 0)
			elem = coefs[offseta];
		else if (offsetb >= 0)
			elem = NegMod(b.coefs[offsetb], p);
		res.coefs[i] = elem;
	}
	delete [] ds;
	return res;
}

SmallPrimeFieldDistributedDenseMultivariateModularPolynomial SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator- () const {
	SmallPrimeFieldDistributedDenseMultivariateModularPolynomial res (var, degs, p);
	std::copy(names, names+var+1, res.names);
	for (int i = 0; i < res.n; ++i)
		res.coefs[i] = NegMod(coefs[i], p);
	return res;
}

void SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::plain_multiplication(SmallPrimeFieldDistributedDenseMultivariateModularPolynomial* c, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& a, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const {
	int* ds = new int[a.var];
	for (int i = 0; i < a.n; ++i) {
		for (int v = 0; v < a.var; ++v)
			ds[v] = (i / a.acc[v]) % (a.degs[v] + 1);
		for (int j = 0; j < b.n; ++j) {
			int k = 0;
			for (int v = 0; v < b.var; ++v) {
				int e = (j / b.acc[v]) % (b.degs[v] + 1);
				if (v < a.var)
					k += (ds[v] + e) * c->acc[v];
				else
					k += e * c->acc[v];
			}
			for (int v = b.var; v < a.var; ++v)
				k += ds[v] * c->acc[v];
			sfixn t = MulMod(a.coefs[i], b.coefs[j], p);
			c->coefs[k] = AddMod(c->coefs[k], t, p);
		}
	}
	delete [] ds;
}

bool SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::fftbased_multiplication(SmallPrimeFieldDistributedDenseMultivariateModularPolynomial* c, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& a, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const {
	bool isComputed = 0;
	int l = logceiling(c->n), s = 1 << l;
	if (s < c->n) { s <<= 1; l++; }
	sfixn* f = new sfixn[s];
	sfixn* g = new sfixn[s];
	f[0] = a.coefs[0], g[0] = b.coefs[0];
	for (int i = 1; i < s; ++i)
		f[i] = g[i] = 0;
	for (int i = 1; i < a.n; ++i) {
		int k = 0;
		for (int j = 0; j < a.var; ++j) {
			int e = (i / a.acc[j]) % (a.degs[j] + 1);
			k += e * c->acc[j];
		}
		f[k] = a.coefs[i];
	}
	for (int i = 1; i < b.n; ++i) {
		int k = 0;
		for (int j = 0; j < b.var; ++j) {
			int e = (i / b.acc[j]) % (b.degs[j] + 1);
			k += e * c->acc[j];
		}
		g[k] = b.coefs[i];
	}

	isComputed = PBPAS::ks_mul(l, s, f, g, p, SDMPBASESIZE);

	for (int i = 0; i < c->n; ++i)
		c->coefs[i] = f[i];

	delete [] f;
	delete [] g;
	return isComputed;
}

SmallPrimeFieldDistributedDenseMultivariateModularPolynomial SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator* (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const {
	if (p != b.p) {
		std::cout << "BPAS: error, trying to multiply between Z/" << p << "Z and Z/" << b.p << "Z from SFDDMMP." << std::endl;
		exit(1);
	}
	if (isConstant()) { return b * coefs[0]; }
	if (b.isConstant()) { return *this * b.coefs[0]; }

	bool isSame = isSameRing(b);
	if (!isSame) {
		std::cout << "BPAS: error, trying to multiply between Z/" << p << "Z[";
		for (int i = 1; i <= var; ++i) {
			std::cout << names[i];
			if (i < var)
				std::cout << ", ";
		}
		std::cout << "] and Z/" << b.p << "Z[";
		for (int i = 1; i <= b.var; ++i) {
			std::cout << b.names[i];
			if (i < b.var)
				std::cout << ", ";
		}
		std::cout << "] from SFDDMMP." << std::endl;
		exit(1);
	}

	int* ds = new int[var];
	for (int i = 0; i < var; ++i)
		ds[i] = degs[i] + b.degs[i];
	SmallPrimeFieldDistributedDenseMultivariateModularPolynomial res (var, ds, p);
	std::copy(names, names+var+1, res.names);
	bool isComputed = 0;
	if (res.n > SDMPBASESIZE) { // SDMPBASESIZE=1024
		// two convolution multiplication only with 2 variables
		isComputed = PBPAS::DMPMul(p, var, res.coefs, res.degs, coefs, degs, b.coefs, b.degs, 0);
	}
	if (!isComputed && res.n > 64) {
		// FFT-based multiplication using KS
		isComputed = fftbased_multiplication(&res, *this, b);
	}
	if (!isComputed) {
		plain_multiplication(&res, *this, b);
	}

	delete [] ds;
	return res;
}

SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator*= (const sfixn& s) {
	sfixn e = s;
	if (e != 0 && e != 1) {
		if (e < 0) { e = e % p + p; }
		for (int i = 0; i < n; ++i)
			coefs[i] = MulMod(coefs[i], e, p);
	}
	else if (e == 0) { zero(); }
	return *this;
}

void SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::print(std::ostream &out) const {
	bool isFirst = 0;
	for (int i = 0; i < this->n; ++i) {
		if (this->coefs[i] != 0) {
			if (isFirst) {
				if (this->coefs[i] >= 0)
					out << "+";
				else if (this->coefs[i] == -1)
					out << "-";
				if (this->coefs[i] != 1 && this->coefs[i] != -1)
					out << this->coefs[i];
				bool isIt = 1;
				for (int j = 0; j < this->var; ++j) {
					int exp = (i / this->acc[j]) % (this->degs[j] + 1);
					if (exp) {
						if ((this->coefs[i] != 1 && this->coefs[i] != -1 && isIt) || !isIt)
							out << "*";
						out << this->names[j+1];
						if (exp > 1)
							out << "^" << exp;
						isIt = 0;
					}
				}
			}
			else { out << this->coefs[i]; }
			isFirst = 1;
		}
	}
	if (!isFirst) { out << "0"; }
}

ExpressionTree SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::convertToExpressionTree() const {
	std::cerr << "BPAS ERROR: SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::convertToExpressionTree() NOT YET IMPLEMENTED" << std::endl;
	//TODO
	exit(1);
	return ExpressionTree();
}
