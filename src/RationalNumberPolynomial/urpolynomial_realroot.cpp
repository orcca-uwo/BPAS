#include "RationalNumberPolynomial/urpolynomial.h"

int DenseUnivariateRationalPolynomial::numberOfSignChanges() {
	int num = 0;
	int prev = 0, cur = 0;
	for (int i = 0; i <= curd; i++) {
		if (coef[i] > 0)
			cur = 1;
		else if (coef[i] < 0)
			cur = -1;
		if (!prev)
			prev = cur;
		if (cur != prev) {
			++num;
			prev = cur;
		}
	}
	return num;
}

lfixz DenseUnivariateRationalPolynomial::positiveRootBound() {
	lfixq* q = new lfixq[curd+1];
	if (coef[curd] < 0) {
		for (int i = 0; i <= curd; i++)
			q[i] = -coef[i];
	}
	else {
		for (int i = 0; i <= curd; i++)
			q[i] = coef[i];
	}

	lfixz maxelem = 2;
	for (int i = 0; i < curd; i++) {
		if (q[i] < 0) {
			lfixz elem = lfixz(-q[i] / q[curd]) + 1;
			int isExact = mpz_root(elem.get_mpz_t(), elem.get_mpz_t(), curd-i);
			if (!isExact) { elem++; }
			elem <<= 1;

			if (maxelem < elem)
				maxelem = elem;
		}
	}

	delete [] q;

	return maxelem;
}

lfixz DenseUnivariateRationalPolynomial::cauchyRootBound() {
	lfixq a0, ad;
	lfixq cur_max, elem;

	if (coef[0] < 0)
		a0 = -coef[0];
	else
		a0 = coef[0];

	if (coef[curd] < 0)
		ad = -coef[curd];
	else
		ad = coef[curd];

	cur_max = a0;
	for (int i = 1; i < curd; i++) {
		if (coef[i] < 0)
			elem = -coef[i];
		else
			elem = coef[i];
		if (cur_max < elem)
			cur_max = elem;
	}
	elem = ad;
	lfixz tmp = lfixz(cur_max / elem) + 2;

	return tmp;
}

lfixz DenseUnivariateRationalPolynomial::complexRootBound() {
	if (curd) {
		DenseUnivariateRationalPolynomial p(curd);
		p.curd = curd - 1;
		for (int i = 0; i < curd; ++i) {
			if (coef[i] < 0)
				p.coef[i] = -coef[i];
			else
				p.coef[i] = coef[i];
		}
		lfixq c = abs(coef[curd]);

		int k = 0;
		lfixz r = 1;
		lfixq val = c * r;
		while (val <= p.evaluate(r)) {
			r <<= 1;
			val = c;
			for (int i = 0; i < curd; ++i)
				val *= r;
			k++;
		}

		k = std::min(k, 3);
		lfixz s = r;
		for (int i = 0; i < k; ++i) {
			s >>= 1;
			lfixz rr = r - s;
			val = c;
			for (int i = 0; i < curd; ++i)
				val *= rr;
				if (p.evaluate(rr) < val)
				r = rr;
		}

		return r;
	}
	return 0;
}

int DenseUnivariateRationalPolynomial::rootBoundExp() {
	lfixz k1 = positiveRootBound();
	lfixz k2 = cauchyRootBound();
	lfixz k3 = complexRootBound();

	lfixz s = (k1 > k2)? k2 : k1;
	if (k3 < s) { s = k3; }

	int k = 0;
	while (s > 0) {
		s >>= 1;
		k++;
	}
	return k;
}

lfixz DenseUnivariateRationalPolynomial::rootBound() {
	lfixz k1 = positiveRootBound();
	lfixz k2 = cauchyRootBound();
	lfixz k3 = complexRootBound();
	lfixz k4 = (k1 > k2)? k2 : k1;
	return (k3 > k4)? k4 : k3;
}
