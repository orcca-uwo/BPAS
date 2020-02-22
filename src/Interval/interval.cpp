#include "Interval/interval.h"

void intervalMultiplication(Interval* res, Interval* a, Interval* b) {
	if (a->left < 0 && b->left < 0) {
		res->left = a->right * b->right;
		res->right = a->left * b->left;
	}
	else if (a->left > 0 && b->left > 0) {
		res->left = a->left * b->left;
		res->right = a->right * b->right;
	}
	else if (a->left > 0 && b->left < 0) {
		res->left = a->right * b->left;
		res->right = a->left * b->right;
	}
	else if (a->left < 0 && b->right > 0) {
		res->left = a->left * b->right;
		res->right = a->right * b->left;
	}
	else {
		res->left = 0;
		res->right = 0;
	}
}

bool Intervals::isExactEqual(Intervals& pIs) {
	if (numberOfIntervals() != pIs.numberOfIntervals())
		return 0;

	for (int i = 0; i < roots.size(); ++i) {
		int k = i / var, j = i % var;
		if (roots[i].left != pIs.interval(k, j)->left || roots[i].right != pIs.interval(k, j)->right)
			return 0;
	}

	return 1;
}

void Intervals::copyFrom(Intervals& pIs, int k) {
	var = pIs.numberOfVariables();
	roots.clear();
	for (int i = 0; i < var; i++) {
		Interval* pI = pIs.interval(k, i);
		roots.push_back(*pI);
	}
}

void Intervals::scale(int k) {
	if (k > 0) {
		int n = roots.size();

		for (int i = 0; i < n; i++) {
			roots[i].left <<= k;
			roots[i].right <<= k;
		}
	}
}

void Intervals::transformLeft() {
	int n = roots.size();

	for (int i = 0; i < n; i++) {
		roots[i].left >>= 1;
		roots[i].right >>= 1;
	}
}

void Intervals::transformRight() {
	int n = roots.size();

	for (int i = 0; i < n; i++) {
		roots[i].left++;
		roots[i].left >>= 1;
		roots[i].right++;
		roots[i].right >>= 1;
	}
}

void Intervals::negate() {
	int n = roots.size();

	lfixq elem;
	for (int i = 0; i < n/2; i++) {
		elem = roots[i].left;
		roots[i].left = -roots[i].right;
		roots[i].right = -elem;

		elem = roots[n-i-1].left;
		roots[n-i-1].left = -roots[n-i-1].right;
		roots[n-i-1].right = -elem;

		swapInterval(&roots[i], &roots[n-i-1]);
	}
	if (n%2) {
		elem = roots[n/2].left;
		roots[n/2].left = -roots[n/2].right;
		roots[n/2].right = -elem;
	}
}

void Intervals::concatenate(Intervals& pIs) {
	int n = pIs.numberOfIntervals();

	if (n) {
		for (int i = 0; i < n; i++) {
			int v = pIs.numberOfVariables();
			for (int j = 0; j < v; j++) {
				Interval* pI = pIs.interval(i, j);
				roots.push_back(*pI);
			}
		}
	}
}

std::ostream& operator<< (std::ostream &out, Intervals& a) {
	for (int i = 1; i <= a.var; ++i) {
		out << a.names[i];
		if (i < a.var) { out << " * "; }
		else { out << " = "; }
	}
	int n = a.roots.size() / a.var;
	if (!n) { out << "[]"; }
	for (int i = 0; i < n; i++) {
		if (i) { out << "\t"; }
		for (int j = 0; j < a.var; ++j) {
			out << "[" << a.roots[i*a.var+j].left << ", " << a.roots[i*a.var+j].right << "]";
			if (j < a.var-1) { out << " * "; }
		}
		if (i < n-1) { out << ",\n"; }
	}
	return out;
}
