#include "RegularChain/rationalregularchain.h"

lfixq bpas_root_width;

int SparseMultivariateRationalPolynomial::isIntervalsMatchable(Intervals* a, Intervals* b, DenseUnivariateRationalPolynomial* lo, DenseUnivariateRationalPolynomial* up, lfixq width) {
        int an = a->numberOfIntervals();
        int bn = b->numberOfIntervals();

        if (!an && !bn)
                return 0;
        if (an != bn)
                return -1;

	// Check the sequence, should be one of the following:
	// 1, 2, 2, 1
	// 1, 2, 2, 1, 1, 2
	// 2, 1, 1, 2
	// 2, 1, 1, 2, 2, 1

	int order = 0;  // 0: 1 2; 1: 2 1
        for (int i = 0; i < an; i++) {
                Interval* aI = a->interval(0, i);
                Interval* bI = b->interval(0, i);

		// Check the width between two intervals
		if (abs(bI->left - aI->left) > width)
			return -1;

		// Checking the signs
		int as = cilk_spawn sign(lo->evaluate(aI->right));
		int bs = sign(up->evaluate(bI->right));
        	cilk_sync;
		if (as * bs < 0)
                	return -1;

		// Checking whether overlaped
		if (i < an-1) {
                        Interval* anext = a->interval(0, i+1);
                        Interval* bnext = b->interval(0, i+1);

                        lfixq max = (aI->right > bI->right)? aI->right : bI->right;
                        lfixq min = (anext->left < bnext->left)? anext->left : bnext->left;
                        if (max >= min)
                                return -2;
                }

		// Checking the order
		if (aI->right < bI->left) {
			if (i && !order) { return -1; }
			order = 0;
		}
		else if (bI->right < aI->left) {
			if (i && order) { return -1; }
			order = 1;
		}
		else { return -2; }
        }

        return 1;
}

int SparseMultivariateRationalPolynomial::refineSleeveUnivariateIntervals(Intervals* pI, Intervals* a, Intervals* b, DenseUnivariateRationalPolynomial* lo, DenseUnivariateRationalPolynomial* up, lfixq width) {
	int isOrder = isIntervalsMatchable(a, b, lo, up, bpas_root_width);
	// Refine the sleeve polynomials
        while (isOrder == -2) {
                width >>= 3;
                Intervals a1(*a), b1(*b);
                *a = lo->refineRoots(a1, width);
                *b = up->refineRoots(b1, width);
                isOrder = isIntervalsMatchable(a, b, lo, up, bpas_root_width);
        }

        if (!isOrder)
                return 0;
        else if (isOrder < 0)
                return -1;

	// Assign the roots based on the two lists
        int n = a->numberOfIntervals();
        for (int i = 0; i < n; i++) {
                Interval* aI = a->interval(0, i);
                Interval* bI = b->interval(0, i);

                Interval elem;
                elem.left = (aI->left > bI->left)? bI->left : aI->left;
                elem.right = (aI->right > bI->right)? aI->right : bI->right;
                pI->pushInterval(elem);
        }

        return 1;
}

void RationalRegularChain::refineMultivariateInterval(Intervals* apI, lfixq width, int ts) {
	p.refineRoot(apI->interval(0, 0), width);

	int step = apI->numberOfVariables() - 1;
        for (int i = 0; i < step; i++) {
                int v = mps[i].numberOfVariables();
		// int n = mps[i].numberOfTerms();
		SparseMultivariateRationalPolynomial A(mps[i]);
		Interval* c = apI->interval(0, i+1);
		Intervals pIs;
		int s = (c->left >= 0)? 0 : 1;
		A.positiveRealRootIsolation(&pIs, *apI, width, ts, s, 0);
		if (s) { pIs.negate(); }

		// Search for the root of variable i+1
		if (c->left != 0) {
			int left = 0, right = pIs.numberOfIntervals()-1;
			int middle = right;
			while (middle > -1) {
				Interval* d = pIs.interval(0, middle);
				if (d->left > c->right) 
					right = middle;
				else if (d->right < c->left) 
					left = middle;
				else {
					c->left = d->left;
					c->right = d->right;
					break;
				}
				middle = (left + right) / 2;
			}
		}
        }
}

int SparseMultivariateRationalPolynomial::positiveRealRootIsolation(Intervals* pIs, Intervals& apIs, lfixq width, int ts, bool s, bool check) {
	int d = leadingVariableDegree().get_si();
	DenseUnivariateRationalPolynomial up(d+1), lo(d+1);
	sleeveBoundURPolynomials(&up, &lo, apIs, 0, s);

	if (check && up.coefficient(d) * lo.coefficient(d) <= 0)
		return 1;

	if (up == lo)
		*pIs = up.positiveRealRootIsolate(width, ts);
	else {
#if __GNUC__ == 4
		Intervals a = cilk_spawn lo.positiveRealRootIsolate(width, ts);
		Intervals b = up.positiveRealRootIsolate(width, ts);
		cilk_sync;
#else
		Intervals a, b;
		cilk_spawn lo.positiveRealRootIsolate(width, a, ts);
		up.positiveRealRootIsolate(width, b, ts);
		cilk_sync;
#endif
		int res = refineSleeveUnivariateIntervals(pIs, &a, &b, &lo, &up, width);
		if (res < 0) { return 1; }
	}

	return 0;
}


void RationalRegularChain::isolatePositiveMultivariateRealRoots(Intervals* pIs, Intervals* apI, int s, int v, lfixq width, int ts) {
	// int n = mps[v-1].numberOfTerms();
	SparseMultivariateRationalPolynomial A(mps[v-1]);

	int posres = A.positiveRealRootIsolation(pIs, *apI, width, ts, s, 1);
	// Need to refine previous root
	while (posres == 1) {
		width >>= 3;
		refineMultivariateInterval(apI, width, ts);
		pIs->clear();
		posres = A.positiveRealRootIsolation(pIs, *apI, width, ts, s, 1);
	}
}

void RationalRegularChain::isolateMultivariateRealRoots(Intervals* mpIs, Intervals* apIs, int k, lfixq width, int ts) {
	int v = apIs->numberOfVariables();
	Intervals ppIs, npIs;

	Intervals apI, bpI;
	apI.copyFrom(*apIs, k);
	bpI.copyFrom(*apIs, k);

	// Isolate positive and negative roots
	cilk_spawn isolatePositiveMultivariateRealRoots(&ppIs, &apI, 0, v, width, ts);
	isolatePositiveMultivariateRealRoots(&npIs, &bpI, 1, v, width, ts);
	cilk_sync;

	// Merge positive and negative roots
	int negres = npIs.numberOfIntervals();
	npIs.negate();
	for (int i = 0; i < negres; i++) {
		for (int j = 0; j < v; j++)
			mpIs->pushInterval(*(bpI.interval(0, j)));
		mpIs->pushInterval(*(npIs.interval(i, 0)));
	}

	if (mps[v-1].isConstantTermZero()) {
		for (int j = 0; j < v; j++)
			mpIs->pushInterval(*(apIs->interval(k, j)));
		mpIs->pushInterval(0/1, 0/1);
	}

	int posres = ppIs.numberOfIntervals();
	for (int i = 0; i < posres; i++) {
		for (int j = 0; j < v; j++)
			mpIs->pushInterval(*(apI.interval(0, j)));
		mpIs->pushInterval(*(ppIs.interval(i, 0)));
	}
}

int RationalRegularChain::rcRealRootIsolation(Intervals* mpIs, lfixq width, int ts) {
        // Univariate real root isolation
	*mpIs = p.realRootIsolate(width, ts);

        int m = mpIs->numberOfIntervals();
        if (m) {
                for (int i = 0; i < var-1; i++) {
                        int v = mps[i].numberOfVariables();
                        if (v - mpIs->numberOfVariables() != 1) {
                                std::cout << "BPAS: Warning, the input is not a regular chain!" << std::endl;
                                return -1;
                        }
                        // Convert to SLP
                        mps[i].straightLineProgram();

                        m = mpIs->numberOfIntervals();
                        Intervals* curpIs = new Intervals[m];
                        // Parallel isolate each root
                        cilk_for (int j = 0; j < m; j++)
                                isolateMultivariateRealRoots(&curpIs[j], mpIs, j, width, ts);

                        // Merge all roots
                        Intervals pIs(v);
                        for (int j = 0; j < m; j++)
                                pIs.concatenate(curpIs[j]);
                        if (!pIs.numberOfIntervals())
                                return i+2;
                        *mpIs = pIs;
                        delete [] curpIs;
                }
                return 0;
        }
        else { return 1; }
}

Intervals RealRootIsolation(std::vector<RationalRegularChain> chains, lfixq width, int ts) {
	int n = chains.size();
	if (n) {
		bpas_root_width = width;
	        Intervals* pIs = new Intervals[n];
	        // Parallel isolate each regular chain
	        cilk_for (int i = 0; i < n; ++i)
	                pIs[i] = chains[i].realRootIsolate(width, ts);

	        // Merge all roots
	        int vars = chains[0].numberOfVariables();
		Intervals mpIs(vars);
	        for (int i = 0; i < n; ++i) {
	                if (vars == pIs[i].numberOfVariables())
	                        mpIs.concatenate(pIs[i]);
	        }
	        delete [] pIs;
		mpIs.setVariableNames(chains[0].variables());
		return mpIs;
	}
	return Intervals();
}
