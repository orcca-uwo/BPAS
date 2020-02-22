#include "RationalNumberPolynomial/urpolynomial.h"

void DenseUnivariateRationalPolynomial::isolateScaledUnivariatePolynomial(Intervals* pIs, DenseUnivariateRationalPolynomial* p, int ts) {
	DenseUnivariateRationalPolynomial A(*p);
	A.reciprocal();
	A.taylorShift(ts);

	int num = A.numberOfSignChanges();
	if (num == 1) {
		pIs->pushInterval((double)0/1, (double)1/1);
		return;
	}
	else if (!num)
		return;

	Intervals ppIs;
	Intervals npIs;
	Intervals zpIs;

	A = *p;
	A.homothetic();
	DenseUnivariateRationalPolynomial B(A);
	B.taylorShift(ts);

	if (B.divideByVariableIfCan()) {
		zpIs.pushInterval((double)1/2, (double)1/2);
	}

	cilk_spawn isolateScaledUnivariatePolynomial(&npIs, &A, ts);
	isolateScaledUnivariatePolynomial(&ppIs, &B, ts);
	cilk_sync;

	npIs.transformLeft();
	ppIs.transformRight();
	pIs->concatenate(npIs);
	pIs->concatenate(zpIs);
	pIs->concatenate(ppIs);
}

void DenseUnivariateRationalPolynomial::isolatePositiveUnivariateRealRoots(Intervals* pIs, DenseUnivariateRationalPolynomial* p, int ts) {
	int k = p->rootBoundExp();
	DenseUnivariateRationalPolynomial A(*p);
	A.scaleTransform(k);

	isolateScaledUnivariatePolynomial(pIs, &A, ts);
	//genDescartes(pIs, &A, ts);
	pIs->scale(k);
}

void DenseUnivariateRationalPolynomial::isolateUnivariateRealRoots(Intervals* pIs, DenseUnivariateRationalPolynomial* p, int ts) {
	DenseUnivariateRationalPolynomial A(*p);

	Intervals zpIs;
	if (A.divideByVariableIfCan())
		zpIs.pushInterval((double)0/1, (double)0/1);

	Intervals ppIs;
	Intervals npIs;

	DenseUnivariateRationalPolynomial B(A);
	B.negativeVariable();

	cilk_spawn isolatePositiveUnivariateRealRoots(&ppIs, &A, ts);
	isolatePositiveUnivariateRealRoots(&npIs, &B, ts);
	cilk_sync;

	npIs.negate();
	pIs->concatenate(npIs);
	pIs->concatenate(zpIs);
	pIs->concatenate(ppIs);
}

void DenseUnivariateRationalPolynomial::refineUnivariateInterval(Interval* pI, lfixq npI, DenseUnivariateRationalPolynomial* p, lfixq width) {
        int ls = -1, rs = 1;
        lfixq left = pI->left, right = pI->right;

         if (left == right)
                return;

        while (sign(left) * sign(right) <= 0 || width < (right - left) || npI <= right) {
                lfixq middle = (left + right) >> 1;

                ls = cilk_spawn sign(p->evaluate(left));
                rs = cilk_spawn sign(p->evaluate(right));
                int ms = sign(p->evaluate(middle));
                cilk_sync;

                if (!ms) {
                        left = middle;
                        right = middle;
                        break;
                }
                else if ((ms > 0 && ls < rs) || (ms < 0 && rs < ls))
                	right = middle;
                else { left = middle; }
        }
        pI->left = left;
        pI->right = right;
}

void DenseUnivariateRationalPolynomial::refineUnivariateIntervals(Intervals* repIs, Intervals* pIs, DenseUnivariateRationalPolynomial* p, lfixq width) {
        int m = pIs->numberOfIntervals();

	*repIs = *pIs;
        cilk_for (int i = 0; i < m; i++) {
		lfixq border = (i < m-1)? pIs->interval(i+1, 0)->left : pIs->interval(i, 0)->right+1;
                refineUnivariateInterval(repIs->interval(i, 0), border, p, width);
        }
}

void DenseUnivariateRationalPolynomial::DenseUnivariateRationalPolynomial::univariatePositiveRealRootIsolation(Intervals* pIs, DenseUnivariateRationalPolynomial* p, lfixq width, int ts) {
        Intervals repIs;
        isolatePositiveUnivariateRealRoots(&repIs, p, ts);
	refineUnivariateIntervals(pIs, &repIs, p, width);
}

void DenseUnivariateRationalPolynomial::univariateRealRootIsolation(Intervals* pIs, DenseUnivariateRationalPolynomial* p, lfixq width, int ts) {
        Intervals repIs;
        isolateUnivariateRealRoots(&repIs, p, ts);
	refineUnivariateIntervals(pIs, &repIs, p, width);
	std::vector<Symbol> xs;
	xs.push_back(Symbol(p->variable()));
	pIs->setVariableNames(xs);
}
