#include "RationalNumberPolynomial/mrpolynomial.h"

/****************
* Multi-Divisor Division
*****************/

/** Normal Form Algorithm in lexicographical polynomial ordering
 * @param[in] superNames The vector of variable names in order
 * @param[in] ts The triangular-set
 * @param[out] r The remainder
 * @param[out] quoSet The quotient-set
 * \brief This algorithm runs by f and ts[s], and NULL quoSet[s] where s is the size of triangular-set,
 * Computes r, quoSet[0], ... and quoSet[s] such that f = quoSet[0]*ts[0] + ... + quoSet[s-1]*ts[s-1] + r.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::lexNormalForm(const std::vector<Symbol>& superNames,
		const std::vector<SparseMultivariateRationalPolynomial>& ts, std::vector<SparseMultivariateRationalPolynomial>* quoSet) const{

	// the type of MDD: (0:heapMDD, 1:triangularSetMDD, 2:primitiveFactorTriangularSetMDD)
	int TYPE = 1;
	int numVar = superNames.size();
	int tsSize = ts.size();
	// TODO: test symbol sets of all input polynomials belong to superNames

	// make a zero polynomial w.r.t. superNames
	SparseMultivariateRationalPolynomial zeroRef(numVar);
	zeroRef.setVariableNames(superNames);
	int zeroRefDegs[numVar] = {0};
	zeroRefDegs[0] = 1;
	zeroRef.setCoefficient(numVar, zeroRefDegs, RationalNumber(1));

	// update the dividend polynomial w.r.t. superNames:  --> divEx
	SparseMultivariateRationalPolynomial divEx(numVar);
	divEx = *this;
	divEx += zeroRef;
	divEx -= zeroRef;

	// update the triangular set w.r.t. superNames:       --> tsEx
	std::vector<SparseMultivariateRationalPolynomial> tsEx;
	tsEx.clear(); // it is not necessary!
	for (int i = 0; i < tsSize; ++i){
		SparseMultivariateRationalPolynomial tmp; // = ts[i];
		tmp = ts[i];
		tmp += zeroRef;
		tmp -= zeroRef;
		tsEx.push_back(tmp);
	}

	// make Node* polynomials for MDD's inputs:
	Node* g[tsSize];
	Node* q[tsSize];
	for (int i = 0; i < tsSize; ++i){
		g[i] = tsEx[tsSize-i-1].poly;
		q[i] = NULL;
	}
	Node* f = divEx.poly;
	Node* r = NULL;

	// compute the MDD:
	multiDivisorDivision(f, g, q, &r, tsSize, numVar, TYPE);

	// convert outputs into SMQP:
	SparseMultivariateRationalPolynomial outR {r, numVar, zeroRef.names};

	if (true){ // TODO:
		quoSet->clear();
		for (int i = 0; i < tsSize; ++i){
			SparseMultivariateRationalPolynomial qi {q[tsSize-i-1], numVar, zeroRef.names};
			quoSet->push_back(qi);
		}
	}

	return outR;
}

/*
* Multi-Divisor Division (MDD)
 * @param[in] dividend The dividend
 * @param[in] divisorSet The divisor-set
 * @param[in] t The type of MDD (type = 0 ? HeapMDD : (type = 1 ? triangularSetMDD : primitiveFactorTriangularSetMDD))
 * @param[out] rem The remainder
 * @param[out] quoSet The quotient-set
 * @param[out] 0-1 Return 1, unless the C version of algorithm doesn't work
 * \brief This algorithm computes the remainder and quoSet.

bool multiDivisorDivision(const SparseMultivariateRationalPolynomial& dividend,
		const std::vector<SparseMultivariateRationalPolynomial>& divisorSet,
		std::vector<SparseMultivariateRationalPolynomial>* quoSet,
		SparseMultivariateRationalPolynomial* rem, int t) {
	bool status = false;
	int numvar = dividend.nvar; // the number of variables
	int s = divisorSet.size(); // the size of divisor set
	int type = t; // 0: heapMDD, 1: tirangularSetMDD, 2: primitiveFactorMDD

	Node* g[s]; // input divisor set
	Node* q[s]; // quotient set (it should be convert to the SparseMultivariateRationalPolynomial at the end!)
	for (int i = 0; i < s; ++i) {
		g[i] = divisorSet[s - i - 1].poly;
		q[i] = NULL;
	}

	for (int i = 0; i < s; ++i) {
		convertExponentPoly(g[i], i + 1, s);
	}

	Node* f = dividend.poly;
	//Node** r = &(rem->poly);
	SparseMultivariateRationalPolynomial tmpR { NULL, numvar, dividend.names };
	Node* r = NULL;

	if (status != true) {
		multiDivisorDivision(f, g, q, &r, s, numvar, type);
		tmpR.poly = r;
	}

	quoSet->clear();
	for (int i = 0; i < s; ++i) {
		SparseMultivariateRationalPolynomial qi { q[s - i - 1], numvar,
				dividend.names };
		quoSet->push_back(qi);
	}

	*rem = tmpR;

	return true;
}
*/
