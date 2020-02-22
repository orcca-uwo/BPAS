#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include "../../include/TriangularSet/triangularset.hpp"

/** 
 * Normal Form (or Multi-Divisor Division (MDD)) 
 * Given the dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}},
 * such that f = q_0*g_0 + ... + q_{s-1}*g_{s-1} + r.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::lexNormalForm(const std::vector<Symbol>& superNames, const std::vector<SparseMultivariateRationalPolynomial>& ts, std::vector<SparseMultivariateRationalPolynomial>* quoSet) const
{    
    // the type of MDD: (0:heapMDD, 1:triangularSetMDD, 2:primitiveFactorTriangularSetMDD)
    int TYPE = 0;
    int numVar = superNames.size();
    int tsSize = ts.size();
    
    // make a zero polynomial w.r.t. superNames
    SparseMultivariateRationalPolynomial zeroRef(numVar);
    zeroRef.setRingVariables(superNames);
    
    int zeroRefDegs[numVar];
    for (int i = 0; i < numVar; ++i) {
		zeroRefDegs[i] = 0;
    }
    zeroRefDegs[0] = 1;
    zeroRef.setCoefficient(numVar, zeroRefDegs, RationalNumber(1));
    
    // update the dividend polynomial w.r.t. superNames:  --> divEx
    SparseMultivariateRationalPolynomial divEx(numVar);
    divEx = *this;
    divEx += zeroRef;
    divEx -= zeroRef;
    
    // update the triangular set w.r.t. superNames:       --> tsEx
    std::vector<SparseMultivariateRationalPolynomial> tsEx;
    tsEx.clear();
    for (int i = 0; i < tsSize; ++i){
		SparseMultivariateRationalPolynomial tmp;
		tmp = ts[i];
		tmp += zeroRef;
		tmp -= zeroRef;
		tsEx.push_back(tmp);
    }
    
    AltArr_t* g[tsSize];
    AltArr_t* q[tsSize];
    for (int i = 0; i < tsSize; ++i){
		g[i] = tsEx[tsSize-i-1].poly;
		q[i] = NULL;
    }
    AltArr_t* f = divEx.poly;
    AltArr_t* r = NULL;
    
    // compute the MDD:
    multiDivisorDivision_AA(f, g, q, &r, tsSize, numVar, TYPE);
        
    // convert outputs into SMQP:
    SparseMultivariateRationalPolynomial outR {r, numVar, zeroRef.names};
    
    quoSet->clear();
    for (int i = 0; i < tsSize; ++i){
		SparseMultivariateRationalPolynomial qi {q[tsSize-i-1], numVar, zeroRef.names};
		quoSet->push_back(qi);
    }
    
    return outR;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::lexNormalizeDim0 (const std::vector<Symbol>& superNames, const std::vector<SparseMultivariateRationalPolynomial>& ts, SparseMultivariateRationalPolynomial* A) const
{    
    int numVar = superNames.size();
    int tsSize = ts.size();
    
    // make a zero polynomial w.r.t. superNames
    SparseMultivariateRationalPolynomial zeroRef(numVar);
    zeroRef.setRingVariables(superNames);
    
    int zeroRefDegs[numVar];
    for (int i = 0; i < numVar; ++i) {
        zeroRefDegs[i] = 0;
    }

    zeroRefDegs[0] = 1;
    zeroRef.setCoefficient(numVar, zeroRefDegs, RationalNumber(1));
    
    // update the dividend polynomial w.r.t. superNames:  --> divEx
    SparseMultivariateRationalPolynomial divEx(numVar);
    divEx = *this;
    divEx += zeroRef;
    divEx -= zeroRef;
    
    // update the triangular set w.r.t. superNames:       --> tsEx
    std::vector<SparseMultivariateRationalPolynomial> tsEx;
    tsEx.clear();
    for (int i = 0; i < tsSize; ++i){
        SparseMultivariateRationalPolynomial tmp;
        tmp = ts[i];
        tmp += zeroRef;
        tmp -= zeroRef;
        tsEx.push_back(tmp);
    }
    
    AltArr_t* g[tsSize];
    for (int i = 0; i < tsSize; ++i){
        g[i] = tsEx[tsSize-i-1].poly;
    }
    AltArr_t* f = divEx.poly;
    AltArr_t* a = NULL;

    AltArr_t* r = normalizePolynomial_AA (f, g, &a, tsSize, numVar);
        
    // convert outputs into SMQP:
    SparseMultivariateRationalPolynomial outR {r, numVar, zeroRef.names};
    
    SparseMultivariateRationalPolynomial aa {a, numVar, zeroRef.names};
    *A = aa;
    
    return outR;
}


/** 
 * Multi-Divisor Division (MDD) where the divisor-set is a Triangular Set
 * Given the  dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s],
 * Return (by reference) the remainder and the quotient set Q[s] = {q_0, ..., q_{s-1}}.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::triangularSetNormalForm(const TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial>& ts, std::vector<SparseMultivariateRationalPolynomial>* quoSet) const
{
    int verf = 0;
    int TYPE = 1; // the type of MDD:
                  // (0:heapMDD, 1:triangularSetMDD, 2:primitiveFactorTriangularSetMDD)
    int nvar = ts.numberOfVariables();
    std::vector<Symbol> tsVars = ts.allVariables();
    std::vector<Symbol> vars = orderPreservingSetDifference (variables(), tsVars);
    vars.reserve (vars.size() + tsVars.size());
    vars.insert (vars.end(), tsVars.begin(), tsVars.end());
    int tnvar = vars.size();
	
    std::vector<SparseMultivariateRationalPolynomial> polys = ts.polynomials();
    int nSet = polys.size();
	
    for (int i = 0; i < polys.size(); ++i){
		polys[i].setRingVariables(vars);
    }
	
    // update the dividend polynomials w.r.t superNames:  -->divEx
    SparseMultivariateRationalPolynomial divEx(nvar);
    divEx = *this;
    divEx.setRingVariables(vars);
    
    // make Node* polynomials for MDD's inputs:
    AltArr_t* g[nSet];
    AltArr_t* q[nSet];
    for (int i = 0; i < nSet; ++i){
		// g[i] = polys[nSet-i-1].poly;
		g[i] = polys[i].poly;
		q[i] = NULL;
    }

    AltArr_t* f = divEx.poly;
    AltArr_t* r = NULL;
    
    // compute the MDD:
    multiDivisorDivision_AA (f, g, q, &r, nSet, tnvar, TYPE);
    
    if (verf != 0){
		AltArr_t* h = NULL;
		verf = multiDivisorDivisionVerification_AA (f, g, q, r, h, nSet, tnvar);
		if (verf != 1){
			fprintf (stderr,
					 "SMQP Error: the correctness of multiDivisorDivision_AA algorithm\t FAILED\n");
			exit(1);
		}
    }
    
    // convert outputs into SMQP:
    SparseMultivariateRationalPolynomial outR {r, tnvar, divEx.names};
	
    quoSet->clear();
    for (int i = 0; i < nSet; ++i){
		SparseMultivariateRationalPolynomial qi {q[nSet-i-1], tnvar, divEx.names};
        // SparseMultivariateRationalPolynomial qi {q[i], tnvar, divEx.names};
        quoSet->push_back(qi);
    }

    return SparseMultivariateRationalPolynomial(outR);
}

/** 
 * Specialized Normal Form where the divisor-set is a Triangular Set
 * Given the  dividend, f, and a divisor-set of polynomials of size s,
 * G[s] = {g_0, ..., g_{s-1}} to compute the reduce polynomial (remainder) r with respect to the G[s].
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::triangularSetOnlyNormalForm(const TriangularSet<RationalNumber, SparseMultivariateRationalPolynomial>& ts) const
{
    int verf = 0;
    int nvar = ts.numberOfVariables();
    std::vector<Symbol> tsVars = ts.allVariables();
    std::vector<Symbol> vars = orderPreservingSetDifference (variables(), tsVars);
    vars.reserve (vars.size() + tsVars.size());
    vars.insert (vars.end(), tsVars.begin(), tsVars.end());
    int tnvar = vars.size();
    
    std::vector<SparseMultivariateRationalPolynomial> polys = ts.polynomials();
    int nSet = polys.size();
    
    for (int i = 0; i < polys.size(); ++i){
        polys[i].setRingVariables(vars);
    }
    
    // update the dividend polynomials w.r.t superNames:  -->divEx
    SparseMultivariateRationalPolynomial divEx(nvar);
    divEx = *this;
    divEx.setRingVariables(vars);
    
    // make Node* polynomials for MDD's inputs:
    AltArr_t* g[nSet];
    for (int i = 0; i < nSet; ++i){
        // g[i] = polys[nSet-i-1].poly;
        g[i] = polys[i].poly;
    }

    AltArr_t* f = divEx.poly;
    AltArr_t* r;
    
    // compute the MDD:
    r = onlyNormalForm_AA (f, g, nSet, tnvar);
    
    // convert outputs into SMQP:
    SparseMultivariateRationalPolynomial outR {r, tnvar, divEx.names};
    
    return SparseMultivariateRationalPolynomial(outR);
}

/** 
 * Do the pseudo division of c by the triangular-set (divisor set) of B in the naive principle
 * such that hPow*f = quoSet_0 * B_0 + ... + quoSet_{nSet-1} * B_{nSet-1} + rem.
 * Quotients are returned in quoSet, and remainder in rem. 
 * If hPow is not NULL then *hPow is set to the initial of b to the power of 
 * the exact number of division steps which occurred..
 * nvar : size of exponent vectors
 * nSet : the size of the divisor set 
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::triangularSetPseudoDivide (const TriangularSet<RationalNumber,SparseMultivariateRationalPolynomial>& ts, std::vector<SparseMultivariateRationalPolynomial>* quoSet, SparseMultivariateRationalPolynomial* h) const
{
    
    int lazy = 0;
    int verf = 0; 
    int nvar = ts.numberOfVariables();
    std::vector<Symbol> tsVars =ts.allVariables();
    std::vector<Symbol> vars = orderPreservingSetDifference (variables(), tsVars);
    vars.reserve (vars.size() + tsVars.size());
    vars.insert (vars.end(), tsVars.begin(), tsVars.end());
    int tnvar = vars.size();
    
    std::vector<SparseMultivariateRationalPolynomial> polys = ts.polynomials();
    int nSet = polys.size(); //ts.numberOfVariables();
    
    for (int i = 0; i < polys.size(); ++i){
	polys[i].setRingVariables(vars);
    }
    
    // update the dividend polynomial w.r.t. superNames:  --> divEx
    SparseMultivariateRationalPolynomial divEx(nvar);
    divEx = *this;
    divEx.setRingVariables(vars);
    
    // make Node* polynomials for MDPD's inputs:
    AltArr_t* g[nSet];
    AltArr_t* q[nSet];
    for (int i = 0; i < nSet; ++i){
	g[i] = polys[nSet-i-1].poly;
	// g[i] = polys[i].poly;
	q[i] = NULL;
    }
    AltArr_t* f = divEx.poly;
    AltArr_t* r = NULL;
    AltArr_t* hPow = NULL;
    
    //TODO: this will break if q or hPow is null... sometimes we only want remainder.
    // compute MDPD:
	multiDivisorPseudoDivide_AA (f, g, q, &r, &hPow, tnvar, lazy, nSet);
    
    if (verf != 0){
	verf = multiDivisorDivisionVerification_AA (f, g, q, r, hPow, nSet, tnvar);
	if (verf != 1){
	    fprintf (stderr, "SMQP Error: the correctness of multiDivisorPseudoDivision algorithm in C side \t FAILED\n");
	    exit(1);
	}
    }
    
    // convert outputs into SMQP:
    SparseMultivariateRationalPolynomial outR {r, tnvar, divEx.names};
    
    quoSet->clear();
    for (int i = 0; i < nSet; ++i){
    	SparseMultivariateRationalPolynomial qi {q[nSet-i-1], tnvar, divEx.names};
	// SparseMultivariateRationalPolynomial qi {q[i], tnvar, divEx.names};
	quoSet->push_back(qi);
	}
    
    if (h != NULL) {
        SparseMultivariateRationalPolynomial h_tmp {hPow,  tnvar, divEx.names};
        *h = h_tmp;
    } else {
        freePolynomial_AA(hPow);
    }
    return SparseMultivariateRationalPolynomial(outR);
}
