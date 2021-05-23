#include "RationalNumberPolynomial/mrpolynomial-altarr.hpp"

////////// Construtors ////////////////////////////////////

/**
 * Private constructor to create SMQP directly from a Node*.
 * Makes a copy of input varNames
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(AltArr_t* aa, int vars, Symbol* varNames) :
    poly(aa),
    nvar(vars)
{
    names = new Symbol[nvar+1];
    std::copy(varNames, varNames+nvar+1, names);
}

/**
 * Construct a multivariate polynomial
 *
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial() :
    poly(NULL),
    nvar(0)
{
    names = new Symbol[1];
    names[0] = "1";
}

/**
 * Construct a multivariate polynomial with specific number of variables
 *
 * @param v: Number of variables
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(int v) :
    poly(NULL),
    nvar(v)
{
    names = new Symbol[nvar+1];
    names[0] = "1";
    for (int i = 1; i <= nvar; ++i) {
        std::ostringstream convert;
        convert << i;
        names[i] = "x_";
        names[i] += convert.str();
    }
}

/**
 * Construct with a variable name such that f(x) = x;
 *
 * @param x: The variable name
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial (const Symbol& x) :
    poly(NULL),
    nvar(1)
{
    names = new Symbol[2];
    names[0] = "9";
    names[1] = x;

    mpq_t coef;
    mpq_init(coef);
    mpq_set_ui(coef, 1ul, 1ul);
    poly = makeConstPolynomial_AA(1, nvar, coef);
    degree_t degsList[nvar] = {0};
    degsList[0] = 1;
    setDegrees_AA_inp(poly, 0, degsList, 1);
    // poly->elems->degs = 1;
    mpq_clear(coef);
}

SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial (const std::string& str) :
poly(NULL),
nvar(0),
names(NULL)
{
    this->fromString(str);
}


/**
 * Copy Constructor.
 *
 * Does not reuse underlying memory allocated by b.
 *
 * @param b: A sparse multivariate polynomial
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(const SparseMultivariateRationalPolynomial& b) :
    nvar(b.nvar)
{
    poly = deepCopyPolynomial_AA(b.poly);
    names = new Symbol[nvar+1];
    std::copy(b.names, b.names+nvar+1, names);
    slp = std::vector<SLPRepresentation>(b.slp);
}

/**
 * Move Constructor.
 *
 * @params b: The r-value reference polynomial.
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(SparseMultivariateRationalPolynomial&& b) {
    nvar = b.nvar;
    poly = b.poly;
    names = new Symbol[nvar+1]; //many functions rely on the fact that moved from polys retain their variable ordering
    std::copy(b.names, b.names+nvar+1, names);
    slp = b.slp;

    b.poly = NULL;
    b.slp.clear();
}


/**
 * Create an SMZP from a SMQP.
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(const SparseMultivariateIntegerPolynomial& b) {
    nvar = b.nvar;
    poly = deepCopyPolynomial_AAFromAAZ(b.poly);
    names = new Symbol[nvar+1];
    std::copy(b.names, b.names+nvar+1, names);
}

/**
 * Create a SMQP from a Integer.
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(const Integer& r, int nv) :
    nvar(nv),
    poly(NULL)
{
    mpq_t coef;
    mpq_init(coef);
    mpz_set(mpq_numref(coef), r.get_mpz_t());
    poly = makeConstPolynomial_AA(1, nvar, coef);
    mpq_clear(coef);

    names = new Symbol[nvar+1];
    names[0] = "1";
    for (int i = 1; i <= nvar; ++i) {
        std::ostringstream convert;
        convert << i;
        names[i] = "x_";
        names[i] += convert.str();
    }
}

/**
 * Create a SMQP from a RationalNumber.
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(const RationalNumber& r, int nv) :
    nvar(nv),
    poly(NULL)
{
    poly = makeConstPolynomial_AA(1, nvar, r.get_mpq_t());

    names = new Symbol[nvar+1];
    names[0] = "1";
    for (int i = 1; i <= nvar; ++i) {
        std::ostringstream convert;
        convert << i;
        names[i] = "x_";
        names[i] += convert.str();
    }

}

/**
 * Create a SMQP from a univariate rational polynomial.
 *
 * @param p: A SUQP polynomial.
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial (const DenseUnivariateRationalPolynomial& p) :
    nvar(1),
    poly(NULL)
{
    names = new Symbol[nvar+1];
    names[0] = "9";
    names[1] = p.variable();

    RationalNumber coef;
    long int pDeg = p.degree().get_si();
    poly = makePolynomial_AA(pDeg+1, nvar);
    int curSize = 0;
    degree_t degsList[nvar] = {0};
    for (long int i = pDeg; i >= 0; --i) {
        coef = p.coefficient(i);
        if (coef != 0) {
            mpq_init(poly->elems[curSize].coef);
            mpq_set(poly->elems[curSize].coef, coef.get_mpq_t());
            degsList[0] = i;
            setDegrees_AA_inp(poly, curSize, degsList, nvar);
            // poly->elems[curSize].degs = i;
            ++curSize;
        }
    }
    poly->size = curSize;
}

/**
 * Construct from a SUP<SMQP> polynomial
 *
 * @param s: The SUP<SMQP> polynomial
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial (const SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>& s) :
    nvar(0),
    poly(NULL)
{
    names = new Symbol[1];
    names[0] = "1";

    long long int d = s.degree().get_si();
    SparseMultivariateRationalPolynomial c = s.coefficient(d);
    while (c.isZero()) {
        --d;
        if (d < 0) {
            return;
        }
        c = s.coefficient(d);
    }
    if (c.nvar == 0) {
        //then we have a bunch of "constant" SMQP's as coefficients in SUP.
        delete[] names;
        names = new Symbol[2];
        names[0] = "9";
        names[1] = s.variable();
        nvar = 1;

        poly = makePolynomial_AA(d+1, nvar);
        int curSize = 0;
        degree_t degsList[nvar] = {0};
        for (int k = d; k >= 0; --k) {
            SparseMultivariateRationalPolynomial coef = s.coefficient(k);
            if (coef.isZero()) {
                continue;
            }
            mpq_init(poly->elems[curSize].coef);
            mpq_set(poly->elems[curSize].coef, coef.poly->elems[0].coef);
            degsList[0] = k;
            setDegrees_AA_inp(poly, curSize, degsList, nvar);
            // poly->elems[curSize].degs = k;
            ++curSize;
        }
        poly->size = curSize;

        return;
    }

    nvar = c.nvar + 1;
    AltArr_t* newPoly = makePolynomial_AA(1, nvar);
    delete[] names;
    names = new Symbol[nvar+1];
    names[0] = "9";
    names[1] = s.variable();
    std::copy(c.names+1, c.names+1+c.nvar, names+2);

    for (int k = 0; k <= d; ++k) {
        c = s.coefficient(k);
        //coefficient is 0 so do nothing
        if (c.isZero()) {
            continue;
        }

        //Get the coef poly, expanding it's num vars and multiply through
        //by x^k.
        expandNumVarsLeft_AA(c.poly, nvar);
        c.poly = mainLShiftPolynomial_AA_inp(c.poly, k);
        // for (int i = 0; i < c.poly->size; ++i) {
        //     c.poly->elems[i].degs |= ((degrees_t) k << mvarDegOffset);
        // }
        newPoly = addPolynomials_AA_inp(newPoly, c.poly, nvar);
    }

    this->poly = newPoly;
}

/**
 * Destroy the polynomial and release underlying node memory.
 **/
SparseMultivariateRationalPolynomial::~SparseMultivariateRationalPolynomial() {
    delete[] names;
    freePolynomial_AA(poly);
    slp.clear();
}



////////// BPASRing /////////////////////////////////

bool SparseMultivariateRationalPolynomial::isZero() const {
    return isZero_AA(poly);
}

void SparseMultivariateRationalPolynomial::zero() {
    freePolynomial_AA(poly);
    poly = NULL;
}

bool SparseMultivariateRationalPolynomial::isOne() const {
    return isOne_AA(poly);
}

void SparseMultivariateRationalPolynomial::one() {
    freePolynomial_AA(poly);
    RationalNumber r(1);
    poly = makeConstPolynomial_AA(1, nvar, r.get_mpq_t());
}

bool SparseMultivariateRationalPolynomial::isNegativeOne() const {
    return isNegativeOne_AA(poly);
}

void SparseMultivariateRationalPolynomial::negativeOne() {
    freePolynomial_AA(poly);
    RationalNumber r(-1);
    poly = makeConstPolynomial_AA(1, nvar, r.get_mpq_t());
}

int SparseMultivariateRationalPolynomial::isConstant() const {
    return isConstant_AA(poly);
}

void SparseMultivariateRationalPolynomial::preparePolysForSRC(const SparseMultivariateRationalPolynomial& q, const Symbol& v, std::vector<Symbol>& superRing, bool sameRing, int* varMap, int& swapIdx, AltArrZ_t** ppZ, AltArrZ_t** qqZ) const {

    // printVariables(this->ringVariables(), "pRingVars");
    // printVariables(q.ringVariables(), "qRingVars");
    // printVariables(superRing,  "superRing");
    int swappedIdx = 0;
    //could we check if v is not the main variable after shrinking the nvars below?
    if (superRing[0] != v) {
        //by previous check ons leadingVariable and superRing[0] iterator must not be invalid;
        std::vector<Symbol>::iterator it = std::find(superRing.begin(), superRing.end(), v);

        swappedIdx = it - superRing.begin();
        Symbol tmpV = superRing[0];
        superRing[0] = superRing[swappedIdx];
        superRing[swappedIdx] = tmpV;
        sameRing = false;
    }

    swapIdx = swappedIdx;


    AltArrZ_t *pZ, *qZ;
    if (!sameRing) {
        //Get a CPP version of the polys for the sake of doing the set ring vars.
        SparseMultivariateIntegerPolynomial P = this->primitivePartSMZP();
        SparseMultivariateIntegerPolynomial Q = q.primitivePartSMZP();
        P.setRingVariables(superRing);
        Q.setRingVariables(superRing);
        pZ = P.poly;
        qZ = Q.poly;
        P.poly = NULL;
        Q.poly = NULL;
    } else {
        mpq_t tmp;
        mpq_init(tmp);
        pZ = primitivePartAndContent_AAZFromAA(this->poly, tmp);
        qZ = primitivePartAndContent_AAZFromAA(q.poly, tmp);
        mpq_clear(tmp);
    }

    if (mainDegree_AAZ(pZ) < mainDegree_AAZ(qZ)) {
        AltArrZ_t* tmp = pZ;
        pZ = qZ;
        qZ = tmp;
    }

    //at this point, all corner cases checked, p and q are in same ring, and deg(p, v) >= deg(q,v);
    //now, we must shrink to smallest space where p or q have positive degree

    // for (int i = 0; i < origNvar; ++i) {
    //     std::cerr << "varmap[" << i << "]: " << varMap[i] << std::endl;
    // }
    int trueNvar = computeShrinkMapPolys_AAZ(pZ, qZ, varMap);

    // std::cerr << "truNvar: " << trueNvar << std::endl;
    // for (int i = 0; i < origNvar; ++i) {
    //     std::cerr << "varmap[" << i << "]: " << varMap[i] << std::endl;
    // }

    if (pZ->nvar != trueNvar) {
        shrinkAndReorderVars_AAZ(pZ, varMap, pZ->nvar);
        shrinkAndReorderVars_AAZ(qZ, varMap, qZ->nvar);
    }

    *ppZ = pZ;
    *qqZ = qZ;
}

/**
 * Given a polynomial q, compute the subresultant chain between this and q, viewing both polynomial
 * recursively with main variable v.
 *
 * @note, if this and q do not exist in the same ambient space, the space of p will define the ordering of the subresultants.
 *
 * @param q, the other polynomial for which to compute the subresultant chain.
 * @param v, the main variable to be used when computing the subresultant chain.
 *
 * @return a vector containing the subresultant chain whereby index 0 is resultant and subresultant degrees increase with index.
 *
 */
std::vector<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::subresultantChain(const SparseMultivariateRationalPolynomial& q, const Symbol& v, bool filled) const {
    if (this->leadingVariable() != v && q.leadingVariable() != v) {
        throw std::invalid_argument("In SMQP::subresultantChain at least one of the input polynomials must have " + v.toString() + " as leading variable.");
    }

    // // Debug Print
    // std::cerr << "[SMQP+] subresultantChain mdeg(p)= " << this->mainDegree() << " mdeg(q)= " << q.mainDegree() << std::endl; 

    std::vector<Symbol> superRing;
    bool sameRing = orderPreservingSetUnion(this->ringVariables(), q.ringVariables(), superRing);

    std::vector<SparseMultivariateRationalPolynomial> chain;

    if (isZero()) {
        chain.emplace_back(q);
        chain[0].setRingVariables(superRing);
        return chain;
    }
    if (q.isZero()) {
        chain.emplace_back(superRing.size());
        chain[0].setRingVariables(superRing);
        return chain;
    }

    bool foundVP = false;
    for (int i = 0; i < nvar && !foundVP; ++i) {
        if (names[i+1] == v) {
            foundVP = true;
        }
    }
    bool foundVQ = false;
    for (int i = 0; i < q.nvar && !foundVQ; ++i) {
        if (q.names[i+1] == v) {
            foundVQ = true;
        }
    }
    foundVP = foundVP || (this->degree(v) == 0);
    foundVQ = foundVQ || (q.degree(v) == 0);


    // std::cerr << "foundVP : " << foundVP << " , foundVQ: " << foundVQ << std::endl;

    if (!foundVP && !foundVQ) {
        chain.emplace_back(*this); //do this so chain exists in proper ambient space
        chain[0].one();
        chain[0].setRingVariables(superRing);
        return chain;
    } else if (!foundVQ) {
        //degree of q in v is 0
        chain.push_back(q ^ this->mainDegree());
        chain[0].setRingVariables(superRing);
        return chain;
    } else if (!foundVP) {
        chain.push_back(*this ^ q.mainDegree());
        chain[0].setRingVariables(superRing);
        return chain;
    }

    AltArrZ_t* pZ = NULL;
    AltArrZ_t* qZ = NULL;
    int origNvar = superRing.size();
    int varMap[origNvar];
    int swappedIdx = 0;
    preparePolysForSRC(q, v, superRing, sameRing, varMap, swappedIdx, &pZ, &qZ);


    AltArrsZ_t* SC = NULL;
    int size = 0;
    SubresultantChainMode mode = AUTO;

#if !(defined(SERIAL) && SERIAL) && (defined(PARALLEL_SUBRES) && PARALLEL_SUBRES)
    long min_mdeg = (long) MIN(mainLeadingDegree_AAZ(pZ), mainLeadingDegree_AAZ(qZ));
    long max_pdeg = (long) MAX(partialDegree_AAZ(pZ, 1), partialDegree_AAZ(qZ, 1));
    if (pZ->nvar == 2 && ( min_mdeg > 20 && max_pdeg > 10) ) {
        DUSP::ParallelSubresultantChainZ(pZ, qZ, min_mdeg, max_pdeg, &SC, &size);
    } else {
        SubresultantChainZ (mode, pZ,qZ, &SC, &size);
    }
#else
    SubresultantChainZ (mode, pZ,qZ, &SC, &size);
#endif

    // // Debug Print
    // std::cerr << "[SMQP+] subresultantChain... done!" << std::endl; 

    freePolynomial_AAZ (pZ);
    freePolynomial_AAZ (qZ);

    AltArr_t* tmp;
    AltArrsZ_t* cur = SC;

    chain.reserve (size);

    int swapMap[origNvar];
    if (swappedIdx != 0) {
        for (int i = 0; i < origNvar; ++i) {
            swapMap[i] = i;
        }
        swapMap[swappedIdx] = 0;
        swapMap[0] = swappedIdx;

        //reset superRing also
        Symbol tmpV = superRing[0];
        superRing[0] = superRing[swappedIdx];
        superRing[swappedIdx] = tmpV;
    }

    //despite the name, it is only used as the zero polynomial in the second loop if filled==1
    SparseMultivariateRationalPolynomial zero;
    zero.zero();
    zero.setRingVariables(superRing);

    for (int i = 0; cur != NULL && i < size; ++i){
        reverseShrinkVariables_AAZ_inp(cur->poly, origNvar, varMap);
        if (swappedIdx != 0) {
            reorderVars_AAZ (cur->poly, swapMap, origNvar);
        }

        tmp = deepCopyPolynomial_AAFromAAZ (cur->poly);
        zero.poly = tmp;
        chain.emplace_back(std::move(zero));
        zero.poly = NULL;
        cur = cur->next;
    }

    freeAltArrsZ (SC); // free

    if (filled && chain.size() > 1) {
        degree_t fullSize = mainLeadingDegree_AA (chain[chain.size()-2].poly) + 2;
        degree_t delta;

        if (chain.size() < fullSize) {
            chain.reserve (fullSize);
            for (int i = chain.size()-2; i > 0; --i){
                if (mainLeadingDegree_AA (chain[i].poly) != mainLeadingDegree_AA (chain[i-1].poly) + 1) {
                    delta = mainLeadingDegree_AA (chain[i].poly) - mainLeadingDegree_AA (chain[i-1].poly);
                    if (i > 1) {
                        i = i-1;
                        for (int j = 0; j < delta-2; ++j)
                            chain.insert (chain.begin()+i,zero);
                    }
                    else {
                        for (int j=0; j<delta-1; ++j)
                            chain.insert(chain.begin()+i,zero);
                    }
                }
            }
            if (mainLeadingDegree_AA (chain[0].poly) != 0){
                for (int j = 0; j < mainLeadingDegree_AA (chain[0].poly); ++j){
                    chain.insert (chain.begin(),zero);
                }
            }
        }
    }

    return chain;
}

/**
 * Subresultant Chain
 * Return the list of subresultants
 *
 * @param q: The other sparse univariate polynomial
 **/
std::vector<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::subresultantChain (const SparseMultivariateRationalPolynomial& q, int filled) const {

    // // Debug Print
    // std::cerr << "[SMQP+] subresultantChain mdeg(p)= " << this->mainDegree() << " mdeg(q)= " << q.mainDegree() << std::endl; 

    if (isZero()) {
    	std::vector<SparseMultivariateRationalPolynomial> zeroChain_p;
    	zeroChain_p.push_back (q);
    	return zeroChain_p;
    }

    if (q.isZero()) {
    	std::vector<SparseMultivariateRationalPolynomial> zeroChain_q;
    	zeroChain_q.push_back (*this);
    	return zeroChain_q;
    }

    bool TYPE = 0; // 0: SMZP Ducos Algorithm, 1: SUP Ducos Algorothm

    Symbol v = q.leadingVariable();
    // commented to support constant polynomials in algebraic extension
    // if (v != leadingVariable()) {
    //     std::cout << "BPAS: error, cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
    //     exit(1);
    // }

    if (TYPE == 0){
    // // Subresultant from SMZP:
        std::vector<int> xs;
        bool isOrdered = isOrderedRing(q, xs);
        SparseMultivariateRationalPolynomial ppP, ppQ;
        std::vector<SparseMultivariateRationalPolynomial> subr;

        if (!isOrdered) {
            std::cout << "BPAS: error, trying to compute subresultant chain between Q[";
            for (int i = 1; i <= nvar; ++i) {
                std::cout << names[i];
                if (i < nvar) { std::cout << ", "; }
            }
            std::cout << "] and Q[";
            for (int i = 1; i <= q.nvar; ++i) {
                std::cout << q.names[i];
                if (i < q.nvar) { std::cout << ", "; }
            }
            std::cout << "]." << std::endl;
            exit(1);
        }

        int superNvar = xs.size() / 2;
        if (superNvar != nvar || superNvar != q.nvar) {

    //map indices to the expanded superset.
            int varmap[nvar];
            int qvarmap[q.nvar];
            for (int i = 0; i < xs.size(); i += 2) {
                if (xs[i] != 0) {
                    varmap[xs[i]-1] = i/2;
                }
                if (xs[i+1] != 0) {
                    qvarmap[xs[i+1]-1] = i/2;
                }
            }

    //create new combined names array
            Symbol newnames[superNvar+1];
    //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
            newnames[0] = this->names[0];
            for (int i = 0; i < xs.size(); i += 2) {
                if (xs[i] != 0) {
                    newnames[(i/2) + 1] = names[xs[i]];
                } else {
                    newnames[(i/2) + 1] = q.names[xs[i+1]];
                }
            }

            ppP = this->expandVariables(superNvar, newnames, varmap);
            ppQ = q.expandVariables(superNvar, newnames, qvarmap);
    // ppP = ppP.primitivePart();
    // ppQ = ppQ.primitivePart();
        } else {
            ppP = *this;
            ppQ = q;

    // ppP = this->primitivePart();
    // ppQ = q.primitivePart();
        }

        AltArr_t* Pq = primitivePart_AA (ppP.poly);
        AltArr_t* Qq = primitivePart_AA (ppQ.poly);

        AltArrZ_t* P = deepCopyPolynomial_AAZFromAA (Pq);
        AltArrZ_t* Q = deepCopyPolynomial_AAZFromAA (Qq);

        freePolynomial_AA (Pq);
        freePolynomial_AA (Qq);

        // AltArrZ_t* P = deepCopyPolynomial_AAZFromAA (ppP.poly);
        // AltArrZ_t* Q = deepCopyPolynomial_AAZFromAA (ppQ.poly);

        int lvarP = leadingVariable_AAZ (P);
        int lvarQ = leadingVariable_AAZ (Q);
        bool isShrinked = 0;

        // commented to support constant polynomials in algebraic extension
        // if (lvarP != lvarQ){
        //     std::cout << "BPAS Error: cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
        //     exit(1);
        // }

        // if (lvarP > 0){
        //     for (int i = 0; i < lvarP; ++i){
        //         shrinkNumVarsAtIdx_AAZ (P, 0);
        //         shrinkNumVarsAtIdx_AAZ (Q, 0);
        //     }
        //     isShrinked = 1;
        // }

		if (lvarP > 0) {
			int varmap[superNvar];
			int j = 0;

			for (int i = 0; i < lvarP; i++) {
				varmap[i] = -1;
			}

			for (int i = lvarP; i < superNvar; i++) {
				varmap[i] = j;
				j++;
			}

			shrinkAndReorderVars_AAZ (P, varmap, superNvar);
			shrinkAndReorderVars_AAZ (Q, varmap, superNvar);
			isShrinked = 1;
		}

        // if (leadingVariable_AAZ (P) != 0 || leadingVariable_AAZ (Q) != 0){
        //     std::cout << "BPAS Error: shrink in subresultantChain does not work properly!" << std::endl;
        //     exit(1);
        // }

        AltArrsZ_t* SC;
        AltArr_t* tmp;
        int size = 0;
        SubresultantChainMode mode = AUTO;

#if !(defined(SERIAL) && SERIAL) && defined(PARALLEL_SUBRES) && PARALLEL_SUBRES
    long min_mdeg = (long) MIN(mainLeadingDegree_AAZ(P), mainLeadingDegree_AAZ(Q));
    long max_pdeg = (long) MAX(partialDegree_AAZ(P, 1), partialDegree_AAZ(Q, 1));
    if (P->nvar == 2 && (min_mdeg > 20 && max_pdeg > 10)) {
            DUSP::ParallelSubresultantChainZ(P, Q, min_mdeg, max_pdeg, &SC, &size);
    } else {
            SubresultantChainZ (mode, P,Q, &SC, &size);
    }
#else
    SubresultantChainZ (mode, P,Q, &SC, &size);
#endif

        freePolynomial_AAZ (P); 
        freePolynomial_AAZ (Q);

        AltArrsZ_t* cur = SC;
        subr.reserve (size);
        for (int i = 0; cur != NULL && i < size; ++i){
        tmp = deepCopyPolynomial_AAFromAAZ (cur->poly);
        if (isShrinked && tmp != NULL && tmp->size != 0){
            expandNumVarsLeft_AA (tmp, ppP.nvar);
            if (tmp->nvar != ppP.nvar){
                std::cout << "BPAS Error: expand in subresultantChain does not work properly!" << std::endl;
                exit(1);
            }
        }
        subr.push_back (SparseMultivariateRationalPolynomial (tmp, ppP.nvar, ppP.names));

        cur = cur->next;
        }

    freeAltArrsZ (SC); // free

    if (filled && subr.size() > 1){

    // std::cout << "[Ali-TEST] In SMQP, subr.size =  " << subr.size() << std::endl;

        std::vector<SparseMultivariateRationalPolynomial> chain = subr;
        SparseMultivariateRationalPolynomial zero;
        zero.zero();

    // std::cout << "[Ali-TEST] In SMQP, poly =  " << chain[chain.size()-2] << std::endl;

        degree_t fullSize = mainLeadingDegree_AA (chain[chain.size()-2].poly) + 2;
        degree_t delta;

    // std::cerr << "chain.size() = " << chain.size() << std::endl;
    // std::cerr << "fullSize = " << fullSize << std::endl;

        if (chain.size() < fullSize){
            chain.reserve (fullSize);
            for (int i = chain.size()-2; i > 0; --i){
                if (mainLeadingDegree_AA (chain[i].poly) != mainLeadingDegree_AA (chain[i-1].poly) + 1) {
                    delta = mainLeadingDegree_AA (chain[i].poly) - mainLeadingDegree_AA (chain[i-1].poly);
                    if (i > 1) {
                        i = i-1;
                        for (int j = 0; j < delta-2; ++j)
                            chain.insert (chain.begin()+i,zero);
                    }
                    else {
                        for (int j=0; j<delta-1; ++j)
                            chain.insert(chain.begin()+i,zero);
                    }
                }
            }
            if (mainLeadingDegree_AA (chain[0].poly) != 0){
                for (int j = 0; j < mainLeadingDegree_AA (chain[0].poly); ++j){
                    chain.insert (chain.begin(),zero);
                }
            }
        }
    // std::cerr << "chain.size() = " << chain.size() << std::endl;

        return chain;
    }

    return subr;
    } else {
    // Subresultant from SUP:
        std::cout << "SUP Subresultant" << std::endl;
        std::vector<SparseMultivariateRationalPolynomial> Out;
        std::vector<SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>> S;
        SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> a, b;
        SparseMultivariateRationalPolynomial temp;
        if (q.leadingVariableDegree() > leadingVariableDegree()) {
            a = q.convertToSUP(v);
            b = convertToSUP(v);
        }
        else {
            a = convertToSUP(v);
            b = q.convertToSUP(v);
        }
        S = a.subresultantChain(b);
        Out.reserve(S.size());
        for (int i=0; i<S.size(); ++i) {
            temp = SparseMultivariateRationalPolynomial(S[i]);
            Out.push_back(temp);
        }
        return Out;
    }
}

/**
 * Given a polynomial q, compute the subresultant of index idx between this and q, viewing both polynomial
 * recursively with main variable v. This function in fact computes the subresultant of index idx and idx+1
 * and so both are returned, unless idx+1 does not exist.
 *
 * @note, if the subresultant of index idx+1 is degenerative then many subresultants will be returned until the first non-degenerative one is found.
 * @note, if this and q do not exist in the same ambient space, the space of p will define the ordering of the subresultants.
 *
 * @param q, the other polynomial for which to compute the subresultant chain.
 * @param v, the main variable to be used when computing the subresultant chain.
 *
 * @return a vector containing the subresultant chain whereby index 0 is requested subresultant idx subresultant degrees increase with index.
 *
 */
std::vector<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::subresultantAtIdx (const SparseMultivariateRationalPolynomial& q, const Symbol& v, int idx, specSRC_AAZ ** lazyInfo) const {
    
    // // Debug Print
    // std::cerr << "[SMQP+] subresultantAtIdx mdeg(p)= " << this->mainDegree() << " mdeg(q)= " << q.mainDegree() << " index= " << idx << std::endl; 
    // std::cerr << "[SMQP+] p = " << this << std::endl;
    // std::cerr << "[SMQP+] q = " << q << std::endl;


    if (this->leadingVariable() != v && q.leadingVariable() != v) {
        throw std::invalid_argument("In SMQP::subresultantChainAtIdx_new at least one of the input polynomials must have " + v.toString() + " as leading variable.");
    }

    std::vector<Symbol> superRing;
    bool sameRing = orderPreservingSetUnion(this->ringVariables(), q.ringVariables(), superRing);

    std::vector<SparseMultivariateRationalPolynomial> chain;

    if (isZero()) {
        chain.emplace_back(q);
        chain[0].setRingVariables(superRing);
        return chain;
    }
    if (q.isZero()) {
        chain.emplace_back(superRing.size());
        chain[0].setRingVariables(superRing);
        return chain;
    }

    bool foundVP = false;
    for (int i = 0; i < nvar && !foundVP; ++i) {
        if (names[i+1] == v) {
            foundVP = true;
        }
    }
    bool foundVQ = false;
    for (int i = 0; i < q.nvar && !foundVQ; ++i) {
        if (q.names[i+1] == v) {
            foundVQ = true;
        }
    }
    foundVP = foundVP || (this->degree(v) == 0);
    foundVQ = foundVQ || (q.degree(v) == 0);


    // std::cerr << "foundVP : " << foundVP << " , foundVQ: " << foundVQ << std::endl;

    if (!foundVP && !foundVQ) {
        chain.emplace_back(*this); //do this so chain exists in proper ambient space
        chain[0].one();
        chain[0].setRingVariables(superRing);
        return chain;
    } else if (!foundVQ) {
        //degree of q in v is 0
        chain.push_back(q ^ this->mainDegree());
        chain[0].setRingVariables(superRing);
        return chain;
    } else if (!foundVP) {
        chain.push_back(*this ^ q.mainDegree());
        chain[0].setRingVariables(superRing);
        return chain;
    }


    AltArrZ_t* pZ = NULL;
    AltArrZ_t* qZ = NULL;
    int origNvar = superRing.size();
    int varMap[origNvar];
    int swappedIdx = 0;
    preparePolysForSRC(q, v, superRing, sameRing, varMap, swappedIdx, &pZ, &qZ);

    AltArrZ_t* SC_idx = NULL;
    AltArrZ_t* SC_idx1 = NULL;
    AltArr_t* tmp;

    SubresultantChainMode mode = AUTO;

#if !(defined(SERIAL) && SERIAL) && (defined(PARALLEL_SUBRES) && PARALLEL_SUBRES)
    // fprintf(stderr, "PARALLEL_HGCD-SUBRES... ON\n");
    long min_mdeg = (long) MIN(mainLeadingDegree_AAZ(pZ), mainLeadingDegree_AAZ(qZ));
    long max_pdeg = (long) MAX(partialDegree_AAZ(pZ, 1), partialDegree_AAZ(qZ, 1));
    if (pZ->nvar == 2 && ( min_mdeg > 20 && max_pdeg > 10) ) {
        AltArrsZ_t* SS = NULL;
        int size = 0;
        int flag = DUSP::hgcdBiModularFFTSubresultantChainZ(pZ, qZ, idx, &SS, &size);
		if (flag) {
			if (size == 2) {
				SC_idx1 = deepCopyPolynomial_AAZ(SS->poly);
				SC_idx = deepCopyPolynomial_AAZ(SS->next->poly);
			} else if (size == 1) {
				SC_idx1 = deepCopyPolynomial_AAZ(SS->poly);
				SC_idx = NULL;
			} else {
				SC_idx1 = NULL;
				SC_idx = NULL;
			}
			freeAltArrsZ (SS);
		} else {
        SubresultantChainAtIdxZ (mode, pZ, qZ, idx, &SC_idx, &SC_idx1, lazyInfo);
        }
    } else {
        SubresultantChainAtIdxZ (mode, pZ, qZ, idx, &SC_idx, &SC_idx1, lazyInfo);
    }
#else
    // fprintf(stderr, "PARALLEL_HGCD-SUBRES... OFF\n");
        SubresultantChainAtIdxZ (mode, pZ, qZ, idx, &SC_idx, &SC_idx1, lazyInfo);
#endif

    // int degQ = mainDegree_AAZ(qZ);
    freePolynomial_AAZ (pZ); // free
    freePolynomial_AAZ (qZ); // free

    // int mdeg = mainDegree_AAZ(SC_idx1);
    chain.reserve (2);

    int swapMap[origNvar];
    if (swappedIdx != 0) {
        for (int i = 0; i < origNvar; ++i) {
            swapMap[i] = i;
        }
        swapMap[swappedIdx] = 0;
        swapMap[0] = swappedIdx;

        //reset superRing also
        Symbol tmpV = superRing[0];
        superRing[0] = superRing[swappedIdx];
        superRing[swappedIdx] = tmpV;
    }

    //despite the name, it is only used as the zero polynomial in the second loop if filled==1
    SparseMultivariateRationalPolynomial zero;
    zero.zero();
    zero.setRingVariables(superRing);

    reverseShrinkVariables_AAZ_inp(SC_idx, origNvar, varMap);
    if (swappedIdx != 0) {
        reorderVars_AAZ(SC_idx, swapMap, origNvar);
    }
    tmp = deepCopyPolynomial_AAFromAAZ (SC_idx);
    freePolynomial_AAZ(SC_idx);
    zero.poly = tmp;
    chain.emplace_back(std::move(zero));
    zero.poly = NULL;

    reverseShrinkVariables_AAZ_inp(SC_idx1, origNvar, varMap);
    if (swappedIdx != 0) {
        reorderVars_AAZ(SC_idx1, swapMap, origNvar);
    }
    tmp = deepCopyPolynomial_AAFromAAZ(SC_idx1);
    freePolynomial_AAZ(SC_idx1);
    zero.poly = tmp;
    chain.emplace_back(std::move(zero));
    zero.poly = NULL;

    return chain;
}

std::vector<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::subresultantInitialAtIdx (const SparseMultivariateRationalPolynomial& q, const Symbol& v, int idx, Integer& mdegIdx, Integer& mdegIdx1, specSRC_AAZ **lazyInfo) const {

    // // Debug Print
    // std::cerr << "[SMQP+] subresultantInitialAtIdx mdeg(p)= " << this->mainDegree() << " mdeg(q)= " << q.mainDegree() << " index= " << idx << std::endl; 

    if (this->leadingVariable() != v && q.leadingVariable() != v) {
        throw std::invalid_argument("In SMQP::PrincipleCoefficientAtIdx_new at least one of the input polynomials must have " + v.toString() + " as leading variable.");
    }

    std::vector<Symbol> superRing;
    bool sameRing = orderPreservingSetUnion(this->ringVariables(), q.ringVariables(), superRing);

    std::vector<SparseMultivariateRationalPolynomial> chain;

    if (isZero()) {
        chain.emplace_back(q);
        chain[0].setRingVariables(superRing);
        return chain;
    }
    if (q.isZero()) {
        chain.emplace_back(superRing.size());
        chain[0].setRingVariables(superRing);
        return chain;
    }

    bool foundVP = false;
    for (int i = 0; i < nvar && !foundVP; ++i) {
        if (names[i+1] == v) {
            foundVP = true;
        }
    }
    bool foundVQ = false;
    for (int i = 0; i < q.nvar && !foundVQ; ++i) {
        if (q.names[i+1] == v) {
            foundVQ = true;
        }
    }
    foundVP = foundVP || (this->degree(v) == 0);
    foundVQ = foundVQ || (q.degree(v) == 0);


    // std::cerr << "foundVP : " << foundVP << " , foundVQ: " << foundVQ << std::endl;

    if (!foundVP && !foundVQ) {
        chain.emplace_back(*this); //do this so chain exists in proper ambient space
        chain[0].one();
        chain[0].setRingVariables(superRing);
        return chain;
    } else if (!foundVQ) {
        //degree of q in v is 0
        chain.push_back(q ^ this->mainDegree());
        chain[0].setRingVariables(superRing);
        return chain;
    } else if (!foundVP) {
        chain.push_back(*this ^ q.mainDegree());
        chain[0].setRingVariables(superRing);
        return chain;
    }


    AltArrZ_t* pZ = NULL;
    AltArrZ_t* qZ = NULL;
    int origNvar = superRing.size();
    int varMap[origNvar];
    int swappedIdx = 0;
    preparePolysForSRC(q, v, superRing, sameRing, varMap, swappedIdx, &pZ, &qZ);

    AltArr_t* tmp;
    PolyPairZ_t *pc_idx=NULL, *pc_idx1=NULL;
    SubresultantChainMode mode = AUTO;

    SubresultantInitAtIdxZ (mode, pZ, qZ, idx, &pc_idx, &pc_idx1, lazyInfo);
    AltArrZ_t* SC_idx = NULL;     AltArrZ_t* SC_idx1 = NULL;
    if (pc_idx != NULL) {SC_idx = pc_idx->poly;}
    if (pc_idx1 != NULL) {SC_idx1 = pc_idx1->poly;}

    // // Debug Print
    // std::cerr << "[SMQP+] subresultantInitialAtIdx... done!"<< std::endl; 

    //return main degrees of subreusltants.
    mdegIdx = pc_idx->d;
    mdegIdx1 = pc_idx1->d;

    if (pc_idx != NULL) {
        free (pc_idx);
    }
    if (pc_idx1 != NULL) {
        free (pc_idx1);
    }

    // int degQ = mainDegree_AAZ(qZ);
    freePolynomial_AAZ (pZ); // free
    freePolynomial_AAZ (qZ); // free

    // int mdeg = mainDegree_AAZ(SC_idx1);
    chain.reserve (2);

    int swapMap[origNvar];
    if (swappedIdx != 0) {
        for (int i = 0; i < origNvar; ++i) {
            swapMap[i] = i;
        }
        swapMap[swappedIdx] = 0;
        swapMap[0] = swappedIdx;

        //reset superRing also
        Symbol tmpV = superRing[0];
        superRing[0] = superRing[swappedIdx];
        superRing[swappedIdx] = tmpV;
    }

    //despite the name, it is only used as the zero polynomial in the second loop if filled==1
    SparseMultivariateRationalPolynomial zero;
    zero.zero();
    zero.setRingVariables(superRing);

    reverseShrinkVariables_AAZ_inp(SC_idx, origNvar, varMap);
    if (swappedIdx != 0) {
        reorderVars_AAZ(SC_idx, swapMap, origNvar);
    }
    tmp = deepCopyPolynomial_AAFromAAZ (SC_idx);
    freePolynomial_AAZ(SC_idx);
    zero.poly = tmp;
    chain.emplace_back(std::move(zero));
    zero.poly = NULL;

    reverseShrinkVariables_AAZ_inp(SC_idx1, origNvar, varMap);
    if (swappedIdx != 0) {
        reorderVars_AAZ(SC_idx1, swapMap, origNvar);
    }
    tmp = deepCopyPolynomial_AAFromAAZ(SC_idx1);
    freePolynomial_AAZ(SC_idx1);
    zero.poly = tmp;
    chain.emplace_back(std::move(zero));
    zero.poly = NULL;

    return chain;
}


// /**
//  * Subresultant Chain At Idx
//  * Return idx-th subresultant chain and next (with higher degree) polynomial in the chain.
//  * @note idx default is 0 (to return last two subresultants)
//  **/
// std::vector<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::subresultantChainAtIdx (const SparseMultivariateRationalPolynomial& q, int idx) const {

//     if (isZero()) {
//     	std::vector<SparseMultivariateRationalPolynomial> zeroChain_p;
//     	zeroChain_p.push_back (q);
//     	return zeroChain_p;
//     }

//     if (q.isZero()) {
//     	std::vector<SparseMultivariateRationalPolynomial> zeroChain_q;
//     	zeroChain_q.push_back (*this);
//     	return zeroChain_q;
//     }

//     // bool TYPE = 0; // 0: SMZP Ducos Algorithm, 1: DUZP modular algorithm
//     // Now, these cases handle in C-size

//     std::cerr << "nvar = " << nvar << std::endl;
//     std::cerr << "q.nvar = " << q.nvar << std::endl;


//     Symbol v = q.leadingVariable();
//     // commented to support constant polynomials in algebraic extension
//     // if (v != leadingVariable()) {
//     //     std::cout << "BPAS: error, cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
//     //     exit(1);
//     // }

//     // // Subresultant from SMZP:
//     std::vector<int> xs;
//     bool isOrdered = isOrderedRing(q, xs);
//     SparseMultivariateRationalPolynomial ppP, ppQ;
//     std::vector<SparseMultivariateRationalPolynomial> subr;

//     if (!isOrdered) {
//         std::cout << "BPAS: error, trying to compute subresultant chain between Q[";
//         for (int i = 1; i <= nvar; ++i) {
//             std::cout << names[i];
//             if (i < nvar) { std::cout << ", "; }
//         }
//         std::cout << "] and Q[";
//         for (int i = 1; i <= q.nvar; ++i) {
//             std::cout << q.names[i];
//             if (i < q.nvar) { std::cout << ", "; }
//         }
//         std::cout << "]." << std::endl;
//         exit(1);
//     }

//     int superNvar = xs.size() / 2;
//     if (superNvar != nvar || superNvar != q.nvar) {

// //map indices to the expanded superset.
//         int varmap[nvar];
//         int qvarmap[q.nvar];
//         for (int i = 0; i < xs.size(); i += 2) {
//             if (xs[i] != 0) {
//                 varmap[xs[i]-1] = i/2;
//             }
//             if (xs[i+1] != 0) {
//                 qvarmap[xs[i+1]-1] = i/2;
//             }
//         }

// //create new combined names array
//         Symbol newnames[superNvar+1];
// //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
//         newnames[0] = this->names[0];
//         for (int i = 0; i < xs.size(); i += 2) {
//             if (xs[i] != 0) {
//                 newnames[(i/2) + 1] = names[xs[i]];
//             } else {
//                 newnames[(i/2) + 1] = q.names[xs[i+1]];
//             }
//         }

//         ppP = this->expandVariables(superNvar, newnames, varmap);
//         ppQ = q.expandVariables(superNvar, newnames, qvarmap);
// // ppP = ppP.primitivePart();
// // ppQ = ppQ.primitivePart();
//     } else {
//         ppP = *this;
//         ppQ = q;

// // ppP = this->primitivePart();
// // ppQ = q.primitivePart();
//     }

//     AltArr_t* Pq = primitivePart_AA (ppP.poly);
//     AltArr_t* Qq = primitivePart_AA (ppQ.poly);

//     AltArrZ_t* P = deepCopyPolynomial_AAZFromAA (Pq);
//     AltArrZ_t* Q = deepCopyPolynomial_AAZFromAA (Qq);

//     freePolynomial_AA (Pq);
//     freePolynomial_AA (Qq);

//     // AltArrZ_t* P = deepCopyPolynomial_AAZFromAA (ppP.poly);
//     // AltArrZ_t* Q = deepCopyPolynomial_AAZFromAA (ppQ.poly);

//     int lvarP = leadingVariable_AAZ (P);
//     int lvarQ = leadingVariable_AAZ (Q);
//     bool isShrinked = 0;

//     //TODO:
//     //if main variables are NOT same, even if one (or both) has degree 0 in that first variable.
//     // ppP.names[1] != ppQ.names[1];




//     // commented to support constant polynomials in algebraic extension
//     // if (lvarP != lvarQ){
//     //     std::cout << "BPAS Error: cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
//     //     exit(1);
//     // }

//     // if (lvarP > 0){
//     //     for (int i = 0; i < lvarP; ++i){
//     //         shrinkNumVarsAtIdx_AAZ (P, 0);
//     //         shrinkNumVarsAtIdx_AAZ (Q, 0);
//     //     }
//     //     isShrinked = 1;
//     // }

//     if (lvarP > 0) {
//         int varmap[superNvar];
//         int j = 0;

//         for (int i = 0; i < lvarP; i++) {
//             varmap[i] = -1;
//         }

//         for (int i = lvarP; i < superNvar; i++) {
//             varmap[i] = j;
//             j++;
//         }

//         shrinkAndReorderVars_AAZ (P, varmap, superNvar);
//         shrinkAndReorderVars_AAZ (Q, varmap, superNvar);
//         isShrinked = 1;
//     }

//     // if (leadingVariable_AAZ (P) != 0 || leadingVariable_AAZ (Q) != 0){
//     //     std::cout << "BPAS Error: shrink in subresultantChain does not work properly!" << std::endl;
//     //     exit(1);
//     // }

//     AltArrZ_t* SC_idx = NULL;
//     AltArrZ_t* SC_idx1 = NULL;
//     AltArr_t* tmp;

//     std::cerr << "P_Z->nvar = " << P->nvar << std::endl;
//     std::cerr << "Q_Z->nvar = " << Q->nvar << std::endl;

//     DucosSubresultantChainAtIdxZ (P, Q, idx, &SC_idx, &SC_idx1);

//     freePolynomial_AAZ (P); // free
//     freePolynomial_AAZ (Q); // free

//     subr.reserve (2);
//     // cerr << "size := " << size << std::endl;                        // TEST

//     tmp = deepCopyPolynomial_AAFromAAZ (SC_idx);
//     if (isShrinked && tmp != NULL && tmp->size != 0){
//         expandNumVarsLeft_AA (tmp, ppP.nvar);
//         if (tmp->nvar != ppP.nvar){
//             std::cout << "BPAS Error: expand in subresultantChain does not work properly!" << std::endl;
//             exit(1);
//         }
//     }
//     subr.push_back (SparseMultivariateRationalPolynomial (tmp, ppP.nvar, ppP.names));
//     tmp = deepCopyPolynomial_AAFromAAZ (SC_idx1);
//     if (isShrinked && tmp != NULL && tmp->size != 0){
//         expandNumVarsLeft_AA (tmp, ppP.nvar);
//         if (tmp->nvar != ppP.nvar){
//             std::cout << "BPAS Error: expand in subresultantChain does not work properly!" << std::endl;
//             exit(1);
//         }
//     }
//     subr.push_back (SparseMultivariateRationalPolynomial (tmp, ppP.nvar, ppP.names));

//     freePolynomial_AAZ (SC_idx);
//     freePolynomial_AAZ (SC_idx1);

//     return subr;
// }


/**
 * Extended Subresultant Chain
 * Return the list of subresultants with Besout Coefficients
 *
 * @param q: The other sparse univariate polynomial
 **/
std::vector<std::vector<SparseMultivariateRationalPolynomial>> SparseMultivariateRationalPolynomial::exSubresultantChain (const SparseMultivariateRationalPolynomial& q, const Symbol& v) const {

    // // Debug Print
    // std::cerr << "[SMQP+] exSubresultantChain mdeg(p)= " << this->mainDegree() << " mdeg(q)= " << q.mainDegree() << std::endl; 

    if (isZero()) {
        std::vector<std::vector<SparseMultivariateRationalPolynomial>> zeroChain_p;
        std::vector<SparseMultivariateRationalPolynomial> zeroElem_p;
        SparseMultivariateRationalPolynomial one_p, zero_p;
        zero_p.zero();
        one_p.one();
        zeroElem_p.push_back (q);      // q = 0*p + 1*q
        zeroElem_p.push_back (zero_p);
        zeroElem_p.push_back (one_p);
        zeroChain_p.push_back (zeroElem_p);
        return zeroChain_p;
    }

    if (q.isZero()) {
        std::vector<std::vector<SparseMultivariateRationalPolynomial>> zeroChain_q;
        std::vector<SparseMultivariateRationalPolynomial> zeroElem_q;
        SparseMultivariateRationalPolynomial one_q, zero_q;
        one_q.one();
        zero_q.zero();
        zeroElem_q.push_back (*this);      // p = 1*p + 0*q
        zeroElem_q.push_back (one_q);
        zeroElem_q.push_back (zero_q);
        zeroChain_q.push_back (zeroElem_q);
        return zeroChain_q;
    }

    if (this->degree(v) != q.degree(v)) {
        std::cout << "BPAS: error, cannot compute (extended) subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
        exit(1);
    }

    std::vector<Symbol> superRing;
    bool sameRing = orderPreservingSetUnion(this->ringVariables(), q.ringVariables(), superRing);

    AltArrZ_t* pZ = NULL;
    AltArrZ_t* qZ = NULL;
    int origNvar = superRing.size();
    int varMap[origNvar];
    int swappedIdx = 0;
    preparePolysForSRC(q, v, superRing, sameRing, varMap, swappedIdx, &pZ, &qZ);

    exgcds_t* SC;
    AltArr_t* tmp;
    AltArr_t* tmpA;
    AltArr_t* tmpB;
    int size = 0;

    AltArr_t* P = deepCopyPolynomial_AAFromAAZ(pZ);
    freePolynomial_AAZ(pZ);
    AltArr_t* Q = deepCopyPolynomial_AAFromAAZ(qZ);
    freePolynomial_AAZ(qZ);

    exDucosSubresultantChain_AA (P,Q, &SC, &size);

    freePolynomial_AA (P); // free
    freePolynomial_AA (Q); // free

    exgcds_t* cur = SC;
    std::vector<std::vector<SparseMultivariateRationalPolynomial>> subr;
    subr.reserve (size);
        // cerr << "size := " << size << std::endl;                        // TEST

    int swapMap[origNvar];
    if (swappedIdx != 0) {
        for (int i = 0; i < origNvar; ++i) {
            swapMap[i] = i;
        }
        swapMap[swappedIdx] = 0;
        swapMap[0] = swappedIdx;

        //reset superRing also
        Symbol tmpV = superRing[0];
        superRing[0] = superRing[swappedIdx];
        superRing[swappedIdx] = tmpV;
    }

    //despite the name, it is only used as the zero polynomial in the second loop if filled==1
    SparseMultivariateRationalPolynomial zero;
    zero.zero();
    zero.setRingVariables(superRing);

    for (int i = 0; cur != NULL && i < size; ++i){

        std::vector<SparseMultivariateRationalPolynomial> vElems;

        tmp = cur->r; // deepCopyPolynomial_AA (cur->r);
        tmpA = cur->a; // deepCopyPolynomial_AA (cur->a);
        tmpB = cur->b; // deepCopyPolynomial_AA (cur->b);

        reverseShrinkVariables_AA_inp(tmp, origNvar, varMap);
        if (swappedIdx != 0) {
            reorderVars_AA(tmp, swapMap, origNvar);
        }
        reverseShrinkVariables_AA_inp(tmpA, origNvar, varMap);
        if (swappedIdx != 0) {
            reorderVars_AA(tmpA, swapMap, origNvar);
        }
        reverseShrinkVariables_AA_inp(tmpB, origNvar, varMap);
        if (swappedIdx != 0) {
            reorderVars_AA(tmpB, swapMap, origNvar);
        }

        zero.poly = tmp;
        vElems.emplace_back(std::move(zero));
        zero.poly = tmpA;
        vElems.emplace_back(std::move(zero));
        zero.poly = tmpB;
        vElems.emplace_back(std::move(zero));

        subr.emplace_back(std::move(vElems));
        cur = cur->next;
    }

    freeExgcds_AA(SC);

    return subr;
}



SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::resultant (const SparseMultivariateRationalPolynomial& q, const Symbol& v) const
{
    // // Debug Print
    // std::cerr << "[SMQP+] resultant mdeg(p)= " << this->mainDegree() << " mdeg(q)= " << q.mainDegree() << std::endl; 
    
    if (this->degree(v) == 0 && q.degree(v) == 0) {
        std::cout << "BPAS: error, cannot compute subresultant chain if both polynomials do not contain variable v." << std::endl;

    }
    if (this->degree(v) == 0) {
    	return q;
    }
    if (q.degree(v) == 0) {
    	return *this;
    }

    std::vector<Symbol> superRing;
    bool sameRing = orderPreservingSetUnion(this->ringVariables(), q.ringVariables(), superRing);

    AltArrZ_t* pZ = NULL;
    AltArrZ_t* qZ = NULL;
    int origNvar = superRing.size();
    int varMap[origNvar];
    int swappedIdx = 0;
    preparePolysForSRC(q, v, superRing, sameRing, varMap, swappedIdx, &pZ, &qZ);

    AltArrZ_t* res = DucosResultantZ (pZ, qZ);

    int swapMap[origNvar];
    if (swappedIdx != 0) {
        for (int i = 0; i < origNvar; ++i) {
            swapMap[i] = i;
        }
        swapMap[swappedIdx] = 0;
        swapMap[0] = swappedIdx;

        //reset superRing also
        Symbol tmpV = superRing[0];
        superRing[0] = superRing[swappedIdx];
        superRing[swappedIdx] = tmpV;
    }


    //despite the name, it is only used as the zero polynomial in the second loop if filled==1
    SparseMultivariateRationalPolynomial tmpSMQP;
    tmpSMQP.zero();
    tmpSMQP.setRingVariables(superRing);

    reverseShrinkVariables_AAZ_inp(res, origNvar, varMap);
    if (swappedIdx != 0) {
        reorderVars_AAZ(res, swapMap, origNvar);
    }
    AltArr_t* tmp = deepCopyPolynomial_AAFromAAZ (res);
    freePolynomial_AAZ(res);
    tmpSMQP.poly = tmp;
    return tmpSMQP;

}

/**
 * Subresultant Chain GCD
 * Return the last non-zero subresultant of the current polynomial and the input polynomial if the resultant is zero and return 1 otherwise
 *
 * @param q: The other sparse univariate polynomial
 **/
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::subresultantGCD (const SparseMultivariateRationalPolynomial& q) const {
    Symbol v = q.leadingVariable();
    SparseMultivariateRationalPolynomial one;
    one.one();
    if (v != leadingVariable()) {
        return one;
        // std::cout << "BPAS: error, cannot compute subresultant gcd if leading variable of input is different from leading variable of the current object." << std::endl;
        // exit(1);
    }

    if (isConstant() || q.isConstant()) {
        return one;
    }

    if (this->numberOfRingVariables() == 1 || q.numberOfRingVariables() == 1) {
        return this->primitiveGCD(q);
    }
    // // //
    // GCD from C++ side:
    //startTimer(&start);
    std::vector<SparseMultivariateRationalPolynomial> src = subresultantChain(q);

    //  for (int i=0; i<src.size(); ++i)
    //      std::cerr << "src[" << i << "] = " << src[i] << std::endl;
    if (!src[0].isZero()) {
    	//      stopTimer(&start,&elapsed);
    	//      std::cerr << "SUP gcd time: " << elapsed << std::endl;
    	return one;
    }
    else {
    	if (src[2].degree(v) == src[1].degree(v)) {
    	    //          stopTimer(&start,&elapsed);
    	    //          std::cerr << "SUP gcd time: " << elapsed << std::endl;
    	    return src[2];
    	}
    	else {
    	    //          stopTimer(&start,&elapsed);
    	    //          std::cerr << "SUP gcd time: " << elapsed << std::endl;
    	    return src[1];
    	}
    }
    // // //

    // // // //
    // // GCD from C side:

    // SparseMultivariateRationalPolynomial ppP = this->primitivePart();
    // SparseMultivariateRationalPolynomial ppQ = q.primitivePart();
    // AltArrZ_t* P = deepCopyPolynomial_AAZFromAA (ppP.poly);
    // AltArrZ_t* Q = deepCopyPolynomial_AAZFromAA (ppQ.poly);
    // AltArr_t* gcd;
    // gcd = deepCopyPolynomial_AAFromAAZ (DucosGCDZ (P, Q));
    // return SparseMultivariateRationalPolynomial (gcd, nvar, names);
    // // // //
}

/**
 * Get the GCD between this and b.
 * If this and b have all integer coefficients, the gcd will have integer coefficients
 * with proper GCD among those coefficients. Otherwise, the returned GCD is monic
 * with rational number coefficients.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::gcd(const SparseMultivariateRationalPolynomial& b) const {

    bool TYPE = 0; // 0: gcd_AAZ, 1: primitiveGCD
    SparseMultivariateRationalPolynomial ret;

    if (TYPE == 0){
	// std::cout << "[Ali-TEST] computing SMQP::gcd ... " << std::endl;
	std::vector<int> xs;
	bool isOrdered = isOrderedRing(b, xs);
	SparseMultivariateRationalPolynomial ppP, ppQ;

	// // TEST
	// std::cout << "BPAS : *this := " << *this << std::endl;
	// std::cout << "BPAS : b := " << b << std::endl;

	if (!isOrdered) {
	    std::cout << "BPAS: error, trying to compute gcd  between Q[";
	    for (int i = 1; i <= nvar; ++i) {
		std::cout << names[i];
		if (i < nvar) { std::cout << ", "; }
	    }
	    std::cout << "] and Q[";
	    for (int i = 1; i <= b.nvar; ++i) {
		std::cout << b.names[i];
		if (i < b.nvar) { std::cout << ", "; }
	    }
	    std::cout << "]." << std::endl;
	    exit(1);
	}

	int superNvar = xs.size() / 2;
	if (superNvar != nvar || superNvar != b.nvar) {

	    //map indices to the expanded superset.
	    int varmap[nvar];
	    int bvarmap[b.nvar];
	    for (int i = 0; i < xs.size(); i += 2) {
		if (xs[i] != 0) {
		    varmap[xs[i]-1] = i/2;
		}
		if (xs[i+1] != 0) {
		    bvarmap[xs[i+1]-1] = i/2;
		}
	    }

	    //create new combined names array
	    Symbol newnames[superNvar+1];
	    //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
	    newnames[0] = this->names[0];
	    for (int i = 0; i < xs.size(); i += 2) {
		if (xs[i] != 0) {
		    newnames[(i/2) + 1] = names[xs[i]];
		} else {
		    newnames[(i/2) + 1] = b.names[xs[i+1]];
		}
	    }

	    ppP = this->expandVariables(superNvar, newnames, varmap);
	    ppQ = b.expandVariables(superNvar, newnames, bvarmap);
	    // ppP = ppP.primitivePart();
	    // ppQ = ppQ.primitivePart();
	} else {
	    ppP = *this;
	    ppQ = b;

	    // ppP = this->primitivePart();
	    // ppQ = b.primitivePart();
	}

	AltArr_t* Pq = primitivePart_AA (ppP.poly);
	AltArr_t* Qq = primitivePart_AA (ppQ.poly);

	AltArrZ_t* P = deepCopyPolynomial_AAZFromAA (Pq);
	AltArrZ_t* Q = deepCopyPolynomial_AAZFromAA (Qq);
    int Pnvar = P == NULL ? 0 : P->nvar;

	// TODO: uncomment :
	freePolynomial_AA (Pq);
	freePolynomial_AA (Qq);

	// std::cout << "\n\n\n\n";
	// std::cout << "ppP := " << ppP << std::endl;
	// std::cout << "ppQ := " << ppQ << std::endl;

	// AltArrZ_t* P = deepCopyPolynomial_AAZFromAA (ppP.poly);
	// AltArrZ_t* Q = deepCopyPolynomial_AAZFromAA (ppQ.poly);

	// std::cout << "[Ali-TEST] nvar(P) = " << P->nvar << std::endl;
	// std::cout << "[Ali-TEST] nvar(Q) = " << Q->nvar << std::endl;

	// SparseMultivariateRationalPolynomial test_P (deepCopyPolynomial_AA (Pq), ppP.nvar, ppP.names);
	// std::cout << "[Ali-TEST] P = "  << test_P << std::endl;

	// SparseMultivariateRationalPolynomial test_Q (deepCopyPolynomial_AA (Qq), ppP.nvar, ppP.names);
	// std::cout << "[Ali-TEST] Q = "  << test_Q << std::endl;


#if defined(WITH_BLAD) && WITH_BLAD

    AltArrZ_t* gz = NULL;
    if (Pnvar > 1) {
        if (P->size < 10 || Q->size < 10) {
            gz = gcd_AAZ (P,Q);
        } else {
            char* c_names[Pnvar];
            for (int i = 0; i < Pnvar; ++i) {
                std::string str = ppP.names[i+1].toString();
                c_names[i] = (char*) malloc(sizeof(char)*str.length()+1);
                strcpy(c_names[i], str.c_str());
            }
            gcdBLAD_AAZ(P, Q, (const char**) c_names, &gz);
            for (int i = 0; i < nvar; ++i) {
                free(c_names[i]);
            }
        }
    } else if (Pnvar == 1) {
        gz = univariateGCD_AAZ(P, Q);
    } else {
        gz = makePolynomial_AAZ(1, 0);
        mpz_init(gz->elems->coef);
        mpz_set_si(gz->elems->coef, 1l);
        gz->size = 1;
    }
#else
    AltArrZ_t* gz = gcd_AAZ (P, Q);
#endif

	// std::cout << "[Ali-TEST] gcd_AAZ from SMQP::gcd is done... " << std::endl;

	// std::cout << "[Ali-TEST] gcd_AAZ from SMQP::gcd is called... " << std::endl;

	AltArr_t* g = deepCopyPolynomial_AAFromAAZ (gz);
	freePolynomial_AAZ (P);
	freePolynomial_AAZ (Q);
	freePolynomial_AAZ (gz);

	SparseMultivariateRationalPolynomial gcd (g, ppP.nvar, ppP.names);
	ret = gcd;
	// gcd += ppP.content() * ppQ.content();
	// return SparseMultivariateRationalPolynomial (gcd);

    } else {
	ret = this->primitiveGCD(b);
    }

    mpz_t g1;
    mpz_t g2;
    mpz_init(g1);
    mpz_init(g2);
    integerPolynomialTestCont_AA(this->poly, g1);
    integerPolynomialTestCont_AA(b.poly, g2);
    if (mpz_cmp_si(g1, 0l) != 0 && mpz_cmp_si(g2, 0l) != 0) {
	//they are both integer polys
	mpz_gcd(g1, g1, g2);
	ret *= RationalNumber(g1);
    } else if (!ret.isConstant()) {
	ret /= ret.leadingCoefficient();
    }

    mpz_clear(g1);
    mpz_clear(g2);

    return ret;

}


/**
 * Get GCD between *this and b as a primitive polynomial.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::primitiveGCD(const SparseMultivariateRationalPolynomial& b) const {
    if (isZero()) {
        return b;
    }
    if (b.isZero()) {
        return *this;
    }

    if (isConstant() || b.isConstant()) {
        SparseMultivariateRationalPolynomial ret;
        ret.one();
        return ret;
    //     mpz_t mpzG;
    //     mpz_init(mpzG);
    //     mpz_set_ui(mpzG, 1l);
    //     if (b.isConstant() && mpz_cmp_si(mpq_denref(b.poly->elems->coef), 1l) == 0) {
    //         std::cerr << "b is constant and integer" << std::endl;
    //         mpz_set(mpzG, mpq_numref(b.poly->elems->coef));
    //         for (int i = 0; i < this->poly->size; ++i) {
    //             gmp_fprintf(stderr, "gcd: %Zd, this coef: %Qd\n", mpzG, this->poly->elems[i].coef);
    //             if (mpz_cmp_si(mpq_denref(this->poly->elems[i].coef), 1l) != 0) {
    //                 mpz_set_ui(mpzG, 1l);
    //                 break;
    //             }
    //             mpz_gcd(mpzG, mpzG, mpq_numref(this->poly->elems[i].coef));
    //             gmp_fprintf(stderr, "gcd: %Zd\n", mpzG);
    //             if (mpz_cmp_si(mpzG, 1l) == 0) {
    //                 break;
    //             }
    //         }
    //     } else if (mpz_cmp_si(mpq_denref(poly->elems->coef), 1l) == 0) {
    //         std::cerr << "this is constant and integer" << std::endl;
    //         mpz_set(mpzG, mpq_numref(this->poly->elems->coef));
    //         for (int i = 0; i < b.poly->size; ++i) {
    //             gmp_fprintf(stderr, "gcd: %Zd, b coef: %Qd\n", mpzG, b.poly->elems[i].coef);
    //             if (mpz_cmp_si(mpq_denref(b.poly->elems[i].coef), 1l) != 0) {
    //                 mpz_set_ui(mpzG, 1l);
    //                 break;
    //             }
    //             mpz_gcd(mpzG, mpzG, mpq_numref(b.poly->elems[i].coef));
    //             gmp_fprintf(stderr, "gcd: %Zd\n", mpzG);
    //             if (mpz_cmp_si(mpzG, 1l) == 0) {
    //                 break;
    //             }
    //         }
    //     }

    //     RationalNumber mpzc(mpzG);
    //     SparseMultivariateRationalPolynomial ret(mpzc);
    //     std::cerr << "mpzc: " << mpzc << " ret: " << ret << std::endl;
    //     mpz_clear(mpzG);
    //     return ret;
    }

    std::vector<int> xs;
    if (!isOrderedRing(b, xs)) {
        SparseMultivariateRationalPolynomial ret(0);
        ret.one();
        return ret;
    }

    if (this->nvar == 1 && b.nvar == 1) {
        if (xs.size() > 2) {
            SparseMultivariateRationalPolynomial ret(0);
            ret.one();
            return ret;
        }

        //otherwise they are the same variables
        AltArr_t* g = univariateGCD_AA(this->poly, b.poly);
        if (isConstant_AA(g)) {
            SparseMultivariateRationalPolynomial ret(g, 0, names);
            return ret;
        }
        SparseMultivariateRationalPolynomial ret(g, nvar, names);
        return ret;
    }
    if (this->nvar == 1) {
        std::vector<Symbol> vars = b.ringVariables();
        int match = -1;
        for (int i = 0; i < b.nvar; ++i) {
            if (vars[i] == names[1]) {
                match = i;
                break;
            }
        }
        if (match < 0) {
            SparseMultivariateRationalPolynomial ret(0);
            ret.one();
            return ret;
        }

        vars.erase(vars.begin()+match);

        SparseMultivariateRationalPolynomial newB = b.content(vars);

        //it is possible that the content that comes back is an integral content.
        if (newB.nvar == 0) {
            //recursive call to handle is isContstant() case above
            return this->primitiveGCD(newB);
        }

        AltArr_t* g = univariateGCD_AA(this->poly, newB.poly);
        if (isConstant_AA(g)) {
            SparseMultivariateRationalPolynomial ret(g, 0, names);
            return ret;
        }
        SparseMultivariateRationalPolynomial ret(g, nvar, names);
        return ret;
    } else if (b.nvar == 1) {
        std::vector<Symbol> vars = this->ringVariables();
        int match = -1;
        for (int i = 0; i < nvar; ++i) {
            if (vars[i] == b.names[1]) {
                match = i;
                break;
            }
        }
        if (match < 0) {
            SparseMultivariateRationalPolynomial ret(0);
            ret.one();
            return ret;
        }

        vars.erase(vars.begin()+match);

        SparseMultivariateRationalPolynomial newThis = this->content(vars);

        //it is possible that the content that comes back is an integral content.
        if (newThis.nvar == 0) {
            //recursive call to handle is isContstant() case above
            return b.primitiveGCD(newThis);
        }

        AltArr_t* g = univariateGCD_AA(newThis.poly, b.poly);
        if (isConstant_AA(g)) {
            SparseMultivariateRationalPolynomial ret(g, 0, names);
            return ret;
        }
        SparseMultivariateRationalPolynomial ret(g, b.nvar, b.names);
        return ret;
    }

    int superNvar = xs.size() / 2;
    if (superNvar != nvar || superNvar != b.nvar) {
        //map indices to the expanded superset.
        int varmap[nvar];
        int bvarmap[b.nvar];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                varmap[xs[i]-1] = i/2;
            }
            if (xs[i+1] != 0) {
                bvarmap[xs[i+1]-1] = i/2;
            }
        }

        //create new combined names array
        Symbol newnames[superNvar+1];
        //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
        newnames[0] = this->names[0];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                newnames[(i/2) + 1] = names[xs[i]];
            } else {
                newnames[(i/2) + 1] = b.names[xs[i+1]];
            }
        }

        SparseMultivariateRationalPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);


        std::vector<Symbol> vars;
        vars.push_back(newnames[1]);
        SparseMultivariateRationalPolynomial tC = tempa.content(vars);
        SparseMultivariateRationalPolynomial bC = tempb.content(vars);
        SparseMultivariateRationalPolynomial g = tC.primitiveGCD(bC);

        SparseMultivariateRationalPolynomial tPP = tempa / tC;
        SparseMultivariateRationalPolynomial bPP = tempb / bC;

        SparseMultivariateRationalPolynomial ppG = tPP.subresultantGCD(bPP);
        SparseMultivariateRationalPolynomial ppGP = ppG.primitivePart(vars);

        SparseMultivariateRationalPolynomial ret = g * ppGP;

        ret.setRingVariables(ret.variables());
        return ret;
    } else {
        //exact same ring, easy.
        std::vector<Symbol> vars;
        vars.push_back(names[1]);
        SparseMultivariateRationalPolynomial tC = this->content(vars);
        SparseMultivariateRationalPolynomial bC = b.content(vars);

        std::vector<Symbol> rVars;

        SparseMultivariateRationalPolynomial g = tC.primitiveGCD(bC);

        SparseMultivariateRationalPolynomial tPP = *this / tC;
        SparseMultivariateRationalPolynomial bPP = b / bC;

        SparseMultivariateRationalPolynomial ppG = tPP.subresultantGCD(bPP);
        SparseMultivariateRationalPolynomial ppGP = ppG.primitivePart(vars);

        SparseMultivariateRationalPolynomial ret = g * ppGP;

        ret.setRingVariables(ret.variables());
        return ret;
    }

}

Factors<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::squareFree() const {
    return this->squareFree(this->ringVariables());
}

Factors<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::squareFree(const std::vector<Symbol>& vars) const {
    Factors<SparseMultivariateRationalPolynomial> sf;
    if (isConstant() || isZero() || vars.size() == 0) {
        sf.setRingElement(*this);
        return sf;
    }

    std::vector<Symbol> nextVars;
    int strIdx = 0;
    for (int j = 0; j < vars.size(); ++j) {
        if (strIdx == 0) {
            for (int i = 0; i < nvar; ++i) {
                if (vars[j] == names[i+1]) {
                    strIdx = i+1;
                    break;
                }
            }
            if (strIdx == 0){
                nextVars.push_back(vars[j]);
            }
        } else {
            nextVars.push_back(vars[j]);
        }
    }

    if (strIdx == 0) {
        sf.setRingElement(*this);
        return sf;
    }


    std::vector<Symbol> curVars;
    curVars.push_back(names[strIdx]);
    SparseMultivariateRationalPolynomial content;
    SparseMultivariateRationalPolynomial primPart = this->primitivePart(curVars, content);

    // std::cout << "primPart = " << primPart << std::endl;

    if (primPart.isConstant()) {
        if (!primPart.isOne()) {
            sf.addFactor(primPart, 1);
        }
    } else if (primPart.degree(names[strIdx]) == 1) {
        // SparseMultivariateRationalPolynomial df = primPart.leadingCoefficient(names[strIdx]);
        // SparseMultivariateRationalPolynomial g = primPart.gcd(df);
        // SparseMultivariateRationalPolynomial next = primPart / g;
        // next *= g;
        // sf.push_back(next);

        sf.addFactor(primPart, 1);
    } else {
        SparseMultivariateRationalPolynomial dx = primPart.derivative(names[strIdx]);

    // std::cout << "dx = " << dx << std::endl;

        SparseMultivariateRationalPolynomial g = (primPart.gcd(dx)).primitivePart();

        SparseMultivariateRationalPolynomial next = primPart / g;
        if (next.leadingCoefficient() < 1) {
            next.negate();
        }

        int k = 1;
        if (g.isOne()) {
            AltArr_t* remFact = NULL;
            AltArr_t* comFact = commonFactor_AA(next.poly, &remFact);
            if (!isConstant_AA(comFact)) {
                //next /= comFact
                AltArr_t* temp = next.poly;
                next.poly = remFact;
                freePolynomial_AA(temp);

                //set comFact exp to 1 to add to sf.
                temp = makeConstPolynomial_AA(1, nvar, RationalNumber(1).get_mpq_t());
                degree_t tempDegs[nvar] = {0};
                degree_t deg = partialDegreeTerm_AA(comFact, 0, 0);
                SparseMultivariateRationalPolynomial fac(temp, nvar, names);
                if (deg != 0) {
                    tempDegs[0] = 1;
                    setDegrees_AA_inp(temp, 0, tempDegs, nvar);
                    sf.addFactor(fac, deg);
                }
                for (int j = 1; j < nvar; ++j){
                    deg = partialDegreeTerm_AA(comFact, 0, j);
                    if (deg == 0) {
                        continue;
                    }
                    tempDegs[j-1] = 0;
                    tempDegs[j] = 1;
                    setDegrees_AA_inp(temp, 0, tempDegs, nvar);
                    //updating temp automatically updates fac
                    sf.addFactor(fac, deg);
                }
            }
        }
        while (g.degree(names[strIdx]) > 0) {
        // std::cout << "next = " << next << std::endl;
        // std::cout << "g = " << g << std::endl;

            SparseMultivariateRationalPolynomial y = (next.gcd(g)).primitivePart();
            // std::cerr << "after primitiveGCD(g)" << std::endl;
            //don't add factors of 1.
            if (next != y) {
                sf.addFactor(next / y, k);
            }
            g /= y;
            next = y;
            if (next.leadingCoefficient() < 1) {
                next.negate();
            }
            ++k;
        }
        sf.addFactor(next, k);
    }

    //if content is constant then this recursive call will set the ringElement
    //to be that constant. Do this instead of setting the ring constant to
    //the integral content of this.
    // if (!content.isConstant()) {
        Factors<SparseMultivariateRationalPolynomial> contSf = content.squareFree(nextVars);
        sf.addFactors(contSf);
        sf.multiplyRingElement(contSf.ringElement());
    // }

    // //Calculate leading multiplier to return;
    // RationalNumber iCont = this->content();
    // sf.setRingElement(iCont);
    // std::cerr << "leaving squareFree(vs)..." << std::endl;
    return sf;
}

/**
 * Computes the square free part of this polynomail. That is, the polynomial of this
 * divided by all square factors. This is with respect to all variables.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::squareFreePart() const {

    int TYPE = 0; // 0: SMZP (C version), 1: SMQP (C++ version)

    if (isZero()) {
        return *this;
    }

    if (!TYPE) {
        mpq_t qCont;
        mpq_init(qCont);
        AltArrZ_t* aaz = primitivePartAndContent_AAZFromAA(this->poly, qCont);
        mpq_clear(qCont);
        // SparseMultivariateRationalPolynomial ppP = this->primitivePart();
		// AltArrZ_t* aaz = deepCopyPolynomial_AAZFromAA (ppP.poly);
        AltArrZ_t* sqrPart_AAZ = squareFreePart_AAZ (aaz, this->nvar);
        freePolynomial_AAZ (aaz);
        AltArr_t* sqrPart_AA = deepCopyPolynomial_AAFromAAZ (sqrPart_AAZ);
        freePolynomial_AAZ (sqrPart_AAZ);
        SparseMultivariateRationalPolynomial result (sqrPart_AA, this->nvar, this->names);
        return result;
    }

    std::vector<Symbol> vars = this->ringVariables();
    return this->squareFreePart(vars);
}

/**
 * Computes the square free part of this polynomail. That is, the polynomial of this
 * divided by all square factors. This is with respect to all variables.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::squareFreePart(std::vector<Symbol>& vars) const {
	if (isZero()) {
		return *this;
    }

    std::vector<Symbol> matchedVars;
    for (int i = 0; i < vars.size(); ++i) {
        for (int j = 0; j <= nvar; ++j) {
            if (vars[i] == names[j+1]) {
                matchedVars.push_back(vars[i]);
                break;
            }
        }
    }

	// std::cerr << "[SMQP-Ali-SFP] aa := " << *this << std::endl;

    SparseMultivariateRationalPolynomial cont;
    SparseMultivariateRationalPolynomial fact;
    fact.one();
    std::vector<Symbol> var;
    SparseMultivariateRationalPolynomial primpart = *this;
    for (int i = 0; i < matchedVars.size(); ++i) {
        var.clear();
        var.push_back(matchedVars[i]);

        primpart = primpart.primitivePart(var, cont);
        SparseMultivariateRationalPolynomial diff = primpart.derivative(var[0]);

		// std::cerr << "[SMQP-Ali] primpart["<<i<<"] := " << primpart << std::endl;
		// std::cerr << "[SMQP-Ali] diff["<<i<<"] := " << diff << std::endl;
		// std::cerr << "[SMQP-Ali] nvar := " << nvar << std::endl;

		SparseMultivariateRationalPolynomial g = (primpart.gcd(diff)).primitivePart();
        SparseMultivariateRationalPolynomial next = primpart / g;
        if (next.leadingCoefficient() < 1) {
            next.negate();
        }
        fact *= next;
        if (cont.isConstant()) {
            break;
        }
        primpart = cont;
    }

    return fact;
}

////////// Private Helpers ////////////////////////////////////

/**
 * Finds the variable superset which contains both variable orderings from *this and b, if possible.
 *
 * xs is filled such that it's length is twice that of the number of variables in
 * the resulting superset and it's indices are such that for variable i in [0, ..., size(xs)/2)
 * in the superset is equal to this.names[xs[2*i]] if xs[2*i] != 0 and equal to
 * b.names[xs[2*i+1]] if xs[2*i+1] != 0.
 * This means that if both xs[2*i] and xs[2*i+1] are non-zero, then the variable is shared
 * between *this and b.
 *
 * For example, [9,x,y] and [9,s,t,x] produces [0,1,0,2,1,3,2,0]. Recall leading "9" indicates
 * user-supplied ordering.
 *
 * returns true iff such an ordering is possible.
 */
bool SparseMultivariateRationalPolynomial::isOrderedRing(const SparseMultivariateRationalPolynomial& b, std::vector<int>& xs) const {
    // if (names[0] != b.names[0]) { return 0; }
    // if (names[0] == "1") {
    //     int k = 1;
    //     for (; k <= nvar && k <= b.nvar; ++k) {
    //         xs.push_back(k);
    //         xs.push_back(k);
    //     }
    //     for (int i = k; i <= nvar; ++i) {
    //         xs.push_back(i);
    //         xs.push_back(0);
    //     }
    //     for (int i = k; i <= b.nvar; ++i) {
    //         xs.push_back(0);
    //         xs.push_back(i);
    //     }
    //     return 1;
    // }
    if (!nvar) {
        for (int i = 1; i <= b.nvar; ++i) {
            xs.push_back(0);
            xs.push_back(i);
        }
        return 1;
    }
    if (!b.nvar) {
        for (int i = 1; i <= nvar; ++i) {
            xs.push_back(i);
            xs.push_back(0);
        }
        return 1;
    }

    bool isFound = 0;
    int* pos = new int[nvar];
    for (int i = 1; i <= nvar; ++i) {
        isFound = 0;
        for (int j = 1; j <= b.nvar; ++j) {
            if (names[i] == b.names[j]) {
                pos[i-1] = j;
                isFound = 1;
                break;
            }
        }
        if (!isFound) { pos[i-1] = 0; }
    }

    isFound = 0;
    int ak = 1, bk = pos[0];
    if (pos[0]) {
        for(int j = 1; j < pos[0]; ++j) {
            xs.push_back(0);
            xs.push_back(j);
        }
        xs.push_back(1);
        xs.push_back(pos[0]);
    }
    for (int i = 1; i < nvar; ++i) {
        if (pos[i] > bk) {
            while (pos[i] - bk > 1 && (i > ak || !bk)) {
                if (bk) { ak++; }
                bk++;
                if (names[ak] < b.names[bk]) {
                    xs.push_back(ak);
                    xs.push_back(0);
                    xs.push_back(0);
                    xs.push_back(bk);
                }
                else if (names[ak] > b.names[bk]) {
                    xs.push_back(0);
                    xs.push_back(bk);
                    xs.push_back(ak);
                    xs.push_back(0);
                }
            }
            for(int j = bk+1; j < pos[i]; ++j) {
                xs.push_back(0);
                xs.push_back(j);
            }
            for(int j = (!bk)? ak : ak+1; j <= i; ++j) {
                xs.push_back(j);
                xs.push_back(0);
            }
            xs.push_back(i+1);
            xs.push_back(pos[i]);
            bk = pos[i];
            ak = i + 1;
        }
        else if (pos[i] && pos[i] < bk) {
            isFound = 1;
            break;
        }
    }
    if (!isFound) {
        for (int i = (bk)? ak+1 : ak, j = bk+1; i <= nvar && j <= b.nvar; ++i, ++j) {
            if (names[i] < b.names[j]) {
                xs.push_back(i);
                xs.push_back(0);
                xs.push_back(0);
                xs.push_back(j);
            }
            else if (names[i] > b.names[j]) {
                xs.push_back(0);
                xs.push_back(j);
                xs.push_back(i);
                xs.push_back(0);
            }
            ak = i, bk = j;
        }
        for (int i = ak+1; i <= nvar; ++i) {
            xs.push_back(i);
            xs.push_back(0);
        }
        for (int j = bk+1; j <= b.nvar; ++j) {
            xs.push_back(0);
            xs.push_back(j);
        }
    }

    delete [] pos;

    if (isFound) { xs.clear(); return 0; }
    else { return 1; }
}

/**
 * Returns a copy of *this under the new variable ordering supplied.
 * varmap is such that this.names[i] = newvars[varmap[i]]
 * Returns an SMQP equal to *this but expended to newvars.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::expandVariables(int vars, Symbol* newvars, int varmap[]) const {
    if (isZero()) {
        SparseMultivariateRationalPolynomial temp(NULL, vars, newvars);
        return temp;
    }
    if (isConstant()) {
        AltArr_t* tempPoly = makeConstPolynomial_AA(1, vars, poly->elems->coef);
        SparseMultivariateRationalPolynomial temp(tempPoly, vars, newvars);
        return temp;
    }

    AltArr_t* tempPoly = deepCopyPolynomial_AA(poly);
    expandNumVars_AA(tempPoly, vars);
    reorderVars_AA(tempPoly, varmap, nvar);

    SparseMultivariateRationalPolynomial temp(tempPoly, vars, newvars);
    return temp;
}

void SparseMultivariateRationalPolynomial::expandVarsInPlace(int vars, Symbol* newvars, int varmap[]) {
    if (isZero() || isConstant()) {
        if (poly != NULL) {
            expandNumVars_AA(poly, vars);
        }
    } else {
        expandNumVars_AA(poly, vars);

        //here we use nvar and vars as # of vars as varmap is of size nvar.
        reorderVars_AA(poly, varmap, nvar);
    }

    nvar = vars;
    delete[] names;
    names = new Symbol[vars+1];
    std::copy(newvars, newvars+vars+1, names);
}

/**
 * Rearrange exponent vectors in place and then re-sort the polynomial.
 */
void SparseMultivariateRationalPolynomial::reorderVarsInPlace(int varmap[]) {
    if (isZero()) {
        return;
    }
    if (isConstant()) {
        return;
    }

    reorderVars_AA(poly, varmap, nvar);
    Symbol* newvars = new Symbol[nvar+1];
    newvars[0] = names[0];
    for (int i = 0; i < nvar; ++i) {
        newvars[varmap[i]+1] = names[i+1];
    }
    delete[] names;
    names = newvars;

}

////////// BPASPolynomial ////////////////////////////////////

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator= (const SparseMultivariateRationalPolynomial& b) {
    if (this != &b) {
        freePolynomial_AA(poly);
        poly = deepCopyPolynomial_AA(b.poly);

        nvar = b.nvar;
        delete[] names;
        names = new Symbol[nvar+1];
        std::copy(b.names, b.names+nvar+1, names);
        slp = b.slp;
    }
    return *this;
}

/**
 * Movement assignment: move b to be this polynomail.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator= (SparseMultivariateRationalPolynomial&& b) {
    if (this != &b) {
        freePolynomial_AA(poly);
        poly = b.poly;
        b.poly = NULL;
        nvar = b.nvar;
        b.nvar = 0;

        delete[] names;
        names = new Symbol[nvar+1];
        std::copy(b.names, b.names+nvar+1, names);

        slp = b.slp;
        b.slp.clear();
    }
    return *this;
}

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator= (const RationalNumber& r) {
    if (poly != NULL) {
        freePolynomial_AA(poly);
        poly = NULL;
    }

    slp.clear();
    poly = makeConstPolynomial_AA(1, nvar, r.get_mpq_t());

    return *this;
}


SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator+ (const SparseMultivariateRationalPolynomial& b) const {
    if (b.isZero()) {
        return *this;
    }
    if (isZero()) {
        return b;
    }
    if (this->isConstant() != 0) {
        return (b + this->poly->elems->coef);
    }
    if (b.isConstant() != 0) {
        return (*this + b.poly->elems->coef);
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);

    if (!isOrdered) {
        std::cout << "BPAS: error, trying to add between Q[";
        for (int i = 1; i <= nvar; ++i) {
            std::cout << names[i];
            if (i < nvar) { std::cout << ", "; }
        }
        std::cout << "] and Q[";
        for (int i = 1; i <= b.nvar; ++i) {
            std::cout << b.names[i];
            if (i < b.nvar) { std::cout << ", "; }
        }
        std::cout << "]." << std::endl;
        exit(1);
    }

    int superNvar = xs.size() / 2;
    if (superNvar != nvar || superNvar != b.nvar) {

        //map indices to the expanded superset.
        int varmap[nvar];
        int bvarmap[b.nvar];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                varmap[xs[i]-1] = i/2;
            }
            if (xs[i+1] != 0) {
                bvarmap[xs[i+1]-1] = i/2;
            }
        }

        //create new combined names array
        Symbol newnames[superNvar+1];
        //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
        newnames[0] = this->names[0];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                newnames[(i/2) + 1] = names[xs[i]];
            } else {
                newnames[(i/2) + 1] = b.names[xs[i+1]];
            }
        }

        SparseMultivariateRationalPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        AltArr_t* sum = addPolynomials_AA(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateRationalPolynomial ret(sum, superNvar, newnames);
        return ret;

    } else {
        AltArr_t* sum = addPolynomials_AA(poly, b.poly, nvar);
        SparseMultivariateRationalPolynomial ret(sum, nvar, names);
        return ret;
    }
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator+ (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret += *this;
    return ret;
}

// SparseMultivariateRationalPolynomial operator+ (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
//     SparseMultivariateRationalPolynomial ret = a;
//     ret += b;
//     return ret;
// }

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator+= (const SparseMultivariateRationalPolynomial& b) {
    if (nvar == b.nvar) {
//    	std::cerr << "nvar: " << nvar << std::endl;
//    	std::cerr << "this: " << *this << std::endl;
//    	std::cerr << "b   : " << b << std::endl;


        std::vector<int> xs;
        bool isOrdered = isOrderedRing(b, xs);
        int superNvar = xs.size() / 2;
        if (isOrdered && superNvar == nvar) {
            this->poly = addPolynomials_AA_inp(this->poly, b.poly, nvar);
            return *this;
        }
    }

    *this = (*this + b);
    return *this;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator- () const {
    SparseMultivariateRationalPolynomial temp = *this;
    temp.negate();
    return temp;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator- (const SparseMultivariateRationalPolynomial& b) const {
    if (b.isZero()) {
        return *this;
    }
    if (isZero()) {
        return -b;
    }
    if (this->isConstant() != 0) {
        SparseMultivariateRationalPolynomial negB = -b;
        return (negB + this->poly->elems->coef);
    }
    if (b.isConstant() != 0) {
        return (*this - b.poly->elems->coef);
    }


    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);

    if (!isOrdered) {
        std::cout << "BPAS: error, trying to subtract between Q[";
        for (int i = 1; i <= nvar; ++i) {
            std::cout << names[i];
            if (i < nvar) { std::cout << ", "; }
        }
        std::cout << "] and Q[";
        for (int i = 1; i <= b.nvar; ++i) {
            std::cout << b.names[i];
            if (i < b.nvar) { std::cout << ", "; }
        }
        std::cout << "]." << std::endl;
        exit(1);
    }

    int superNvar = xs.size() / 2;
    if (superNvar != nvar || superNvar != b.nvar) {

        //map indices to the expanded superset.
        int varmap[nvar];
        int bvarmap[b.nvar];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                varmap[xs[i]-1] = i/2;
            }
            if (xs[i+1] != 0) {
                bvarmap[xs[i+1]-1] = i/2;
            }
        }

        //create new combined names array
        Symbol newnames[superNvar+1];
        //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
        newnames[0] = this->names[0];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                newnames[(i/2) + 1] = names[xs[i]];
            } else {
                newnames[(i/2) + 1] = b.names[xs[i+1]];
            }
        }

        SparseMultivariateRationalPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        AltArr_t* sum = subPolynomials_AA(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateRationalPolynomial ret(sum, superNvar, newnames);
        return ret;

    } else {
        AltArr_t* sum = subPolynomials_AA(poly, b.poly, nvar);
        SparseMultivariateRationalPolynomial ret(sum, nvar, names);
        return ret;
    }
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator- (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret.negate();
    ret += *this;
    return ret;
}

// SparseMultivariateRationalPolynomial operator- (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
//     SparseMultivariateRationalPolynomial ret = a;
//     ret -= b;
//     return ret;
// }

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator-= (const SparseMultivariateRationalPolynomial& b) {
    if (nvar == b.nvar) {
        std::vector<int> xs;
        bool isOrdered = isOrderedRing(b, xs);
        int superNvar = xs.size() / 2;
        if (isOrdered && superNvar == nvar) {
            this->poly = subPolynomials_AA_inp(this->poly, b.poly, nvar);
            return *this;
        }
    }

    *this = (*this - b);
    return *this;
}

/**
 * Multiply *this by the specified polynomail.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator* (const SparseMultivariateRationalPolynomial& b) const {
    if (b.isZero() || this->isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }
    if (this->isConstant() != 0) {
        return (b * this->poly->elems->coef);
    }
    if (b.isConstant() != 0) {
        return (*this * b.poly->elems->coef);
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);

    if (!isOrdered) {
        std::cout << "BPAS: error, trying to multiply between Q[";
        for (int i = 1; i <= nvar; ++i) {
            std::cout << names[i];
            if (i < nvar) { std::cout << ", "; }
        }
        std::cout << "] and Q[";
        for (int i = 1; i <= b.nvar; ++i) {
            std::cout << b.names[i];
            if (i < b.nvar) { std::cout << ", "; }
        }
        std::cout << "]." << std::endl;
        exit(1);
    }

    int superNvar = xs.size() / 2;
    if (superNvar != nvar || superNvar != b.nvar) {

        //map indices to the expanded superset.
        int varmap[nvar];
        int bvarmap[b.nvar];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                varmap[xs[i]-1] = i/2;
            }
            if (xs[i+1] != 0) {
                bvarmap[xs[i+1]-1] = i/2;
            }
        }

         //create new combined names array
        Symbol newnames[superNvar+1];
        //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
        newnames[0] = this->names[0];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                newnames[(i/2) + 1] = names[xs[i]];
            } else {
                newnames[(i/2) + 1] = b.names[xs[i+1]];
            }
        }

        SparseMultivariateRationalPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);
        AltArr_t* prod = multiplyPolynomials_AA(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateRationalPolynomial ret(prod, superNvar, newnames);
        return ret;

    } else {
        AltArr_t* prod = multiplyPolynomials_AA(poly, b.poly, nvar);
        SparseMultivariateRationalPolynomial ret(prod, nvar, names);
        return ret;
    }
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator* (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret *= *this;
    return ret;
}

// SparseMultivariateRationalPolynomial operator* (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
//     SparseMultivariateRationalPolynomial ret = a;
//     ret *= b;
//     return ret;
// }

/**
 * Update this by multiplying by the specified polynomail.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator*= (const SparseMultivariateRationalPolynomial& b) {
    *this = (*this * b);
    return *this;
}

/**
 * Divide *this by the specified polynomial.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator/ (const SparseMultivariateRationalPolynomial& b) const {
    std::vector<Symbol> vars = this->ringVariables();
    // std::cerr << "in operator/" << std::endl;
    // std::cerr << "t nvar: " << this->nvar << std::endl;
    // std::cerr << "this v: " << this->names[0] << std::endl;
    // for (auto v : vars) {
    //     std::cerr << "this v: " << v << std::endl;
    // }
    // vars = b.ringVariables();
    // std::cerr << "b nvar: " << b.nvar << std::endl;
    // std::cerr << "b    v: " << b.names[0] << std::endl;
    // for (auto v : vars) {
    //     std::cerr << "b    v: " << v << std::endl;
    // }

    SparseMultivariateRationalPolynomial q,r;
    // std::cerr << "calling divide" << std::endl;

    this->divide(b, q, r);
    if (!r.isZero()) {
        std::cerr << "BPAS ERROR: SMQP non-exact division." << std::endl;
        std::cerr << "dividend: " << *this << std::endl;
        std::cerr << "divisor: " << b << std::endl;
        std::cerr << "quoteint: " << q << std::endl;
        std::cerr << "remainder: " << r << std::endl;
        exit(1);
    }
    // std::cerr << "quoteint: " << q << std::endl;
    // std::cerr << "remainder: " << r << std::endl;

    return q;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator/ (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret /= *this;
    return ret;
}

// SparseMultivariateRationalPolynomial operator/ (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
//     SparseMultivariateRationalPolynomial ret = a;
//     ret /= b;
//     return ret;
// }

/**
 * Update *this by dividing by the specified polynomial.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator/= (const SparseMultivariateRationalPolynomial& b){
    *this = (*this / b);
    return *this;
}

/**
 * Exponentiate *this by the input exponent integer.
 * Treats negative exponents as positive.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator^ (long long int e) const {
    if (e == 0) {
        SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
        ret.one();
        return ret;
    }

    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    if (e == 1) {
        return SparseMultivariateRationalPolynomial(*this);
    }

    e = (e < 0) ? -e : e;
    AltArr_t* retPoly = exponentiatePoly_AA(poly, e, nvar);
    SparseMultivariateRationalPolynomial ret(retPoly, nvar, names);
    return ret;
}

/**
 * Update *this by exponentiating this to the input integer.
 * Treats negative exponents as positive.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator^= (long long int e) {
    *this = (*this ^ e);
    return *this;
}

/**
 * Determine if *this is equal to the specified polynomial.
 * This takes into account the variable ordering on both poylnomials
 * in such a way that the same polynomial under different variable orderings
 * are NOT equal.
 */
bool SparseMultivariateRationalPolynomial::operator== (const SparseMultivariateRationalPolynomial& b) const {
    return this->isEqual(b);
}

/**
 * Determine if *this is not equal to the specified polynomial.
 * This takes into account the variable ordering on both poylnomials
 * in such a way that the same polynomial under different variable orderings
 * are NOT equal.
 */
bool SparseMultivariateRationalPolynomial::operator!= (const SparseMultivariateRationalPolynomial& b) const {
    return !this->isEqual(b);
}

/**
 * Output the string representation of *this to the input ostream.
 */
void SparseMultivariateRationalPolynomial::print(std::ostream& os) const {
    if (this->poly != NULL) {
        std::string tempVars[this->nvar];
        for (int i = 0 ; i < this->nvar; ++i) {
            tempVars[i] = this->names[i+1].toString();
        }
        // if (this->poly->unpacked) {
        //     os << "(unpacked) ";
        // } else {
        //     os << "( packed ) ";
        // }
        os << polyToString_AA(poly, tempVars);
    } else {
        os << "0";
    }
}

std::istream& operator>>(std::istream& in, SparseMultivariateRationalPolynomial& p) {

    std::string str;
    getline(in, str);

    p.fromString(str);
    return in;
}


void SparseMultivariateRationalPolynomial::fromString(const std::string& str) {

    altarr_pack* pack = generate_altarr_pack(str.c_str());

    freePolynomial_AA(this->poly);
    this->poly = pack->altarr_t_data;

    delete[] this->names;
    this->names = new Symbol[pack->numVars + 1];
    for (int i = 0; i < pack->numVars; ++i) {
        names[i+1] = Symbol(std::string(pack->vars[i]));
        free(pack->vars[i]);
    }

    this->nvar = pack->numVars;
    free(pack->vars);
    pack->numVars = 0;
    free(pack);
}


RationalNumber SparseMultivariateRationalPolynomial::content() const {
    if (isZero()) {
        return RationalNumber(0);
    }
    if (isConstant()) {
        RationalNumber ret(poly->elems->coef);
        return ret;
    }

    mpq_t ret;
    mpq_init(ret);
    integralContent_AA(poly, ret);
    RationalNumber rn(ret);
    mpq_clear(ret);
    return rn;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::content(const std::vector<Symbol>& v) const {
    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    bool allFound = 1;
    std::vector<Symbol> matchingSyms;
    std::vector<Symbol> actualSyms = this->variables();
    for (int i = 0; i < actualSyms.size(); ++i) {
        bool found = 0;
        for (int j = 0; j < v.size(); ++j) {
            if (actualSyms[i] == v[j]) {
                found = 1;
                matchingSyms.push_back(v[j]);
                break;
            }
        }
        if (!found) {
            allFound = 0;
        }
    }

    // std::cerr << "\n\nGetting content: " << std::endl;
    // for (auto v : matchingSyms) {
    //     std::cerr << "matching v:" << v << std::endl;
    // }
    // for (int i = 0; i < nvar; ++i) {
    //     std::cerr << "this vars: " << names[i+1] << std::endl;
    // }
    // std::cerr << "all found: " << allFound << std::endl;

    if (allFound) {
        return this->content();
    }

    if (matchingSyms.size() == 0) {
        return *this;
    }


    SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> sup = this->convertToSUP(matchingSyms[0]);
    if (sup.numberOfTerms() <= 1 && matchingSyms.size() == 1) {
        //in this case, don't multiply by numerical content.
        //it is sufficient to return just the (recursive) coefficient in this
        //case where the sup has only one term;
        SparseMultivariateRationalPolynomial ret = sup.leadingCoefficient();
        return ret;
    }

    SparseMultivariateRationalPolynomial content = sup.content().primitivePart();
    if (content.isOne()) {
        return this->content();
    }

    for (int i = 1; i < matchingSyms.size(); ++i) {
        sup = this->convertToSUP(matchingSyms[i]);
        SparseMultivariateRationalPolynomial temp = sup.content().primitivePart();
        content = content.gcd(temp);
        if (content.isOne()) {
            return this->content();
        }
    }

    content *= this->content();
    return content;
}


SparseMultivariateIntegerPolynomial SparseMultivariateRationalPolynomial::primitivePartSMZP() const {
    mpq_t content;
    mpq_init(content);
    AltArrZ_t* aaz = primitivePartAndContent_AAZFromAA(poly, content);
    mpq_clear(content);

    return SparseMultivariateIntegerPolynomial(aaz, nvar, names);
}

SparseMultivariateIntegerPolynomial SparseMultivariateRationalPolynomial::primitivePartSMZP(RationalNumber& content) const {
    mpq_t cont;
    mpq_init(cont);
    AltArrZ_t* aaz = primitivePartAndContent_AAZFromAA(poly, cont);
    content = RationalNumber(cont);
    mpq_clear(cont);

    return SparseMultivariateIntegerPolynomial(aaz, nvar, names);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::primitivePart() const {
    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    // std::cout << "[Ali-TEST] computing primitivePart_AA... " << std::endl;
    AltArr_t* pp = primitivePart_AA(this->poly);
    // std::cout << "[Ali-TEST] primitivePart_AA is done... " << std::endl;


    return SparseMultivariateRationalPolynomial(pp, nvar, names);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::primitivePart(RationalNumber& content) const {
    if (isZero()) {
        content.zero();
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    mpq_t cont;
    mpq_init(cont);
    AltArr_t* pp = primitivePartAndContent_AA(this->poly, cont);
    content = RationalNumber(cont);
    mpq_clear(cont);

    return SparseMultivariateRationalPolynomial(pp, nvar, names);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::primitivePart(const Symbol& s) const {
    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

	std::vector<Symbol> v;
	v.push_back(s);
    return this->primitivePart(v);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::primitivePart(const std::vector<Symbol>& v) const {
    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    SparseMultivariateRationalPolynomial cont = this->content(v);
    return (*this / cont);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::primitivePart(const std::vector<Symbol>& v, SparseMultivariateRationalPolynomial& content) const {
    if (isZero()) {
        content.zero();
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    // std::cout << "[Ali-TEST] computing content(v) in PrimitivePart ... " << std::endl;
    content = this->content(v);
    // std::cout << "[Ali-TEST] content(v) in PrimitivePart is done! ... " << std::endl;

    return (*this / content);
}


SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::mainPrimitivePart() const {
	SparseMultivariateRationalPolynomial content;
	return mainPrimitivePart(content);
}


SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::mainPrimitivePart(SparseMultivariateRationalPolynomial& content) const {
    if (isZero()) {
        content.one();
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }



    mpq_t qcont;
    mpq_init(qcont);
    AltArrZ_t* pPoly = primitivePartAndContent_AAZFromAA(this->poly, qcont);
    AltArrZ_t* cont = NULL;
    AltArrZ_t* res = mainPrimitiveFactorization_AAZ (pPoly, &cont);

    AltArr_t* retCont = deepCopyPolynomial_AAFromAAZ(cont);
    multiplyByRational_AA_inp(retCont, qcont);
    mpq_clear(qcont);

    SparseMultivariateRationalPolynomial ref (retCont, nvar, names);

    freePolynomial_AAZ (cont);
    freePolynomial_AAZ (pPoly);

    content = ref;
    return SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ (res), nvar, names);
}



////////// BPASMultivariatePolynomial ////////////////////////////////////

/**
 * Get the number of variables in this polynomial.
 */
int SparseMultivariateRationalPolynomial::numberOfVariables() const {
    int foundVars[nvar] = {0};
    nonZeroVariables_AA(poly, foundVars);
    int res = 0;
    for (int i = 0; i < nvar; ++i) {
        if (foundVars[i]) {
            ++res;
        }
    }
    return res;
}


/**
 * Get the number of non-zero terms
 */
Integer SparseMultivariateRationalPolynomial::numberOfTerms() const {
    if (poly != NULL) {
        return poly->size;
    }
    return 0;
}

/**
 * Total degree.
 */
Integer SparseMultivariateRationalPolynomial::degree() const {
    if (isZero()) {
        return -1;
    }
    return totalDegree_AA(poly);
}

/**
 * Get the degree of a variable
 */
Integer SparseMultivariateRationalPolynomial::degree(const Symbol& str) const {
    if (isZero()) {
        return -1;
    }

    int strIdx = -1;
    for (int i = 1; i <= nvar; ++i) {
        if (names[i] == str) {
            strIdx = i-1;
            break;
        }
    }
    if (strIdx == -1) {
        return 0;
    }

    if (strIdx == 0) {
        //then str is the leading variable and we can just take the first
        //as we assume lexicographical order.
        return partialDegreeTerm_AA(poly, 0, 0);
    } else {
        return partialDegree_AA(poly, strIdx);
    }
}

/**
 * Get the leading coefficient
 */
RationalNumber SparseMultivariateRationalPolynomial::leadingCoefficient() const {
    if (isZero()) {
        return RationalNumber(0);
    }

    return RationalNumber(poly->elems->coef);
}

RationalNumber SparseMultivariateRationalPolynomial::trailingCoefficient() const {
    if (isZero()) {
        return RationalNumber(0);
    }

    return RationalNumber(poly->elems[poly->size-1].coef);
}


/**
 * Get a coefficient, given the exponent of each variable in d.
 * v is the number of variables in d. It is assumed that the first this.nvar
 * variables of d match the variables of this polynomial
 */
RationalNumber SparseMultivariateRationalPolynomial::coefficient(int v, const int* d) const {
    if (v < nvar) {
        std::cerr << "BPAS ERROR: SMQP calling coefficient without enough variables." << std::endl;
        exit(1);
    }

    RationalNumber ret;
    coefficient_AA(poly, d, v, ret.get_mpq_t());
    return ret;
}

/**
 * Set a coefficient, given the exponent of each variable
 */
void SparseMultivariateRationalPolynomial::setCoefficient(int v, const int* d, const RationalNumber& rn) {
    if (v != nvar) {
        std::cout << "BPAS: error, SMQP(" << nvar << "), but trying to setCoefficient with " << v << " variables." << std::endl;
        exit(1);
    }

    if (poly == NULL) {
        poly = makeConstPolynomial_AA(1, nvar, rn.get_mpq_t());
        if (nvar != 0) {
            setDegrees_AA_inp(poly, 0, d, nvar);
        }
        return;
    }

    setCoefficient_AA(poly, d, v, rn.get_mpq_t());
}

/**
 * Set variables' name.
 */
void SparseMultivariateRationalPolynomial::setRingVariables (const std::vector<Symbol>& xs) {
    int ns = xs.size();

    //can easily expand
    if (nvar == 0) {
        if (ns > 0) {
            delete[] names;
            names = new Symbol[ns+1];
            names[0] = "9";
            std::copy(xs.begin(), xs.end(), names+1);
            nvar = ns;
            if (poly !=NULL) {
                expandNumVars_AA(poly, ns);
            }
        }
        return;
    }

    if (ns == 0 && !(isZero() || isConstant())) {
        std::cerr << "BPAS ERROR: SMQP: Trying to remove all variables from a non-constant polynomial." << std::endl;
        exit(1);
    }

    //check for duplicate variable names.
    for (int i = 0; i < ns; ++i) {
        for (int j = i+1; j < ns; ++j) {
            if (xs[i] == xs[j]) {
                std::cerr << "BPAS ERROR: SMQP: Duplicate variable name in input variables for setVariableNames: " << xs[i] << " " << xs[j] << std::endl;
                exit(1);
            }
        }
    }

    if (ns >= nvar) {
        //just a reorder or rename of vars.

        Symbol newnames[ns+1];
        newnames[0] = "9";
        std::copy(xs.begin(), xs.end(), newnames+1);

        bool onlyRename = (ns == nvar);
        bool foundArr[nvar];
        bool usedXs[ns];
        int varmap[nvar];
        for (int i = 0; i < nvar; ++i) {
            foundArr[i] = 0;
            varmap[i] = -1;
        }
        for (int i = 0; i < ns; ++i) {
            usedXs[i] = 0;
        }

        for (int i = 0; i < nvar; ++i) {
            bool found = 0;
            for (int j = 0; j < ns; ++j) {
                if (names[i+1] == xs[j]) {
                    found = 1;
                    varmap[i] = j;
                    found = 1;
                    foundArr[i] = 1;
                    usedXs[j] = 1;
                    // newnames[j+1] = names[i+1];
                    break;
                }
            }
            if (found) {
                onlyRename = 0;
            }
        }

        if (onlyRename) {
            names[0] = "9";
            for (int i = 0; i < nvar; ++i) {
                names[i+1] = xs[i];
            }
            return;
        }

        //at this point, foundArr holds all variable in both this and xs.
        //we need to now pair up our vars with those in xs to rename the rest.
        for (int i = 0; i < nvar; ++i) {
            if (foundArr[i]) {
                continue;
            }
            for (int j = 0; j < ns; ++j) {
                if (!usedXs[j]) {
                    varmap[i] = j;
                    usedXs[j] = 1;
                    // newnames[j+1] = xs[j];
                    break;
                }
            }
        }

        this->expandVarsInPlace(ns, newnames, varmap);
        return;

    }

    //otherwise, we are decreasing the variables of the ring
    //we must check that this is a valid operation.
    std::vector<Symbol> nonZeros = this->variables();
    bool foundArr[nvar];
    bool usedXs[nvar];
    int varmap[nvar];
    for (int i = 0; i < nvar; ++i) {
        foundArr[i] = 0;
        varmap[i] = -1;
    }
    for (int i = 0; i < ns; ++i) {
        usedXs[i] = 0;
    }


    //search first for variables that much in this.names and xs.
    for (int i = 0; i < nvar; ++i) {
        for (int j = 0; j < ns; ++j) {
            if (names[i+1] == xs[j]) {
                varmap[i] = j;
                foundArr[i] = 1;
                usedXs[j] = 1;
                break;
            }
        }
    }

    //now assign, right to left, renames from this to xs.
    for (int i = 0; i < nvar; ++i) {
        if (!foundArr[i]) {
            bool set = 0;
            for (int j = 0; j < ns; ++j) {
                if (!usedXs[j]) {
                    varmap[i] = j;
                    usedXs[j] = 1;

                    foundArr[i] = 1;
                    break;
                }
            }
        }
    }

    //after the previous two loops we now know which variables we are to remove.
    for (int i = 0; i < nvar; ++i) {
        if (!foundArr[i]) {
            for (int j = 0; j < nonZeros.size(); ++j) {
                if (names[i+1] == nonZeros[j]) {
                    std::cerr << "BPAS ERROR: SMQP trying to remove variable " << names[i+1] << " which is non-zero in this polynomial." << std::endl;
                    exit(1);
                }
            }
        }
    }

    shrinkAndReorderVars_AA(poly, varmap, nvar);
    nvar = ns;
    delete[] names;
    names = new Symbol[nvar+1];
    names[0] = "9";
    std::copy(xs.begin(), xs.end(), names+1);
}

/**
 * Get variable names of variables with non-zero degree;
 */
std::vector<Symbol> SparseMultivariateRationalPolynomial::variables() const {
    if (nvar == 0 || isConstant()) {
        return std::vector<Symbol>();
    }

    int foundVar[nvar] = {0};
    nonZeroVariables_AA(poly, foundVar);

    std::vector<Symbol> retSyms;
    for (int i = 0; i < nvar; ++i) {
        if (foundVar[i]) {
            retSyms.push_back(names[i+1]);
        }
    }

    return retSyms;
}

/**
 * Get variable names of all variables available to this polynomial,
 * even those that have zero degree.
 */
std::vector<Symbol> SparseMultivariateRationalPolynomial::ringVariables() const {
    std::vector<Symbol> varNames;
    for (int i = 0; i < nvar; ++i) {
        varNames.push_back(names[i+1]);
    }
    return varNames;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::derivative(const Symbol& s, int k) const {
    if (k <= 0) {
        return *this;
    }

    if (isZero() || isConstant()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);;
    }

    int idx = -1;
    for (int i = 0; i < nvar; ++i) {
        if (s == names[i+1]) {
            idx = i;
            break;
        }
    }

    if (idx == -1) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    AltArr_t* temp = derivative_AA(poly, idx, k);
    return SparseMultivariateRationalPolynomial(temp, nvar, names);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::integral(const Symbol& s, int k) const {
    if (k <= 0) {
        return *this;
    }

    int idx = -1;
    for (int i = 0; i < nvar; ++i) {
        if (s == names[i+1]) {
            idx = i;
            break;
        }
    }

    if (idx == -1) {
        Symbol newnames[nvar+2];
        newnames[0] = "9";
        newnames[1] = s;
        for (int j = 1; j < nvar+1; ++j) {
            newnames[j+1] = names[j];
        }
        AltArr_t* temp = integral_AA(poly, idx, k);
        return SparseMultivariateRationalPolynomial(temp, nvar+1, newnames);
    }

    AltArr_t* temp = integral_AA(poly, idx, k);
    return SparseMultivariateRationalPolynomial(temp, nvar, names);
}


////////// RecursivelyViewedPolynomial  ////////////////////////////////////

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::initial() const {
    if (nvar == 0 || isConstant() || isZero()) {
        return *this;
    }
    return leadingCoefficientInVariable(leadingVariable());
}

int SparseMultivariateRationalPolynomial::mainDegree() const {
    if (isZero()) {
        return -1;
    }
    if (isConstant()) {
        return 0;
    }

    return mainDegree_AA(poly);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::rank() const {
    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    if (isConstant()) {
        SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
        ret.one();

        return ret;
    }

    int lvIdx = mainVariable_AA(poly);
    Symbol lv = names[lvIdx + 1];

    SparseMultivariateRationalPolynomial ret(lv);
    degree_t deg = this->mainDegree();
    setDegrees_AA_inp(ret.poly, 0, &deg, 1);

    return ret;
}


SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::head() const {
    if (nvar == 0 || isZero() || isConstant()) {
        return *this;
    }

    SparseMultivariateRationalPolynomial tmp = leadingCoefficientInVariable(leadingVariable()) * this->rank();
    tmp.setRingVariables(this->ringVariables());
    return tmp;
}


SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::tail() const {
    if (nvar == 0 || isZero() || isConstant()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    return *this - head();
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::separant() const {
    std::cerr << "BPAS ERROR: SparseMultivariateRationalPolynomial::separant NOT YET IMPLEMENTED" << std::endl;
    return *this;
}



////////// SMQP-Specific ////////////////////////////////////

bool SparseMultivariateRationalPolynomial::isEqual(const SparseMultivariateRationalPolynomial& b) const {
    if (isZero()) {
        if (b.isZero()) {
            return 1;
        } else {
            return 0;
        }
    }
    if (b.isZero()) {
        return 0;
    }
    if (isConstant()) {
        if (b.isConstant()) {
            return mpq_cmp(poly->elems->coef, b.poly->elems->coef) == 0;
        } else {
            return 0;
        }
    }
    if (poly->size != b.poly->size) {
        return 0;
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);
    if (!isOrdered) { return 0; }

    int ret = isEqualWithVariableOrdering_AA(poly, b.poly, xs.data(), xs.size());

    return ret;
}

/**
 * Evaluate *this polynomial given the variables and their values in vars, and values.
 * vars must be a list of variables which are a (not necessarily proper) subset.
 * vars and values should have matching indices. i.e. the values of vars[i] is values[i].
 *
 * returns a new SMQP where all variables in vars have been evaluated using the values given.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::evaluate(const std::vector<Symbol>& vars, const std::vector<RationalNumber>& values) const {

    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }
    if (isConstant()) {
        return SparseMultivariateRationalPolynomial(*this);
    }

    if (vars.size() != values.size()) {
        std::cerr << "BPAS: SMQP: ERROR number of variables does not much number of points in evaluate." << std::endl;
        exit(1);
    }

    if (nvar == 1) {
        if (names[1] != vars[0]) {
            std::cerr << "BPAS: SMQP: error: trying to evaluate variable that does not exist in this SMQP: " << vars[0] << std::endl;
            exit(1);
        }

        mpq_t res;
        mpq_init(res);
        univarEvaluate_AA(poly, values[0].get_mpq_t(), res);
        RationalNumber r(res);
        SparseMultivariateRationalPolynomial ret(r);
        mpq_clear(res);
        return ret;
    }

    int active[nvar];
    mpq_t vals[nvar];
    for (int i = 0; i < nvar; ++i) {
        active[i] = 0;
        mpq_init(vals[i]);
    }

    int newNvar = nvar;
    for (int i = 0; i < vars.size(); ++i) {
        bool found = 0;
        for (int j = 0; j < nvar; ++j) {
            if (names[j+1] == vars[i]) {
                found = 1;

                --newNvar;
                active[j] = 1;
                mpq_set(vals[j], values[i].get_mpq_ref().get_mpq_t());
            }
        }
        if (!found) {
            std::cerr << "BPAS: SMQP: error: trying to evaluate variable that does not exist in this SMQP: " << vars[i] << std::endl;
            exit(1);
        }
    }

    Symbol newnames[newNvar+1];
    if (newNvar == 0) {
        newnames[0] = "1";
    } else {
        newnames[0] = names[0];
        int idx = 1;
        for(int i = 0; i < nvar; ++i) {
            if (!active[i]) {
                newnames[idx] = names[i+1];
                ++idx;
            }
        }
    }

    AltArr_t* evalPoly = evaluatePoly_AA(poly, active, vals, nvar);
    return SparseMultivariateRationalPolynomial(evalPoly, newNvar, newnames);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::interpolate(const std::vector<std::vector<RationalNumber>>& points, const std::vector<RationalNumber>& vals) {
    if (points.size() != vals.size()) {
        std::cerr << "BPAS ERROR: SMQP: In interpolation, number of points does not match number of values." << std::endl;
        exit(1);
    }
    if (points.size() == 0) {
        return 0;
    }

    if (points[0].size() == 1) {
        mpq_t mpqPoints[points.size()];
        mpq_t mpqVals[points.size()];
        for (int i = 0; i < points.size(); ++i) {
            mpq_init(mpqPoints[i]);
            mpq_init(mpqVals[i]);
            mpq_set(mpqPoints[i], points[i][0].get_mpq_t());
            mpq_set(mpqVals[i], vals[i].get_mpq_t());
        }

        AltArr_t* interp = univarInterpolate_AA(mpqPoints, mpqVals, points.size());

        for (int i = 0; i < points.size(); ++i) {
            mpq_clear(mpqPoints[i]);
            mpq_clear(mpqVals[i]);
        }

        SparseMultivariateRationalPolynomial ret(1);
        ret.poly = interp;
        return ret;
    }

    std::cerr << "BPAS SMQP: interpolate(points, vals) NOT YET IMPLEMENTED" << std::endl;
    return 0;
}


/**
 * Divide this by polynomial b, returning the quotient and remainder in q and r, respectively.
 *
 * returns a boolean indicating if the division was exact.
 */
bool SparseMultivariateRationalPolynomial::divide(const SparseMultivariateRationalPolynomial& b, SparseMultivariateRationalPolynomial& q, SparseMultivariateRationalPolynomial& r) const {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
        exit(1);
    }

    if (b.isConstant()) {
        q = *this / b.poly->elems->coef;
        r = SparseMultivariateRationalPolynomial(NULL, nvar, names);
        return true;
    }
    if (isConstant()) {
        q = SparseMultivariateRationalPolynomial(NULL, nvar, names);
        r = *this;
        return false;
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);

    if (!isOrdered) {
        std::cout << "BPAS: error, trying to divide between Q[";
        for (int i = 1; i <= nvar; ++i) {
            std::cout << names[i];
            if (i < nvar) { std::cout << ", "; }
        }
        std::cout << "] and Q[";
        for (int i = 1; i <= b.nvar; ++i) {
            std::cout << b.names[i];
            if (i < b.nvar) { std::cout << ", "; }
        }
        std::cout << "]." << std::endl;
        exit(1);
    }

    int superNvar = xs.size() / 2;
    if (superNvar != nvar || superNvar != b.nvar) {
        //map indices to the expanded superset.
        int varmap[nvar];
        int bvarmap[b.nvar];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                varmap[xs[i]-1] = i/2;
            }
            if (xs[i+1] != 0) {
                bvarmap[xs[i+1]-1] = i/2;
            }
        }

        //create new combined names array
        Symbol newnames[superNvar+1];
        //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
        newnames[0] = this->names[0];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                newnames[(i/2) + 1] = names[xs[i]];
            } else {
                newnames[(i/2) + 1] = b.names[xs[i+1]];
            }
        }

        // std::cerr << "newnames: "  << std::endl;
        // for (auto v : newnames) {
        //     std::cerr << "v: " << v << std::endl;
        // }

        // for (int i = 0; i < superNvar; ++i) {
        //     std::cerr << "avarmap[" << i << "]: " << varmap[i] << std::endl;
        // }
        // for (int i = 0; i < superNvar; ++i) {
        //     std::cerr << "bvarmap[" << i << "]: " << bvarmap[i] << std::endl;
        // }
        // std::cerr << std::endl;

        SparseMultivariateRationalPolynomial tempc = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        AltArr_t* qAA = NULL, * rAA = NULL;
        dividePolynomials_AA(tempc.poly, tempb.poly, &qAA, &rAA, superNvar);

        q = SparseMultivariateRationalPolynomial(qAA, superNvar, newnames);
        r = SparseMultivariateRationalPolynomial(rAA, superNvar, newnames);

        return (rAA == NULL || rAA->size == 0);

    } else {

        // std::cerr << "divide on same ring" << std::endl;

        AltArr_t* qAA = NULL, * rAA = NULL;
        dividePolynomials_AA(this->poly, b.poly, &qAA, &rAA, nvar);
        q = SparseMultivariateRationalPolynomial(qAA, nvar, names);
        r = SparseMultivariateRationalPolynomial(rAA, nvar, names);
        return (rAA == NULL || rAA->size == 0);
    }
}

/**
 * Get the remainder of *this divided by b.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator% (const SparseMultivariateRationalPolynomial& b) const {
    SparseMultivariateRationalPolynomial q,r;
    this->divide(b, q, r);
    return r;
}

/**
 * Update *this by setting it to the remainder of *this divided by b.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator%= (const SparseMultivariateRationalPolynomial& b) {
    *this = (*this % b);
    return *this;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::pseudoDivideBySMZP(const SparseMultivariateRationalPolynomial& b, SparseMultivariateRationalPolynomial* quo, SparseMultivariateRationalPolynomial* mult, bool lazy) const {
    if (b.isZero()) {
        std::cerr << "BPAS: error, pseudo-dividend is zero from SMQP." << std::endl;
        exit(1);
    }
    if (b.isConstant()) {
        std::cerr << "BPAS: error, pseudo-dividend is a constant in SMQP." << std::endl;
        exit(1);
    }

    if (isZero()) {
        if (quo != NULL) {
            (*quo).zero();
        }
        if (mult != NULL) {
            (*mult).one();
        }
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    //This will also capture case when nvar = 0 for both.
    if (b.isConstant() && isConstant()) {
        if (quo != NULL) {
            *quo = *this;
        }
        if (mult != NULL) {
            *mult = b;
        }
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);

    if (!isOrdered) {
        std::cerr << "p: " << *this << std::endl;
        std::cerr << "q: " << b << std::endl;
        std::cout << "BPAS: error, trying to pseudo divide between Q[";
        for (int i = 1; i <= nvar; ++i) {
            std::cout << names[i];
            if (i < nvar) { std::cout << ", "; }
        }
        std::cout << "] and Q[";
        for (int i = 1; i <= b.nvar; ++i) {
            std::cout << b.names[i];
            if (i < b.nvar) { std::cout << ", "; }
        }
        std::cout << "]." << std::endl;
        exit(1);
    }

    SparseMultivariateRationalPolynomial q, r;

    int superNvar = xs.size() / 2;
    if (superNvar != nvar || superNvar != b.nvar) {
        //map indices to the expanded superset.
        int varmap[nvar];
        int bvarmap[b.nvar];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                varmap[xs[i]-1] = i/2;
            }
            if (xs[i+1] != 0) {
                bvarmap[xs[i+1]-1] = i/2;
            }
        }

        //create new combined names array
        Symbol newnames[superNvar+1];
        //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
        newnames[0] = this->names[0];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                newnames[(i/2) + 1] = names[xs[i]];
            } else {
                newnames[(i/2) + 1] = b.names[xs[i+1]];
            }
        }

        SparseMultivariateRationalPolynomial tempc = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        //Now that we have expanded this and B to live in the same ring, we need
        //to if the main variable of b is the same as c.

        Symbol leadVar = tempb.leadingVariable();
        if (leadVar != tempb.names[1]) {
            //need to rearrange *this to have same leading var and then convert to rec form.
            Symbol newvars[superNvar+1];
            newvars[0] = tempb.names[0];
            int reorderVarmap[superNvar];
            int before = 1;
            for (int i = 0; i < superNvar; ++i) {
                if (tempb.names[i+1] == leadVar) {
                    before = 0;
                    newvars[1] = leadVar;
                    reorderVarmap[i] = 0;
                    //leadVar already in newvars
                    continue;
                }

                //direct map for vars after leadVar, otherwise shift up by 1.
                newvars[i+1+before] = tempb.names[i+1];
                reorderVarmap[i] = i+before;
            }

            int invVarmap[superNvar];
            for (int i = 0; i < superNvar; ++i) {
                invVarmap[reorderVarmap[i]] = i;
            }


            tempc.reorderVarsInPlace(reorderVarmap);
            tempb.reorderVarsInPlace(reorderVarmap);

            AltArrZ_t* cAA = deepCopyPolynomial_AAZFromAA(tempc.poly);
            AltArrZ_t* bAA = deepCopyPolynomial_AAZFromAA(tempb.poly);

            RecArrZ_t* recC = convertToRecursiveArrayZ(cAA);
            RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

            AltArrZ_t* qAA = NULL, * rAA = NULL;
            AltArrZ_t* hPow = NULL;
            int e = 0;
            pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);

            freeRecArrayZAndCoef(recC); //this frees cAA and bAA too
            freeRecArrayZAndCoef(recB);

            r = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(rAA), superNvar, newvars);
            r.reorderVarsInPlace(invVarmap);
            freePolynomial_AAZ(rAA);

            if (quo != NULL) {
                q = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(qAA), superNvar, newvars);
                q.reorderVarsInPlace(invVarmap);
                *quo = q;
            }
            if (mult != NULL) {
                SparseMultivariateRationalPolynomial smqpHPow(deepCopyPolynomial_AAFromAAZ(hPow), superNvar, newvars);
                smqpHPow.reorderVarsInPlace(invVarmap);
                *mult = smqpHPow;
            }

            freePolynomial_AAZ(qAA);
            freePolynomial_AAZ(hPow);
            return r;
        }

        AltArrZ_t* cAA = deepCopyPolynomial_AAZFromAA(tempc.poly);
        AltArrZ_t* bAA = deepCopyPolynomial_AAZFromAA(tempb.poly);

        RecArrZ_t* recC = convertToRecursiveArrayZ(cAA);
        RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

        AltArrZ_t* qAA = NULL, *rAA = NULL;
        AltArrZ_t* hPow = NULL;
        int e = 0;
        pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);

        freeRecArrayZAndCoef(recC); //this frees cAA and bAA too
        freeRecArrayZAndCoef(recB);

        r = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(rAA), superNvar, newnames);
        freePolynomial_AAZ(rAA);

        if (quo != NULL) {
            q = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(qAA), superNvar, newnames);
            *quo = q;
        }
        if (mult != NULL) {
            SparseMultivariateRationalPolynomial smqpHPow(deepCopyPolynomial_AAFromAAZ(hPow), superNvar, newnames);
            *mult = smqpHPow;
        }

        freePolynomial_AAZ(qAA);
        freePolynomial_AAZ(hPow);
        return r;

    } else {
        //in this case, this and b have same variables in same order.

        Symbol leadVar = b.leadingVariable();
        if (leadVar != b.names[1]) {

            //need to rearrange *this to have same leading var and then convert to rec form.
            Symbol newvars[nvar+1];
            newvars[0] = b.names[0];
            int varmap[nvar];
            int before = 1;
            for (int i = 0; i < nvar; ++i) {
                if (b.names[i+1] == leadVar) {
                    before = 0;
                    newvars[1] = leadVar;
                    varmap[i] = 0;
                    //leadVar already in newvars
                    continue;
                }

                //direct map for vars after leadVar, otherwise shift up by 1.
                newvars[i+1+before] = b.names[i+1];
                varmap[i] = i+before;
            }

            int invVarmap[nvar];
            for (int i = 0; i < nvar; ++i) {
                invVarmap[varmap[i]] = i;
            }
            SparseMultivariateRationalPolynomial cReordered = this->expandVariables(nvar, newvars, varmap);
            SparseMultivariateRationalPolynomial bReordered = b.expandVariables(nvar, newvars, varmap);

            AltArrZ_t* cAA = deepCopyPolynomial_AAZFromAA(cReordered.poly);
            AltArrZ_t* bAA = deepCopyPolynomial_AAZFromAA(bReordered.poly);

            RecArrZ_t* recC = convertToRecursiveArrayZ(cAA);
            RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

            AltArrZ_t* qAA = NULL, * rAA = NULL;
            AltArrZ_t* hPow = NULL;
            int e = 0;
            pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);

            freeRecArrayZAndCoef(recC); //this frees cAA and bAA too
            freeRecArrayZAndCoef(recB);

            r = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(rAA), nvar, newvars);
            r.reorderVarsInPlace(invVarmap);
            freePolynomial_AAZ(rAA);

            if (quo != NULL) {
                q = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(qAA), nvar, newvars);
                q.reorderVarsInPlace(invVarmap);
                *quo = q;
            }
            if (mult != NULL) {
                SparseMultivariateRationalPolynomial smqpHPow(deepCopyPolynomial_AAFromAAZ(hPow), superNvar, newvars);
                smqpHPow.reorderVarsInPlace(invVarmap);
                *mult = smqpHPow;
            }

            freePolynomial_AAZ(qAA);
            freePolynomial_AAZ(hPow);
            return r;
       }

       if (nvar == 1 && b.nvar == 1) {
           AltArr_t* qAA = NULL, *rAA = NULL;
           int e = 0;
           univariatePseudoDividePolynomials_AA(this->poly, b.poly, &qAA, &rAA, &e, lazy);

           if (quo != NULL) {
               *quo = SparseMultivariateRationalPolynomial(qAA, nvar, names);
           }

           if (mult != NULL) {
                RationalNumber h(b.poly->elems->coef);
                h ^= e;
                *mult = h;
            }
            return SparseMultivariateRationalPolynomial(rAA, nvar, names);
        }


        AltArrZ_t* cAA = deepCopyPolynomial_AAZFromAA(this->poly);
        AltArrZ_t* bAA = deepCopyPolynomial_AAZFromAA(b.poly);

        RecArrZ_t* recC = convertToRecursiveArrayZ(cAA);
        RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

        AltArrZ_t* qAA = NULL, *rAA = NULL;
        AltArrZ_t* hPow = NULL;
        int e = 0;
        pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);

        freeRecArrayZAndCoef(recC); //this frees cAA and bAA too
        freeRecArrayZAndCoef(recB);

        r = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(rAA), superNvar, names);
        freePolynomial_AAZ(rAA);

        if (quo != NULL) {
            q = SparseMultivariateRationalPolynomial(deepCopyPolynomial_AAFromAAZ(qAA), superNvar, names);
            *quo = q;
        }
        if (mult != NULL) {
            SparseMultivariateRationalPolynomial smqpHPow(deepCopyPolynomial_AAFromAAZ(hPow), superNvar, names);
            *mult = smqpHPow;
        }

        freePolynomial_AAZ(qAA);
        freePolynomial_AAZ(hPow);

        return r;
    }
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::pseudoDivide(const SparseMultivariateRationalPolynomial& b, SparseMultivariateRationalPolynomial* quo, SparseMultivariateRationalPolynomial* mult, bool lazy) const {
    if (b.isZero()) {
        std::cerr << "BPAS: error, pseudo-dividend is zero from SMQP." << std::endl;
        exit(1);
    }
    if (b.isConstant()) {
        std::cerr << "BPAS: error, pseudo-dividend is a constant in SMQP." << std::endl;
        exit(1);
    }

    if (isZero()) {
        if (quo != NULL) {
            (*quo).zero();
        }
        if (mult != NULL) {
            (*mult).one();
        }
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    //This will also capture case when nvar = 0 for both.
    if (b.isConstant() && isConstant()) {
        if (quo != NULL) {
            *quo = *this;
        }
        if (mult != NULL) {
            *mult = b;
        }
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    if (integerPolynomialTest_AA(this->poly) && integerPolynomialTest_AA(b.poly)) {
        return this->pseudoDivideBySMZP(b, quo, mult, lazy);
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);

    if (!isOrdered) {
        std::cerr << "p: " << *this << std::endl;
        std::cerr << "q: " << b << std::endl;
        std::cout << "BPAS: error, trying to pseudo divide between Q[";
        for (int i = 1; i <= nvar; ++i) {
            std::cout << names[i];
            if (i < nvar) { std::cout << ", "; }
        }
        std::cout << "] and Q[";
        for (int i = 1; i <= b.nvar; ++i) {
            std::cout << b.names[i];
            if (i < b.nvar) { std::cout << ", "; }
        }
        std::cout << "]." << std::endl;
        exit(1);
    }

    SparseMultivariateRationalPolynomial q, r;

    int superNvar = xs.size() / 2;
    if (superNvar != nvar || superNvar != b.nvar) {
        //map indices to the expanded superset.
        int varmap[nvar];
        int bvarmap[b.nvar];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                varmap[xs[i]-1] = i/2;
            }
            if (xs[i+1] != 0) {
                bvarmap[xs[i+1]-1] = i/2;
            }
        }

        //create new combined names array
        Symbol newnames[superNvar+1];
        //isOrderedRing returns false if this.names[0] != b.names[0], therefore safe
        newnames[0] = this->names[0];
        for (int i = 0; i < xs.size(); i += 2) {
            if (xs[i] != 0) {
                newnames[(i/2) + 1] = names[xs[i]];
            } else {
                newnames[(i/2) + 1] = b.names[xs[i+1]];
            }
        }

        SparseMultivariateRationalPolynomial tempc = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        //Now that we have expanded this and B to live in the same ring, we need
        //to if the main variable of b is the same as c.

        Symbol leadVar = tempb.leadingVariable();
        if (leadVar != tempb.names[1]) {
            //need to rearrange *this to have same leading var and then convert to rec form.
            Symbol newvars[superNvar+1];
            newvars[0] = tempb.names[0];
            int reorderVarmap[superNvar];
            int before = 1;
            for (int i = 0; i < superNvar; ++i) {
                if (tempb.names[i+1] == leadVar) {
                    before = 0;
                    newvars[1] = leadVar;
                    reorderVarmap[i] = 0;
                    //leadVar already in newvars
                    continue;
                }

                //direct map for vars after leadVar, otherwise shift up by 1.
                newvars[i+1+before] = tempb.names[i+1];
                reorderVarmap[i] = i+before;
            }

            int invVarmap[superNvar];
            for (int i = 0; i < superNvar; ++i) {
                invVarmap[reorderVarmap[i]] = i;
            }

            tempc.reorderVarsInPlace(reorderVarmap);
            tempb.reorderVarsInPlace(reorderVarmap);

            AltArr_t* cAA = tempc.poly;
            AltArr_t* bAA = tempb.poly;

            RecArr_t* recC = convertToRecursiveArray(cAA);
            RecArr_t* recB = convertToRecursiveArray(bAA);

            AltArr_t* qAA = NULL, * rAA = NULL;
            AltArr_t* hPow = NULL;
            int e = 0;
            pesudoDivide_RecArray(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);
            cAA = convertFromRecursiveArray(recC, superNvar);
            bAA = convertFromRecursiveArray(recB, superNvar);
            tempc.poly = cAA;
            tempb.poly = bAA;

            q = SparseMultivariateRationalPolynomial(qAA, superNvar, newvars);
            r = SparseMultivariateRationalPolynomial(rAA, superNvar, newvars);
            q.reorderVarsInPlace(invVarmap);
            r.reorderVarsInPlace(invVarmap);

            SparseMultivariateRationalPolynomial smqpHPow(hPow, superNvar, newvars);
            smqpHPow.reorderVarsInPlace(invVarmap);
            if (quo != NULL) {
                *quo = q;
            }
            if (mult != NULL) {
                *mult = smqpHPow;
            }
            return r;
        }

        AltArr_t* cAA = tempc.poly;
        AltArr_t* bAA = tempb.poly;
        RecArr_t* recC = convertToRecursiveArray(cAA);
        RecArr_t* recB = convertToRecursiveArray(bAA);

        AltArr_t* qAA = NULL, * rAA = NULL;
        AltArr_t* hPow = NULL;
        int e = 0;
        pesudoDivide_RecArray(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);
        cAA = convertFromRecursiveArray(recC, superNvar);
        bAA = convertFromRecursiveArray(recB, superNvar);
        tempc.poly = cAA;
        tempb.poly = bAA;

        q = SparseMultivariateRationalPolynomial(qAA, superNvar, newnames);
        r = SparseMultivariateRationalPolynomial(rAA, superNvar, newnames);
        SparseMultivariateRationalPolynomial smqpHPow(hPow, superNvar, newnames);
        if (quo != NULL) {
            *quo = q;
        }
        if (mult != NULL) {
            *mult = smqpHPow;
        }
        return r;

    } else {
        //in this case, this and b have same variables in same order.

        Symbol leadVar = b.leadingVariable();
        if (leadVar != b.names[1]) {

            //need to rearrange *this to have same leading var and then convert to rec form.
            Symbol newvars[nvar+1];
            newvars[0] = b.names[0];
            int varmap[nvar];
            int before = 1;
            for (int i = 0; i < nvar; ++i) {
                if (b.names[i+1] == leadVar) {
                    before = 0;
                    newvars[1] = leadVar;
                    varmap[i] = 0;
                    //leadVar already in newvars
                    continue;
                }

                //direct map for vars after leadVar, otherwise shift up by 1.
                newvars[i+1+before] = b.names[i+1];
                varmap[i] = i+before;
            }

            int invVarmap[nvar];
            for (int i = 0; i < nvar; ++i) {
                invVarmap[varmap[i]] = i;
            }
            SparseMultivariateRationalPolynomial cReordered = this->expandVariables(nvar, newvars, varmap);
            SparseMultivariateRationalPolynomial bReordered = b.expandVariables(nvar, newvars, varmap);

            AltArr_t* cAA = cReordered.poly;
            AltArr_t* bAA = bReordered.poly;
            RecArr_t* recC = convertToRecursiveArray(cAA);
            RecArr_t* recB = convertToRecursiveArray(bAA);

            AltArr_t* qAA = NULL, * rAA = NULL;
            AltArr_t* hPow = NULL;
            int e = 0;
            pesudoDivide_RecArray(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);
            cReordered.poly = cAA = convertFromRecursiveArray(recC, nvar);
            bReordered.poly = bAA = convertFromRecursiveArray(recB, nvar);

            q = SparseMultivariateRationalPolynomial(qAA, nvar, newvars);
            r = SparseMultivariateRationalPolynomial(rAA, nvar, newvars);

            q.reorderVarsInPlace(invVarmap);
            r.reorderVarsInPlace(invVarmap);

            SparseMultivariateRationalPolynomial smqpHPow(hPow, superNvar, newvars);
            smqpHPow.reorderVarsInPlace(invVarmap);
            if (quo != NULL) {
                *quo = q;
            }
            if (mult != NULL) {
               *mult = smqpHPow;
           }
           return r;
       }

       if (nvar == 1 && b.nvar == 1) {
           AltArr_t* qAA = NULL, *rAA = NULL;
           int e = 0;
           univariatePseudoDividePolynomials_AA(this->poly, b.poly, &qAA, &rAA, &e, lazy);

           if (quo != NULL) {
               *quo = SparseMultivariateRationalPolynomial(qAA, nvar, names);
           }

           if (mult != NULL) {
                RationalNumber h(b.poly->elems->coef);
                h ^= e;
                *mult = h;
            }
            return SparseMultivariateRationalPolynomial(rAA, nvar, names);
        }

        AltArr_t* cAA = NULL; //(this->poly);
        AltArr_t* bAA = b.poly;
        RecArr_t* recC = convertToRecursiveArray(poly);
        RecArr_t* recB = convertToRecursiveArray(bAA);

        AltArr_t* qAA = NULL, *rAA = NULL;
        AltArr_t* hPow = NULL;
        int e = 0;
        pesudoDivide_RecArray(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);

        poly = convertFromRecursiveArray(recC, superNvar);
        bAA = convertFromRecursiveArray(recB, superNvar);

        q = SparseMultivariateRationalPolynomial(qAA, superNvar, names);
        r = SparseMultivariateRationalPolynomial(rAA, superNvar, names);

        SparseMultivariateRationalPolynomial smqpHPow(hPow, superNvar, names);
        if (quo != NULL) {
            *quo = q;
        }
        if (mult != NULL) {
            *mult = smqpHPow;
        }

        return r;
    }
}



/**
 * Add *this and a ratNum_t.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator+ (const RationalNumber& c) const {
    SparseMultivariateRationalPolynomial ret = *this;

    if (isZero()) {
        if (ret.poly != NULL) {
            freePolynomial_AA(ret.poly);
        }
        ret.poly = makeConstPolynomial_AA(1, ret.nvar, c.get_mpq_t());
        return ret;
    }

    addRationalNumber_AA_inp(ret.poly, c.get_mpq_t());
    return ret;
}


/**
 * Update *this by adding r
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator+= (const RationalNumber& c) {
    slp.clear();

    if (isZero()) {
        if (poly != NULL) {
            freePolynomial_AA(poly);
        }
        poly = makeConstPolynomial_AA(1, nvar, c.get_mpq_t());
        return *this;
    }

    addRationalNumber_AA_inp(poly, c.get_mpq_t());
    return *this;
}

/**
 * Subtract the ratNum_t r from *this.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator- (const RationalNumber& c) const {
    RationalNumber negC = -c;
    return *this + negC;
}

/**
 * Update *this by subtracting ratNum_t r.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator-= (const RationalNumber& c) {
    RationalNumber negC = -c;
    return *this += negC;
}

/**
 * Multiply *this by ratNum_t r.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator* (const RationalNumber& c) const {
    if (c == 0 || isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    SparseMultivariateRationalPolynomial ret = *this;

    if (c == 1) {
        return ret;
    }

    mpq_t mult;
    mpq_init(mult);
    mpq_set(mult, c.get_mpq_t());
    for(int i = 0; i < ret.poly->size; ++i) {
        mpq_mul(ret.poly->elems[i].coef, ret.poly->elems[i].coef, mult);
    }
    mpq_clear(mult);

    return ret;
}

/**
 * Update *this by multiplying by ratNum_t r.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator*= (const RationalNumber& c) {
    slp.clear();
    if (isZero()) {
        return *this;
    }
    if (c == 0) {
        freePolynomial_AA(poly);
        poly = NULL;
        return *this;
    }
    if (c == 1) {
        return *this;
    }

    mpq_t mult;
    mpq_init(mult);
    mpq_set(mult, c.get_mpq_t());
    for(int i = 0; i < poly->size; ++i) {
        mpq_mul(poly->elems[i].coef, poly->elems[i].coef, mult);
    }
    mpq_clear(mult);

    return *this;
}

/**
 * Divide *this by ratNum_t r.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator/ (const RationalNumber& c) const {
    if (c == 0) {
        std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
        exit(1);
    }

    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    if (c == 1) {
        return SparseMultivariateRationalPolynomial(*this);
    }

    ratNum_t rInv;
    mpq_init(rInv);
    mpq_inv(rInv, c.get_mpq_t());
    SparseMultivariateRationalPolynomial ret = (*this * rInv);
    mpq_clear(rInv);
    return ret;
}

/**
 * Divide ratNum_t r by SMQP b.
 */
SparseMultivariateRationalPolynomial operator/ (const ratNum_t& r, const SparseMultivariateRationalPolynomial& b) {
    if (b.isConstant()) {
        if (mpq_sgn(b.poly->elems->coef) == 0) {
            std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
            exit(1);
        }
        AltArr_t* rPoly = makeConstPolynomial_AA(1, b.nvar, r);
        mpq_div(rPoly->elems->coef, r, b.poly->elems->coef);
        return SparseMultivariateRationalPolynomial(rPoly, b.nvar, b.names);
    } else {
        std::cerr << "BPAS ERROR: SMQP non-exact division in mpq_t / SparseMultivariateRationalPolynomial" << std::endl;
        exit(1);
    }
}

/**
 * Update *this by dividing by ratNum_t r.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator/= (const RationalNumber& c) {
    *this = *this / c;
    return *this;
}


/**
 * Get the polynomial term at index. Returns 0 if index is beyond the
 * number of terms in this polynomial.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator[] (int index) const {
    AltArr_t* term = termAtIdx_AA(poly, index);
    return SparseMultivariateRationalPolynomial(term, nvar, names);
}


/**
 * Get the leading variable, that is, the highest-order variable with positive degree
 * of this polynomial.
 * returns the leading variable or the empty string if this polynomial has zero variables.
 */
Symbol SparseMultivariateRationalPolynomial::leadingVariable() const {
    if (nvar == 0 || isConstant()) {
        return Symbol("");
    }

    int lvIdx = mainVariable_AA(poly);
    return names[lvIdx + 1];
}

/**
 * Get the degree of this polynomial w.r.t the leading variable.
 */
Integer SparseMultivariateRationalPolynomial::leadingVariableDegree() const {
    if (nvar == 0 || isConstant()) {
        return 0;
    }

    return mainDegree_AA(poly);
}

/**
 * Is the contant term zero.
 */
bool SparseMultivariateRationalPolynomial::isConstantTermZero() const {
    if (isZero()) {
        return 1;
    }

    return isConstantTermZero_AA(poly);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::leadingCoefficientInVariable (const Symbol& x, int* e) const {
    int deg = 0;
    if (e != NULL) {
        *e = 0;
    }
    if (isZero()){
        return *this;
    }

    int k = 0;
    for (int i = 1; i <= nvar; ++i) {
        if (names[i] == x) {
            k = i;
            break;
        }
    }
    if (k == 0) {
        //If supplied variable is not one of our variables, then really all of *this
        //is leading coeff w.r.t input Symbol x.
        return *this;
    }

    int v = nvar - 1;
    SparseMultivariateRationalPolynomial r(v);
    for (int i = 0; i < k; ++i) {
        r.names[i] = names[i];
    }
    for (int i = k; i < v+1; ++i) {
        r.names[i] = names[i+1];
    }

    if (isConstant()) {
        r.poly = makeConstPolynomial_AA(1, v, poly->elems->coef);
        return r;
    }

    --k; //recall names index is +1 of degs index

    int varmap[nvar];
    varmap[k] = 0;
    for (int i = 0; i < k; ++i) {
        varmap[i] = i+1;
    }
    for (int i = k+1; i < nvar; ++i) {
        varmap[i] = i;
    }

    AltArr_t* lc = NULL;
    if (k != 0) {
        AltArr_t* tempPoly = deepCopyPolynomial_AA(poly);
        reorderVars_AA(tempPoly, varmap, nvar);
        deg = mainLeadingDegree_AA(tempPoly);
        lc = mainLeadingCoefficient_AA(tempPoly);
        freePolynomial_AA(tempPoly);
    } else {
        deg = mainLeadingDegree_AA(poly);
        lc = mainLeadingCoefficient_AA(poly);
    }


    //we have consturcted a polynomial in v+1 vars but we want v vars, so re-pack each degs.
    shrinkNumVarsAtIdx_AA(lc, 0);

    r.poly = lc;

    if (e != NULL) {
        *e = deg;
    }
    return r;
}

/**
 * Convert to a SUP<SMQP> given the variable 'x'.
 *
 * returns the SUP<SMQP>.
 */
SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::convertToSUP(const Symbol& x) const {
    SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> r;
    r.setVariableName(x);

    int k = 0;
    for (int i = 1; i <= nvar; ++i) {
        if (names[i] == x) {
            k = i;
            break;
        }
    }

    if (k == 0) {
        r.setCoefficient(0, *this);
        return r;
    }

    if (nvar == 1) {
        degree_t d;
        for (int i = 0; i < poly->size; ++i) {
            SparseMultivariateRationalPolynomial coef(RationalNumber(poly->elems[i].coef), 0);
            d = partialDegreeTerm_AA(poly, i, 0);
            r.setCoefficient(d, coef);
        }
        return r;
    }

    int v = nvar - 1;
    SparseMultivariateRationalPolynomial tempSMQP(v);
    for (int i = 0; i < k; ++i) {
        tempSMQP.names[i] = names[i];
    }
    for (int i = k; i < v+1; ++i) {
        tempSMQP.names[i] = names[i+1];
    }

    --k; //recall names index is +1 of degs index

    int varmap[nvar];
    varmap[k] = 0;
    for (int i = 0; i < k; ++i) {
        varmap[i] = i+1;
    }
    for (int i = k+1; i < nvar; ++i) {
        varmap[i] = i;
    }

    AltArr_t* tempPoly;
    if (k != 0) {
        tempPoly = deepCopyPolynomial_AA(poly);
        reorderVars_AA(tempPoly, varmap, nvar);
    } else {
        tempPoly = poly;
    }

    RecArr_t* recArr = convertToRecursiveArray(tempPoly);
    AltArr_t localAA;
    for (int i = 0; i < recArr->size; ++i) {
        degree_t deg = recArr->elems[i].exp;
        localAA.nvar = nvar;
        localAA.size = localAA.alloc = recArr->elems[i].coefSize;
        localAA.unpacked = recArr->unpacked;
        localAA.elems = recArr->elems[i].coef;

        AltArr_t* toShrinkAA = deepCopyPolynomial_AA(&localAA);
        shrinkNumVarsAtIdx_AA(toShrinkAA, 0);
        tempSMQP.poly = toShrinkAA;

        r.setCoefficient(deg, tempSMQP);

        freePolynomial_AA(tempSMQP.poly);
        tempSMQP.poly = NULL;
    }

    if (k != 0) {
        freeRecArrayAndCoef(recArr);
    } else {
        convertFromRecursiveArray(recArr, nvar);
    }

    // std::vector<SparseMultivariateRationalPolynomial> mpolys(d+1, SparseMultivariateRationalPolynomial(v));

    // Symbol newnames[v+1];
    // for (int i = 0; i < k; ++i) {
    //     newnames[i] = names[i];
    // }
    // for (int i = k; i < v+1; ++i) {
    //     newnames[i] = names[i+1];
    // }
    // //setup names and nvar of each of the polys
    // for (int i = 0; i <= d; ++i) {
    //     std::copy(newnames, newnames+v+1, mpolys[i].names);
    // }

    // unsigned long long int* masks = getExpMaskArray(nvar);
    // int* sizes = getExpOffsetArray(nvar);

    // bool isLeadVar = (this->leadingVariable() == names[k]);
    // --k;

    // degrees_t curDegs;
    // degree_t curD;
    // for (int i = 0; i < poly->size; ++i) {
    //     curD = GET_NTH_EXP(poly->elems[i].degs, masks[k], sizes[k]);
    //     if (mpolys[curD].poly == NULL) {
    //         mpolys[curD].poly = makePolynomial_AA(10, nvar);
    //     }
    //     if (mpolys[curD].poly->alloc <= mpolys[curD].poly->size) {
    //         resizePolynomial_AA(mpolys[curD].poly, mpolys[curD].poly->size * 2);
    //     }
    //     mpq_init(mpolys[curD].poly->elems[mpolys[curD].poly->size].coef);
    //     mpq_set(mpolys[curD].poly->elems[mpolys[curD].poly->size].coef, poly->elems[i].coef);
    //     mpolys[curD].poly->elems[mpolys[curD].poly->size].degs = poly->elems[i].degs & (~masks[k]);
    //     ++(mpolys[curD].poly->size);
    // }

    // for (int i = 0; i <= d; ++i) {
    //     if (mpolys[i].isZero()) {
    //         continue;
    //     }
    //     shrinkNumVarsAtIdx_AA(mpolys[i].poly, k);
    //     if (!(k == 0 || isLeadVar)) {
    //         mergeSortPolynomial_AA(mpolys[i].poly);
    //     }
    //     r.setCoefficient(i, mpolys[i]);
    // }

    // free(masks);
    // free(sizes);

    return r;
}


/**
 * Negate all the coefficients of *this. Note, that due to the
 * sharing nature of underling nodes, this may alter the Nodes of
 * other SMQP.
 */
void SparseMultivariateRationalPolynomial::negate() {
    if (isZero()) {
        return;
    }
    negatePolynomial_AA(poly);
}

/**
 * Get a copy of this polynomial.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::deepCopy() const {
    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    AltArr_t* copy = deepCopyPolynomial_AA(poly);
    return SparseMultivariateRationalPolynomial(copy, nvar, names);
}

Factors<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::factor() const {
    Factors<SparseMultivariateIntegerPolynomial> ZFacts;

    SparseMultivariateIntegerPolynomial f;
    RationalNumber content;
    f = this->primitivePartSMZP(content);
    ZFacts = f.factor();
//    cout << "ZFacts[0] = " << ZFacts.ringElement() << endl;
//for (int i=0; i<ZFacts.size(); ++i)
//	cout << "ZFacts[" << i+1 << "] = " << ZFacts[i] << endl;
    content *= RationalNumber(ZFacts.ringElement().leadingCoefficient());

    std::vector<SparseMultivariateRationalPolynomial> QFacts;
    std::vector<int> QExps;
    QFacts.reserve(ZFacts.size());
    QExps.reserve(ZFacts.size());
    for (int i=0; i<ZFacts.size(); ++i) {
    	QFacts.emplace_back(ZFacts[i].first);
    	QExps.push_back(ZFacts[i].second);
    }

    return Factors<SparseMultivariateRationalPolynomial>(QFacts,QExps,content);
}

static int checkDegs(degree_t* a, degree_t* b, int vars) {
    for (int i = 0; i < vars; ++i) {
        if ( (a[i]) != (b[i]) ) {
            return i;
        }
    }

    return -1;
}

static int firstNonZero(degree_t* d, int vars) {
    for (int i = vars-1; i >= 0; --i) {
        if (d[i]) {
            return i;
        }
    }

    return -1;
}

/**
 * SLP representation of the polynomial
 **/
void SparseMultivariateRationalPolynomial::straightLineProgram() {
    slp.clear();
    if (isZero()) {
        SLPRepresentation r;
        r.type = 4;
        r.b = -1;
        r.a.c = new RationalNumber(0);
        slp.push_back(r);
        return;
    }

    if (isConstant()) {
        SLPRepresentation r;
        r.type = 4;
        r.b = -1;
        r.a.c = new RationalNumber(poly->elems->coef);
        slp.push_back(r);
        return;
    }

    ///////
    // So basically what is happening in this below code is sort of like
    // Horner's method. We start with the highest ordered term and begin
    // at the center of this kind-of horner's method, expanding outwards
    // and adding coefficients of terms as needed, multiplying by variables
    // as well to achieve the proper degree for each term.
    //
    // the variable d sort of keeps track of how many multiplications remain
    // for each variable, for every preceeding term in the SLP.
    // you can see that d[i] is decremented when we add a new multiplicative
    // element to the slp.
    //
    // The stack here allows for multi "expansions" of horner's method for
    // different terms. When a term appears that has incompatible degree
    // with what is left to be multiplied in the preceeding horner's expansion
    // then we push that horner's expansion onto the stack and start another.
    // then finally, popping from the stack, and adding two horner's expansions
    // together.
    //


    std::vector<AAElem_DegList_t*> polyStack;
    std::vector<int> indexStack;
    int prev = 0;
    bool isOp = 1;
    SLPRepresentation elem;

    AltArrDegList_t* listPoly = deepCopyPolynomial_AADegListFromAA(poly);

    degree_t* d = (degree_t*) malloc(sizeof(degree_t)*nvar);
    for (int i = 0; i < nvar; ++i) {
        d[i] = listPoly->elems->degs[i];
    }

    int polySize = listPoly->size;
    for (int i = 0; i < polySize; ++i) {
        degree_t* curd = listPoly->elems[i].degs;
        int k = checkDegs(d, curd, nvar);
        if (k > -1) {
            bool isStart = 0;

            // Mutiply smaller variate than the current one
            for (int j = nvar-1; j > k; --j) {
                if (d[j] > 0 && !polyStack.empty()) {
                    degree_t* prevd = polyStack[polyStack.size()-1]->degs;
                    while(checkDegs(prevd, d, nvar) == j) {
                        elem.op = 0;
                        elem.type = 3;
                        elem.a.i = indexStack[indexStack.size()-1];
                        elem.b = slp.size() - 1;
                        slp.push_back(elem);

                        polyStack.pop_back();
                        indexStack.pop_back();
                        if (!polyStack.empty()) {
                            prevd = polyStack[polyStack.size()-1]->degs;
                        } else {
                            break;
                        }
                    }
                }

                while (d[j] > 0) {
                    elem.op = 1;
                    elem.type = 2;
                    elem.a.i = j;
                    // elem.a = j;
                    elem.b = slp.size() - 1;
                    slp.push_back(elem);
                    d[j]--;
                }
                if (curd[j]) {
                    isStart = 1;
                    d[j] = curd[j];
                }
            }

            // Add previous result if bigger variate doesn't end
            if (!polyStack.empty()) {
                degree_t* prevd = polyStack[polyStack.size()-1]->degs;
                while(checkDegs(prevd, curd, nvar) == k) {
                    elem.op = 0;
                    elem.type = 3;
                    elem.a.i = indexStack[indexStack.size()-1];
                    elem.b = slp.size() - 1;
                    slp.push_back(elem);

                    polyStack.pop_back();
                    indexStack.pop_back();
                    if (!polyStack.empty()) {
                        prevd = polyStack[polyStack.size()-1]->degs;
                    } else {
                        break;
                    }
                }
            }

            if (d[k] < curd[k]) {
                isStart = 1;
                d[k] = curd[k];

            }
            // Mutiply this variate until current degree
            while (d[k] > curd[k]) {
                elem.op = 1;
                elem.type = 2;
                elem.a.i = k;
                elem.b = slp.size() - 1;
                slp.push_back(elem);
                d[k]--;

            }

            // Needed to add this result in future
            if (isStart) {
                isOp = 1;
                polyStack.push_back(&(listPoly->elems[prev]));
                indexStack.push_back(slp.size()-1);
            }
        }

        // Multiply or add this coefficient
        k = firstNonZero(d, nvar);
        if (isOp) {
            elem.op = 1;
            elem.type = 0;
            elem.a.c = new RationalNumber(listPoly->elems[i].coef);
            // elem.a  = i;
            elem.b = k;
            slp.push_back(elem);

            d[k]--;
        }
        else {
            elem.op = 0;
            elem.type = 1;
            elem.a.c = new RationalNumber(listPoly->elems[i].coef);
            // elem.a = i;
            elem.b = slp.size() - 1;
            slp.push_back(elem);
        }
        isOp = 0;

        prev = i;
    }

    int k = firstNonZero(d, nvar);
    if (k > -1) {
        if (!polyStack.empty()) {
            degree_t* prevd = polyStack[polyStack.size()-1]->degs;
            while (checkDegs(prevd, d, nvar) == k) {
                elem.op = 0;
                elem.type = 3;
                elem.a.i = indexStack[indexStack.size()-1];
                elem.b = slp.size() - 1;
                slp.push_back(elem);

                polyStack.pop_back();
                indexStack.pop_back();
                if (!polyStack.empty()) {
                    prevd = polyStack[polyStack.size()-1]->degs;
                } else {
                    break;
                }
            }
        }

        for (int i = k; i >= 0; --i) {
            if (polyStack.size() > 0) {
                degree_t* prevd = polyStack[polyStack.size()-1]->degs;
                while(checkDegs(prevd, d, nvar) == i) {
                    elem.op = 0;
                    elem.type = 3;
                    elem.a.i = indexStack[indexStack.size()-1];
                    elem.b = slp.size() - 1;
                    slp.push_back(elem);

                    polyStack.pop_back();
                    indexStack.pop_back();
                    if (!polyStack.empty()) {
                        prevd = polyStack[polyStack.size()-1]->degs;
                    } else {
                        break;
                    }
                }
            }

            while (d[i] > 0) {
                elem.op = 1;
                elem.type = 2;
                elem.a.i = i;
                elem.b = slp.size() - 1;
                slp.push_back(elem);
                d[i]--;
            }
        }
    }
    free(d);
    freePolynomial_AADL(listPoly);

    // Add the rest results
    if (!indexStack.empty()) {
        for (int i = indexStack.size()-1; i >= 0; --i) {
            elem.op = 0;
            elem.type = 3;
            elem.a.i = indexStack[i];
            elem.b = slp.size() - 1;
            slp.push_back(elem);
        }
    }
}

void SparseMultivariateRationalPolynomial::printSLP(std::ostream& out) const {
    int m = slp.size();
    if (m == 0) {
        out << "BPAS: No form for straight-line program." << std::endl;
        return;
    }

    if (isConstant()) {
        //special form of SLP;
        out << "r_0:\t" << *(slp[0].a.c) << std::endl;
        return;
    }

    std::string opStr;
    for (int i = 0; i < m; ++i) {
        if (slp[i].op == 0) {
            opStr = " + ";
        } else {
            opStr = " * ";
        }

        switch (slp[i].type) {
            //coef & variate
            case 0: {
                out << "r_" << i << ":=\t" << *(slp[i].a.c) << opStr << names[slp[i].b + 1] << ":" << std::endl;
                break;
            }
            //coef & result
            case 1: {
                out << "r_" << i << ":=\t" << *(slp[i].a.c);
                if (slp[i].b >= 0) {
                    out << opStr << "r_" << slp[i].b;
                }
                out << ":" <<  std::endl;
                break;
            }
            //variate & result
            case 2: {
                out << "r_" << i << ":=\t" << names[slp[i].a.i + 1];
                if (slp[i].b >= 0) {
                    out << opStr << "r_" << slp[i].b;
                }
                out << ":" <<  std::endl;
                break;
            }
            //result & result
            case 3: {
                out << "r_" << i << ":=\tr_" << slp[i].a.i << opStr << "r_" << slp[i].b << ":" <<  std::endl;
                break;
            }
            //coef (constant);
            case 4: {
                out << "r_" << i << ":=\t" << *(slp[i].a.c) << ":" <<  std::endl;
            }
        }
    }
}

void SparseMultivariateRationalPolynomial::sleeveBoundURPolynomials(DenseUnivariateRationalPolynomial* up, DenseUnivariateRationalPolynomial* lo, Intervals& pIs, int k, int s) {
    int m = slp.size();
    if (m == 0 || isConstant()) {
        up->setCoefficient(0,0);
        lo->setCoefficient(0,0);
        return;
    }

    int d = 0;
    if (poly != NULL) {
        d = mainLeadingDegree_AA(poly);
    }
    for (int i = 0; i < m; i++) {
        Interval a;
        if (slp[i].type == 0) {
            a.left = (*slp[i].a.c).get_mpq();
            a.right = (*slp[i].a.c).get_mpq();

            // If it is the last variate,
            // set the coefficient of bound polynomials.
            if (slp[i].b > 0) {
                if (slp[i].op)
                    intervalMultiplication(&(slp[i].res), &a, pIs.interval(k, slp[i].b));
                else
                    intervalAddition(&(slp[i].res), &a, pIs.interval(k, slp[i].b));
            }
            else {
                if (s && d%2) {
                    up->setCoefficient(d, RationalNumber(-a.left));
                    lo->setCoefficient(d, RationalNumber(-a.right));
                }
                else {
                    up->setCoefficient(d, RationalNumber(a.right));
                    lo->setCoefficient(d, RationalNumber(a.left));
                }
                d--;
                slp[i].res.left = 0;
                slp[i].res.right = 0;
            }
        }
        else if (slp[i].type == 1) {
            a.left = (*slp[i].a.c).get_mpq();
            a.right = (*slp[i].a.c).get_mpq();
            if (slp[i].op)
                intervalMultiplication(&(slp[i].res), &a, &(slp[slp[i].b].res));
            else
                intervalAddition(&(slp[i].res), &a, &(slp[slp[i].b].res));
        }
        else if (slp[i].type == 2) {
            // If it is the last variate,
            // set the coefficient of bound polynomials.
            if (slp[i].a.i > 0) {
                if (slp[i].op)
                    intervalMultiplication(&(slp[i].res), &(slp[slp[i].b].res), pIs.interval(k, slp[i].a.i));
                else
                    intervalAddition(&(slp[i].res), &(slp[slp[i].b].res), pIs.interval(k, slp[i].a.i));
            }
            else {
                if (s && d%2) {
                    up->setCoefficient(d, RationalNumber(-slp[slp[i].b].res.left));
                    lo->setCoefficient(d, RationalNumber(-slp[slp[i].b].res.right));
                }
                else {
                    up->setCoefficient(d, RationalNumber(slp[slp[i].b].res.right));
                    lo->setCoefficient(d, RationalNumber(slp[slp[i].b].res.left));
                }
                d--;
                slp[i].res.left = 0;
                slp[i].res.right = 0;
            }
        }
        else {      // slp[i].type == 3
            if (slp[i].op)
                intervalMultiplication(&(slp[i].res), &(slp[slp[i].a.i].res), &(slp[slp[i].b].res));
            else
                intervalAddition(&(slp[i].res), &(slp[slp[i].a.i].res), &(slp[slp[i].b].res));
        }
    }

    up->setCoefficient(0, slp[m-1].res.right);
    lo->setCoefficient(0, slp[m-1].res.left);
}

void SparseMultivariateRationalPolynomial::randomPolynomial(int numvar, int nterms, unsigned long int coefBound, degree_t sparsity, bool includeNeg) {

    *this = SparseMultivariateRationalPolynomial(numvar);

    if (numvar == 0) {
        mpq_t tmp;
        mpq_init(tmp);
        poly = makeConstPolynomial_AA(1, nvar, tmp);
        mpq_set_si(poly->elems->coef, (rand() % ((int) pow(2, coefBound-1))) + 1, 1);
        mpq_clear(tmp);

    } else {
        Node* node = buildRandomPoly(numvar, nterms, coefBound, sparsity, includeNeg);
        poly = deepCopyPolynomial_AAFromNode(node, nvar);
        freePolynomial(node);
    }

}

void SparseMultivariateRationalPolynomial::randomPolynomial(std::vector<int> maxDegs, unsigned long int coefBound, float sparsity, bool includeNeg) {

    AltArr_t* aa = buildRandomPolyFromMax(maxDegs.size(), maxDegs.data(), coefBound, sparsity, includeNeg);

    *this = SparseMultivariateRationalPolynomial(maxDegs.size());
    this->poly = aa;
}



/**
 * Construct an ExprTreeNode (well, multiple) which represents a single
 * term. That is, a coefficient and a monomial.
 */
ExprTreeNode* exprTreeNodeFromAAElem(const AltArr_t* n, int idx, int nvar, const Symbol* vars) {
    if (n == NULL) {
        return new ExprTreeNode(0l);
    }

    ExprTreeNode* t = new ExprTreeNode(mpq_class(n->elems[idx].coef));

    if (nvar > 0) {
        degree_t degsList[nvar];
        partialDegreesTerm_AA(n, idx, degsList);
        for (int i = 0; i < nvar; ++i) {
            degree_t deg = degsList[i];

            if (deg > 1) {
                ExprTreeNode* var = new ExprTreeNode(vars[i]);
                ExprTreeNode* num = new ExprTreeNode(deg);
                ExprTreeNode* exp = ExprTreeNode::combineExprTreeNodes(var, num, EXPR_EXP);
                t = ExprTreeNode::combineExprTreeNodes(t, exp, EXPR_MULT);
            } else if (deg == 1) {
                ExprTreeNode* var = new ExprTreeNode(vars[i]);
                t = ExprTreeNode::combineExprTreeNodes(t, var, EXPR_MULT);
            }
        }
    }

    return t;
}

ExpressionTree SparseMultivariateRationalPolynomial::convertToExpressionTree() const {
    if (isZero()) {
        ExprTreeNode* r = new ExprTreeNode(0l);
        ExpressionTree t(r);
        return t;
    }

    ExprTreeNode* prev = exprTreeNodeFromAAElem(poly, 0, nvar, &names[1]);
    for (int i = 1; i < poly->size; ++i) {
        ExprTreeNode* thisNode = exprTreeNodeFromAAElem(poly, i, nvar, &names[1]);
        prev = ExprTreeNode::combineExprTreeNodes(prev, thisNode, EXPR_ADD);
    }

    return ExpressionTree(prev);
}
