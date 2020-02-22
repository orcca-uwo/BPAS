#include "IntegerPolynomial/mzpolynomial.hpp"
#include "IntegerPolynomial/SMZP_Support_Factoring.hpp"

////////// Construtors ////////////////////////////////////

/**
 * Private constructor to create SMQP directly from a Node*.
 * Makes a copy of input varNames
 */
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(AltArrZ_t* aa, int vars, Symbol* varNames) : 
    poly(aa),
    nvar(vars) 
{
    names = new Symbol[nvar+1];
    std::copy(varNames, varNames+nvar+1, names);
}

SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(const RationalNumber& r, AltArr_t* aa, int vars, Symbol* varNames) : 
    poly(NULL),
    nvar(vars)
{
    poly = deepCopyPolynomial_AAZFromAA(aa);
    freePolynomial_AA(aa);

    names = new Symbol[nvar+1];
    std::copy(varNames, varNames+nvar+1, names);
}


/**
 * Construct a multivariate polynomial
 *
 **/
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial() : 
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
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(int v) :
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
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial (const Symbol& x) :
    poly(NULL),
    nvar(1)
{
    names = new Symbol[2];
    names[0] = "9";
    names[1] = x;

    mpz_t coef;
    mpz_init(coef);
    mpz_set_ui(coef, 1ul);
    poly = makeConstPolynomial_AAZ(1, nvar, coef);
    degree_t degsList[nvar] = {0};
    degsList[0] = 1;
    setDegrees_AAZ_inp(poly, 0, degsList, 1);
    // poly->elems->degs = 1;
    mpz_clear(coef);
}


SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial (const std::string& str) :
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
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(const SparseMultivariateIntegerPolynomial& b) :
    nvar(b.nvar)
{
    poly = deepCopyPolynomial_AAZ(b.poly);
    names = new Symbol[nvar+1];
    std::copy(b.names, b.names+nvar+1, names);
    slp = std::vector<SLPZRepresentation>(b.slp);
}

/**
 * Move Constructor.
 *
 * @params b: The r-value reference polynomial.
 */
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(SparseMultivariateIntegerPolynomial&& b) {
    nvar = b.nvar;
    poly = b.poly;
    names = new Symbol[nvar+1];
    std::copy(b.names, b.names+nvar+1, names);
    slp = b.slp;

    b.poly = NULL;
    b.slp.clear();
}

/**
 * Create an SMZP from a SMQP.
 */
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(const SparseMultivariateRationalPolynomial& b) {
    nvar = b.nvar;
    poly = deepCopyPolynomial_AAZFromAA(b.poly);
    names = new Symbol[nvar+1];
    std::copy(b.names, b.names+nvar+1, names);
}

/**
 * Create a SMQP from a Integer. 
 */
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(const Integer& r, int nv) :
    poly(NULL),
    nvar(nv)
{   
    poly = makeConstPolynomial_AAZ(1, nvar, r.get_mpz_t());

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
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial(const RationalNumber& r, int nv) :
    poly(NULL),
    nvar(nv)
{
    if (mpz_cmp_si(mpq_denref(r.get_mpq_t()), 1l) != 0) {
        std::cerr << "Invalid conversion from RationalNumber to SparseMultivariateIntegerPolynomial" << std::endl;
        std::cerr << "RationalNumber was: " << r << std::endl;
        exit(1);
    }

    poly = makeConstPolynomial_AAZ(1, nvar, mpq_numref(r.get_mpq_t()));

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
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial (const DenseUnivariateIntegerPolynomial& p) : 
    poly(NULL), 
    nvar(1)
{
    names = new Symbol[nvar+1];
    names[0] = "9";
    names[1] = p.variable();

    Integer coef;
    long int pDeg = p.degree().get_si();
    poly = makePolynomial_AAZ(pDeg+1, nvar);
    int curSize = 0;
    degree_t degsList[nvar] = {0};
    for (long int i = pDeg; i >= 0; --i) {
        coef = p.coefficient(i);
        if (coef != 0) {
            mpz_init(poly->elems[curSize].coef);
            mpz_set(poly->elems[curSize].coef, coef.get_mpz_t());
            degsList[0] = i;
            setDegrees_AAZ_inp(poly, curSize, degsList, nvar);
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
SparseMultivariateIntegerPolynomial::SparseMultivariateIntegerPolynomial (const SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial>& s) :
    poly(NULL),
    nvar(0)
{
    names = new Symbol[1];
    names[0] = "1";

    long long int d = s.degree().get_si();
    SparseMultivariateIntegerPolynomial c = s.coefficient(d);
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

        poly = makePolynomial_AAZ(d+1, nvar);
        int curSize = 0;
        degree_t degsList[nvar] = {0};
        for (int k = d; k >= 0; --k) {
            SparseMultivariateIntegerPolynomial coef = s.coefficient(k);
            if (coef.isZero()) {
                continue;
            }
            mpz_init(poly->elems[curSize].coef);
            mpz_set(poly->elems[curSize].coef, coef.poly->elems[0].coef);
            degsList[0] = k;
            setDegrees_AAZ_inp(poly, curSize, degsList, nvar);
            // poly->elems[curSize].degs = k;
            ++curSize;
        }
        poly->size = curSize;

        return;
    }

    nvar = c.nvar + 1;
    AltArrZ_t* newPoly = makePolynomial_AAZ(1, nvar);
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
        expandNumVarsLeft_AAZ(c.poly, nvar);
        c.poly = mainLShiftPolynomial_AAZ_inp(c.poly, k);
        // for (int i = 0; i < c.poly->size; ++i) {
            // c.poly->elems[i].degs |= ((degrees_t) k << mvarDegOffset);
        // }
        newPoly = addPolynomials_AAZ_inp(newPoly, c.poly, nvar);
    }

    this->poly = newPoly;
}

/**
 * Destroy the polynomial and release underlying node memory.
 **/
SparseMultivariateIntegerPolynomial::~SparseMultivariateIntegerPolynomial() {
    delete[] names;
    freePolynomial_AAZ(poly);
    slp.clear();
}



////////// BPASRing /////////////////////////////////

bool SparseMultivariateIntegerPolynomial::isZero() const {
    return isZero_AAZ(poly);
}

void SparseMultivariateIntegerPolynomial::zero() {
    freePolynomial_AAZ(poly);
    poly = NULL;
}

bool SparseMultivariateIntegerPolynomial::isOne() const {
    return isOne_AAZ(poly);
}

void SparseMultivariateIntegerPolynomial::one() {
    freePolynomial_AAZ(poly);
    Integer r(1);
    poly = makeConstPolynomial_AAZ(1, nvar, r.get_mpz_t());
}

bool SparseMultivariateIntegerPolynomial::isNegativeOne() const {
    return isNegativeOne_AAZ(poly);
}

void SparseMultivariateIntegerPolynomial::negativeOne() {
    freePolynomial_AAZ(poly);
    Integer r(-1);
    poly = makeConstPolynomial_AAZ(1, nvar, r.get_mpz_t());
}

int SparseMultivariateIntegerPolynomial::isConstant() const {
    return isConstant_AAZ(poly);
}


/**
 * Subresultant Chain
 * Return the list of subresultants
 *
 * @param q: The other sparse univariate polynomial
 **/
std::vector<SparseMultivariateIntegerPolynomial> SparseMultivariateIntegerPolynomial::subresultantChain (const SparseMultivariateIntegerPolynomial& q, int filled) const
{
    bool TYPE = 0; // 0: SMZP Ducos Algorithm, 1: SUP Ducos Algorothm
    Symbol v = q.leadingVariable();
    if (v != leadingVariable()) {
        std::cout << "BPAS: error, cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
        exit(1);
    }

    if (TYPE == 0) {
        // // Subresultant from SMZP:
        std::vector<int> xs;
        bool isOrdered = isOrderedRing(q, xs);
        SparseMultivariateIntegerPolynomial ppP, ppQ;
        std::vector<SparseMultivariateIntegerPolynomial> subr;

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
            ppP = ppP.primitivePart();
            ppQ = ppQ.primitivePart();
        } else {
            ppP = this->primitivePart();
            ppQ = q.primitivePart();
        }    

        // std::cout << "\n\n\n\n";
        // std::cout << "ppP := " << ppP << std::endl;
        // std::cout << "ppQ := " << ppQ << std::endl;

        AltArrZ_t* P = deepCopyPolynomial_AAZ (ppP.poly); //deepCopyPolynomial_AAZFromAA (ppP.poly);
        AltArrZ_t* Q = deepCopyPolynomial_AAZ (ppQ.poly); //deepCopyPolynomial_AAZFromAA (ppQ.poly);

        int lvarP = leadingVariable_AAZ (P);
        int lvarQ = leadingVariable_AAZ (Q);
        bool isShrinked = 0;

        if (lvarP != lvarQ){
            std::cout << "BPAS Error: cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
            exit(1);
        }

        if (lvarP > 0){
            for (int i = 0; i < lvarP; ++i){
                shrinkNumVarsAtIdx_AAZ (P, 0);
                shrinkNumVarsAtIdx_AAZ (Q, 0);	    
            }
            isShrinked = 1;
        }

        // if (leadingVariable_AAZ (P) != 0 || leadingVariable_AAZ (Q) != 0){
        //     std::cout << "BPAS Error: shrink in subresultantChain does not work properly!" << std::endl;
        //     exit(1);
        // }

        AltArrsZ_t* SC;
        AltArrZ_t* tmp;
        int size = 0;
        DucosSubresultantChainZ (P,Q, &SC, &size);

        freePolynomial_AAZ (P); // free
        freePolynomial_AAZ (Q); // free

        AltArrsZ_t* cur = SC;
        subr.reserve (size);
        // cerr << "size := " << size << std::endl;                        // TEST

        for (int i = 0; cur != NULL && i < size; ++i){
            tmp = deepCopyPolynomial_AAZ (cur->poly);// deepCopyPolynomial_AAFromAAZ (cur->poly);
            if (isShrinked && tmp != NULL && tmp->size != 0){
                expandNumVarsLeft_AAZ (tmp, ppP.nvar);
                if (tmp->nvar != ppP.nvar){
                    std::cout << "BPAS Error: expand in subresultantChain does not work properly!" << std::endl;
                    exit(1);
                }
            }
            subr.push_back (SparseMultivariateIntegerPolynomial (tmp, ppP.nvar, ppP.names));

            cur = cur->next;
        }

        freeAltArrsZ (SC); // free

        if (filled){
            std::vector<SparseMultivariateIntegerPolynomial> chain = subr;
            SparseMultivariateIntegerPolynomial zero;
            zero.zero();

            degree_t fullSize = mainLeadingDegree_AAZ (chain[chain.size()-2].poly) + 2;
            degree_t delta;

            // std::cerr << "chain.size() = " << chain.size() << std::endl;
            // std::cerr << "fullSize = " << fullSize << std::endl;

            if (chain.size() < fullSize){
                chain.reserve (fullSize);
                for (int i = chain.size()-2; i > 0; --i){
                    if (mainLeadingDegree_AAZ (chain[i].poly) != mainLeadingDegree_AAZ (chain[i-1].poly) + 1) {
                        delta = mainLeadingDegree_AAZ (chain[i].poly) - mainLeadingDegree_AAZ (chain[i-1].poly);
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
                if (mainLeadingDegree_AAZ (chain[0].poly) != 0){
                    for (int j = 0; j < mainLeadingDegree_AAZ (chain[0].poly); ++j){
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
        std::vector<SparseMultivariateIntegerPolynomial> Out;
        std::vector<SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial>> S;
        SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial> a, b;
        SparseMultivariateIntegerPolynomial temp;
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
            temp = SparseMultivariateIntegerPolynomial(S[i]);
            Out.push_back(temp);
        }
        return Out;
    }
}


SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::resultant (const SparseMultivariateIntegerPolynomial& q) const
{    
    Symbol v = q.leadingVariable();
    if (v != leadingVariable()) {
        std::cout << "BPAS: error, cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
        exit(1);
    }
    
    // // resultant from SMZP:
    
    std::vector<int> xs;
    bool isOrdered = isOrderedRing(q, xs);
    SparseMultivariateIntegerPolynomial ppP, ppQ;
    std::vector<SparseMultivariateIntegerPolynomial> subr;
    
    if (!isOrdered) {
        std::cout << "BPAS: error, trying to add between Q[";
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
    	// ppP = this->primitivePart();
    	// ppQ = q.primitivePart();
    }    

    // // std::cout << "\n\n\n\n";
    // // std::cout << "ppP := " << ppP << std::endl;
    // // std::cout << "ppQ := " << ppQ << std::endl;
    
    AltArrZ_t* P = deepCopyPolynomial_AAZ (ppP.poly); //deepCopyPolynomial_AAZFromAA (ppP.poly);
    AltArrZ_t* Q = deepCopyPolynomial_AAZ (ppQ.poly); //deepCopyPolynomial_AAZFromAA (ppQ.poly);
    
    AltArrZ_t* res = DucosResultantZ (P,Q);
    freePolynomial_AAZ (P);
    freePolynomial_AAZ (Q);
    return SparseMultivariateIntegerPolynomial (res, ppP.nvar, ppP.names);
}

/**
 * Subresultant Chain GCD
 * Return the last non-zero subresultant of the current polynomial and the input polynomial if the resultant is zero and return 1 otherwise
 *
 * @param q: The other sparse univariate polynomial
 **/
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::subresultantGCD (const SparseMultivariateIntegerPolynomial& q) const {
    Symbol v = q.leadingVariable();
    SparseMultivariateIntegerPolynomial one;
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
//  startTimer(&start);
    std::vector<SparseMultivariateIntegerPolynomial> src = subresultantChain(q);
    
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
    
}

/**
 * Get the GCD between this and b.
 * If this and b have all integer coefficients, the gcd will have integer coefficients
 * with proper GCD among those coefficients. Otherwise, the returned GCD is monic
 * with rational number coefficients. 
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::gcd(const SparseMultivariateIntegerPolynomial& b) const {
    
    bool TYPE = 0; // 0: gcd_AAZ, 1: primitiveGCD
    SparseMultivariateIntegerPolynomial ret;
    
    if (TYPE == 0){
    	std::vector<int> xs;
    	bool isOrdered = isOrderedRing(b, xs);
        SparseMultivariateIntegerPolynomial ppP, ppQ;

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
    	    ppP = ppP.primitivePart();
    	    ppQ = ppQ.primitivePart();
    	} else {
    	    ppP = this->primitivePart();
    	    ppQ = b.primitivePart();
    	}    

    	// std::cout << "\n\n\n\n";
    	// std::cout << "ppP := " << ppP << std::endl;
    	// std::cout << "ppQ := " << ppQ << std::endl;
    	
    	AltArrZ_t* P =  deepCopyPolynomial_AAZ (ppP.poly);// deepCopyPolynomial_AAZFromAA (ppP.poly);
    	AltArrZ_t* Q = deepCopyPolynomial_AAZ (ppQ.poly); //deepCopyPolynomial_AAZFromAA (ppQ.poly);
        int Pnvar = P == NULL ? 0 : P->nvar;

#if defined(WITH_MAPLE) && WITH_MAPLE
        AltArrZ_t* gz = gcd_AAZ(P,Q);
        // AltArrZ_t* gz = NULL;
        // if (Pnvar > 1) {
        //     if (P->size < 10 || Q->size < 10) {
        //         gz = gcd_AAZ (P,Q);
        //     } else {
        //         char* c_names[Pnvar];
        //         for (int i = 0; i < Pnvar; ++i) {
        //             std::string str = ppP.names[i+1].toString();
        //             c_names[i] = (char*) malloc(sizeof(char)*str.length()+1);
        //             strcpy(c_names[i], str.c_str());
        //         } 
        //         gz = gcd_MplInt_AAZ(P, Q, (const char**) c_names);
        //         for (int i = 0; i < nvar; ++i) {
        //             free(c_names[i]);
        //         }
        //     }
        // } else if (Pnvar == 1) {
        //     gz = univariateGCD_AAZ(P, Q);
        // } else {
        //     gz = makePolynomial_AAZ(1, 0);
        //     mpz_init(gz->elems->coef);
        //     mpz_set_si(gz->elems->coef, 1l);
        //     gz->size = 1;
        // }
#elif defined(WITH_BLAD) && WITH_BLAD
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
    	freePolynomial_AAZ (P);
    	freePolynomial_AAZ (Q);
    	
    	SparseMultivariateIntegerPolynomial gcd (gz, ppP.nvar, ppP.names);
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
    integralContent_AAZ(this->poly, g1);
    integralContent_AAZ(b.poly, g2);
    mpz_gcd(g1, g1, g2);
    ret *= Integer(g1);
    
    mpz_clear(g1);
    mpz_clear(g2);
    return ret;
}


/**
 * Get GCD between *this and b as a primitive polynomial.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::primitiveGCD(const SparseMultivariateIntegerPolynomial& b) const {
    if (isZero()) {
        return b;
    }
    if (b.isZero()) {
        return *this;
    }

    if (isConstant() || b.isConstant()) {
        SparseMultivariateIntegerPolynomial ret;
        ret.one();
        return ret;
    }

    std::vector<int> xs;
    if (!isOrderedRing(b, xs)) {
        SparseMultivariateIntegerPolynomial ret(0);
        ret.one();
        return ret;
    }    

    if (this->nvar == 1 && b.nvar == 1) {
        if (xs.size() > 2) {
            SparseMultivariateIntegerPolynomial ret(0);
            ret.one();
            return ret;
        }

        //otherwise they are the same variables
        AltArrZ_t* g = univariateGCD_AAZ(this->poly, b.poly);
        if (isConstant_AAZ(g)) {
            SparseMultivariateIntegerPolynomial ret(g, 0, names);
            return ret;
        }
        SparseMultivariateIntegerPolynomial ret(g, nvar, names);
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
            SparseMultivariateIntegerPolynomial ret(0);
            ret.one();
            return ret;
        }

        vars.erase(vars.begin()+match);

        SparseMultivariateIntegerPolynomial newB = b.content(vars);
        
        //it is possible that the content that comes back is an integral content.
        if (newB.nvar == 0) {
            //recursive call to handle is isContstant() case above
            return this->primitiveGCD(newB);
        }

        AltArrZ_t* g = univariateGCD_AAZ(this->poly, newB.poly);
        if (isConstant_AAZ(g)) {
            SparseMultivariateIntegerPolynomial ret(g, 0, names);
            return ret;
        }
        SparseMultivariateIntegerPolynomial ret(g, nvar, names);
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
            SparseMultivariateIntegerPolynomial ret(0);
            ret.one();
            return ret;
        }

        vars.erase(vars.begin()+match);

        SparseMultivariateIntegerPolynomial newThis = this->content(vars);

        //it is possible that the content that comes back is an integral content.
        if (newThis.nvar == 0) {
            //recursive call to handle is isContstant() case above
            return b.primitiveGCD(newThis);
        }
        AltArrZ_t* g = univariateGCD_AAZ(newThis.poly, b.poly);
        if (isConstant_AAZ(g)) {
            SparseMultivariateIntegerPolynomial ret(g, 0, names);
            return ret;
        }
        SparseMultivariateIntegerPolynomial ret(g, b.nvar, b.names);
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

        SparseMultivariateIntegerPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateIntegerPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);
        
        std::vector<Symbol> vars;
        vars.push_back(newnames[1]);
        SparseMultivariateIntegerPolynomial tC = tempa.content(vars);
        SparseMultivariateIntegerPolynomial bC = tempb.content(vars);
        SparseMultivariateIntegerPolynomial g = tC.primitiveGCD(bC);

        SparseMultivariateIntegerPolynomial tPP = tempa / tC;
        SparseMultivariateIntegerPolynomial bPP = tempb / bC;

        SparseMultivariateIntegerPolynomial ppG = tPP.subresultantGCD(bPP);
        SparseMultivariateIntegerPolynomial ppGP = ppG.primitivePart(vars);

        SparseMultivariateIntegerPolynomial ret = g * ppGP;

        ret.setRingVariables(ret.variables());
        return ret;    
    } else {

        //exact same ring, easy.
        std::vector<Symbol> vars;
        vars.push_back(names[1]);
        SparseMultivariateIntegerPolynomial tC = this->content(vars);
        SparseMultivariateIntegerPolynomial bC = b.content(vars);

        std::vector<Symbol> rVars;

        SparseMultivariateIntegerPolynomial g = tC.primitiveGCD(bC);

        SparseMultivariateIntegerPolynomial tPP = *this / tC;
        SparseMultivariateIntegerPolynomial bPP = b / bC;

        SparseMultivariateIntegerPolynomial ppG = tPP.subresultantGCD(bPP);
        SparseMultivariateIntegerPolynomial ppGP = ppG.primitivePart(vars);

        SparseMultivariateIntegerPolynomial ret = g * ppGP;

        ret.setRingVariables(ret.variables());
        return ret;
    }

}

Factors<SparseMultivariateIntegerPolynomial> SparseMultivariateIntegerPolynomial::squareFree() const {
    AltArrZ_t** facts = NULL;
    int* exps = NULL;
    mpz_t cont;
    mpz_init(cont);
    int nfacts = 0;
    
    squareFree_AAZ(this->poly, cont, &facts, &exps, &nfacts);

    Factors<SparseMultivariateIntegerPolynomial> sf;
    sf.setRingElement(SparseMultivariateIntegerPolynomial(cont));
    mpz_clear(cont);

    for (int i = 0; i < nfacts; ++i) {
        sf.addFactor(std::move(SparseMultivariateIntegerPolynomial(facts[i], nvar, names)), exps[i]);
    }

    return sf;

    // return this->squareFree(this->ringVariables());
}

Factors<SparseMultivariateIntegerPolynomial> SparseMultivariateIntegerPolynomial::squareFree(const std::vector<Symbol>& vars) const {
    Factors<SparseMultivariateIntegerPolynomial> sf;
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
    SparseMultivariateIntegerPolynomial content;
    SparseMultivariateIntegerPolynomial primPart = this->primitivePart(curVars, content);

    if (primPart.isConstant()) {
        if (!primPart.isOne()) {
            sf.addFactor(primPart, 1);
        }
    } else if (primPart.degree(names[strIdx]) == 1) {
        // SparseMultivariateIntegerPolynomial df = primPart.leadingCoefficient(names[strIdx]);
        // SparseMultivariateIntegerPolynomial g = primPart.gcd(df);
        // SparseMultivariateIntegerPolynomial next = primPart / g;
        // next *= g;
        // sf.push_back(next); 

        sf.addFactor(primPart, 1); 
    } else {
        SparseMultivariateIntegerPolynomial dx = primPart.derivative(names[strIdx]);
        SparseMultivariateIntegerPolynomial g = (primPart.gcd(dx)).primitivePart();
        SparseMultivariateIntegerPolynomial next = primPart / g;

        if (next.leadingCoefficient() < 1) {
            next.negate();
        }

        int k = 1;
        if (g.isOne()) {
            AltArrZ_t* remFact = NULL;
            AltArrZ_t* comFact = commonFactor_AAZ(next.poly, &remFact);

            if (!isConstant_AAZ(comFact)) {
                //next /= comFact
                AltArrZ_t* temp = next.poly;
                next.poly = remFact;
                freePolynomial_AAZ(temp);

                //set comFact exp to 1 to add to sf.
                temp = makeConstPolynomial_AAZ(1, nvar, Integer(1).get_mpz_t());
                degree_t tempDegs[nvar] = {0};
                degree_t deg = partialDegreeTerm_AAZ(comFact, 0, 0);
                SparseMultivariateIntegerPolynomial fac(temp, nvar, names);
                if (deg != 0) {
                    tempDegs[0] = 1;
                    setDegrees_AAZ_inp(temp, 0, tempDegs, nvar);
                    sf.addFactor(fac, deg);
                }
                for (int j = 1; j < nvar; ++j){
                    deg = partialDegreeTerm_AAZ(comFact, 0, j);
                    if (deg == 0) {
                        continue;
                    }
                    tempDegs[j-1] = 0;
                    tempDegs[j] = 1;
                    setDegrees_AAZ_inp(temp, 0, tempDegs, nvar);
                    sf.addFactor(fac, deg);
                }
            }
        }
        while (g.degree(names[strIdx]) > 0) {
            SparseMultivariateIntegerPolynomial y = next.gcd(g);
            y = y.primitivePart();
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
        Factors<SparseMultivariateIntegerPolynomial> contSf = content.squareFree(nextVars);
        sf.addFactors(contSf);
        sf.multiplyRingElement(contSf.ringElement());
    // }

    return sf;
}

/**
 * Computes the square free part of this polynomail. That is, the polynomial of this
 * divided by all square factors. This is with respect to all variables.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::squareFreePart() const {
    std::vector<Symbol> vars = this->ringVariables();
    return this->squareFreePart(vars);
}

/**
 * Computes the square free part of this polynomail. That is, the polynomial of this
 * divided by all square factors. This is with respect to all variables.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::squareFreePart(std::vector<Symbol>& vars) const {
    std::vector<Symbol> matchedVars;
    for (int i = 0; i < vars.size(); ++i) {
        for (int j = 0; j <= nvar; ++j) {
            if (vars[i] == names[j+1]) {
                matchedVars.push_back(vars[i]);
                break;
            }
        }
    }

    SparseMultivariateIntegerPolynomial cont;
    SparseMultivariateIntegerPolynomial fact;
    fact.one();
    std::vector<Symbol> var;
    SparseMultivariateIntegerPolynomial primpart = *this;
    for (int i = 0; i < matchedVars.size(); ++i) {
        var.clear();
        var.push_back(matchedVars[i]);
        primpart = primpart.primitivePart(var, cont);
        SparseMultivariateIntegerPolynomial diff = primpart.derivative(var[0]);
        SparseMultivariateIntegerPolynomial g = (primpart.gcd(diff)).primitivePart();
        SparseMultivariateIntegerPolynomial next = primpart / g;
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
bool SparseMultivariateIntegerPolynomial::isOrderedRing(const SparseMultivariateIntegerPolynomial& b, std::vector<int>& xs) const {
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
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::expandVariables(int vars, Symbol* newvars, int varmap[]) const {
    if (isZero()) {
        SparseMultivariateIntegerPolynomial temp(NULL, vars, newvars);
        return temp;
    }
    if (isConstant()) {
        AltArrZ_t* tempPoly = makeConstPolynomial_AAZ(1, vars, poly->elems->coef);
        SparseMultivariateIntegerPolynomial temp(tempPoly, vars, newvars);
        return temp;
    }


    AltArrZ_t* tempPoly = deepCopyPolynomial_AAZ(poly);
    expandNumVars_AAZ(tempPoly, vars);
    reorderVars_AAZ(tempPoly, varmap, nvar);

    SparseMultivariateIntegerPolynomial temp(tempPoly, vars, newvars);       
    return temp;
}

void SparseMultivariateIntegerPolynomial::expandVarsInPlace(int vars, Symbol* newvars, int varmap[]) {
    if (isZero() || isConstant()) {
        if (poly != NULL) {
            poly->nvar = vars;
        }
    } else {
        expandNumVars_AAZ(poly, vars);
        
        //here we use nvar and vars as # of vars as varmap is of size nvar.
        reorderVars_AAZ(poly, varmap, nvar);
    }

    nvar = vars;
    delete[] names;
    names = new Symbol[vars+1];
    std::copy(newvars, newvars+vars+1, names);
}

/** 
 * Rearrange exponent vectors in place and then re-sort the polynomial.
 */
void SparseMultivariateIntegerPolynomial::reorderVarsInPlace(int varmap[]) {
    if (isZero()) {
        return;
    }
    if (isConstant()) {
        return;
    }
    
    reorderVars_AAZ(poly, varmap, nvar);
    Symbol* newvars = new Symbol[nvar+1];
    newvars[0] = names[0];
    for (int i = 0; i < nvar; ++i) {
        newvars[varmap[i]+1] = names[i+1]; 
    }
    delete[] names;
    names = newvars;

} 

////////// BPASPolynomial ////////////////////////////////////

SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator= (const SparseMultivariateIntegerPolynomial& b) {
    if (this != &b) {
        freePolynomial_AAZ(poly);
        poly = deepCopyPolynomial_AAZ(b.poly);
        
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
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator= (SparseMultivariateIntegerPolynomial&& b) {
    if (this != &b) {
        freePolynomial_AAZ(poly);
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

SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator= (const Integer& r) {
    if (poly != NULL) {
        freePolynomial_AAZ(poly);
        poly = NULL;
    }
    
    slp.clear();
    poly = makeConstPolynomial_AAZ(1, nvar, r.get_mpz_t());

    return *this;
}


SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator+ (const SparseMultivariateIntegerPolynomial& b) const {
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

        SparseMultivariateIntegerPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateIntegerPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        AltArrZ_t* sum = addPolynomials_AAZ(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateIntegerPolynomial ret(sum, superNvar, newnames);
        return ret;
    
    } else {
        AltArrZ_t* sum = addPolynomials_AAZ(poly, b.poly, nvar);
        SparseMultivariateIntegerPolynomial ret(sum, nvar, names);
        return ret;
    }
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator+ (SparseMultivariateIntegerPolynomial&& b) const {
    SparseMultivariateIntegerPolynomial ret = b;
    ret += *this;
    return ret;
}

// SparseMultivariateIntegerPolynomial operator+ (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b) {
//     SparseMultivariateIntegerPolynomial ret = a;
//     ret += b;
//     return ret;
// }

SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator+= (const SparseMultivariateIntegerPolynomial& b) {
    if (nvar == b.nvar) {
        std::vector<int> xs;
        bool isOrdered = isOrderedRing(b, xs);
        int superNvar = xs.size() / 2;
        if (isOrdered && superNvar == nvar) {
            this->poly = addPolynomials_AAZ_inp(this->poly, b.poly, nvar); 
            return *this;
        }
    }

    *this = (*this + b);
    return *this;
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator- () const {
    SparseMultivariateIntegerPolynomial temp = *this;
    temp.negate();
    return temp;
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator- (const SparseMultivariateIntegerPolynomial& b) const {
    if (b.isZero()) {
        return *this;
    } 
    if (isZero()) {
        return -b;
    }
    if (this->isConstant() != 0) {
        SparseMultivariateIntegerPolynomial negB = -b;
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

        SparseMultivariateIntegerPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateIntegerPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        AltArrZ_t* sum = subPolynomials_AAZ(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateIntegerPolynomial ret(sum, superNvar, newnames);
        return ret;
    
    } else {
        AltArrZ_t* sum = subPolynomials_AAZ(poly, b.poly, nvar);
        SparseMultivariateIntegerPolynomial ret(sum, nvar, names);
        return ret;
    }
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator- (SparseMultivariateIntegerPolynomial&& b) const {
    SparseMultivariateIntegerPolynomial ret = b;
    ret.negate();
    ret += *this;
    return ret;   
}

// SparseMultivariateIntegerPolynomial operator- (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b) {
//     SparseMultivariateIntegerPolynomial ret = a;
//     ret -= b;
//     return ret;
// }

SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator-= (const SparseMultivariateIntegerPolynomial& b) {
    if (nvar == b.nvar) {
        std::vector<int> xs;
        bool isOrdered = isOrderedRing(b, xs);
        int superNvar = xs.size() / 2;
        if (isOrdered && superNvar == nvar) {
            this->poly = subPolynomials_AAZ_inp(this->poly, b.poly, nvar); 
            return *this;
        }
    }

    *this = (*this - b);
    return *this;
}

/**
 * Multiply *this by the specified polynomail.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator* (const SparseMultivariateIntegerPolynomial& b) const {
    if (b.isZero() || this->isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
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

        SparseMultivariateIntegerPolynomial tempa = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateIntegerPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);
        AltArrZ_t* prod = multiplyPolynomials_AAZ(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateIntegerPolynomial ret(prod, superNvar, newnames);
        return ret;
    
    } else {
        AltArrZ_t* prod = multiplyPolynomials_AAZ(poly, b.poly, nvar);
        SparseMultivariateIntegerPolynomial ret(prod, nvar, names);    
        return ret;
    }
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator* (SparseMultivariateIntegerPolynomial&& b) const {
    SparseMultivariateIntegerPolynomial ret = b;
    ret *= *this;
    return ret;   
}

// SparseMultivariateIntegerPolynomial operator* (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b) {
//     SparseMultivariateIntegerPolynomial ret = a;
//     ret *= b;
//     return ret;
// }

/**
 * Update this by multiplying by the specified polynomail.
 */
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator*= (const SparseMultivariateIntegerPolynomial& b) {
    *this = (*this * b);
    return *this;
}

/**
 * Divide *this by the specified polynomial.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator/ (const SparseMultivariateIntegerPolynomial& b) const {
    std::vector<Symbol> vars = this->ringVariables();
    SparseMultivariateIntegerPolynomial q,r;

    this->divide(b, q, r);
    if (!r.isZero()) {
        std::cerr << "BPAS ERROR: SMQP non-exact division." << std::endl;
        std::cerr << "dividend: " << *this << std::endl;
        std::cerr << "divisor: " << b << std::endl;
        std::cerr << "quoteint: " << q << std::endl;
        std::cerr << "remainder: " << r << std::endl;
        exit(1);
    }

    return q;
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator/ (SparseMultivariateIntegerPolynomial&& b) const {
    SparseMultivariateIntegerPolynomial ret = b;
    ret /= *this;
    return ret;   
}

// SparseMultivariateIntegerPolynomial operator/ (SparseMultivariateIntegerPolynomial&& a, const SparseMultivariateIntegerPolynomial& b) {
//     SparseMultivariateIntegerPolynomial ret = a;
//     ret /= b;
//     return ret;
// }

/**
 * Update *this by dividing by the specified polynomial.
 */
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator/= (const SparseMultivariateIntegerPolynomial& b){
    *this = (*this / b);
    return *this;
}

/**
 * Exponentiate *this by the input exponent integer.
 * Treats negative exponents as positive.
 */ 
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator^ (long long int e) const {
    if (e == 0) {
        SparseMultivariateIntegerPolynomial ret(NULL, nvar, names);
        ret.one();
        return ret;
    } 

    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    if (e == 1) {
        return SparseMultivariateIntegerPolynomial(*this);
    }

    e = (e < 0) ? -e : e;
    AltArrZ_t* retPoly = exponentiatePoly_AAZ(poly, e, nvar);
    SparseMultivariateIntegerPolynomial ret(retPoly, nvar, names);
    return ret;
}

/**
 * Update *this by exponentiating this to the input integer.
 * Treats negative exponents as positive.
 */
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator^= (long long int e) {
    *this = (*this ^ e);
    return *this;
}

/**
 * Determine if *this is equal to the specified polynomial.
 * This takes into account the variable ordering on both poylnomials 
 * in such a way that the same polynomial under different variable orderings
 * are NOT equal.
 */ 
bool SparseMultivariateIntegerPolynomial::operator== (const SparseMultivariateIntegerPolynomial& b) const {
    return this->isEqual(b);
}

/**
 * Determine if *this is not equal to the specified polynomial.
 * This takes into account the variable ordering on both poylnomials 
 * in such a way that the same polynomial under different variable orderings
 * are NOT equal.
 */
bool SparseMultivariateIntegerPolynomial::operator!= (const SparseMultivariateIntegerPolynomial& b) const {
    return !this->isEqual(b);
}

/**
 * Output the string representation of *this to the input ostream.
 */
void SparseMultivariateIntegerPolynomial::print(std::ostream& os) const {
    if (this->poly != NULL) {
        std::string tempVars[this->nvar];
        for (int i = 0 ; i < this->nvar; ++i) {
            tempVars[i] = this->names[i+1].toString();
        }
        os << polyToString_AAZ(poly, tempVars);
    } else {
        os << "0";
    }

}

std::istream& operator>>(std::istream& in, SparseMultivariateIntegerPolynomial& p) {
    std::string str;
    getline(in, str);

    p.fromString(str);
    return in;
}

void SparseMultivariateIntegerPolynomial::fromString(const std::string& str) {

    altarr_pack* pack = generate_altarr_pack(str.c_str());

    freePolynomial_AAZ(this->poly);
    this->poly = deepCopyPolynomial_AAZFromAA(pack->altarr_t_data);
    
    delete[] this->names;
    this->names = new Symbol[pack->numVars + 1];
    for (int i = 0; i < pack->numVars; ++i) {
        names[i+1] = Symbol(std::string(pack->vars[i]));
    }

    this->nvar = pack->numVars;

    freePolynomial_AA(pack->altarr_t_data);
    free(pack->vars);
    free(pack);
}

Integer SparseMultivariateIntegerPolynomial::content() const {
    if (isZero()) {
        return Integer(0);
    }
    if (isConstant()) {
        Integer ret(poly->elems->coef);
        return ret;
    }

    mpz_t ret;
    mpz_init(ret);
    integralContent_AAZ(poly, ret);
    Integer in(ret);
    mpz_clear(ret);
    return in;
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::content(const std::vector<Symbol>& v) const {
    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    bool allFound = 1;
    std::vector<Symbol> matchingSyms;
    for (int i = 0; i < nvar; ++i) {
        bool found = 0;
        for (int j = 0; j < v.size(); ++j) {
            if (names[i+1] == v[j]) {
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


    SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial> sup = this->convertToSUP(matchingSyms[0]);
    if (sup.numberOfTerms() <= 1 && matchingSyms.size() == 1) {
        //in this case, don't multiply by numerical content.
        //it is sufficient to return just the (recursive) coefficient in this
        //case where the sup has only one term;
        SparseMultivariateIntegerPolynomial ret = sup.leadingCoefficient();
        return ret;
    }

    SparseMultivariateIntegerPolynomial content = sup.content().primitivePart();
    if (content.isOne()) {
        return this->content();
    }

    for (int i = 1; i < matchingSyms.size(); ++i) {  
        sup = this->convertToSUP(matchingSyms[i]);
        SparseMultivariateIntegerPolynomial temp = sup.content().primitivePart();
        content = content.gcd(temp);
        content = content.primitivePart();
        if (content.isOne()) {
            return this->content();
        }
    }

    content *= this->content();
    return content;
}


SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::primitivePart() const {
    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    AltArrZ_t* pp = primitivePart_AAZ(this->poly);
    return SparseMultivariateIntegerPolynomial(pp, nvar, names);
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::primitivePart(Integer& content) const {
    if (isZero()) {
        content.zero();
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    mpz_t cont;
    mpz_init(cont);
    AltArrZ_t* pp = primitivePartAndContent_AAZ(this->poly, cont);
    content = Integer(cont);
    mpz_clear(cont);

    return SparseMultivariateIntegerPolynomial(pp, nvar, names);
}


SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::primitivePart(const std::vector<Symbol>& v) const {
    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    SparseMultivariateIntegerPolynomial cont = this->content(v);
    return (*this / cont);
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::primitivePart(const std::vector<Symbol>& v, SparseMultivariateIntegerPolynomial& content) const {
    if (isZero()) {
        content.zero();
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    content = this->content(v);
    return (*this / content);
}



////////// BPASMultivariatePolynomial ////////////////////////////////////

/**
 * Get the number of variables in this polynomial.
 */
int SparseMultivariateIntegerPolynomial::numberOfVariables() const {
    int foundVars[nvar] = {0};
    nonZeroVariables_AAZ(poly, foundVars);
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
Integer SparseMultivariateIntegerPolynomial::numberOfTerms() const {
    if (poly != NULL) {
        return poly->size;
    }
    return 0;
}

/** 
 * Total degree.
 */
Integer SparseMultivariateIntegerPolynomial::degree() const {
    if (isZero()) {
        return -1;
    }
    return totalDegree_AAZ(poly);
}

/**
 * Get the degree of a variable 
 */
Integer SparseMultivariateIntegerPolynomial::degree(const Symbol& str) const {
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
        return partialDegreeTerm_AAZ(poly, 0, 0);
    } else {
        return partialDegree_AAZ(poly, strIdx);
    }
}

/**
 * Get the leading coefficient 
 */
Integer SparseMultivariateIntegerPolynomial::leadingCoefficient() const {
    if (isZero()) {
        return Integer(0);
    }

    return Integer(poly->elems->coef);
}

void SparseMultivariateIntegerPolynomial::setLeadingCoefficient(const Integer& it) {
    if (poly == NULL) {
        poly = makeConstPolynomial_AAZ(1, nvar, it.get_mpz_t());
    } else {
        mpz_set(poly->elems->coef, it.get_mpz_t());
    }
}


Integer SparseMultivariateIntegerPolynomial::trailingCoefficient() const {
    if (isZero()) {
        return Integer(0);
    }

    return Integer(poly->elems[poly->size-1].coef);
}


/**
 * Get a coefficient, given the exponent of each variable in d.
 * v is the number of variables in d. It is assumed that the first this.nvar 
 * variables of d match the variables of this polynomial
 */
Integer SparseMultivariateIntegerPolynomial::coefficient(int v, const int* d) const {
    if (v < nvar) {
        std::cerr << "BPAS ERROR: SMQP calling coefficient without enough variables." << std::endl;
        exit(1);
    }

    Integer ret;
    coefficient_AAZ(poly, d, v, ret.get_mpz_t());
    return ret;
}

/**
 * Set a coefficient, given the exponent of each variable
 */
void SparseMultivariateIntegerPolynomial::setCoefficient(int v, const int* d, const Integer& it) {
    if (v != nvar) {
        std::cout << "BPAS: error, SMQP(" << nvar << "), but trying to setCoefficient with " << v << " variables." << std::endl;
        exit(1);
    }

    if (poly == NULL) {
        poly = makeConstPolynomial_AAZ(1, nvar, it.get_mpz_t());
        if (nvar != 0) {
            setDegrees_AAZ_inp(poly, 0, d, nvar);
        }
        return;
    }

    setCoefficient_AAZ(poly, d, v, it.get_mpz_t());
} 

/**
 * Set variables' name.
 */
void SparseMultivariateIntegerPolynomial::setRingVariables (const std::vector<Symbol>& xs) {
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
                poly->nvar = ns;
            }
        }
        return;
    }

    if (ns == 0 && !(isZero() || isConstant())) {
        std::cerr << "BPAS ERROR: SMPQ: Trying to remove all variables from a non-constant polynomial." << std::endl;
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
                    std::cerr << "BPAS ERROR: SMZP trying to remove variable " << names[i+1] << " which is non-zero in this polynomial." << std::endl;
                    exit(1);
                }
            }
        }
    }

    shrinkAndReorderVars_AAZ(poly, varmap, nvar);

    nvar = ns;
    delete[] names;
    names = new Symbol[nvar+1];
    names[0] = "9";
    std::copy(xs.begin(), xs.end(), names+1);
}

/**
 * Get variable names of variables with non-zero degree;
 */
std::vector<Symbol> SparseMultivariateIntegerPolynomial::variables() const {
    if (nvar == 0 || isConstant()) {
        return std::vector<Symbol>();
    }

    int foundVar[nvar] = {0};
    nonZeroVariables_AAZ(poly, foundVar);

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
std::vector<Symbol> SparseMultivariateIntegerPolynomial::ringVariables() const {
    std::vector<Symbol> varNames;
    for (int i = 0; i < nvar; ++i) {
        varNames.push_back(names[i+1]);
    }
    return varNames;   
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::derivative(const Symbol& s, int k) const {
    if (k <= 0) {
        return *this;
    }

    if (isZero() || isConstant()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);;
    }

    int idx = -1;
    for (int i = 0; i < nvar; ++i) {
        if (s == names[i+1]) {
            idx = i;
            break;
        }
    }

    if (idx == -1) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    AltArrZ_t* temp = derivative_AAZ(poly, idx, k);
    return SparseMultivariateIntegerPolynomial(temp, nvar, names);
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::integral(const Symbol& s, int k) const {
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
        AltArr_t* temp = integral_AAZ(poly, idx, k);
        return SparseMultivariateIntegerPolynomial(RationalNumber(1), temp, nvar+1, newnames);
    }

    AltArr_t* temp = integral_AAZ(poly, idx, k);
    return SparseMultivariateIntegerPolynomial(RationalNumber(1), temp, nvar, names);
}

SparseMultivariateRationalPolynomial SparseMultivariateIntegerPolynomial::rationalIntegral(const Symbol& s) const {
    return this->rationalIntegral(s, 1);
}

SparseMultivariateRationalPolynomial SparseMultivariateIntegerPolynomial::rationalIntegral(const Symbol& s, int k) const {
    if (k <= 0) {
        AltArr_t* aa = deepCopyPolynomial_AAFromAAZ(poly);
        return SparseMultivariateRationalPolynomial(aa, nvar, names);
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
        AltArr_t* temp = integral_AAZ(poly, idx, k);
        return SparseMultivariateRationalPolynomial(temp, nvar+1, newnames);
    }

    AltArr_t* temp = integral_AAZ(poly, idx, k);
    return SparseMultivariateRationalPolynomial(temp, nvar, names);
}



////////// RecursivelyViewedPolynomial  ////////////////////////////////////

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::initial() const {
    if (nvar == 0 || isConstant() || isZero()) {
        return *this;
    }  
    return leadingCoefficientInVariable(leadingVariable());
}

int SparseMultivariateIntegerPolynomial::mainDegree() const {
    if (isZero()) {
        return -1;
    }
    if (isConstant()) {
        return 0;
    }

    return mainDegree_AAZ(poly);
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::rank() const {
    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    } 

    if (isConstant()) {
        SparseMultivariateIntegerPolynomial ret(NULL, nvar, names);
        ret.one();

        return ret;
    }

    int lvIdx = mainVariable_AAZ(poly);
    Symbol lv = names[lvIdx + 1];

    SparseMultivariateIntegerPolynomial ret(lv);
    degree_t deg = this->mainDegree();
    setDegrees_AAZ_inp(ret.poly, 0, &deg, 1);

    return ret;
}


SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::head() const {
    if (nvar == 0 || isZero() || isConstant()) {
        return *this;
    }

    SparseMultivariateIntegerPolynomial tmp = leadingCoefficientInVariable(leadingVariable());
    SparseMultivariateIntegerPolynomial rank = this->rank();
    tmp *= rank;
    tmp.setRingVariables(this->ringVariables());
    return tmp;
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::tail() const {
     if (nvar == 0 || isZero() || isConstant()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    return *this - head();
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::separant() const {
    std::cerr << "BPAS ERROR: SparseMultivariateIntegerPolynomial::separant NOT YET IMPLEMENTED" << std::endl;
    return *this;
}



////////// SMQP-Specific ////////////////////////////////////

bool SparseMultivariateIntegerPolynomial::isEqual(const SparseMultivariateIntegerPolynomial& b) const {
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
            return mpz_cmp(poly->elems->coef, b.poly->elems->coef) == 0;
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

    int ret = isEqualWithVariableOrdering_AAZ(poly, b.poly, xs.data(), xs.size());

    return ret;
}

/**
 * Evaluate *this polynomial given the variables and their values in vars, and values.
 * vars must be a list of variables which are a (not necessarily proper) subset.
 * vars and values should have matching indices. i.e. the values of vars[i] is values[i].
 *
 * returns a new SMQP where all variables in vars have been evaluated using the values given. 
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::evaluate(const std::vector<Symbol>& vars, const std::vector<Integer>& values) const {

    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }
    if (isConstant()) {
        return SparseMultivariateIntegerPolynomial(*this);
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

        mpz_t res;
        mpz_init(res);
        univarEvaluate_AAZ(poly, values[0].get_mpz_t(), res);
        Integer r(res);
        SparseMultivariateIntegerPolynomial ret(r);
        mpz_clear(res);
        return ret;
    }

    int active[nvar];
    mpz_t vals[nvar];
    for (int i = 0; i < nvar; ++i) {
        active[i] = 0;
        mpz_init(vals[i]);
    }

    int newNvar = nvar;
    for (int i = 0; i < vars.size(); ++i) {
        bool found = 0;
        for (int j = 0; j < nvar; ++j) {
            if (names[j+1] == vars[i]) {
                found = 1;

                --newNvar;    
                active[j] = 1;
                mpz_set(vals[j], values[i].get_mpz_ref().get_mpz_t());
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

    AltArrZ_t* evalPoly = evaluatePoly_AAZ(poly, active, vals, nvar);
    return SparseMultivariateIntegerPolynomial(evalPoly, newNvar, newnames);
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::interpolate(const std::vector<std::vector<Integer>>& points, const std::vector<Integer>& vals) {
    std::vector<std::vector<RationalNumber>> qPoints;
    std::vector<RationalNumber> qVals;

    for (auto v : points) {
        std::vector<RationalNumber> singlePoint;
        for (auto w : v) {
            singlePoint.emplace_back(w);
        }
        qPoints.push_back(singlePoint);
    }
    for (auto v : vals) {
        qVals.emplace_back(v);
    }
    return SparseMultivariateIntegerPolynomial::interpolate(qPoints, qVals);
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::interpolate(const std::vector<std::vector<RationalNumber>>& points, const std::vector<RationalNumber>& vals) {
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

        SparseMultivariateIntegerPolynomial ret(1);
        ret.poly = deepCopyPolynomial_AAZFromAA(interp);
        freePolynomial_AA(interp);

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
bool SparseMultivariateIntegerPolynomial::divide(const SparseMultivariateIntegerPolynomial& b, SparseMultivariateIntegerPolynomial& q, SparseMultivariateIntegerPolynomial& r) const {
    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
        exit(1);
    }

    if (isZero()) {
        q = SparseMultivariateIntegerPolynomial(NULL, nvar, names);
        r = q;
        return true;           
    }

    if (b.isConstant()) {
        if (isConstant()) {
            if (mpz_divisible_p(this->poly->elems->coef, b.poly->elems->coef)) {
                AltArrZ_t* ret = deepCopyPolynomial_AAZ(this->poly);
                mpz_divexact(ret->elems->coef, ret->elems->coef, b.poly->elems->coef);
                q = SparseMultivariateIntegerPolynomial(ret, nvar, names);
                r = SparseMultivariateIntegerPolynomial(NULL, nvar, names);
                return true;
            } else {
                q = SparseMultivariateIntegerPolynomial(NULL, nvar, names);
                r = *this; 
                return false;
            }
        }

        if (nvar == b.nvar) {
            AltArrZ_t* qAA = NULL, * rAA = NULL;
            dividePolynomials_AAZ(this->poly, b.poly, &qAA, &rAA, nvar);
            q = SparseMultivariateIntegerPolynomial(qAA, nvar, names);
            r = SparseMultivariateIntegerPolynomial(rAA, nvar, names);
            return (rAA == NULL || rAA->size == 0);
        } else {
            AltArrZ_t* tempB = makeConstPolynomial_AAZ(1, nvar, b.poly->elems->coef);
            AltArrZ_t* qAA = NULL, * rAA = NULL;
            dividePolynomials_AAZ(this->poly, tempB, &qAA, &rAA, nvar);
            q = SparseMultivariateIntegerPolynomial(qAA, nvar, names);
            r = SparseMultivariateIntegerPolynomial(rAA, nvar, names);
            freePolynomial_AAZ(tempB);
            return (rAA == NULL || rAA->size == 0);
        }
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

        SparseMultivariateIntegerPolynomial tempc = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateIntegerPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        AltArrZ_t* qAA = NULL, * rAA = NULL;
        dividePolynomials_AAZ(tempc.poly, tempb.poly, &qAA, &rAA, superNvar);

        q = SparseMultivariateIntegerPolynomial(qAA, superNvar, newnames);
        r = SparseMultivariateIntegerPolynomial(rAA, superNvar, newnames);

        return (rAA == NULL || rAA->size == 0);

    } else {

        // std::cerr << "divide on same ring" << std::endl;

        AltArrZ_t* qAA = NULL, * rAA = NULL;
        dividePolynomials_AAZ(this->poly, b.poly, &qAA, &rAA, nvar);
        q = SparseMultivariateIntegerPolynomial(qAA, nvar, names);
        r = SparseMultivariateIntegerPolynomial(rAA, nvar, names);
        return (rAA == NULL || rAA->size == 0);
    }    
}

bool SparseMultivariateIntegerPolynomial::rationalDivide(const SparseMultivariateIntegerPolynomial& b, SparseMultivariateRationalPolynomial& q, SparseMultivariateRationalPolynomial& r) const {
    if (this->isZero()) {
        q.zero();
        r.zero();
        return true;        
    }

    if (b.isZero()) {
        std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
        exit(1);
    }

    AltArr_t* tempc = deepCopyPolynomial_AAFromAAZ(this->poly);
    AltArr_t* tempb = deepCopyPolynomial_AAFromAAZ(b.poly);

    SparseMultivariateRationalPolynomial rC(tempc, nvar, names);
    SparseMultivariateRationalPolynomial rB(tempb, b.nvar, b.names);
    return rC.divide(rB, q, r);
}


/**
 * Get the remainder of *this divided by b.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator% (const SparseMultivariateIntegerPolynomial& b) const {
    SparseMultivariateIntegerPolynomial q,r;
    this->divide(b, q, r);
    return r;
}

/**
 * Update *this by setting it to the remainder of *this divided by b.
 */ 
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator%= (const SparseMultivariateIntegerPolynomial& b) {
    *this = (*this % b);
    return *this;
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::pseudoDivide(const SparseMultivariateIntegerPolynomial& b, SparseMultivariateIntegerPolynomial* quo, SparseMultivariateIntegerPolynomial* mult, bool lazy) const {
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
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    //This will also capture case when nvar = 0 for both.
    if (b.isConstant() && isConstant()) {
        if (quo != NULL) {
            *quo = *this;
        }
        if (mult != NULL) {
            *mult = b;
        }
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);

    if (!isOrdered) {
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

    SparseMultivariateIntegerPolynomial q, r;

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

        SparseMultivariateIntegerPolynomial tempc = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateIntegerPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

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

            AltArrZ_t* cAA = tempc.poly;
            AltArrZ_t* bAA = tempb.poly;

            RecArrZ_t* recC = convertToRecursiveArrayZ(cAA);
            RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

            AltArrZ_t* qAA = NULL, * rAA = NULL;
            AltArrZ_t* hPow = NULL;
            int e = 0;
            pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);
            cAA = convertFromRecursiveArrayZ(recC, superNvar);
            bAA = convertFromRecursiveArrayZ(recB, superNvar);
            tempc.poly = cAA;
            tempb.poly = bAA;

            q = SparseMultivariateIntegerPolynomial(qAA, superNvar, newvars);
            r = SparseMultivariateIntegerPolynomial(rAA, superNvar, newvars);
            q.reorderVarsInPlace(invVarmap);
            r.reorderVarsInPlace(invVarmap);

            SparseMultivariateIntegerPolynomial smqpHPow(hPow, superNvar, newvars);
            smqpHPow.reorderVarsInPlace(invVarmap);
            if (quo != NULL) {
                *quo = q;
            }
            if (mult != NULL) {
                *mult = smqpHPow;
            }
            return r;
        }

        AltArrZ_t* cAA = tempc.poly;
        AltArrZ_t* bAA = tempb.poly;
        RecArrZ_t* recC = convertToRecursiveArrayZ(cAA);
        RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

        AltArrZ_t* qAA = NULL, * rAA = NULL;
        AltArrZ_t* hPow = NULL;
        int e = 0;
        pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);
        cAA = convertFromRecursiveArrayZ(recC, superNvar);
        bAA = convertFromRecursiveArrayZ(recB, superNvar);
        tempc.poly = cAA;
        tempb.poly = bAA;

        q = SparseMultivariateIntegerPolynomial(qAA, superNvar, newnames);
        r = SparseMultivariateIntegerPolynomial(rAA, superNvar, newnames);
        SparseMultivariateIntegerPolynomial smqpHPow(hPow, superNvar, newnames);
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
            SparseMultivariateIntegerPolynomial cReordered = this->expandVariables(nvar, newvars, varmap);
            SparseMultivariateIntegerPolynomial bReordered = b.expandVariables(nvar, newvars, varmap);

            AltArrZ_t* cAA = cReordered.poly;
            AltArrZ_t* bAA = bReordered.poly;
            RecArrZ_t* recC = convertToRecursiveArrayZ(cAA);
            RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

            AltArrZ_t* qAA = NULL, * rAA = NULL;
            AltArrZ_t* hPow = NULL;
            int e = 0;
            pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);
            cReordered.poly = cAA = convertFromRecursiveArrayZ(recC, nvar);
            bReordered.poly = bAA = convertFromRecursiveArrayZ(recB, nvar);
            
            q = SparseMultivariateIntegerPolynomial(qAA, nvar, newvars);
            r = SparseMultivariateIntegerPolynomial(rAA, nvar, newvars);

            q.reorderVarsInPlace(invVarmap);
            r.reorderVarsInPlace(invVarmap);

            SparseMultivariateIntegerPolynomial smqpHPow(hPow, superNvar, newvars);
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
            AltArrZ_t* qAA = NULL, *rAA = NULL;
            int e = 0;
            univariatePseudoDividePolynomials_AAZ(this->poly, b.poly, &qAA, &rAA, &e, lazy);
            
            if (quo != NULL) {
                *quo = SparseMultivariateIntegerPolynomial(qAA, nvar, names);
            }

            if (mult != NULL) {
                Integer h(b.poly->elems->coef);
                h ^= e;
                *mult = h;                
            }
            return SparseMultivariateIntegerPolynomial(rAA, nvar, names);
        }

        int e = 0;
        AltArrZ_t* bAA = b.poly;
        RecArrZ_t* recC = convertToRecursiveArrayZ(poly);
        RecArrZ_t* recB = convertToRecursiveArrayZ(bAA);

        AltArrZ_t* qAA = NULL, *rAA = NULL;
        AltArrZ_t* hPow = NULL;
        pesudoDivide_RecArrayZ(recC, recB, &qAA, &rAA, &e, &hPow, superNvar, lazy);

        poly = convertFromRecursiveArrayZ(recC, superNvar);
        bAA = convertFromRecursiveArrayZ(recB, superNvar);

        q = SparseMultivariateIntegerPolynomial(qAA, superNvar, names);
        r = SparseMultivariateIntegerPolynomial(rAA, superNvar, names);

        SparseMultivariateIntegerPolynomial smqpHPow(hPow, superNvar, names);
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
 * Add *this and an mpz_class.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator+ (const mpz_class& c) const {
    SparseMultivariateIntegerPolynomial ret = *this;

    if (isZero()) {
        if (ret.poly != NULL) {
            freePolynomial_AAZ(ret.poly);
        }
        ret.poly = makeConstPolynomial_AAZ(1, ret.nvar, c.get_mpz_t());
        return ret;
    }

    addInteger_AAZ_inp(ret.poly, c.get_mpz_t());
    return ret;
}


/**
 * Update *this by adding r
 */
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator+= (const mpz_class& c) {    
    slp.clear();

    if (isZero()) {
        if (poly != NULL) {
            freePolynomial_AAZ(poly);
        }
        poly = makeConstPolynomial_AAZ(1, nvar, c.get_mpz_t());
        return *this;
    }

    addInteger_AAZ_inp(poly, c.get_mpz_t());
    return *this;
}

/**
 * Subtract the mpz_class c from *this.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator- (const mpz_class& c) const {
    mpz_class negC = -c;
    return *this + negC;
}

/**
 * Update *this by subtracting mpz_class c.
 */
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator-= (const mpz_class& c) {
    mpz_class negC = -c;
    return *this += negC;
}

/**
 * Multiply *this by mpz_class c.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator* (const mpz_class& c) const {
    if (c == 0 || isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    SparseMultivariateIntegerPolynomial ret = *this;

    if (c == 1) {
        return ret;
    }

    mpz_t mult;
    mpz_init(mult);
    mpz_set(mult, c.get_mpz_t());
    for(int i = 0; i < ret.poly->size; ++i) {
        mpz_mul(ret.poly->elems[i].coef, ret.poly->elems[i].coef, mult);
    }
    mpz_clear(mult);

    return ret;
}

/**
 * Update *this by multiplying by mpz_class c.
 */
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator*= (const mpz_class& c) {
    slp.clear();
    if (isZero()) {
        return *this;
    }
    if (c == 0) {
        freePolynomial_AAZ(poly);
        poly = NULL;
        return *this;
    }
    if (c == 1) {
        return *this;
    }

    mpz_t mult;
    mpz_init(mult);
    mpz_set(mult, c.get_mpz_t());
    for(int i = 0; i < poly->size; ++i) {
        mpz_mul(poly->elems[i].coef, poly->elems[i].coef, mult);
    }
    mpz_clear(mult);

    return *this;
}

/**
 * Divide *this by mpz_class c
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator/ (const mpz_class& c) const {
    if (c == 0) {
        std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
        exit(1);
    }

    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    if (c == 1) {
        return SparseMultivariateIntegerPolynomial(*this);
    }
    

    SparseMultivariateIntegerPolynomial ret = *this;

    mpz_t r;
    mpz_init(r);
    mpz_set(r, c.get_mpz_t());
    for(int i = 0; i < ret.poly->size; ++i) {
        if (mpz_divisible_p(ret.poly->elems[i].coef, r) == 0) {
            std::cerr << "BPAS ERROR: SMZP: Non-exact division dividing by mpz_class: " << c << std::endl;
            exit(1);
        }
        mpz_divexact(ret.poly->elems[i].coef, ret.poly->elems[i].coef, r);
    }

    mpz_clear(r);
    return ret;
}

/**
 * Divide mpz_t r by SMQP b.
 */
SparseMultivariateIntegerPolynomial operator/ (const mpz_t r, const SparseMultivariateIntegerPolynomial& b) {
    if (b.isConstant()) {
        if (mpz_sgn(b.poly->elems->coef) == 0) {
            std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
            exit(1);
        }
        if (mpz_divisible_p(r, b.poly->elems->coef) == 0) {
            std::cerr << "BPAS ERROR: SMZP non-exact division in mpz_t / SparseMultivariateIntegerPolynomial" << std::endl;
            exit(1);
        }

        AltArrZ_t* rPoly = makeConstPolynomial_AAZ(1, b.nvar, r);
        mpz_divexact(rPoly->elems->coef, rPoly->elems->coef, b.poly->elems->coef);
        return SparseMultivariateIntegerPolynomial(rPoly, b.nvar, b.names);
    } else {
        std::cerr << "BPAS ERROR: SMZP non-exact division in mpz_t / SparseMultivariateIntegerPolynomial" << std::endl;
        exit(1);
    }
}

/** 
 * Update *this by dividing by mpz_class c.
 */
SparseMultivariateIntegerPolynomial& SparseMultivariateIntegerPolynomial::operator/= (const mpz_class& c) {
    if (c == 0) {
        std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
        exit(1);
    }

    if (isZero()) {
        return *this;
    }

    if (c == 1) {
        return *this;
    }
    
    mpz_t r;
    mpz_init(r);
    mpz_set(r, c.get_mpz_t());
    for(int i = 0; i < poly->size; ++i) {
        if (mpz_divisible_p(poly->elems[i].coef, r) == 0) {
            std::cerr << "BPAS ERROR: SMZP: Non-exact division dividing by mpz_class: " << c << std::endl;
            exit(1);
        }
        mpz_divexact(poly->elems[i].coef, poly->elems[i].coef, r);
    }

    mpz_clear(r);
    return *this;
}


/**
 * Get the polynomial term at index. Returns 0 if index is beyond the 
 * number of terms in this polynomial.
 */
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::operator[] (int index) const {
    AltArrZ_t* term = termAtIdx_AAZ(poly, index);
    return SparseMultivariateIntegerPolynomial(term, nvar, names);
}


/**
 * Get the leading variable, that is, the highest-order variable with positive degree
 * of this polynomial. 
 * returns the leading variable or the empty string if this polynomial has zero variables.
 */
Symbol SparseMultivariateIntegerPolynomial::leadingVariable() const {
    if (nvar == 0 || isConstant()) {
        return Symbol("");
    }

    int lvIdx = mainVariable_AAZ(poly);
    return names[lvIdx + 1];
}

/**
 * Get the degree of this polynomial w.r.t the leading variable.
 */
Integer SparseMultivariateIntegerPolynomial::leadingVariableDegree() const {
    if (nvar == 0 || isConstant()) {
        return 0;
    }

    return mainDegree_AAZ(poly);
}

/**
 * Is the contant term zero.
 */
bool SparseMultivariateIntegerPolynomial::isConstantTermZero() const {
    if (isZero()) {
        return 1;
    }

    return isConstantTermZero_AAZ(poly);
}

SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::leadingCoefficientInVariable (const Symbol& x, int* e) const {
    int deg = 0;
    if (e != NULL) {
        *e = 0;
    }
    if (isZero()) {
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
    SparseMultivariateIntegerPolynomial r(v);
    for (int i = 0; i < k; ++i) {
        r.names[i] = names[i];
    }
    for (int i = k; i < v+1; ++i) {
        r.names[i] = names[i+1];
    }
    
    if (isConstant()) {
        r.poly = makeConstPolynomial_AAZ(1, v, poly->elems->coef);
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

    AltArrZ_t* tempPoly = deepCopyPolynomial_AAZ(poly);
    if (k != 0) {
        reorderVars_AAZ(tempPoly, varmap, nvar);
    }

    deg = mainLeadingDegree_AAZ(tempPoly);
    AltArrZ_t* lc = mainLeadingCoefficient_AAZ(tempPoly);
    std::vector<const char*> syms;
    for (int i = 0; i < nvar; ++i) {
        syms.push_back(names[i+1].toString().c_str());
    }

    freePolynomial_AAZ(tempPoly);
    shrinkNumVarsAtIdx_AAZ(lc, 0);
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
SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial> SparseMultivariateIntegerPolynomial::convertToSUP(const Symbol& x) const {
    SparseUnivariatePolynomial<SparseMultivariateIntegerPolynomial> r;
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
            SparseMultivariateIntegerPolynomial coef(Integer(poly->elems[i].coef), 0);
            d = partialDegreeTerm_AAZ(poly, i, 0);
            r.setCoefficient(d, coef);
        }
        return r;
    }

    int v = nvar - 1;
    SparseMultivariateIntegerPolynomial tempSMZP(v);
    for (int i = 0; i < k; ++i) {
        tempSMZP.names[i] = names[i];
    }
    for (int i = k; i < v+1; ++i) {
        tempSMZP.names[i] = names[i+1];
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

    AltArrZ_t* tempPoly;
    if (k != 0) {
        tempPoly = deepCopyPolynomial_AAZ(poly);
        reorderVars_AAZ(tempPoly, varmap, nvar);
    } else {
        tempPoly = poly;
    }

    RecArrZ_t* recArr = convertToRecursiveArrayZ(tempPoly);
    AltArrZ_t localAA;
    for (int i = 0; i < recArr->size; ++i) {
        degree_t deg = recArr->elems[i].exp;
        localAA.nvar = nvar;
        localAA.size = localAA.alloc = recArr->elems[i].coefSize;
        localAA.unpacked = recArr->unpacked;
        localAA.elems = recArr->elems[i].coef;

        AltArrZ_t* toShrinkAA = deepCopyPolynomial_AAZ(&localAA);
        shrinkNumVarsAtIdx_AAZ(toShrinkAA, 0);
        tempSMZP.poly = toShrinkAA;

        r.setCoefficient(deg, tempSMZP);

        freePolynomial_AAZ(tempSMZP.poly);
        tempSMZP.poly = NULL;
    }

    if (k != 0) {
        freeRecArrayZAndCoef(recArr);
    } else {
        convertFromRecursiveArrayZ(recArr, nvar);
    }

    return r;

    // int d = this->degree(x).get_si();
    // std::vector<SparseMultivariateIntegerPolynomial> mpolys(d+1, SparseMultivariateIntegerPolynomial(v));
    
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
    // degree_t curD;
    // for (int i = 0; i < poly->size; ++i) {
    //     curD = GET_NTH_EXP(poly->elems[i].degs, masks[k], sizes[k]);
    //     if (mpolys[curD].poly == NULL) {
    //         mpolys[curD].poly = makePolynomial_AAZ(10, nvar);
    //     }
    //     if (mpolys[curD].poly->alloc <= mpolys[curD].poly->size) {
    //         resizePolynomial_AAZ(mpolys[curD].poly, mpolys[curD].poly->size * 2);
    //     }
    //     mpz_init(mpolys[curD].poly->elems[mpolys[curD].poly->size].coef);
    //     mpz_set(mpolys[curD].poly->elems[mpolys[curD].poly->size].coef, poly->elems[i].coef);
    //     mpolys[curD].poly->elems[mpolys[curD].poly->size].degs = poly->elems[i].degs & (~masks[k]);
    //     ++(mpolys[curD].poly->size);
    // }
    
    // for (int i = 0; i <= d; ++i) {
    //     if (mpolys[i].isZero()) {
    //         continue;
    //     }
    //     shrinkNumVarsAtIdx_AAZ(mpolys[i].poly, k);
    //     if (!(k == 0 || isLeadVar)) {
    //         mergeSortPolynomial_AAZ(mpolys[i].poly);
    //     }
    //     r.setCoefficient(i, mpolys[i]);
    // }

    // free(masks);
    // free(sizes);

    // return r;
}


/**
 * Negate all the coefficients of *this. Note, that due to the 
 * sharing nature of underling nodes, this may alter the Nodes of
 * other SMQP.
 */
void SparseMultivariateIntegerPolynomial::negate() {
    if (isZero()) {
        return;
    }
    negatePolynomial_AAZ(poly);
}

/** 
 * Get a copy of this polynomial. 
 */ 
SparseMultivariateIntegerPolynomial SparseMultivariateIntegerPolynomial::deepCopy() const {
    if (isZero()) {
        return SparseMultivariateIntegerPolynomial(NULL, nvar, names);
    }

    AltArrZ_t* copy = deepCopyPolynomial_AAZ(poly);
    return SparseMultivariateIntegerPolynomial(copy, nvar, names);
}

Factors<SparseMultivariateIntegerPolynomial> SparseMultivariateIntegerPolynomial::factor() const {
    Factors<SparseMultivariateIntegerPolynomial> ret;
    if (isConstant() || isZero() || nvar == 0) {
        ret.setRingElement(*this);
        return ret;
    }

    degree_t degsList[nvar];
    partialDegrees_AAZ(poly, degsList);
    bool tryFact = 0;
    for (int i = 0; i < nvar; ++i) {
        if (degsList[i] > 1) {
            tryFact = 1;
            break;
        }
    }
    if (!tryFact) {
        Integer cont;
        SparseMultivariateIntegerPolynomial prim = this->primitivePart(cont);
        ret.setRingElement(SparseMultivariateIntegerPolynomial(cont));
        ret.addFactor(std::move(prim), 1);
        return ret;
    }

//NTL factoring is currently buggy. Jan 2020.
// #if defined(WITH_NTL) && WITH_NTL

//     AltArrZ_t** factors;
//     int* exps;
//     int nfacts = 0;
//     mpz_t cont;
//     mpz_init(cont);

//     SMZP::Factoring::factor(this->poly, &factors, &exps, &nfacts, cont);

//     SparseMultivariateIntegerPolynomial coefPoly(cont);
//     ret.setRingElement(std::move(coefPoly));

//     for (int i = 0; i < nfacts; ++i) {
//         SparseMultivariateIntegerPolynomial factPoly(factors[i], nvar, names);
//         ret.addFactor(std::move(factPoly), exps[i]);
//     }
//     free(exps);
//     free(factors);
// #elif defined(WITH_MAPLE) && WITH_MAPLE

#if defined(WITH_MAPLE) && WITH_MAPLE

    MapleInterfaceStream& mis = MapleInterfaceStream::instance();
    ret = mis.factor(*this);

    // char* c_names[nvar];
    // for (int i = 0; i < nvar; ++i) {
    //     std::string str = names[i+1].toString();
    //     c_names[i] = (char*) malloc(sizeof(char)*str.length()+1);
    //     strcpy(c_names[i], str.c_str());
    // } 

    // mpz_t numericFact;
    // mpz_init(numericFact);
    // AltArrZ_t** facts = NULL;
    // int* exps = NULL;
    // int numFacts;
    // factorPolynomial_MplInt_AAZ(poly, (const char**) c_names, &numFacts, &facts, &exps, numericFact);
    
    // SparseMultivariateIntegerPolynomial coefPoly(numericFact);
    // ret.setRingElement(coefPoly);

    // for (int i = 0; i < numFacts; ++i) {
    //     SparseMultivariateIntegerPolynomial factPoly(facts[i], nvar, names);
    //     ret.addFactor(factPoly, exps[i]);
    // }
    // free(exps);
    // free(facts);
    
    // for (int i = 0; i < nvar; ++i) {
    //     free(c_names[i]);
    // }

#elif defined(WITH_BLAD) && WITH_BLAD
    char* c_names[nvar];
    for (int i = 0; i < nvar; ++i) {
        std::string str = names[i+1].toString();
        c_names[i] = (char*) malloc(sizeof(char)*str.length()+1);
        strcpy(c_names[i], str.c_str());
    } 

    mpz_t numericFact;
    mpz_init(numericFact);
    AltArrZ_t** facts = NULL;
    int* exps;
    int numFacts;
    factorPolynomialBLAD_AAZ(poly, (const char**) c_names, &numFacts, &facts, &exps, numericFact);
   
    SparseMultivariateIntegerPolynomial coefPoly(numericFact);
    ret.setRingElement(coefPoly);

    for (int i = 0; i < numFacts; ++i) {
        SparseMultivariateIntegerPolynomial factPoly(facts[i], nvar, names);
        ret.addFactor(factPoly, exps[i]);
    }

    free(exps);
    free(facts);
    for (int i = 0; i < nvar; ++i) {
        free(c_names[i]);
    }
#else
    ret.addFactor(*this, 1);
#endif

    return ret;
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
void SparseMultivariateIntegerPolynomial::straightLineProgram() {
    slp.clear();
    if (isZero()) {
        SLPZRepresentation r; 
        r.type = 4;
        r.b = -1;
        r.a.c = new Integer(0);
        slp.push_back(r);
        return; 
    }

    if (isConstant()) {
        SLPZRepresentation r;
        r.type = 4;
        r.b = -1;
        r.a.c = new Integer(poly->elems->coef);
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


    std::vector<AAZElem_DegList_t*> polyStack;
    std::vector<int> indexStack;
    int prev = 0;
    bool isOp = 1;
    SLPZRepresentation elem;

    AltArrZDegList_t* listPoly = deepCopyPolynomial_AAZDegListFromAA(poly);

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
            elem.a.c = new Integer(listPoly->elems[i].coef);
            // elem.a  = i;
            elem.b = k;
            slp.push_back(elem);

            d[k]--;
        }
        else {
            elem.op = 0;
            elem.type = 1;
            elem.a.c = new Integer(listPoly->elems[i].coef);
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
    freePolynomial_AAZDL(listPoly);

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

void SparseMultivariateIntegerPolynomial::printSLP(std::ostream& out) const {
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

void SparseMultivariateIntegerPolynomial::randomPolynomial(int numvar, int nterms, unsigned long int coefBound, degree_t sparsity, bool includeNeg) {

    *this = SparseMultivariateIntegerPolynomial(numvar);

    if (numvar == 0) {
        poly = makePolynomial_AAZ(1, nvar);
        mpz_init(poly->elems->coef);
        mpz_set_si(poly->elems->coef, (rand() % ((int) pow(2, coefBound-1))) + 1);
        poly->elems->degs = 0;
        poly->size = 1;

    } else {
        Node* node = buildRandomZPoly(numvar, nterms, coefBound, sparsity, includeNeg);
        poly = deepCopyPolynomial_AAZFromNode(node, nvar);
        freePolynomial(node);
    }

}

void SparseMultivariateIntegerPolynomial::randomPolynomial(std::vector<int> maxDegs, unsigned long int coefBound, float sparsity, bool includeNeg) {

    AltArrZ_t* aa = buildRandomZPolyFromMax(maxDegs.size(), maxDegs.data(), coefBound, sparsity, includeNeg);

    *this = SparseMultivariateIntegerPolynomial(maxDegs.size());
    this->poly = aa;
}



/**
 * Construct an ExprTreeNode (well, multiple) which represents a single
 * term. That is, a coefficient and a monomial.
 */
ExprTreeNode* exprTreeNodeFromAAZElem(const AltArrZ_t* n, int idx, int nvar, const Symbol* vars) {
    if (n == NULL) {
        return new ExprTreeNode(0l);
    }

    ExprTreeNode* t = new ExprTreeNode(mpz_class(n->elems[idx].coef));
    
    if (nvar > 0) {
        degree_t degsList[nvar] = {0};
        int othernvar = n->nvar;
        partialDegreesTerm_AAZ(n, idx, degsList);
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

ExpressionTree SparseMultivariateIntegerPolynomial::convertToExpressionTree() const {
    if (isZero()) {
        ExprTreeNode* r = new ExprTreeNode(0l);
        ExpressionTree t(r);
        return t; 
    } 
   
    ExprTreeNode* prev = exprTreeNodeFromAAZElem(poly, 0, nvar, &names[1]);
    for (int i = 1; i < poly->size; ++i) {
        ExprTreeNode* thisNode = exprTreeNodeFromAAZElem(poly, i, nvar, &names[1]);
        prev = ExprTreeNode::combineExprTreeNodes(prev, thisNode, EXPR_ADD);
    }

    return ExpressionTree(prev);
}


