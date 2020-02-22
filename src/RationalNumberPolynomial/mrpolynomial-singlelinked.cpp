#include "RationalNumberPolynomial/mrpolynomial.h"

////////// Construtors ////////////////////////////////////

/**
 * Private constructor to create SMQP directly from a Node*.
 * Makes a copy of input varNames
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(Node* head, int vars, Symbol* varNames) : 
    poly(head),
    tailNode(NULL),
    constTerm(NULL),
    nvar(vars) 
{
    poly = head;
    tailNode = NULL;
    nvar = vars;
    names = new Symbol[nvar+1];
    std::copy(varNames, varNames+nvar+1, names);
    constTerm = NULL;
}

/**
 * Construct a multivariate polynomial
 *
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial() : 
    nvar(0),
    poly(NULL),
    tailNode(NULL),
    constTerm(NULL)
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
    nvar(v),
    poly(NULL),
    tailNode(NULL),
    constTerm(NULL)
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
    nvar(1),
    poly(NULL),
    tailNode(NULL),
    constTerm(NULL)
{
    names = new Symbol[2];
    names[0] = "9";
    names[1] = x;

    degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
    degs[0] = 1;
    ratNum_t coef;
    mpq_init(coef);
    mpq_set_ui(coef, 1ul, 1ul);
    poly = addTerm(NULL, degs, coef);
    mpq_clear(coef);
}

/**
 * Copy Constructor.
 * 
 * Does not reuse underlying memory allocated by b. 
 *
 * @param b: A sparse multivariate polynomial
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(const SparseMultivariateRationalPolynomial& b) :
    nvar(b.nvar),
    tailNode(NULL),
    constTerm(NULL)
{
    poly = deepCopyPolynomial(b.poly, b.nvar);
    names = new Symbol[nvar+1];
    std::copy(b.names, b.names+nvar+1, names);
    slp = std::vector<SLPRepresentation>(b.slp);

    if (b.constTerm != NULL) {
        constTerm = new RationalNumber(*(b.constTerm));
    }
}

/**
 * Move Constructor.
 *
 * @params b: The r-value reference polynomial.
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(SparseMultivariateRationalPolynomial&& b) {
    nvar = b.nvar;
    poly = b.poly;
    tailNode = b.tailNode;
    names = new Symbol[nvar+1];
    std::copy(b.names, b.names+nvar+1, names);
    slp = b.slp;
    constTerm = b.constTerm;
    b.poly = NULL;
    b.tailNode = NULL;
    b.slp.clear();
    b.constTerm = NULL;
}

/**
 * Create a SMQP from a Integer. 
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(const Integer& r, int nv) :
    nvar(nv),
    poly(NULL),
    tailNode(NULL),
    constTerm(NULL)
{
    if (nv == 0) {
        constTerm = new RationalNumber(r);   
        names = new Symbol[1];
        names[0] = "1";
    } else {
        names = new Symbol[nvar+1];
        names[0] = "1";
        for (int i = 1; i <= nvar; ++i) {
            std::ostringstream convert;
            convert << i;
            names[i] = "x_";
            names[i] += convert.str();
        }
        degrees_t degs = (degrees_t) calloc(nv, sizeof(degree_t));
        poly = addZeroTerm(NULL, degs);
    }
}

/**
 * Create a SMQP from a RationalNumber. 
 */
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial(const RationalNumber& r, int nv) :
    nvar(nv),
    poly(NULL),
    tailNode(NULL),
    constTerm(NULL)
{
    if (nv == 0) {
        names = new Symbol[1];
        names[0] = "1";
        constTerm = new RationalNumber(r);   
    } else {
        names = new Symbol[nvar+1];
        names[0] = "1";
        for (int i = 1; i <= nvar; ++i) {
            std::ostringstream convert;
            convert << i;
            names[i] = "x_";
            names[i] += convert.str();
        }
        degrees_t degs = (degrees_t) calloc(nv, sizeof(degree_t));
        poly = addZeroTerm(NULL, degs);
    }
}

/**
 * Create a SMQP from a univariate rational polynomial. 
 * 
 * @param p: A SUQP polynomial. 
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial (const DenseUnivariateRationalPolynomial& p) {
    nvar = 1;
    tailNode = NULL;
    names = new Symbol[nvar+1];
    names[0] = "9";
    names[1] = p.variable();
    constTerm = NULL;

    mpq_class coef;
    poly = NULL;
    Node* trailingTerm = NULL;
    size_t pDeg = p.degree().get_si();
    for (int i = pDeg; i >= 0; --i) {
        coef = p.coefficient(i).get_mpq();
        if (coef != 0) {
            degrees_t degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
            degs[0] = i;
            mpq_t mpqCoef;
            trailingTerm = addTerm(trailingTerm, degs, coef.get_mpq_t());
            if (poly == NULL) {
                poly = trailingTerm;
            }
        }
    }
}

/**
 * Construct from a SUP<SMQP> polynomial
 *
 * @param s: The SUP<SMQP> polynomial
 **/
SparseMultivariateRationalPolynomial::SparseMultivariateRationalPolynomial (const SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial>& s) {
    nvar = 0;
    poly = NULL;
    tailNode = NULL;
    names = new Symbol[1];
    names[0] = "1";
    constTerm = NULL;

    long int d = s.degree().get_si();
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

        Node* trail = NULL; 
        degrees_t degs;
        for (int k = d; k >= 0; --k) {  
            SparseMultivariateRationalPolynomial coef = s.coefficient(k);
            if (coef.isZero()) {
                continue;
            }
            degs = (degrees_t) malloc(sizeof(degree_t));
            degs[0] = k;
            mpq_class mpqCoef = coef.constTerm->get_mpq();
            trail = addTerm(trail, degs, mpqCoef.get_mpq_t());
            if (poly == NULL) {
                poly = trail;
            }
        }

        return;
    }

    for (int k = 0; k <= d; ++k) {
        SparseMultivariateRationalPolynomial t;
        c = s.coefficient(k);
        
        //coefficient is 0 so do nothing
        if (c.isZero()) {
            continue;
        }

        t.nvar = c.nvar + 1; 
        t.names = new Symbol[t.nvar+1];
        t.names[0] = "9";
        t.names[1] = s.variable();
        if (c.nvar > 0) {
            std::copy(c.names+1, c.names+c.nvar+1, t.names+2);
        }


        t.poly = NULL;
        Node* curC = c.poly;
        Node* trailT = NULL;
        while (curC != NULL) {
            degrees_t degs = (degrees_t) malloc(sizeof(degree_t)*t.nvar);
            degs[0] = k;
            for (int j = 0; j < c.nvar; ++j) {
                degs[j+1] = curC->degs[j];
            }
            trailT = addTerm(trailT, degs, curC->coef);
            if (t.poly == NULL) {
                t.poly = trailT;
            }
            curC = curC->next;
        }

       *this += t;
    }

}

/**
 * Destroy the polynomial and release underlying node memory.
 **/
SparseMultivariateRationalPolynomial::~SparseMultivariateRationalPolynomial() {
    delete[] names;
    freePolynomial(poly);
    slp.clear();
    delete constTerm;
}



////////// BPASRing /////////////////////////////////

bool SparseMultivariateRationalPolynomial::isZero() const {

    if (nvar > 0 && poly == NULL) {
        return 1;
    }
    if (nvar == 0 && (constTerm == NULL || *constTerm == 0)) {
        return 1;
    }
    return 0;
}

void SparseMultivariateRationalPolynomial::zero() {
    freePolynomial(poly);
    poly = NULL;
    tailNode = NULL;
    delete constTerm;
    constTerm = NULL;
}

bool SparseMultivariateRationalPolynomial::isOne() const {
    if (poly != NULL) {
        if (isZeroExponentVector(poly->degs, nvar) && mpq_cmp_ui(poly->coef, 1ul, 1ul) == 0) { 
            return 1;
        }
    }
    if (constTerm != NULL && *constTerm == 1) {
        return 1;
    }
    return 0;
}

void SparseMultivariateRationalPolynomial::one() {
    freePolynomial(poly);
    poly = NULL;
    tailNode = NULL;

    if (nvar > 0) {
        degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
        ratNum_t c;
        mpq_init(c);
        mpq_set_ui(c, 1ul, 1ul);
        poly = addTerm(NULL, degs, c);
        mpq_clear(c);    
    } else {
        constTerm = new RationalNumber(1);
    }
}

bool SparseMultivariateRationalPolynomial::isNegativeOne() const {
    if (poly != NULL) {
        if (isZeroExponentVector(poly->degs, nvar) && mpq_cmp_si(poly->coef, -1l, 1l) == 0) {
            return 1;
        }
    }
    if (constTerm != NULL && *constTerm == -1) {
        return 1;
    }
    return 0;
}

void SparseMultivariateRationalPolynomial::negativeOne() {
    freePolynomial(poly);
    poly = NULL;
    tailNode = NULL;

    if (nvar > 0) {
        degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
        ratNum_t c;
        mpq_init(c);
        mpq_set_si(c, -1l, 1l);
        poly = addTerm(NULL, degs, c);
        mpq_clear(c);
    } else {
        constTerm = new RationalNumber(-1);
    }
}

int SparseMultivariateRationalPolynomial::isConstant() const {
    if (poly == NULL) {
        if (constTerm == NULL) {
            return 1;
        }
        if (*constTerm >= 0) {
            return 1;
        } else {
            return -1;
        }
    }
    if (isZeroExponentVector(poly->degs, nvar)) {
        if (mpq_cmp_ui(poly->coef, 0ul, 1ul) >= 0) {
            return 1;
        } else {
            return -1;
        }
    }
    return 0;
}

/**
 * Get GCD between *this and b.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::gcd(const SparseMultivariateRationalPolynomial& b) const {
    std::cerr << "BPAS ERROR: SparseMultivariateRationalPolynomial::gcd NOT YET IMPLEMENTED" << std::endl;
    return *this;
}



////////// Private Helpers ////////////////////////////////////

/**
 * Get the tail Node.
 */
Node* SparseMultivariateRationalPolynomial::getTailNode() const {
    if (tailNode != NULL) {
        return tailNode;
    }

    if (poly == NULL) {
        return NULL;
    }

    Node* node = poly;
    while (node->next != NULL) {
        node = node->next;
    }
    tailNode = node;
    return node;
}

void SparseMultivariateRationalPolynomial::removeTailNode() {
    if (poly == NULL) {
        return;
    }

    Node* node = poly;
    Node* prev = poly;
    while (node->next != NULL) {
        prev = node;
        node = node->next;
    }

    //node is now tail;
    if (node == poly) {
        poly = NULL;
        tailNode = NULL;
        freePolynomial(node);
    } else {
        prev->next = NULL;
        tailNode = prev;
        freePolynomial(node);
    }
}


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
    if (names[0] != b.names[0]) { return 0; }
    if (names[0] == "1") {
        int k = 1;
        for (; k <= nvar && k <= b.nvar; ++k) {
            xs.push_back(k);
            xs.push_back(k);
        }
        for (int i = k; i <= nvar; ++i) {
            xs.push_back(i);
            xs.push_back(0);
        }
        for (int i = k; i <= b.nvar; ++i) {
            xs.push_back(0);
            xs.push_back(i);
        }
        return 1;
    }
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
    SparseMultivariateRationalPolynomial temp(NULL, vars, newvars);

    Node* trailT = NULL;
    Node* curNode = poly;
    int i;
    degrees_t curDegs, degs;
    while (curNode != NULL) {
        curDegs = curNode->degs;
        degs = (degrees_t) calloc(vars, sizeof(degree_t));
        for (i = 0; i < nvar; ++i) {
            degs[varmap[i]] = curDegs[i];
        }

        trailT = addTerm(trailT, degs, curNode->coef);
        if (temp.poly == NULL) {
            temp.poly = trailT;
        }
        curNode = curNode->next;
    }

    //reorder polynomial so that it maintains lexicographical in the new variable ordering
    sortPolynomial(&temp.poly, vars); 

    return temp;
}

/** 
 * Rearrange exponent vectors in place and then re-sort the polynomial.
 */
void SparseMultivariateRationalPolynomial::reorderVarsInPlace(int varmap[]) {

    Node* curNode = poly;
    degrees_t degs;
    degrees_t tempDegs = (degrees_t) malloc(nvar*sizeof(degree_t));
    while(curNode != NULL) {
        degs = curNode->degs;
        for (int i = 0; i < nvar; ++i) {
            tempDegs[i] = degs[i];
        }
        for (int i = 0; i < nvar; ++i) {
            degs[varmap[i]] = tempDegs[i];
        }
        curNode = curNode->next;
    }
    free(tempDegs);

    Symbol* newvars = new Symbol[nvar+1];
    newvars[0] = names[0];
    for (int i = 0; i < nvar; ++i) {
        newvars[varmap[i]+1] = names[i+1]; 
    }
    delete[] names;
    names = newvars;

    sortPolynomial(&poly, nvar);
} 

////////// BPASPolynomial ////////////////////////////////////

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator= (const SparseMultivariateRationalPolynomial& b) {
    if (this != &b) {
        freePolynomial(poly);
        poly = deepCopyPolynomial(b.poly, b.nvar);
        tailNode = NULL;
        nvar = b.nvar;
        delete[] names; 
        names = new Symbol[nvar+1];
        std::copy(b.names, b.names+nvar+1, names);
        slp = b.slp;
        delete constTerm;
        constTerm = NULL;
        if (b.constTerm != NULL) {
            constTerm = new RationalNumber(*(b.constTerm));
        }
    }
    return *this;
}

/**
 * Movement assignment: move b to be this polynomail.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator= (SparseMultivariateRationalPolynomial&& b) {
    if (this != &b) {
        freePolynomial(poly);
        poly = b.poly;
        b.poly = NULL;
        tailNode = b.tailNode;
        b.tailNode = NULL;

        nvar = b.nvar;
        b.nvar = 0;

        delete[] names;
        names = new Symbol[nvar+1];
        std::copy(b.names, b.names+nvar+1, names);

        slp = b.slp;
        b.slp.clear();

        delete constTerm;
        constTerm = b.constTerm;
        b.constTerm = NULL;
    }
    return *this;
}

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator= (const RationalNumber& r) {
    if (poly != NULL) {
        freePolynomial(poly);
        poly = NULL;
    }
    tailNode = NULL;

    if (constTerm != NULL) {
        *constTerm = r;
        return *this;
    }

    if (nvar == 0) {
        constTerm = new RationalNumber(r);
    } else {
        degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
        poly = addZeroTerm(NULL, degs);
    }
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
        if (nvar > 0) {
            return (b + this->poly->coef);
        } else {
            return (b + *constTerm);
        }
    }
    if (b.isConstant() != 0) {
        if (b.nvar > 0){
            return (*this + b.poly->coef);
        } else {
            return (*this + *(b.constTerm));
        }
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
        Node* sum = addPolynomials(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateRationalPolynomial ret(sum, superNvar, newnames);
        return ret;
    
    } else {
        Node* sum = addPolynomials(poly, b.poly, nvar);
        SparseMultivariateRationalPolynomial ret(sum, nvar, names);
        return ret;
    }
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator+ (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret += *this;
    return ret;
}

SparseMultivariateRationalPolynomial operator+ (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
    SparseMultivariateRationalPolynomial ret = a;
    ret += b;
    return ret;
}

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator+= (const SparseMultivariateRationalPolynomial& b) {
    slp.clear();
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
        if (nvar > 0) {
            return (negB + this->poly->coef);
        } else {
            return (negB + (*constTerm));
        }
    }
    if (b.isConstant() != 0) {
        if (b.nvar > 0) {
            return (*this - b.poly->coef);
        } else {
            return (*this - *(b.constTerm));
        }
    }

    SparseMultivariateRationalPolynomial negB = -b;
    return (*this + negB);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator- (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret -= *this;
    return ret;   
}

SparseMultivariateRationalPolynomial operator- (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
    SparseMultivariateRationalPolynomial ret = a;
    ret -= b;
    return ret;
}

SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator-= (const SparseMultivariateRationalPolynomial& b) {
    slp.clear();
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
        if (nvar > 0) {
            return (b * this->poly->coef);
        } else {
            return (b * *constTerm);
        }
    }
    if (b.isConstant() != 0) {
        if (b.nvar > 0) {
            return (*this * b.poly->coef);
        } else {
            return (*this * *(b.constTerm));
        }
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
        Node* prod = multiplyPolynomials(tempa.poly, tempb.poly, superNvar);
        SparseMultivariateRationalPolynomial ret(prod, superNvar, newnames);
        return ret;
    
    } else {
        Node* prod = multiplyPolynomials(poly, b.poly, nvar);
        SparseMultivariateRationalPolynomial ret(prod, nvar, names);
        return ret;
    }
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator* (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret *= *this;
    return ret;   
}

SparseMultivariateRationalPolynomial operator* (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
    SparseMultivariateRationalPolynomial ret = a;
    ret *= b;
    return ret;
}

/**
 * Update this by multiplying by the specified polynomail.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator*= (const SparseMultivariateRationalPolynomial& b) {
    slp.clear();
    *this = (*this * b);
    return *this;
}

/**
 * Divide *this by the specified polynomial.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator/ (const SparseMultivariateRationalPolynomial& b) const {
    SparseMultivariateRationalPolynomial q,r;
    this->divide(b, q, r);
    return q;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator/ (SparseMultivariateRationalPolynomial&& b) const {
    SparseMultivariateRationalPolynomial ret = b;
    ret /= *this;
    return ret;   
}

SparseMultivariateRationalPolynomial operator/ (SparseMultivariateRationalPolynomial&& a, const SparseMultivariateRationalPolynomial& b) {
    SparseMultivariateRationalPolynomial ret = a;
    ret /= b;
    return ret;
}

/**
 * Update *this by dividing by the specified polynomial.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator/= (const SparseMultivariateRationalPolynomial& b){
    slp.clear();
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

    if (nvar == 0) {
        SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
        ret.constTerm = new RationalNumber(*constTerm);
        for (int i = 2; i <= e; ++i) {
            *(ret.constTerm) *= *constTerm;
        }
        return ret;
    } 


    e = (e < 0) ? -e : e;
    Node* retNode = exponentiatePoly(poly, e, nvar);
    SparseMultivariateRationalPolynomial ret(retNode, nvar, names);
    return ret;
}

/**
 * Update *this by exponentiating this to the input integer.
 * Treats negative exponents as positive.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator^= (long long int e) {
    slp.clear();
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
        os << polyToString(this->poly, this->nvar, tempVars);
    } else if (this->constTerm != NULL) {
        os << *(this->constTerm);
    } else {
        os << "0";
    }
}


/**
 * Parse an input string into the current polynomial.
 */
std::istream& operator>> (std::istream& in, const SparseMultivariateRationalPolynomial& p) {
    std::cerr << "BPAS: SMQP: input stream not yet implemented" << std::endl;
}


RationalNumber SparseMultivariateRationalPolynomial::content() const {
    //TODO 
    std::cerr << "BPAS ERROR: SparseMultivariateRationalPolynomial::content() NOT YET IMPLEMENTED" << std::endl;
    return RationalNumber(1); 
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::primitivePart() const {
    //TODO
    std::cerr << "BPAS ERROR: SparseMultivariateRationalPolynomial::primitivePart() NOT YET IMPLEMENTED" << std::endl;
    return *this;
}



////////// BPASMultivariatePolynomial ////////////////////////////////////

/**
 * Get the number of variables in this polynomial.
 */
int SparseMultivariateRationalPolynomial::numberOfVariables() const {
    return nvar;
}

/**
 * Get the number of variables in this polynomial with non-zero degree.
 */
int SparseMultivariateRationalPolynomial::numberOfNonZeroVariables() const {
    if (nvar == 0 || isConstant()) {
        return 0;
    }

    Node* node = poly;
    degrees_t degs = node->degs;

    bool foundVar[nvar];
    foundVar[0] = (degs[0] > 0);
    int searchIdx = 1;
    while(node != NULL && searchIdx < nvar) {
        degs = node->degs;

        if (degs[searchIdx] > 0) {
            foundVar[searchIdx] = 1;
            while(searchIdx < nvar && foundVar[searchIdx] == 1) {
                ++searchIdx;
            }
        }

        for (int i = searchIdx; i < nvar; ++i) {
            if (degs[i] != 0) {
                foundVar[i] = 1;
            }
        }

        node = node->next;
    }

    int res = 0;
    for (int i = 0; i < nvar; ++i) {
        if (foundVar[i]) {
            ++res;
        }
    }
    return res;
}


/**
 * Get the number of non-zero terms 
 */
Integer SparseMultivariateRationalPolynomial::numberOfTerms() const {
    if (constTerm != NULL) {
        return 1;
    }
    return numberOfTermsNode(poly);
}

/** 
 * Total degree.
 */
Integer SparseMultivariateRationalPolynomial::degree() const {
    if (isConstant()) {
        return 0;
    }
    degree_t total = 0, totalMax = 0;
    Node* node = poly;
    degrees_t degs;
    while(node != NULL) {
        total = 0;
        degs = node->degs;
        for (int i = 0; i < nvar; ++i) {
            total += degs[i];
        }
        if (total > totalMax) {
            totalMax = total;
        }
        node = node->next;
    }

    return totalMax;
}

/**
 * Get the degree of a variable 
 */
Integer SparseMultivariateRationalPolynomial::degree(const Symbol& str) const {
    if (poly == NULL || nvar == 0) {
        return 0;
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
        return poly->degs[0];
    } else {
        Node* node = poly;
        int d = 0;
        while (node != NULL) {
            if (node->degs[strIdx] > d) {
                d = node->degs[strIdx];
            }
            node = node->next;
        }
        return d;
    }
}

/**
 * Get the leading coefficient 
 */
RationalNumber SparseMultivariateRationalPolynomial::leadingCoefficient() const {
    if (isZero()) {
        return RationalNumber(0);
    }

    if (nvar == 0) {
        return RationalNumber(*constTerm);
    }

    return RationalNumber(poly->coef);
}

RationalNumber SparseMultivariateRationalPolynomial::trailingCoefficient() const {
    if (isZero()) {
        return RationalNumber(0);
    }

    if (nvar == 0) {
        return RationalNumber(*constTerm);
    }

    Node* tailN = getTailNode();
    return RationalNumber(tailN->coef);
}


/**
 * Get a coefficient, given the exponent of each variable in d.
 * v is the number of variables in d. It is assumed that the first this.nvar 
 * variables of d match the variables of this polynomial
 */
RationalNumber SparseMultivariateRationalPolynomial::coefficient(int v, int* d) const {
    if (nvar == 0 && v == 0) {
        if (constTerm == NULL) {
            return 0;
        } else {
            return *constTerm;
        }
    }

    Node* node = poly;

    degrees_t degs = (degrees_t) calloc(nvar, sizeof(degrees_t));
    for(int i = 0; i < nvar; i++) {
        degs[i] = d[i];
    }

    while (node != NULL) {
        int cmp = compareExponentVectors(node->degs, degs, nvar);
        if (cmp == 0) {
            return RationalNumber(node->coef);
        } else if (cmp < 0) {
            break; //term with exponent vector d doesn't exist in this poly
        }
        node = node->next;
    }
    return RationalNumber(0, 1);
}

/**
 * Set a coefficient, given the exponent of each variable
 */
void SparseMultivariateRationalPolynomial::setCoefficient(int v, int* d, const RationalNumber& rn) {
    if (v != nvar) {
        std::cout << "BPAS: error, SMQP(" << nvar << "), but trying to setCoefficient with " << v << " variables." << std::endl;
        exit(1);
    }

    if (v == 0) {
        delete constTerm;
        constTerm = new RationalNumber(rn);
        return;
    }

    mpq_class r = rn.get_mpq();
    degrees_t degs = (degrees_t) malloc(sizeof(degree_t) * nvar);
    for (int i = 0; i < nvar; ++i) {
        degs[i] = d[i];
    }

    if (poly == NULL) {
        poly = addTerm(NULL, degs, r.get_mpq_t());
        return;
    }

    Node* prev = NULL;
    Node* node = poly;
    while (node != NULL) {
        int cmp = compareExponentVectors(node->degs, degs, nvar);
        if (cmp == 0) {
            mpq_set(node->coef, r.get_mpq_t());
            free(degs);
            return;
        }
        if (cmp < 0) {
            //then we need to insert d just before the current Node.
            if (prev == NULL) {
                //we insert as the head. 
                poly = addTerm(NULL, degs, r.get_mpq_t());
                poly->next = node;
    
            } else {
                Node* newNode = addTerm(NULL, degs, r.get_mpq_t());
                newNode->next = node;
                prev->next = newNode;
            }
            return;
        }

        prev = node;
        node = node->next;
    }

    //if we get here then we have went through the entire polynomail without finding
    //the matching term, or one which is smaller than d. Therefore insert d at the end. 
    //prev is now tail of the linked list.
    addTerm(prev, degs, r.get_mpq_t());
} 

/**
 * Set variables' name.
 */
void SparseMultivariateRationalPolynomial::setVariableNames (const std::vector<Symbol>& xs) {
    int ns = xs.size();
    if (ns != nvar) {
        std::cout << "BPAS: error, SMQP(" << nvar << "), but trying to setVariableNames with " << ns << " variables." << std::endl;
        exit(1);
    }

    if (ns == 0) {
        return;
    }

    if (names[0] == "1") {
        names[0] = "9";
        for (int i = 1; i < nvar+1 ; ++i) {
            names[i] = xs[i-1];
        }
        return;
    }

    bool isReorder = 1;
    int varmap[nvar];
    for (int i = 0; i < nvar; ++i) {
        bool found = 0;
        for (int j = 0; j < ns; ++j) {
            if (names[i+1] == xs[j]) {
                varmap[i] = j;
                found = 1;
                break;
            }
        }
        if (!found) {
            isReorder = 0;
        }
    }


    if (isReorder) {
        Symbol newvars[nvar+1];
        newvars[0] = "9";
        std::copy(xs.begin(), xs.end(), newvars+1);
        *this = this->expandVariables(nvar, newvars, varmap);
    } else {
        names[0] = "9";
        for (int i = 0; i < nvar; ++i) {
            names[i+1] = xs[i];
        }
    }
}

/**
 * Get variable names of variables with non-zero degree;
 */
std::vector<Symbol> SparseMultivariateRationalPolynomial::nonZeroVariables() const {
    if (nvar == 0 || isConstant()) {
        return std::vector<Symbol>();
    }

    Node* node = poly;
    degrees_t degs = node->degs;

    bool foundVar[nvar];
    for (int i = 0; i < nvar; ++i) {
        foundVar[i] = 0;
    }

    foundVar[0] = (degs[0] > 0);
    int searchIdx = 1;
    while(node != NULL && searchIdx < nvar) {
        degs = node->degs;

        if (degs[searchIdx] > 0) {
            foundVar[searchIdx] = 1;
            while(searchIdx < nvar && foundVar[searchIdx] == 1) {
                ++searchIdx;
            }
        }

        for (int i = searchIdx; i < nvar; ++i) {
            if (degs[i] > 0) {
                foundVar[i] = 1;
            }
        }

        node = node->next;
    }

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
std::vector<Symbol> SparseMultivariateRationalPolynomial::variables() const {
    std::vector<Symbol> varNames;
    for (int i = 0; i < nvar; ++i) {
        varNames.push_back(names[i+1]);
    }
    return varNames;   
}



////////// RecursivelyViewedPolynomial  ////////////////////////////////////

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::initial() const {
    if (nvar == 0) {
        return SparseMultivariateRationalPolynomial(0);
    }  
    return leadingCoefficientInVariable(leadingVariable());
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::mainDegree() const {
    if (isConstant()) {
        return 0;
    }

    for (int i = 0; i < nvar; ++i) {
        if (poly->degs[i] != 0) {
            return poly->degs[i];
        }
    }

    return 0;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::head() const {

    Symbol leadVar;
    int k = 0;
    int kdeg = 0;
    for (int i = 0; i < nvar; ++i) {
        if (poly->degs[i] != 0) {
            leadVar = names[i+1];
            k = i;
            kdeg = poly->degs[i];
            break;
        }
    }

    SparseMultivariateRationalPolynomial initial = leadingCoefficientInVariable(leadVar);
    Node* node = initial.poly;
    while (node != NULL) {
        degrees_t newDegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
        memcpy(newDegs+1, node->degs, sizeof(degree_t)*(nvar-1));
        newDegs[0] = 0; //incase k != 0
        newDegs[k] = kdeg;
        free(node->degs);
        node->degs = newDegs;

        node = node->next;
    }
    SparseMultivariateRationalPolynomial actualInitial(initial.poly, nvar, names);
    initial.poly = NULL;

    return actualInitial;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::tail() const {
    SparseMultivariateRationalPolynomial head = this->head();
    return (*this - head);
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
            RationalNumber tempa, tempb;
            if (nvar == 0) {
                tempa = *constTerm;
            } else {
                tempa = mpq_class(poly->coef);
            }
            if (b.nvar == 0) {
                tempb = *(b.constTerm);
            } else {
                tempb = mpq_class((b.poly)->coef);
            }
            return tempa == tempb;
        } else {
            return 0;
        }
    }

    std::vector<int> xs;
    bool isOrdered = isOrderedRing(b, xs);
    if (!isOrdered) { return 0; }

    int v = xs.size() / 2;
    Node* curA = poly, * curB = b.poly;
    while (curA != NULL && curB != NULL) {
        if (mpq_cmp(curA->coef, curB->coef) != 0) {
            return 0;
        }

        for (int i = 0; i < v; ++i) {
            if (xs[2*i] && xs[2*i+1] && (curA->degs[xs[2*i]-1] != curB->degs[xs[2*i+1]-1]) ) {
                return 0;
            } else if (!xs[2*i] && xs[2*i+1] && (curB->degs[xs[2*i+1]-1] != 0) ) {
                return 0;
            } else if (xs[2*i] && !xs[2*i+1] && (curA->degs[xs[2*i]-1] != 0) ) {
                return 0;
            }
        }

        curA = curA->next;
        curB = curB->next;
    }

    if (curA != NULL || curB != NULL) {
        return 0;
    }

    return 1;
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

    
    int* active = (int*) malloc(nvar*sizeof(int));
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

    Node* evalNode = evaluatePoly(poly, active, vals, nvar);
    return SparseMultivariateRationalPolynomial(evalNode, newNvar, newnames);
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
        if (b.nvar == 0) {
            q = *this / *(b.constTerm);
        } else {
            q = *this / b.poly->coef;
        }
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

        SparseMultivariateRationalPolynomial tempc = this->expandVariables(superNvar, newnames, varmap);
        SparseMultivariateRationalPolynomial tempb = b.expandVariables(superNvar, newnames, bvarmap);

        Node* qNode = NULL, * rNode = NULL;
        dividePolynomials(tempc.poly, tempb.poly, &qNode, &rNode, superNvar);

        q = SparseMultivariateRationalPolynomial(qNode, superNvar, newnames);
        r = SparseMultivariateRationalPolynomial(rNode, superNvar, newnames);

        return (rNode == NULL);

    } else {
        Node* qNode = NULL,* rNode = NULL;
        dividePolynomials(this->poly, b.poly, &qNode, &rNode, nvar);
        q = SparseMultivariateRationalPolynomial(qNode, nvar, names);
        r = SparseMultivariateRationalPolynomial(rNode, nvar, names);
        return (rNode == NULL);
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
    slp.clear();
    *this = (*this % b);
    return *this;
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::pseudoDivide(const SparseMultivariateRationalPolynomial& b, SparseMultivariateRationalPolynomial* quo, SparseMultivariateRationalPolynomial* mult, bool lazy) const {
    if (b.isZero()) {
        std::cerr << "BPAS: error, dividend is zero from SMQP." << std::endl;
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

            Node* cNode = tempc.poly;
            Node* bNode = tempb.poly;
            RecNode_t* recC = convertToRecursiveNode(cNode);
            RecNode_t* recB = convertToRecursiveNode(bNode);

            RecNode_t* recQ, *recR;
            Node* hPow;
            int e;
            pesudoDivide(recC, recB, &recQ, &recR, &e, &hPow, superNvar, lazy);
            cNode = convertFromRecursiveNode(recC);
            bNode = convertFromRecursiveNode(recB);
            Node* qNode = convertFromRecursiveNode(recQ);
            Node* rNode = convertFromRecursiveNode(recR);

            q = SparseMultivariateRationalPolynomial(qNode, superNvar, newvars);
            r = SparseMultivariateRationalPolynomial(rNode, superNvar, newvars);
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

        Node* cNode = tempc.poly;
        Node* bNode = tempb.poly;
        RecNode_t* recC = convertToRecursiveNode(cNode);
        RecNode_t* recB = convertToRecursiveNode(bNode);
        RecNode_t* recQ, *recR;
        Node* hPow;
        int e;
        pesudoDivide(recC, recB, &recQ, &recR, &e, &hPow, superNvar, lazy);
        cNode = convertFromRecursiveNode(recC);
        bNode = convertFromRecursiveNode(recB);
        Node* qNode = convertFromRecursiveNode(recQ);
        Node* rNode = convertFromRecursiveNode(recR);

        q = SparseMultivariateRationalPolynomial(qNode, superNvar, newnames);
        r = SparseMultivariateRationalPolynomial(rNode, superNvar, newnames);
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

            Node* cNode = cReordered.poly;
            Node* bNode = bReordered.poly;
            RecNode_t* recC = convertToRecursiveNode(cNode);
            RecNode_t* recB = convertToRecursiveNode(bNode);

            RecNode_t* recQ, *recR;
            Node* hPow;
            int e;
            pesudoDivide(recC, recB, &recQ, &recR, &e, &hPow, superNvar, lazy);
            cNode = convertFromRecursiveNode(recC);
            bNode = convertFromRecursiveNode(recB);
            Node* qNode = convertFromRecursiveNode(recQ);
            Node* rNode = convertFromRecursiveNode(recR);

            q = SparseMultivariateRationalPolynomial(qNode, nvar, newvars);
            r = SparseMultivariateRationalPolynomial(rNode, nvar, newvars);

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

        Node* cNode = this->poly;
        Node* bNode = b.poly;
        RecNode_t* recC = convertToRecursiveNode(cNode);
        RecNode_t* recB = convertToRecursiveNode(bNode);
        RecNode_t* recQ, *recR;
        Node* hPow;
        int e;
        pesudoDivide(recC, recB, &recQ, &recR, &e, &hPow, superNvar, lazy);
        cNode = convertFromRecursiveNode(recC);
        bNode = convertFromRecursiveNode(recB);
        Node* qNode = convertFromRecursiveNode(recQ);
        Node* rNode = convertFromRecursiveNode(recR);

        q = SparseMultivariateRationalPolynomial(qNode, superNvar, names);
        r = SparseMultivariateRationalPolynomial(rNode, superNvar, names);

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
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator+ (const mpq_class& c) const {
    if (isZero()) {
        if (nvar == 0) {
            SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
            ret.constTerm = new RationalNumber(c);
            return ret;
        } else {
            degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
            Node* node = addTerm(NULL, d, c.get_mpq_t());
            return SparseMultivariateRationalPolynomial(node, nvar, names);
        }
    }

    if (isConstant() && nvar == 0) {
        SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
        ret.constTerm = new RationalNumber(c);
        *(ret.constTerm) += *constTerm;
        return ret;
    }

    SparseMultivariateRationalPolynomial ret = *this;
    Node* rTail = ret.getTailNode();
    if (isZeroExponentVector(rTail->degs, nvar)) {
        mpq_add(rTail->coef, rTail->coef, c.get_mpq_t());
        if (mpq_sgn(rTail->coef) == 0) {
            ret.removeTailNode();
        }
    } else {
        degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
        ret.tailNode = addTerm(rTail, d, c.get_mpq_t());
    }

    return ret;
}


/**
 * Update *this by adding r
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator+= (const mpq_class& c) {    
    slp.clear();
    if (isZero()) {
        if (nvar == 0) {
            if (constTerm == NULL) {
                constTerm = new RationalNumber(0);
            }
            *constTerm += mpq_class(c);
        } else {
            degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
            poly = addTerm(NULL, d, c.get_mpq_t());
            tailNode = poly;
        }
        return *this;
    }
    if (isConstant() && nvar == 0) {
        *constTerm += mpq_class(c);
        return *this;
    }

    Node* t = getTailNode();
    if (isZeroExponentVector(t->degs, nvar)) {
        mpq_add(t->coef, t->coef, c.get_mpq_t());
        if (mpq_sgn(t->coef) == 0) {
            removeTailNode();
        }
    } else {
        degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
        tailNode = addTerm(t, d, c.get_mpq_t());
    }
    return *this;
}

/**
 * Subtract the ratNum_t r from *this.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator- (const mpq_class& c) const {
    if (isZero()) {
        ratNum_t negR;
        mpq_init(negR);
        mpq_neg(negR, c.get_mpq_t());
        if (nvar == 0) {
            SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
            ret.constTerm = new RationalNumber(negR);
            mpq_clear(negR);
            return ret;
        }
        degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
        Node* node = addTerm(NULL, d, negR);
        mpq_clear(negR);
        return SparseMultivariateRationalPolynomial(node, nvar, names);
    }
    if (isConstant() && nvar == 0) {
        SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
        ret.constTerm = new RationalNumber(*constTerm);
        *(ret.constTerm) -= c;
        return ret;
    }

    SparseMultivariateRationalPolynomial ret = *this;
    Node* t = ret.getTailNode();
    if (isZeroExponentVector(t->degs, nvar)) {
        mpq_sub(t->coef, t->coef, c.get_mpq_t());
        if (mpq_sgn(t->coef) == 0) {
            ret.removeTailNode();
        }
    } else {
        ratNum_t negR;
        mpq_init(negR);
        mpq_neg(negR, c.get_mpq_t());
        degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
        ret.tailNode = addTerm(t, d, negR);
        mpq_clear(negR);        
    }
    return ret;
}

/**
 * Update *this by subtracting ratNum_t r.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator-= (const mpq_class& c) {
    slp.clear();
    if (isZero()) {
        if (nvar == 0) {
            if (constTerm == NULL) {
                constTerm = new RationalNumber(0);
            }
            *constTerm -= c;
        } else {
            ratNum_t negR;
            mpq_init(negR);
            mpq_neg(negR, c.get_mpq_t());
            degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
            poly = addTerm(NULL, d, negR);
            tailNode = poly;
            mpq_clear(negR);
        }
        return *this;
    }
    if (isConstant() && nvar == 0) {
        *constTerm -= c;
        return *this;
    }

    Node* t = getTailNode();
    if (isZeroExponentVector(t->degs, nvar)) {
        mpq_sub(t->coef, t->coef, c.get_mpq_t());
        if (mpq_sgn(t->coef) == 0) {
            removeTailNode();
        }
    } else {
        ratNum_t negR;
        mpq_init(negR);
        mpq_neg(negR, c.get_mpq_t());
        degrees_t d = (degrees_t) calloc(nvar, sizeof(degree_t));
        tailNode = addTerm(t, d, negR);
        mpq_clear(negR);
    }
    return *this;
}

/**
 * Multiply *this by ratNum_t r.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator* (const mpq_class& c) const {
    if (c == 0 || isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }

    SparseMultivariateRationalPolynomial ret = *this;

    if (c == 1) {
        return ret;
    }
    if (isConstant() && nvar == 0) {
        *(ret.constTerm) *= c;
        return ret;
    }

    Node* node = ret.poly;
    while (node != NULL) {
        mpq_mul(node->coef, node->coef, c.get_mpq_t());
        node = node->next;
    }
    return ret;
}

/**
 * Update *this by multiplying by ratNum_t r.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator*= (const mpq_class& c) {
    slp.clear();

    if (isZero()) {
        return *this;
    }

    if (c == 0) {
        freePolynomial(poly);
        poly = NULL;
        tailNode = NULL;
        return *this;
    }

    if (c == 1) {
        return *this;
    }
    if (isConstant() && nvar == 0) {
        *constTerm *= c;
        return *this;
    }

    Node* node = poly;
    while (node != NULL) {
        mpq_mul(node->coef, node->coef, c.get_mpq_t());
        node = node->next;
    }
    return *this;
}

/**
 * Divide *this by ratNum_t r.
 */
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::operator/ (const mpq_class& c) const {
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
    if (isConstant() && nvar == 0) {
        SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
        ret.constTerm = new RationalNumber(*constTerm);
        *(ret.constTerm) /= c;
        return ret;       
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
        if (b.constTerm != NULL) {
            RationalNumber quo = RationalNumber(r) / *(b.constTerm);
            SparseMultivariateRationalPolynomial ret(b.nvar);
            ret.constTerm = new RationalNumber(quo);
            return ret;
        }

        if (mpq_sgn(b.poly->coef) == 0) {
            std::cout << "BPAS: error, dividend is zero from SMQP." << std::endl;
            exit(1);
        }
        degrees_t d = (degrees_t) calloc(b.nvar, sizeof(degree_t));
        Node* rNode = addTerm(NULL, d, r);
        mpq_div(rNode->coef, rNode->coef, b.poly->coef);
        return SparseMultivariateRationalPolynomial(rNode, b.nvar, b.names);
    } else {
        SparseMultivariateRationalPolynomial ret(NULL, b.nvar, b.names);
        return ret;
    }
}

/** 
 * Update *this by dividing by ratNum_t r.
 */
SparseMultivariateRationalPolynomial& SparseMultivariateRationalPolynomial::operator/= (const mpq_class& c) {

    slp.clear();
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

    if (isConstant() && nvar == 0) {
        *constTerm /= c;
        return *this;
    }

    ratNum_t rInv;
    mpq_init(rInv);
    mpq_inv(rInv, c.get_mpq_t());
    *this *= rInv;
    mpq_clear(rInv);
    return *this;
}

/**
 * Get the leading variable, that is, the highest-order variable with positive degree
 * of this polynomial. 
 * returns the leading variable or the empty string if this polynomial has zero variables.
 */
Symbol SparseMultivariateRationalPolynomial::leadingVariable() const {
    if (nvar == 0 || isConstant()) {
        return "";
    }

    for (int i = 0; i < nvar; ++i) {
        if (poly->degs[i] != 0) {
            return names[i+1];
        }
    }

    return "";
}

/**
 * Get the degree of this polynomial w.r.t the leading variable.
 */
Integer SparseMultivariateRationalPolynomial::leadingVariableDegree() const {
    if (nvar == 0 || isConstant()) {
        return 0;
    }

    for (int i = 0; i < nvar; ++i) {
        if (poly->degs[i] != 0) {
            return poly->degs[i];
        }
    }

    return 0;
}

/**
 * Is the contant term zero.
 */
bool SparseMultivariateRationalPolynomial::isConstantTermZero() const {
    if (isZero()) {
        return 1;
    }

    if (constTerm != NULL) {
        return (*constTerm == 0);
    }

    Node* last = poly;
    while (last->next != NULL) {
        last = last->next;
    }
    return !isZeroExponentVector(last->degs, nvar);
}

SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::leadingCoefficientInVariable (const Symbol& x, int* e) const {
    int deg = 0;
    if (e != NULL) {
        *e = 0;
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
        if (constTerm != NULL) {
            r.constTerm = new RationalNumber(*constTerm);
        } else if (v == 0) {
            r.constTerm = new RationalNumber(poly->coef);
        } else {
            degrees_t degs = (degrees_t) calloc(v, sizeof(degree_t));
            r.poly = addTerm(NULL, degs, poly->coef);
        }
        return r;
    }
    
    --k; //recall names index is +1 of degs index
    Node* cur = poly;
    Node* rTail = NULL;
    while (cur != NULL) {
        if (cur->degs[k] >= deg) {
            if (cur->degs[k] > deg) {
                deg = cur->degs[k];
                freePolynomial(r.poly);
                r.poly = NULL;
                rTail = NULL;
            }
            degrees_t degs = (degrees_t) malloc(sizeof(degree_t)*v);
            for (int i = 0; i < k; ++i) {
                degs[i] = cur->degs[i];
            }
            for (int i = k; i < v; ++i) {
                degs[i] = cur->degs[i+1];
            }
            rTail = addTerm(rTail, degs, cur->coef);
            if (r.poly == NULL) {
                r.poly = rTail;
            }
        }
        cur = cur->next;
    }

    //probably not needed
    sortPolynomial(&r.poly, v);

	if (v == 0) {
		r.constTerm = new RationalNumber(r.poly->coef);
		freePolynomial(r.poly);
		r.poly = NULL;	
	}

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
    r.setVariableName(x.toString()); 

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
        Node* node = poly;
        while (node != NULL) {
            SparseMultivariateRationalPolynomial coef(RationalNumber(node->coef), 0);
            r.setCoefficient(node->degs[0], coef);            
            node = node->next;
        }        
        return r;
    }

    int v = nvar - 1;
    int d = 0;
    Node* curNode = poly;
    while (curNode != NULL) {
        d = curNode->degs[k-1] > d ? curNode->degs[k-1] : d;
        curNode = curNode->next;
    }

    SparseMultivariateRationalPolynomial mpolys[d+1];
    Node* tails[d+1];
    for (int i = 0; i < d+1; ++i) {
        tails[i] = NULL;
    }

    Symbol newnames[v+1];
    for (int i = 0; i < k; ++i) {
        newnames[i] = names[i];
    }
    for (int i = k; i < v+1; ++i) {
        newnames[i] = names[i+1];
    }
    //setup names and nvar of each of the polys
    for (int i = 0; i <= d; ++i) {
        mpolys[i].nvar = v;
        delete[] mpolys[i].names;
        mpolys[i].names = new Symbol[v+1];
        std::copy(newnames, newnames+v+1, mpolys[i].names);
    }

    --k;
    curNode = poly;
    int curD;
    while (curNode != NULL) {
        curD = curNode->degs[k];
        degrees_t degs = (degrees_t) malloc(sizeof(degree_t)*v);
        for (int i = 0; i < k; ++i) {
            degs[i] = curNode->degs[i];
        }
        for (int i = k; i < v; ++i) {
            degs[i] = curNode->degs[i+1];
        }
        tails[curD] = addTerm(tails[curD], degs, curNode->coef);
        if (mpolys[curD].poly == NULL) {
            mpolys[curD].poly = tails[curD];
        }
        curNode = curNode->next;
    }

    for (int i = 0; i <= d; ++i) {
        if (mpolys[i].isZero()) {
            continue;
        }
        sortPolynomial(&mpolys[i].poly, mpolys[i].nvar);
        r.setCoefficient(i, mpolys[i]);
    }

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
    if (isConstant() && nvar == 0) {
        *constTerm *= -1;
    }
    negatePolynomial(poly);
}

/** 
 * Get a copy of this polynomial. 
 */ 
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::deepCopy() const {
    if (isZero()) {
        return SparseMultivariateRationalPolynomial(NULL, nvar, names);
    }
    if (isConstant() && nvar == 0) {
        SparseMultivariateRationalPolynomial ret(NULL, nvar, names);
        ret.constTerm = new RationalNumber(*constTerm);
        return ret;
    }

    Node* copy = deepCopyPolynomial(poly, nvar);
    return SparseMultivariateRationalPolynomial(copy, nvar, names);
}

int checkDegs(degrees_t a, degrees_t b, int vars) {
    for (int i = 0; i < vars; ++i) {
        if (a[i] != b[i])
            return i;
    }
    return -1;
}

int firstNonZero(degrees_t d, int vars) {
    for (int i = vars-1; i >= 0; --i) {
        if (d[i])
            return i;
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
        if (nvar == 0) {
            r.a.c = new RationalNumber(*constTerm);
        } else {
            r.a.c = new RationalNumber(poly->coef);
        }
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

    Node* node = poly;
    Node* prev = poly;
    if (poly != NULL) {
        degrees_t d = (degrees_t) malloc(sizeof(degree_t)*nvar);
        for (int i = 0; i < nvar; ++i) {
            d[i] = poly->degs[i];
        }
        std::vector<Node*> nodeStack;
        std::vector<int> indexStack;

        bool isOp = 1;
        SLPRepresentation elem;
        while (node != NULL) {
            degrees_t curd = node->degs;
            int k = checkDegs(d, curd, nvar);
            if (k > -1) {
                bool isStart = 0;

                // Mutiply smaller variate than the current one
                for (int j = nvar-1; j > k; --j) {
                    if (d[j] > 0 && !nodeStack.empty()) {
                        degrees_t prevd = nodeStack[nodeStack.size()-1]->degs;
                        while(checkDegs(prevd, d, nvar) == j) {
                            elem.op = 0;
                            elem.type = 3;
                            elem.a.i = indexStack[indexStack.size()-1];
                            elem.b = slp.size() - 1;
                            slp.push_back(elem);

                            nodeStack.pop_back();
                            indexStack.pop_back();
                            if (!nodeStack.empty()) {
                                prevd = nodeStack[nodeStack.size()-1]->degs;
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
                if (!nodeStack.empty()) {
                    degrees_t prevd = nodeStack[nodeStack.size()-1]->degs;
                    while(checkDegs(prevd, curd, nvar) == k) {
                        elem.op = 0;
                        elem.type = 3;
                        elem.a.i = indexStack[indexStack.size()-1];
                        elem.b = slp.size() - 1;
                        slp.push_back(elem);

                        nodeStack.pop_back();
                        indexStack.pop_back();
                        if (!nodeStack.empty()) {
                            prevd = nodeStack[nodeStack.size()-1]->degs;
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
                    nodeStack.push_back(prev);
                    indexStack.push_back(slp.size()-1);
                }
            }

            // Multiply or add this coefficient
            k = firstNonZero(d, nvar);
            if (isOp) {
                elem.op = 1;
                elem.type = 0;
                elem.a.c = new RationalNumber(node->coef);
                // elem.a  = i;
                elem.b = k;
                slp.push_back(elem);

                d[k]--;
            }
            else {
                elem.op = 0;
                elem.type = 1;
                elem.a.c = new RationalNumber(node->coef);
                // elem.a = i;
                elem.b = slp.size() - 1;
                slp.push_back(elem);
            }
            isOp = 0;

            prev = node;
            node = node->next;
        }

        int k = firstNonZero(d, nvar);
        if (k > -1) {
            if (!nodeStack.empty()) {
                degrees_t prevd = nodeStack[nodeStack.size()-1]->degs;
                while (checkDegs(prevd, d, nvar) == k) {
                    elem.op = 0;
                    elem.type = 3;
                    elem.a.i = indexStack[indexStack.size()-1];
                    elem.b = slp.size() - 1;
                    slp.push_back(elem);

                    nodeStack.pop_back();
                    indexStack.pop_back();
                    if (!nodeStack.empty()) {
                        prevd = nodeStack[nodeStack.size()-1]->degs;
                    } else {
                        break;
                    }
                }
            }

            for (int i = k; i >= 0; --i) {
                if (nodeStack.size() > 0) {
                    degrees_t prevd = nodeStack[nodeStack.size()-1]->degs;
                    while(checkDegs(prevd, d, nvar) == i) {
                        elem.op = 0;
                        elem.type = 3;
                        elem.a.i = indexStack[indexStack.size()-1];
                        elem.b = slp.size() - 1;
                        slp.push_back(elem);

                        nodeStack.pop_back();
                        indexStack.pop_back();
                        if (!nodeStack.empty()) {
                            prevd = nodeStack[nodeStack.size()-1]->degs;
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
        d = poly->degs[0];
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
        constTerm = new RationalNumber((rand() % (coefBound-1)) + 1);
    } else {
        poly = buildRandomPoly(numvar, nterms, coefBound, sparsity, includeNeg);
    }
}


/**
 * Construct an ExprTreeNode (well, multiple) which represents a single
 * term. That is, a coefficient and a monomial.
 */
ExprTreeNode* exprTreeNodeFromNode(Node* n, int nvar, Symbol* vars) {
    if (n == NULL) {
        return new ExprTreeNode(0l);
    }

    degrees_t degs = n->degs;
    ExprTreeNode* t = new ExprTreeNode(mpq_class(n->coef));
    for (int i = 0; i < nvar; ++i) {
        if (degs[i] > 1) {
            //TODO pass symbol directly to expression tree
            ExprTreeNode* var = new ExprTreeNode(vars[i]);
            ExprTreeNode* num = new ExprTreeNode(degs[i]);
            ExprTreeNode* exp = ExprTreeNode::combineExprTreeNodes(var, num, EXPR_EXP);
            t = ExprTreeNode::combineExprTreeNodes(t, exp, EXPR_MULT);
        } else if (degs[i] == 1) {
            ExprTreeNode* var = new ExprTreeNode(vars[i]);
            t = ExprTreeNode::combineExprTreeNodes(t, var, EXPR_MULT);
        }
    }

    return t;
}

ExpressionTree SparseMultivariateRationalPolynomial::convertToExpressionTree() const {
    if (constTerm != NULL) {
        ExprTreeNode* r = new ExprTreeNode(constTerm->get_mpq());
        ExpressionTree t(r);
        return t;
    }
    if (poly == NULL) {
        ExprTreeNode* r = new ExprTreeNode(0l);
        ExpressionTree t(r);
        return t; 
    }

    Node* node = poly;
    ExprTreeNode* prev = exprTreeNodeFromNode(node, nvar, &names[1]);
    node = node->next;
    while (node != NULL) {
        ExprTreeNode* thisNode = exprTreeNodeFromNode(node, nvar, &names[1]);
        prev = ExprTreeNode::combineExprTreeNodes(prev, thisNode, EXPR_ADD);
        node = node->next;
    }


    return ExpressionTree(prev);
}

/****************
* SubResultantChain
*****************/

		

/**
 * Subresultant Chain
 * Return the list of subresultants
 *
 * @param q: The other sparse univariate polynomial
 **/
std::vector<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::subresultantChain (const SparseMultivariateRationalPolynomial& q) const {
	Symbol v = q.leadingVariable();
	if (v != leadingVariable()) {
		std::cout << "BPAS: error, cannot compute subresultant chain if leading variable of input is different from leading variable of the current object." << std::endl;
		exit(1);
	}

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


/**
 * Subresultant Chain GCD
 * Return the last non-zero subresultant of the current polynomial and the input polynomial if the resultant is zero and return 1 otherwise
 *
 * @param q: The other sparse univariate polynomial
 **/
SparseMultivariateRationalPolynomial SparseMultivariateRationalPolynomial::subresultantGCD (const SparseMultivariateRationalPolynomial& q) const {
	Symbol v = q.leadingVariable();
	if (v != leadingVariable()) {
		std::cout << "BPAS: error, cannot compute subresultant gcd if leading variable of input is different from leading variable of the current object." << std::endl;
		exit(1);
	}

	SparseMultivariateRationalPolynomial one;
	one.one();
	if (isConstant() != 0 || q.isConstant() != 0)
		return one;
//	unsigned long long start;
//	float elapsed;
//	startTimer(&start);
	// This code needs an efficient direct conversion to DUQP to be competitive
//	if (numberOfVariables() == 1) {
//		SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> a,b;
//		a = convertToSUP(leadingVariable());
//		b = q.convertToSUP(leadingVariable());
//		DenseUnivariateRationalPolynomial f = a.convertToDUQP();
//		DenseUnivariateRationalPolynomial g = b.convertToDUQP();
//		std::cerr << "f = " << f << std::endl;
//		std::cerr << "g = " << g << std::endl;
//		DenseUnivariateRationalPolynomial z = f.gcd(g);
//		std::cerr << "gcd(f,g) = " << z << std::endl;
//		SparseUnivariatePolynomial<SparseMultivariateRationalPolynomial> r(z);
//		SparseMultivariateRationalPolynomial s(r);
////		return s;
//	}
//	stopTimer(&start,&elapsed);
//	std::cerr << "DUQP gcd time: " << elapsed << std::endl;

//	startTimer(&start);
	std::vector<SparseMultivariateRationalPolynomial> src = subresultantChain(q);
//	for (int i=0; i<src.size(); ++i)
//		std::cerr << "src[" << i << "] = " << src[i] << std::endl;
	if (!src[0].isZero()) {
//		stopTimer(&start,&elapsed);
//		std::cerr << "SUP gcd time: " << elapsed << std::endl;
		return one;
	}
	else {
		if (src[2].degree() == src[1].degree()) {
//			stopTimer(&start,&elapsed);
//			std::cerr << "SUP gcd time: " << elapsed << std::endl;
			return src[2];
		}
		else {
//			stopTimer(&start,&elapsed);
//			std::cerr << "SUP gcd time: " << elapsed << std::endl;
			return src[1];
		}
	}
	
}
