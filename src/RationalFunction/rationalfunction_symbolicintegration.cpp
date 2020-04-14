#include "../../include/ring.h"
#include "../../include/RationalNumberPolynomial/urpolynomial.h"
#include "../../include/RingPolynomial/upolynomial.h"
#include "../../include/RationalFunction/rationalfunction_symbolicintegration.h"
#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"

template <class UnivariatePolynomialOverField, class Field>
void _hermiteReduce(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, std::vector<UnivariatePolynomialOverField> *g, std::vector<UnivariatePolynomialOverField> *h){
	/* hermiteReduce(A, D, g, h):                         */
	/* Hermite reduction of a rational function           */
	/*                                                    */
	/* Given a UPoF A and D in Q[x] with D nozero and     */
	/* coprime with A, return g, h in Q(x), i.e. rational */
	/* functions in x with rational coefficients, such    */
	/* that A/D = dg/dx + h and h has a squarefree        */
	/* denominator.                                       */
	/*                                                    */
	/* The rational functions g and h are given in parts  */
	/* in the following way:                              */
	/*  1) g, a sum of quotients of polynomials, is given */
	/*     as a UPoF vector, where the numerators and     */
	/*     denominators are listed pairwise with          */
	/*     numerator preceding denominator;               */
	/*  2) h, a single quotient of polynomials, is given  */
	/*     as a UPoF vector, where the numerator and      */
	/*     denominator are the first and second elements  */

	assert(!D.isZero());
	UnivariatePolynomialOverField a,b,c,dm,dmstar,dm2,dstar,deriv;
	a.setVariableName(A.variable());
	b.setVariableName(A.variable());
	c.setVariableName(A.variable());
	dm.setVariableName(A.variable());
	dmstar.setVariableName(A.variable());
	dm2.setVariableName(A.variable());
	dstar.setVariableName(A.variable());
	deriv.setVariableName(A.variable());

	a = A;
	g->clear();
	deriv = D;
	deriv.differentiate(1);
	//std::cout << "deriv = " << deriv << std::endl;
	dm = D.gcd(deriv); // deriv = dD/dx
	//std::cout << "dm = " << dm << std::endl;
	dstar = D;
	dstar /= dm;
	//std::cout << "dstar = " << dstar << std::endl;

	while (dm.degree() > 0){
		deriv = dm;
		deriv.differentiate(1);
		//std::cout << "deriv = " << deriv << std::endl;
		dm2 = dm.gcd(deriv); // deriv = d(dm)/dx
		//std::cout << "dm2 = " << dm2 << std::endl;
		dmstar = dm;
		dmstar /= dm2;
		//std::cout << "dmstar = " << dmstar << std::endl;
		deriv *= -dstar;
		deriv /= dm;
		//std::cout << "-deriv*dstar/dm = " << deriv << std::endl;
		_extendedEuclidean<UnivariatePolynomialOverField,Field>(deriv,dmstar,a,&b,&c); // deriv = d(dm)/dx(-dmstar/dm)
		a = -b;
		//std::cout << "b = " << b << ", c = " << c << std::endl;
		a.differentiate(1);
		a *= dstar;
		a /= dmstar;
		a += c;
		//std::cout << "a*dstar/dmstar + c = " << a << std::endl;
		if (!b.isZero()){
			// canonicalize must be called for these rational functions
			g->push_back(b);
			g->push_back(dm);
		}
		dm = dm2;
		//std::cout << "dm = " << dm << std::endl;
	}
	// canonicalize must be called for this rational function
	h->push_back(a);
	h->push_back(dstar);
}

template <class UnivariatePolynomialOverField, class Field>
SparseUnivariatePolynomial<UnivariatePolynomialOverField> _UPoFToSUP(UnivariatePolynomialOverField& P, Symbol variable){
	UnivariatePolynomialOverField upofTemp;
	upofTemp.setVariableName(variable);
	SparseUnivariatePolynomial<UnivariatePolynomialOverField> A;
	A.setVariableName(P.variable());
	Field fTemp(RationalNumber(0));

	for (int i=0; i <= P.degree().get_si(); i++){
		fTemp = P.coefficient(i);
		if (fTemp != 0){
			upofTemp.setCoefficient(0,fTemp);
			A.setCoefficient(i,upofTemp);
		}
	}
	return A;
}

template <class UnivariatePolynomialOverField, class Field>
void _prepRothsteinTragerResultant(SparseUnivariatePolynomial<UnivariatePolynomialOverField> *a, SparseUnivariatePolynomial<UnivariatePolynomialOverField> *d, SparseUnivariatePolynomial<UnivariatePolynomialOverField> *amtdd, UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, Symbol indeterminate){
	/* prepRothsteinTragerResultant:                  */
	/*                                                */
	/* Converts the numerator A and denominator D of  */
	/* the log part of the integral from UPoF to SUP  */
	/* and returns A-t*dD/dx as an SUP variable,      */
	/* where t is an indeterminate.                   */

	int i;
	//Field n(D.degree()); uncomment this to use Davenport divisor
	Field fTemp(RationalNumber(0));
	UnivariatePolynomialOverField upofTemp;
	Symbol variable(A.variable());
	if (variable == indeterminate) {
		if (variable == "%")
			indeterminate = "t";
		else
			indeterminate = "%";
	}
	upofTemp.setVariableName(indeterminate);

	// Create SUP version of A
	*a = _UPoFToSUP<UnivariatePolynomialOverField,Field>(A,indeterminate);

	// Create SUP version of D
	*d = _UPoFToSUP<UnivariatePolynomialOverField,Field>(D,indeterminate);

	// Create second argument for Rothstein-Trager resultant A-t*dD/dx
	// To fix exponential coefficient growth for 1/(1+x^n), use (n*A-t*dD/dx)/n
	// with the "Davenport divisor" n to prevent expression swell
	amtdd->setVariableName(variable);
	UnivariatePolynomialOverField dD(D);
	dD.differentiate(1);
	for (i=0; i <= dD.degree().get_si(); i++){
		fTemp = dD.coefficient(i);
		if (fTemp != 0){
			upofTemp.zero();
			upofTemp = (upofTemp + fTemp << 1);
			amtdd->setCoefficient(i,upofTemp);
		}
	}
	*amtdd = -*amtdd;

	//upofTemp.zero();
	//upofTemp.setCoefficient(0,n);
	//supTemp *= upofTemp; // uncomment these three lines to use Davenport divisor

	*amtdd += *a;

	//*amtdd /= upofTemp; //  uncomment this to use Davenport divisor
}

template <class UnivariatePolynomialOverField, class Field>
void _computeLogArguments(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, std::vector<UnivariatePolynomialOverField> *U, SparseUnivariatePolynomial<UnivariatePolynomialOverField> &amtdd, SparseUnivariatePolynomial<UnivariatePolynomialOverField> &d, UnivariatePolynomialOverField &D, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > &R, std::vector<UnivariatePolynomialOverField> &Q, Symbol indeterminate, bool PROFILING){

	int i,j,m;
	Field fTemp;
	UnivariatePolynomialOverField upofTemp,upofTemp2;
	upofTemp.setVariableName(D.variable());
	upofTemp2.setVariableName(D.variable());
	std::vector<UnivariatePolynomialOverField> B;
	Factors<UnivariatePolynomialOverField> Fs;
	UnivariatePolynomialOverField onePoly;
	onePoly.setVariableName(Q[0].variable());
	onePoly.one();

	// Profiling variables
	unsigned long long start;
	float elapsed = 0;

	for (i=1; i < Q.size(); i++){
		if (Q[i].degree() > 0){
			if (D.degree() == i){
				S->push_back(d);
				U->push_back(Q.at(i));
			}
			else {
				if (i < R.size()){
					// check that there is a subresultant of degree i
					if (R.at(i).degree() == i){
						m = i;
						upofTemp = R.at(m).leadingCoefficient();
						S->push_back(R.at(m));
						U->push_back(Q.at(i));
						Fs = upofTemp.squareFree();
						B.clear();
						B.push_back(Fs.ringElement());
						int curExp = 1;
						for (int i = 0; i < Fs.size(); ++i) {
							for (int k = curExp; k < Fs[i].second; ++k) {
								B.push_back(onePoly);
							}
							curExp = Fs[i].second;
							B.push_back(Fs[i].first);
						}
						upofTemp.one();
						for (j=1; j<B.size(); j++){

							if (PROFILING) {
								startTimer(&start);
							}

							upofTemp2 = B.at(j).gcd(Q.at(i));

							if (PROFILING){
								stopTimerAddElapsed(&start,&elapsed);
							}

							upofTemp2 = upofTemp2^j;
							upofTemp *= upofTemp2;
						}
						S->back() /= upofTemp;
					}
				}
				else if (R.size() == 1){
					upofTemp = amtdd.leadingCoefficient();
					S->push_back(amtdd);
					U->push_back(Q.at(i));
					Fs = upofTemp.squareFree();
					B.clear();
					B.push_back(Fs.ringElement());
					int curExp = 1;
					for (int i = 0; i < Fs.size(); ++i) {
						for (int k = curExp; k < Fs[i].second; ++k) {
							B.push_back(onePoly);
						}
						curExp = Fs[i].second;
						B.push_back(Fs[i].first);
					}
					upofTemp.one();
					for (j=1; j<B.size(); j++){

						if (PROFILING) {
							startTimer(&start);
						}

						upofTemp2 = B.at(j).gcd(Q.at(i));

						if (PROFILING){
							stopTimerAddElapsed(&start,&elapsed);
						}

						upofTemp2 = upofTemp2^j;
						upofTemp *= upofTemp2;
					}
					S->back() /= upofTemp;
				}
			}
		}
	}
	if (PROFILING)
		std::cout << "\t\tcomputeLogArguments gcd\t\t" << elapsed << std::endl;

	for (i=0; i < U->size(); i++){
		for (j=0; (S->at(i).degree()+1) > j ; j++){
			upofTemp = S->at(i).coefficient(j);
			_remainder<UnivariatePolynomialOverField,Field>(upofTemp,U->at(i),&upofTemp2);
			upofTemp2.setVariableName(indeterminate);
			S->at(i).setCoefficient(j,upofTemp2);
		}
		//fTemp = commonZFactor(S->at(i));
		//upofTemp.zero();
		//upofTemp.setCoefficient(0,fTemp);
		//S->at(i) /= upofTemp;
	}
}

template <class UnivariatePolynomialOverField, class Field>
void _simplifyLogArguments(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, std::vector< UnivariatePolynomialOverField > *U, Symbol indeterminate){
	//Field fTemp;
	DenseUnivariateRationalPolynomial upofTemp;
	DenseUnivariateRationalPolynomial upofTemp2;
	if (U->size() > 0) {
		upofTemp.setVariableName(S->at(0).variable());
		for (int i=0; i < U->size(); i++){
			for (int j=0; (S->at(i).degree()+1) > j; j++){
				upofTemp = S->at(i).coefficient(j);
				_remainder<UnivariatePolynomialOverField,Field>(upofTemp,U->at(i),&upofTemp2);
				upofTemp2.setVariableName(indeterminate);
				S->at(i).setCoefficient(j,upofTemp2);
			}
			//fTemp = commonZFactor(S->at(i));
			//upofTemp.zero();
			//upofTemp.setCoefficient(0,fTemp);
			//S->at(i) /= upofTemp;
		}
	}
}

template <class UnivariatePolynomialOverField, class Field>
void _integrateRationalFunctionLogPart(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, std::vector<UnivariatePolynomialOverField> *U, UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, bool PROFILING){
	/* S = integrateRationalLogPart(A, D, U):             */
	/* Integration of log part of a rational function     */
	/*                                                    */
	/* Given UPoF A,D (in Q[x]) with deg(A) <= deg(D),    */
	/* D nonzero, squarefree and coprime with A, return   */
	/* the integral of A/D in the following form:         */
	/*  1) an SUP vector S containing the arguments of    */
	/*     the logarithms in the integral of A/D, where   */
	/*     each argument is a polynomials in the main     */
	/*     variable of A and an indeterminate; and        */
	/*  2) a corresponding UPoF vector U, where the roots */
	/*     each polynomial determine the values at which  */
	/*     the indeterminates in the corresponding        */
	/*     entries of S must be evaluated to evaluate the */
	/*     integral                                       */

	assert(A.degree()<=D.degree());

	Symbol indeterminate("c");
	Symbol variable(A.variable());
	D.setVariableName(variable);
	int i;
	UnivariatePolynomialOverField resultant;
	resultant.setVariableName(variable);
	std::vector<UnivariatePolynomialOverField> Q;
	SparseUnivariatePolynomial<UnivariatePolynomialOverField> a,d,amtdd;
	a.setVariableName(variable);
	d.setVariableName(variable);
	amtdd.setVariableName(variable);
	std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > R;

	// Profiling variables
	unsigned long long start;
	float elapsed = 0;

	if (PROFILING){
		startTimer(&start);
	}

	_prepRothsteinTragerResultant<UnivariatePolynomialOverField,Field>(&a, &d, &amtdd, A, D, indeterminate);

	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t\tPrep Rothstein\t\t\t" << elapsed << std::endl;
	}

	// Compute subresultant chain
	d.setVariableName(Symbol("%"));
	amtdd.setVariableName(Symbol("%"));
	// handle case of trivial resultant
	if (amtdd.degree() >= 1) {
		if (PROFILING){
			startTimer(&start);
		}

		R = d.subresultantChain(amtdd);

		if (PROFILING){
			stopTimer(&start,&elapsed);
			std::cout << "\t\tSubresultant\t\t\t" << elapsed << std::endl;
		}
	}
	else
		R.push_back(amtdd);
	d.setVariableName(variable);
	amtdd.setVariableName(variable);
	for (i=0; i<R.size(); i++){
		R.at(i).setVariableName(variable);
	}

	resultant = R.at(0).coefficient(0);
	//std::cout << "resultant = " << resultant << std::endl;
	SparseUnivariatePolynomial<UnivariatePolynomialOverField> supTemp;
	supTemp.setCoefficient(0,resultant);
	//Field fTemp = commonZFactor(supTemp);
	//std::cout << "content(resultant) = " << fTemp << std::endl;

	if (PROFILING){
		startTimer(&start);
	}

	Factors<UnivariatePolynomialOverField> f = resultant.squareFree();

	Q.clear();
	Q.push_back(f.ringElement());
	int curExp = 1;
	UnivariatePolynomialOverField onePoly;
	onePoly.one();
	for (int i = 0; i < f.size(); ++i) {
		for (int k = curExp; k < f[i].second; ++k) {
			Q.push_back(onePoly);
		}
		curExp = f[i].second;
		Q.push_back(f[i].first);
	}
	// for (int i = 0; i < f.size(); ++i) {
	// 	Q.push_back(f[i].first);
	// }

	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t\tSquare free fact\t\t" << elapsed << std::endl;
		startTimer(&start);
	}

	_computeLogArguments<UnivariatePolynomialOverField,Field>(S, U, amtdd, d, D, R, Q, indeterminate, PROFILING);

	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t\tLog arguments\t\t\t" << elapsed << std::endl;
	}
}

template <class UnivariatePolynomialOverField, class Field>
void _integrateRationalFunctionRationalPart(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, UnivariatePolynomialOverField *P, std::vector<UnivariatePolynomialOverField> *G, UnivariatePolynomialOverField *R, UnivariatePolynomialOverField *H, bool PROFILING){
	/* integrateRationalFunctionRationalPart(A, D, P, G, U, S):  */
	/* Integrate rational part of a rational function            */
	/*                                                           */
	/* Given UPoF A,D (in F[x]) with deg(A) < deg(D),            */
	/* return the rational part of the integral of A/D in the    */
	/* and the integrand of the log part in thefollowing form:   */
	/*  1) a UPoF P, the polynomial term in the integral,        */
	/*     if it exists;                                         */
	/*  2) a UPoF vector G containing the rational               */
	/*     function terms in the integral, where the             */
	/*     numerators and denominators are listed                */
	/*     pairwise with numerator preceding denominator;        */
	/*  3) the integrand R/H of the log part, where R and H are  */
	/*     polynomials and H is squarefree.                      */

	//assert(A.degree()<D.degree());

	std::vector<UnivariatePolynomialOverField> K;
							// this vector is needed to store the rational function part
							// of the Hermite reduction of A and D, which will be processed
							// to compute the logarithmic part of the integral of A/D

	// set main variable
	Symbol mainVariable;
	mainVariable = A.variable();

	// initialize helper and temporary variables
	int i,j;
	UnivariatePolynomialOverField upofTemp;
	upofTemp.setVariableName(mainVariable);
	SparseUnivariatePolynomial<UnivariatePolynomialOverField> supTemp;
	supTemp.setVariableName(mainVariable);

	// clear vector container
	G->clear();

	// Profiling variables
	unsigned long long start;
	float elapsed = 0;

	// extract polynomial part from the integrand

	_euclideanDivide<UnivariatePolynomialOverField,Field>(A,D,P,R);

	// compute the rational function part of the integral of R/D by Hermite reduction,
	// which is returned in G, along with the derivative of the rest of the integral,
	// viz., the derivatives of the polynomial and logarithmic parts, which is returned
	// as a single rational function as a pair of polynomials in H

	if (PROFILING) {
		startTimer(&start);
	}

	_hermiteReduce<UnivariatePolynomialOverField,Field>(*R,D,G,&K);

	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t\tHermite reduce\t" << elapsed << std::endl;
	}

	// process the derivative of the rest of the integral (polynomial and logarithmic parts),
	// contained in H=(num(H),den(H)), to extract the derivative of the polynomial part Q,
	// and the remainder R/den(H), the derivative of the logarithmic part, which is integrated by
	// integrateRationalLogPart

	UnivariatePolynomialOverField Q;
	Q.setVariableName(mainVariable);

	if (PROFILING){
		startTimer(&start);
	}

	// compute num(H)/den(H) = Q + R/den(H)
	_euclideanDivide<UnivariatePolynomialOverField,Field>(K.at(0),K.at(1),&Q,R);

	// add quotient to polynomial part
	*P += Q;

	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t\tEuclidean div\t" << elapsed << std::endl;
	}

	// integrate polynomial part, if any
	if (!(P->isZero()));
		*P = P->integral();

	// return the log part of the integral as R/H (R already set by second call to _euclideanDivide)
	*H = K.at(1);

}

template <class UnivariatePolynomialOverField, class Field>
void _initializeRationalIntegration(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, bool PROFILING){

	// set main variable
	if (A.variable() != D.variable()){
		std::cout << "BPAS: error, trying to integrate Field[" << A.variable() << "]/Field[" << D.variable() << "]." << std::endl;
		    exit(1);
	}
	Symbol mainVariable;
	mainVariable = A.variable();

	UnivariatePolynomialOverField upofTemp;
	upofTemp.setVariableName(mainVariable);

	// Profiling variables
	unsigned long long start;
	float elapsed = 0;

	if (PROFILING){
		startTimer(&start);
	}

	// remove any common gcd of A and D
	upofTemp = A.gcd(D);
	A /= upofTemp;
	D /= upofTemp;

	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\tCanonicalize\t" << elapsed << std::endl;
	}
}

template <class UnivariatePolynomialOverField, class Field>
void _integrateRationalFunction(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, UnivariatePolynomialOverField *P, std::vector<UnivariatePolynomialOverField> *G, std::vector<UnivariatePolynomialOverField> *U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, bool PROFILING){
	/* integrateRationalFunction(A, D, P, G, U, S):       */
	/* Integration of a rational function                 */
	/*                                                    */
	/* Given UPoF A,D (in F[x]) with deg(A) < deg(D),     */
	/* return the integral of A/D in the following form:  */
	/*  1) a UPoF P, the polynomial term in the integral, */
	/*     if it exists;                                  */
	/*  2) a UPoF vector G containing the rational        */
	/*     function terms in the integral, where the      */
	/*     numerators and denominators are listed         */
	/*     pairwise with numerator preceding denominator; */
	/*  3) an SUP vector S containing the arguments of    */
	/*     the logarithms in the integral of A/D, where   */
	/*     each argument is a polynomials in the main     */
	/*     variable of A and an indeterminate; and        */
	/*  4) a corresponding UPoF vector U, where the roots */
	/*     each polynomial determine the values at which  */
	/*     the indeterminates in the corresponding        */
	/*     entries of S must be evaluated to evaluate the */
	/*     integral                                       */

	//assert(A.degree()<D.degree());

	// variable checking and canonicalization
	_initializeRationalIntegration<UnivariatePolynomialOverField,Field>(A,D,PROFILING);

	// set main variable
	Symbol mainVariable;
	mainVariable = A.variable();

	// initialize helper and temporary variables
	//int i,j;
	//SparseUnivariatePolynomial<UnivariatePolynomialOverField> supTemp;
	//supTemp.setVariableName(mainVariable);

	// clear vector containers
	U->clear();
	S->clear();

	// Profiling variables
	unsigned long long start;
	float elapsed = 0;

	// compute the rational part of the integral of R/D by Hermite reduction, which is
	// returned as G.at(0)/G.at(1), and the polynomial part, returned in P, along with the
	// derivative of the rest of the integral, viz., the integrand of the logarithmic
	// part, which is returned as a single rational function as a pair of polynomials
	// R and H

	UnivariatePolynomialOverField R; // stores numerator of log part of integral of A/D
	UnivariatePolynomialOverField H; // stores squarefee denominator of the log part
	R.setVariableName(mainVariable);
	H.setVariableName(mainVariable);

	if (PROFILING) {
		std::cout << "\tRational function rational part" << std::endl;
		std::cout << "\t------------------------------" << std::endl;
		startTimer(&start);
	}

	_integrateRationalFunctionRationalPart<UnivariatePolynomialOverField,Field>(A,D,P,G,&R,&H,PROFILING);

	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}

	// compute integral of R/H

	if (!R.isZero()) {

		if (PROFILING){
			std::cout << "\tRational function log part" << std::endl;
			std::cout << "\t------------------------------" << std::endl;
			startTimer(&start);
		}

		_integrateRationalFunctionLogPart<UnivariatePolynomialOverField,Field>(S,U,R,H,PROFILING);

		if (PROFILING){
			stopTimer(&start,&elapsed);
			std::cout << "\t------------------------------" << std::endl;
			std::cout << "\t\t\t" << elapsed << std::endl;
		}
	}
}

// to avoid linking errors
template void _hermiteReduce<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, std::vector<DenseUnivariateRationalPolynomial> *g, std::vector<DenseUnivariateRationalPolynomial> *h);
template void _hermiteReduce<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, std::vector< SparseUnivariatePolynomial<RationalNumber> > *g, std::vector< SparseUnivariatePolynomial<RationalNumber> > *h);
//template void _hermiteReduce<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(SparseUnivariatePolynomial<ComplexRationalNumber> &A, SparseUnivariatePolynomial<ComplexRationalNumber> &D, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *g, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *h);

template void _prepRothsteinTragerResultant<DenseUnivariateRationalPolynomial,RationalNumber>(SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> *a, SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> *d, SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> *amtdd, DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, Symbol indeterminate);

//template void _computeLogArguments<DenseUnivariateRationalPolynomial,RationalNumber>(std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > *S, std::vector<DenseUnivariateRationalPolynomial> *U, SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> &amtdd, SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> &d, DenseUnivariateRationalPolynomial &D, std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > &R, std::vector<DenseUnivariateRationalPolynomial> &Q, std::string indeterminate, bool PROFILING);

template void _integrateRationalFunctionLogPart<DenseUnivariateRationalPolynomial,RationalNumber>(std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > *S, std::vector<DenseUnivariateRationalPolynomial> *U, DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, bool PROFILING);
template void _integrateRationalFunctionLogPart<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<RationalNumber> > > *S, std::vector< SparseUnivariatePolynomial<RationalNumber> > *U, SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, bool PROFILING);
//template void _integrateRationalFunctionLogPart<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<ComplexRationalNumber> > > *S, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *U, SparseUnivariatePolynomial<ComplexRationalNumber> &A, SparseUnivariatePolynomial<ComplexRationalNumber> &D, bool PROFILING);


//template void _integrateRationalFunctionRationalPart<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, DenseUnivariateRationalPolynomial *P, std::vector<DenseUnivariateRationalPolynomial> *G, DenseUnivariateRationalPolynomial *R, DenseUnivariateRationalPolynomial *H, bool PROFILING);
//template void _integrateRationalFunctionRationalPart<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, SparseUnivariatePolynomial<RationalNumber> *P, std::vector< SparseUnivariatePolynomial<RationalNumber> > *G, SparseUnivariatePolynomial<RationalNumber> *R, SparseUnivariatePolynomial<RationalNumber> *H, bool PROFILING);
//template void _integrateRationalFunctionRationalPart<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(SparseUnivariatePolynomial<ComplexRationalNumber> &A, SparseUnivariatePolynomial<ComplexRationalNumber> &D, SparseUnivariatePolynomial<ComplexRationalNumber> *P, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *G, SparseUnivariatePolynomial<ComplexRationalNumber> *R, SparseUnivariatePolynomial<ComplexRationalNumber> *H, bool PROFILING);

template void _integrateRationalFunction<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, DenseUnivariateRationalPolynomial *P, std::vector<DenseUnivariateRationalPolynomial> *G, std::vector<DenseUnivariateRationalPolynomial> *U, std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > *S, bool PROFILING);
template void _integrateRationalFunction<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, SparseUnivariatePolynomial<RationalNumber> *P, std::vector< SparseUnivariatePolynomial<RationalNumber> > *G, std::vector< SparseUnivariatePolynomial<RationalNumber> > *U, std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<RationalNumber> > > *S, bool PROFILING);
//template void _integrateRationalFunction<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(SparseUnivariatePolynomial<ComplexRationalNumber> &A, SparseUnivariatePolynomial<ComplexRationalNumber> &D, SparseUnivariatePolynomial<ComplexRationalNumber> *P, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *G, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *U, std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<ComplexRationalNumber> > > *S, bool PROFILING);
