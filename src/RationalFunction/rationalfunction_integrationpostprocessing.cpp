#include "../../include/polynomial.h"
#include "../../include/ring.h"
#include "../../include/RationalNumberPolynomial/urpolynomial.h"
#include "../../include/RingPolynomial/upolynomial.h"
#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"
#include "../../include/RationalFunction/rationalfunction_integrationpostprocessing.h"
#include "../../include/RationalFunction/rationalfunction_complexrationalnumberordering.h"


template <class UnivariatePolynomialOverRealField>
ComplexRationalNumber _complexEvaluate(UnivariatePolynomialOverRealField &P, ComplexRationalNumber &c) {
	ComplexRationalNumber ComplexRationalNumberTemp;
	ComplexRationalNumber output;
	output.set(P.coefficient(P.degree().get_si()),0);
	for (int i=P.degree().get_si()-1; i>-1; i--){
		ComplexRationalNumberTemp.setRealPart(P.coefficient(i));
		output *= c;
		output += ComplexRationalNumberTemp;
	}
	return output;
}

template <class UnivariatePolynomialOverRealField, class RealField>
UnivariatePolynomialOverRealField _realEvaluateCoefficients(SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> &P, RationalNumber c) {
	int i = 0;
	RealField fTemp;
	UnivariatePolynomialOverRealField upofTemp,output;
	output.setVariableName(P.variable());
	
	upofTemp = P.coefficient(P.degree().get_si());
	fTemp = upofTemp.evaluate(c);
	output = output + fTemp << P.degree().get_si();
	for (i=0;i<P.degree().get_si();i++) {
		if (!P.coefficient(i).isZero()) {
			upofTemp = P.coefficient(i);
			fTemp = upofTemp.evaluate(c);
			output.setCoefficient(i,fTemp);
		}
	}
	return output;
}

/*template <class UnivariatePolynomialOverField>
SparseUnivariatePolynomial<RationalNumber> _realEvaluateCoefficients(SparseUnivariatePolynomial<UnivariatePolynomialOverField> &P, RationalNumber c) {
	int i = 0;
	RationalNumber RationalNumberTemp = 0;
	UnivariatePolynomialOverField upofTemp;
	SparseUnivariatePolynomial<RationalNumber> output;
	output.setVariableName(P.variable());
	
	for (i=0;i<=P.degree();i++) {
		if (!P.coefficient(i).isZero()) {
			upofTemp = P.coefficient(i);
			RationalNumberTemp = upofTemp.evaluate(c);
			output.setCoefficient(i,RationalNumberTemp);
		}
	}
	return output;
}*/

template <class UnivariatePolynomialOverRealField>
SparseUnivariatePolynomial<ComplexRationalNumber> _complexEvaluateCoefficients(SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> &P, ComplexRationalNumber &c) {
	int i = 0;
	ComplexRationalNumber ComplexRationalNumberTemp;
	UnivariatePolynomialOverRealField upofTemp;
	SparseUnivariatePolynomial<ComplexRationalNumber> output;
	output.setVariableName(P.variable());
	
	for (i=0;i<=P.degree().get_si();i++) {
		if (!P.coefficient(i).isZero()) {
			upofTemp = P.coefficient(i);
			ComplexRationalNumberTemp = _complexEvaluate<UnivariatePolynomialOverRealField>(upofTemp,c);
			output.setCoefficient(i,ComplexRationalNumberTemp);
		}
	}
	return output;
}

template <class UnivariatePolynomialOverRealField,class RealField>
void _fixZeros(UnivariatePolynomialOverRealField *A, double precision){
	RealField zero;
	double val = 0;
	double prec = precision*precision;
	for (int i=0; i<A->degree().get_si()+1; i++){
		val = A->coefficient(i).get_d();
		if (val < 0)
			val = -val;
		if (val < prec)
			A->setCoefficient(i,zero);
	}
}

template <class RealField>
void _fixZero(RealField c, double precision){
	RealField zero;
	double val = 0;
	val = c.get_d();
	if (val < 0)
		val = -val;
	if (val < precision*precision)
		c = zero;
}

SparseUnivariatePolynomial<ComplexRationalNumber> RNPolyToCRNPoly(SparseUnivariatePolynomial<RationalNumber> &A){
	SparseUnivariatePolynomial<ComplexRationalNumber> output;
	output.setVariableName(A.variable());
	RationalNumber rTemp;
	ComplexRationalNumber cTemp;
	for (int i=0; i<A.degree().get_si()+1; i++){
		rTemp = A.coefficient(i);
		if (!rTemp.isZero()){
			cTemp.setRealPart(rTemp);
			output.setCoefficient(i,cTemp);
		}
	}
	return output;
}

SparseUnivariatePolynomial<RationalNumber> CRNPolyRealPart(SparseUnivariatePolynomial<ComplexRationalNumber> &A){
	SparseUnivariatePolynomial<RationalNumber> output;
	output.setVariableName(A.variable());
	RationalNumber rTemp;
	for (int i=0; i<A.degree().get_si()+1; i++){
		rTemp = A.coefficient(i).realPart();
		if (!rTemp.isZero())
			output.setCoefficient(i,rTemp);
	}
	return output;
}

SparseUnivariatePolynomial<RationalNumber> CRNPolyImaginaryPart(SparseUnivariatePolynomial<ComplexRationalNumber> &A){
	SparseUnivariatePolynomial<RationalNumber> output;
	output.setVariableName(A.variable());
	RationalNumber rTemp;
	for (int i=0; i<A.degree().get_si()+1; i++){
		rTemp = A.coefficient(i).imaginaryPart();
		if (!rTemp.isZero())
			output.setCoefficient(i,rTemp);
	}
	return output;
}

SparseUnivariatePolynomial<ComplexRationalNumber> conjugate(SparseUnivariatePolynomial<ComplexRationalNumber> &A){
	ComplexRationalNumber c;
	SparseUnivariatePolynomial<ComplexRationalNumber> B;
	B.setVariableName(A.variable());
	for (int i=0; i<A.degree().get_si()+1; i++) {
		c = A.coefficient(i);
		if (!c.isZero()) {
			c = c.conjugate();
			B.setCoefficient(i,c);
		}
	}
	return B;
}

template <class UnivariatePolynomialOverRealField, class RealField>
UnivariatePolynomialOverRealField _realSUPToUPoF(SparseUnivariatePolynomial<RationalNumber> &A) {
	RealField fTemp,zero;
	UnivariatePolynomialOverRealField output;
	output.setVariableName(A.variable());
	
	fTemp = A.coefficient(A.degree().get_si());
	output = output + fTemp << A.degree().get_si();
	for (int i=0; i<A.degree().get_si(); i++){
		fTemp = A.coefficient(i);
		if (fTemp != zero)
			output.setCoefficient(i,fTemp);
	}
	return output;
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _arctan2ToArctan(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &B, RealField coef, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn){
	/* _arctan2ToArctan(A,B,coef,atn,Atn1):                      */
	/* Convert two argument arctangent coef*arctan(A,B) to a     */
	/* series of one argument arctangent terms of the form       */
	/* atn[i]*arctan(Atn[i]) using Rioboo's recursive algorithm. */
	
	UnivariatePolynomialOverRealField upofTemp,upofTemp2,C,D,G;
	upofTemp.setVariableName(A.variable());
	C.setVariableName(A.variable());
	D.setVariableName(A.variable());
	G.setVariableName(A.variable());
	
	_remainder<UnivariatePolynomialOverRealField,RealField>(A,B,&upofTemp);
	if (upofTemp.isZero()){
		if (A.degree() > B.degree()){
			upofTemp = A;
			upofTemp /= B;
			atn->push_back(coef);
			Atn->push_back(upofTemp);
		}
		return;
	}
	if (A.degree() < B.degree()){
		upofTemp = -B;
		_arctan2ToArctan<UnivariatePolynomialOverRealField,RealField>(upofTemp,A,coef,atn,Atn);
		return;
	}
	upofTemp = -A;
	_extendedEuclidean<UnivariatePolynomialOverRealField,RealField>(B,upofTemp,&D,&C,&G);
	upofTemp = D;
	upofTemp *= A;
	upofTemp2 = C;
	upofTemp2 *= B;
	if (upofTemp2.degree() >= upofTemp.degree()){ // avoid segmentation fault from DenseUnivariateRationalPolynomial
		upofTemp2 += upofTemp;
		upofTemp = upofTemp2;
	}
	else{
		upofTemp += upofTemp2;
	}
	upofTemp /= G; // upofTemp = (AD + BC)/G
	if (upofTemp.degree() > 0) {
		atn->push_back(coef);
		Atn->push_back(upofTemp);
	}
	_arctan2ToArctan<UnivariatePolynomialOverRealField,RealField>(D,C,coef,atn,Atn);
	return;
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _logToReal(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &S, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<
RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn, int prec, bool PROFILING){
	/* logToReal(E,S,lg,Lg,atn,Atn,prec):                       */
	/* Convert formal sums to log and one argument arctan terms */
	/*                                                          */
	
	double precision = pow(2.0,-prec);
	RealField a,b,fTemp;
	UnivariatePolynomialOverRealField A,B;
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> upofTemp;
	SparseUnivariatePolynomial<ComplexRationalNumber> complexSupTemp;
	SparseUnivariatePolynomial<RationalNumber> realSupTemp;
	std::vector<ComplexRationalNumber> r;
	ComplexRationalNumber ComplexRationalNumberTemp;
	
	// Profiling variables
	unsigned long long start;
	float elapsed_a(0), elapsed_b(0), elapsed_c(0), elapsed_d(0);
	
	for (int k=0; k<E.size(); k++){
		r = E.at(k);
		sort(r.begin(),r.end(),ComplexRationalNumberOrdering(prec));
		for (int j=0; j<E.at(k).size(); j++){
			upofTemp = S.at(k);
			a = r.at(j).realPart();
			b = r.at(j).imaginaryPart();
			if (b == 0){
				if (a != 0){
					A = _realEvaluateCoefficients<UnivariatePolynomialOverRealField,RealField>(upofTemp,a);	
									// note that this code may have to be edited
									// to avoid generating large coefficients
					fTemp = A.leadingCoefficient();
					
					if (PROFILING){
						startTimer(&start);
					}
						
					A /= fTemp;
					
					if (PROFILING){
						stopTimerAddElapsed(&start,&elapsed_d);
					}
					
					//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&A,precision);
					lg->push_back(a);
					Lg->push_back(A);
				}
			}
			else if (a == 0){
				if (PROFILING){
					startTimer(&start);
				}
					
					complexSupTemp = _complexEvaluateCoefficients<UnivariatePolynomialOverRealField>(upofTemp,r.at(j));
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_a);
					startTimer(&start);
				}
				
				realSupTemp = CRNPolyRealPart(complexSupTemp);
				A = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				realSupTemp = CRNPolyImaginaryPart(complexSupTemp);
				B = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				fTemp = (abs(A.leadingCoefficient())+abs(B.leadingCoefficient()))/2;
				//cout << "fTemp = " << fTemp << endl;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_b);
					startTimer(&start);
				}
					
				A /= fTemp;
				B /= fTemp;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_d);
				}
				
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&A,precision);
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&B,precision);
				fTemp = 2*b;
				
				if (PROFILING){
					startTimer(&start);
				}
					
				_arctan2ToArctan<UnivariatePolynomialOverRealField,RealField>(A,B,fTemp,atn,Atn);
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_c);
				}
				
				j++;
			}
			else {
				
				if (PROFILING){
					startTimer(&start);
				}
					
					complexSupTemp = _complexEvaluateCoefficients<UnivariatePolynomialOverRealField>(upofTemp,r.at(j));
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_a);
					startTimer(&start);
				}
				
				realSupTemp = CRNPolyRealPart(complexSupTemp);
				A = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				realSupTemp = CRNPolyImaginaryPart(complexSupTemp);
				B = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				fTemp = (abs(A.leadingCoefficient())+abs(B.leadingCoefficient()))/2;
				//cout << "fTemp = " << fTemp << endl;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_b);
					startTimer(&start);
				}
					
				A /= fTemp;
				B /= fTemp;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_d);
				}
				
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&A,precision);
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&B,precision);
				fTemp = 2*b;
				
				if (PROFILING){
					startTimer(&start);
				}
					
				_arctan2ToArctan<UnivariatePolynomialOverRealField,RealField>(A,B,fTemp,atn,Atn);
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_c);
					startTimer(&start);
				}
				
				A *= A;
				B *= B;
				A += B;
				A /= A.leadingCoefficient();
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_d);
				}
				
				lg->push_back(a);
				Lg->push_back(A);
				j++;
			}
		}
	}
	
	if (PROFILING){
		std::cout << "\t\tComplex subs\t" << elapsed_a << std::endl;
		std::cout << "\t\tRe and Im parts\t" << elapsed_b << std::endl;
		std::cout << "\t\tarctan2ToArctan\t" << elapsed_c << std::endl;
		std::cout << "\t\tpoly arith\t" << elapsed_d << std::endl;
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _logToReal(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &S, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<
RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn1, std::vector<UnivariatePolynomialOverRealField> *Atn2, int prec, bool PROFILING){
	/* logToReal(E,S,lg,Lg,atn,Atn,prec):                       */
	/* Convert formal sums to log and one argument arctan terms */
	/*                                                          */
	
	double precision = pow(2.0,-prec);
	RealField a,b,fTemp;
	UnivariatePolynomialOverRealField A,B;
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> upofTemp;
	SparseUnivariatePolynomial<ComplexRationalNumber> complexSupTemp;
	SparseUnivariatePolynomial<RationalNumber> realSupTemp;
	std::vector<ComplexRationalNumber> r;
	ComplexRationalNumber ComplexRationalNumberTemp;
	
	// Profiling variables
	unsigned long long start;
	float elapsed_a = 0, elapsed_b = 0, elapsed_d = 0;
	
	for (int k=0; k<E.size(); k++){
		r = E.at(k);
		sort(r.begin(),r.end(),ComplexRationalNumberOrdering(prec));
		for (int j=0; j<E.at(k).size(); j++){
			upofTemp = S.at(k);
			a = r.at(j).realPart();
			b = r.at(j).imaginaryPart();
			if (b == 0){
				if (a != 0){
					A = _realEvaluateCoefficients<UnivariatePolynomialOverRealField,RealField>(upofTemp,a);	
									// note that this code may have to be edited
									// to avoid generating large coefficients
					fTemp = A.leadingCoefficient();
					
					if (PROFILING){
						startTimer(&start);
					}
						
					A /= fTemp;
					
					if (PROFILING){
						stopTimerAddElapsed(&start,&elapsed_d);
					}
					
					//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&A,precision);
					lg->push_back(a);
					Lg->push_back(A);
				}
			}
			else if (a == 0){
				if (PROFILING){
					startTimer(&start);
				}
					
					complexSupTemp = _complexEvaluateCoefficients<UnivariatePolynomialOverRealField>(upofTemp,r.at(j));
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_a);
					startTimer(&start);
				}
				
				realSupTemp = CRNPolyRealPart(complexSupTemp);
				A = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				realSupTemp = CRNPolyImaginaryPart(complexSupTemp);
				B = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				fTemp = (abs(A.leadingCoefficient())+abs(B.leadingCoefficient()))/2;
				//cout << "fTemp = " << fTemp << endl;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_b);
				}
				
				if (PROFILING){
					startTimer(&start);
				}
					
				A /= fTemp;
				B /= fTemp;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_d);
				}
				
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&A,precision);
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&B,precision);
				fTemp = 2*b;
				atn->push_back(fTemp);
				Atn1->push_back(A);
				Atn2->push_back(B);
				
				j++;
			}
			else {
				
				if (PROFILING){
					startTimer(&start);
				}
					
					complexSupTemp = _complexEvaluateCoefficients<UnivariatePolynomialOverRealField>(upofTemp,r.at(j));
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_a);
					startTimer(&start);
				}
				
				realSupTemp = CRNPolyRealPart(complexSupTemp);
				A = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				realSupTemp = CRNPolyImaginaryPart(complexSupTemp);
				B = _realSUPToUPoF<UnivariatePolynomialOverRealField,RealField>(realSupTemp);
				fTemp = (abs(A.leadingCoefficient())+abs(B.leadingCoefficient()))/2;
				//cout << "fTemp = " << fTemp << endl;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_b);
					startTimer(&start);
				}
					
				A /= fTemp;
				B /= fTemp;
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_d);
				}
				
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&A,precision);
				//_fixZeros<UnivariatePolynomialOverRealField,RealField>(&B,precision);
				fTemp = 2*b;
				atn->push_back(fTemp);
				Atn1->push_back(A);
				Atn2->push_back(B);
				
				if (PROFILING){
					startTimer(&start);
				}
				
				A *= A;
				B *= B;
				A += B;
				A /= A.leadingCoefficient();
				
				if (PROFILING){
					stopTimerAddElapsed(&start,&elapsed_d);
				}
				
				lg->push_back(a);
				Lg->push_back(A);
				j++;
			}
		}
	}
	
	if (PROFILING){
		std::cout << "\tComplex subs\t" << elapsed_a << std::endl;
		std::cout << "\tRe and Im parts\t" << elapsed_b << std::endl;
		std::cout << "\tpoly arith\t" << elapsed_d << std::endl;
	}
}

// to avoid linking errors
template void _arctan2ToArctan<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &B, RationalNumber coef, std::vector<RationalNumber> *atn, std::vector<DenseUnivariateRationalPolynomial> *Atn);
template void _arctan2ToArctan<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &B, RationalNumber coef, std::vector<RationalNumber> *atn, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn);
//template void _arctan2ToArctan<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(SparseUnivariatePolynomial<ComplexRationalNumber> &A, SparseUnivariatePolynomial<ComplexRationalNumber> &B, ComplexRationalNumber coef, std::vector<ComplexRationalNumber> *atn, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *Atn);

template void _logToReal<DenseUnivariateRationalPolynomial,RationalNumber>(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > &S, std::vector<RationalNumber> *lg, std::vector<DenseUnivariateRationalPolynomial> *Lg, std::vector<RationalNumber> *atn, std::vector<DenseUnivariateRationalPolynomial> *Atn, int prec, bool PROFILING);
template void _logToReal<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<RationalNumber> > > &S, std::vector<RationalNumber> *lg, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Lg, std::vector<RationalNumber> *atn, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn, int prec, bool PROFILING);

template void _logToReal<DenseUnivariateRationalPolynomial,RationalNumber>(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > &S, std::vector<RationalNumber> *lg, std::vector<DenseUnivariateRationalPolynomial> *Lg, std::vector<RationalNumber> *atn, std::vector<DenseUnivariateRationalPolynomial> *Atn1, std::vector<DenseUnivariateRationalPolynomial> *Atn2, int prec, bool PROFILING);
template void _logToReal<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<RationalNumber> > > &S, std::vector<RationalNumber> *lg, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Lg, std::vector<RationalNumber> *atn, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn1, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn2, int prec, bool PROFILING);
