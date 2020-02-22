#ifndef _RATIONALFUNCTION_INTEGRATIONPOSTPROCESSING_H
#define _RATIONALFUNCTION_INTEGRATIONPOSTPROCESSING_H

// Complex Polynomial Routines
extern SparseUnivariatePolynomial<ComplexRationalNumber> RNPolyToCRNPoly(SparseUnivariatePolynomial<RationalNumber> &A);
extern SparseUnivariatePolynomial<RationalNumber> CRNPolyRealPart(SparseUnivariatePolynomial<ComplexRationalNumber> &A);
extern SparseUnivariatePolynomial<RationalNumber> CRNPolyImaginaryPart(SparseUnivariatePolynomial<ComplexRationalNumber> &A);
extern SparseUnivariatePolynomial<ComplexRationalNumber> conjugate(SparseUnivariatePolynomial<ComplexRationalNumber> &A);

//template <class UnivariatePolynomialOverField>
//ComplexRationalNumber _complexEvaluate(UnivariatePolynomialOverField &P, ComplexRationalNumber &c);
//template <class UnivariatePolynomialOverField>
//SparseUnivariatePolynomial<RationalNumber> _realEvaluateCoefficients(SparseUnivariatePolynomial<UnivariatePolynomialOverField> &P, mpq_class c);
//template <class UnivariatePolynomialOverField>
//SparseUnivariatePolynomial<ComplexRationalNumber> _complexEvaluateCoefficients(SparseUnivariatePolynomial<UnivariatePolynomialOverField> *P, ComplexRationalNumber &c);
template <class UnivariatePolynomialOverRealField, class RealField>
extern void _arctan2ToArctan(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &B, RealField coef, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn);
template <class UnivariatePolynomialOverRealField, class RealField>
extern void _logToReal(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &S, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<
RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn, int prec, bool PROFILING);
template <class UnivariatePolynomialOverRealField, class RealField>
extern void _logToReal(std::vector< std::vector<ComplexRationalNumber> > &E, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &S, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<
RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn1, std::vector<UnivariatePolynomialOverRealField> *Atn2, int prec, bool PROFILING);

#endif
