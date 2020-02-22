#ifndef _RATIONALFUNCTION_INTEGRATIONPRINTING_H
#define _RATIONALFUNCTION_INTEGRATIONPRINTING_H


//template <class RealField>
//extern std::string fieldToFPS(RealField a);
template <class UnivariatePolynomialOverField, class Field>
extern void _printFormalIntegral(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, UnivariatePolynomialOverField &P, std::vector<UnivariatePolynomialOverField> &G, std::vector<UnivariatePolynomialOverField> &U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > &S, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting);
template <class UnivariatePolynomialOverField, class Field>
extern void _printIntegral(UnivariatePolynomialOverField A, UnivariatePolynomialOverField D, UnivariatePolynomialOverField P, std::vector<UnivariatePolynomialOverField> G, std::vector<Field> h, std::vector<UnivariatePolynomialOverField> H, std::vector<Field> k, std::vector<UnivariatePolynomialOverField> K1, std::vector<UnivariatePolynomialOverField> K2, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting, std::string outputFormatting);

#endif
