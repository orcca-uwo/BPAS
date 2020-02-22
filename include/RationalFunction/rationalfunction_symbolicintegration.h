#ifndef _RATIONALFUNCTION_SYMBOLICINTEGRATION_H
#define _RATIONALFUNCTION_SYMBOLICINTEGRATION_H

template <class UnivariatePolynomialOverField, class Field>
extern void _hermiteReduce(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, std::vector<UnivariatePolynomialOverField> *g, std::vector<UnivariatePolynomialOverField> *h);
template <class UnivariatePolynomialOverField, class Field>
extern void _prepRothsteinTragerResultant(SparseUnivariatePolynomial<UnivariatePolynomialOverField> *a, SparseUnivariatePolynomial<UnivariatePolynomialOverField> *d, SparseUnivariatePolynomial<UnivariatePolynomialOverField> *amtdd, UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, Symbol indeterminate);
template <class UnivariatePolynomialOverField, class Field>
extern void _computeLogArguments(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, std::vector<UnivariatePolynomialOverField> *U, SparseUnivariatePolynomial<UnivariatePolynomialOverField> &amtdd, SparseUnivariatePolynomial<UnivariatePolynomialOverField> &d, UnivariatePolynomialOverField &D, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > &R, std::vector<UnivariatePolynomialOverField> &Q, Symbol indeterminate, bool PROFILING);
template <class UnivariatePolynomialOverField, class Field>
extern void _simplifyLogArguments(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, std::vector< UnivariatePolynomialOverField > *U, Symbol indeterminate);
template <class UnivariatePolynomialOverField, class Field>
void _initializeRationalIntegration(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, bool PROFILING);
template <class UnivariatePolynomialOverField, class Field>
void _integrateRationalFunctionRationalPart(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, UnivariatePolynomialOverField *P, std::vector<UnivariatePolynomialOverField> *G, UnivariatePolynomialOverField *R, UnivariatePolynomialOverField *H, bool PROFILING);
template <class UnivariatePolynomialOverField, class Field>
extern void _integrateRationalFunctionLogPart(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, std::vector<UnivariatePolynomialOverField> *U, UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, bool PROFILING);
template <class UnivariatePolynomialOverField, class Field>
extern void _integrateRationalFunction(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, UnivariatePolynomialOverField *P, std::vector<UnivariatePolynomialOverField> *G, std::vector<UnivariatePolynomialOverField> *U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > *S, bool PROFILING);

#endif
