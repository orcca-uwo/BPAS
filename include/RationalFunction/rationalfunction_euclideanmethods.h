#ifndef _RATIONALFUNCTION_METHODS_H
#define _RATIONALFUNCTION_METHODS_H


template <class UnivariatePolynomialOverField, class Field>
extern void _euclideanDivide(UnivariatePolynomialOverField &a, UnivariatePolynomialOverField &b, UnivariatePolynomialOverField *q, UnivariatePolynomialOverField *r);
template <class UnivariatePolynomialOverField, class Field>
extern void _quotient(UnivariatePolynomialOverField &a, UnivariatePolynomialOverField &b, UnivariatePolynomialOverField *q);
template <class UnivariatePolynomialOverField, class Field>
extern void _remainder(UnivariatePolynomialOverField &a, UnivariatePolynomialOverField &b, UnivariatePolynomialOverField *r);
template <class UnivariatePolynomialOverField, class Field>
extern void _halfExtendedEuclidean(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &B, UnivariatePolynomialOverField *s, UnivariatePolynomialOverField *g);
template <class UnivariatePolynomialOverField, class Field>
extern void _halfExtendedEuclidean(UnivariatePolynomialOverField &a, UnivariatePolynomialOverField &b, UnivariatePolynomialOverField &c, UnivariatePolynomialOverField *s);
template <class UnivariatePolynomialOverField, class Field>
extern void _extendedEuclidean(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &B, UnivariatePolynomialOverField *s, UnivariatePolynomialOverField *t, UnivariatePolynomialOverField *g);
template <class UnivariatePolynomialOverField, class Field>
void _extendedEuclidean(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &B, UnivariatePolynomialOverField &C, UnivariatePolynomialOverField *s, UnivariatePolynomialOverField *t);

#endif
