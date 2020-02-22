#ifndef _MULTIPRECISION_ROOTFINDING_H
#define _MULTIPRECISION_ROOTFINDING_H

template <class UnivariatePolynomialOverField, class Field>
extern std::vector< std::vector<ComplexRationalNumber> > _rootsMultiprecision(std::vector<UnivariatePolynomialOverField> &U, int prec);
template <class UnivariatePolynomialOverField, class Field>
extern std::vector< std::vector<ComplexRationalNumber> > _rootsDoublePrecision(std::vector<UnivariatePolynomialOverField> &U);

#endif
