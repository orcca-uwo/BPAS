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
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


