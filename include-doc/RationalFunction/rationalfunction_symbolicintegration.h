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


