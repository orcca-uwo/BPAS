#ifndef __CUSTOM_DFT_H 
#define __CUSTOM_DFT_H 
#include "../../../include/FiniteFields/SmallPrimeField_Support.h"
#include "../../../include/FiniteFields/BigPrimeField_Support.h"
#include "../../../include/FiniteFields/SmallPrimeField.hpp"
#include "../../../include/FiniteFields/BigPrimeField.hpp"
#include "../../../include/FiniteFields/GeneralizedFermatPrimeField.hpp"

template<class FiniteField>
inline void swap(FiniteField* a,FiniteField* b){
	FiniteField tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
	}
template void swap<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* b);
template void swap<BigPrimeField>(BigPrimeField*a,BigPrimeField*b);
template void swap<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField*a,GeneralizedFermatPrimeField*b);


template<class FiniteField>
inline void DFT2 (FiniteField* a0,FiniteField* a1){
	FiniteField sum;
	sum=*a0+*a1;
	*a1=*a0-*a1;
	*a0=sum;
}
template void DFT2<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* b);
template void DFT2<BigPrimeField>(BigPrimeField* a,BigPrimeField* b);
template void DFT2<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* b);


template <class FiniteField>
extern  FiniteField* DFT_8(FiniteField* a,FiniteField* omegas);

inline GeneralizedFermatPrimeField* DFT_8(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas);

template <class FiniteField>
extern  FiniteField* DFT_16(FiniteField* a,FiniteField* omegas);

inline GeneralizedFermatPrimeField* DFT_16(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas);

template <class FiniteField>
extern  FiniteField* DFT_32(FiniteField* a,FiniteField* omegas);

inline GeneralizedFermatPrimeField* DFT_32(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas);

template <class FiniteField>
extern  FiniteField* DFT_64(FiniteField* a,FiniteField* omegas);

inline GeneralizedFermatPrimeField* DFT_64(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas);

template<class FiniteField>
void precomputed_DFT_general(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner);


template<class FiniteField>
void precomputed_DFT_general_new(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner);


inline void swap(long int* a,long int* b){
	long int tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}


inline void DFT2 (long int* a0,long int* a1,long int prime){
long int sum;
sum=smallprimefield_add(*a0,*a1,prime);
*a1=smallprimefield_sub(*a0,*a1,prime);
*a0=sum;
}


inline long int host_pow(long int x, long int e, long int prime, long long int R, long long int pP){
	long int m=1;
	for(long int i=0; i<e;i++)
		m=smallprimefield_multi(x, m, prime, R, pP);
	return m;
}

long int* DFT_8(long int* a,long int* omegas, long int prime,long long int R,long long int pP);

long int* DFT_16(long int* a,long int* omegas, long int prime,long long int R,long long int pP);

long int* DFT_32(long int* a,long int* omegas, long int prime,long long int R,long long int pP);

long int* DFT_64(long int* a,long int* omegas, long int prime,long long int R,long long int pP);

long int* precomputeInnerPowersOfOmega(long int omega,int basecase,long int prime,long long int R,long long int pP);

long int** precomputeOuterPowersOfOmega(int K, int e, long int omega, long int prime, long long int R);

void precomputed_DFT_general(long int* vector, int K, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner);

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


