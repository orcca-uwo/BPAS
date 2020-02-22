#ifndef __DFT_UTILS_H
#define __DFT_UTILS_H

#include "../../bpas.h"

//#include "../../FiniteFields/SmallPrimeField_Support.h"
// long int covert_in(long int _a, long int prime, long long int R){
// 	while (_a < 0){
//           _a = _a + prime;
//       }
//     _a = ((_a%prime)*(R%prime))%prime;
//     return _a;
// }
//
// long int covert_out(long int a, long int prime, long long int R, long long int Pp){
//     // long long int y = a + prime*((a*Pp) & (R-1));
// 		long long int y = a + prime*((a*Pp)%R);
// 		y>>=32;
//     //long int z =  y >> 32 ;
// 		long int z=y;
//     if(z >= prime){
//       z = z - prime;
//     }
//     return z;
// }
//
// long long int getPp(long int prime, long long int R){
// 		long long int p1 = R-prime;
// 		for(int i = 0; i < R; i ++){
//       if ((i * p1)%R == 1){
//         return i;
//       }
//     }
// }
//
// long int add(long int a, long int b, long int prime){
// 	return (a + b)%prime;
// }
//
// long int sub( long int a, long int b, long int prime){
// 	if ((a - b)<0){
//         return (prime + (a - b));
//       }
//       else{
//         return a - b;
//       }
// }
//
// long int multi(long int a, long int b, long int prime, long long int R, long long int Pp){
// 	    long long int x = (b * a);
//       long long int w = x*Pp;
//       long long int y = x + prime*(w&(R-1));
//       long long int z =  y>>32 ;
//       if(z >= prime){
//         z = z - prime;
//       }
//       return z;
// }
//
// long int host_pow(long int x,long int e,long int prime,long long int R,long long int pP){
// 	long int m=1;
// 	for(long int i=0; i<e;i++)
// 		m=multi(m,x,prime,R,pP);
// 	return m;
// }
//
// inline void DFT2 (long int* a0,long int* a1,long int prime){
//   long int sum;
//   sum=add(*a0,*a1,prime);
//   *a1=sub(*a0,*a1,prime);
//   *a0=sum;
// }
//
// inline void DFT_2(long int* vector, long int omega_w,long int prime,long long int R,long long int pP){
//
//     long int a0,a1;
//     a0=add(vector[0],vector[1],prime);
//     a1=sub(vector[0],vector[1],prime);
//     vector[0]=a0;
//     vector[1]=a1;
// }

inline void swap(long int* a,long int* b){
  long int tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

template<class FiniteField>
inline void swap(FiniteField* a,FiniteField* b){
  FiniteField tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

template void swap<SmallPrimeField>(SmallPrimeField* a,SmallPrimeField* b);
template void swap<BigPrimeField>(BigPrimeField*a,BigPrimeField*b);

//template void swap<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField*a,GeneralizedFermatPrimeField*b);

// inline void swap(SmallPrimeField* a,SmallPrimeField* b){
//   SmallPrimeField tmp;
//   tmp = *a;
//   *a = *b;
//   *b = tmp;
// }
//
// inline void swap(BigPrimeField* a,BigPrimeField* b){
//   BigPrimeField tmp;
//   tmp = *a;
//   *a = *b;
//   *b = tmp;
// }
//

inline void swap(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* b){
   GeneralizedFermatPrimeField tmp;
   tmp = *a;
   *a = *b;
   *b = tmp;
}



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


// inline void DFT2 (SmallPrimeField* a0,SmallPrimeField* a1){
//   SmallPrimeField sum;
//   sum=*a0+*a1;
//   *a1=*a0-*a1);
//   *a0=sum;
// }
//
// inline void DFT2 (BigPrimeField* a0,BigPrimeField* a1){
//   BigPrimeField sum;
//   sum=*a0+*a1;
//   *a1=*a0-*a1);
//   *a0=sum;
// }
//

inline void DFT2 (GeneralizedFermatPrimeField* a0,GeneralizedFermatPrimeField* a1){
   GeneralizedFermatPrimeField sum;
   sum=*a0+*a1;
   *a1=*a0-*a1;
   *a0=sum;
}



template<class FiniteField>
inline void DFT_2(FiniteField* vector, FiniteField omega_w){
    FiniteField a0,a1;
    a0=vector[0]+vector[1];
    a1=vector[0]-vector[1];
    vector[0]=a0;
    vector[1]=a1;
}

template void DFT_2<SmallPrimeField>(SmallPrimeField* A,SmallPrimeField omega);
template void DFT_2<BigPrimeField>(BigPrimeField* A,BigPrimeField omega);
template void DFT_2<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* A,GeneralizedFermatPrimeField omega);



// precomputes the powers of omega for the provided basecase
// long int* precomputeInnerPowersOfOmega(long int omega,int basecase,long int prime,long long int R,long long int pP){
// 	long int* omegas;
//   int arr_size;
//   arr_size = basecase >> 1;
//   omegas = new long int[arr_size];
//   omegas[0] = 1;
//   omegas[1] = omega;
//   for (int i=2;i<arr_size;++i){
//     omegas[i]=multi(omegas[i-1],omega,prime,R,pP);
//   }
//   return omegas;
// }

// precomputes the powers of omega for the provided basecase
template<class FiniteField>
FiniteField* precomputeInnerPowersOfOmega(FiniteField omega,int basecase){
  FiniteField* omegas;
  int arr_size;
  arr_size = basecase >> 1;
  omegas = new FiniteField[arr_size];
  omegas[0] = 1;
  omegas[1] = omega;
  for (int i=2;i<arr_size;++i){
    omegas[i]=omegas[i-1]*omega;
  }
  return omegas;
}
template SmallPrimeField* precomputeInnerPowersOfOmega<SmallPrimeField>(SmallPrimeField omega,int basecase);
template BigPrimeField* precomputeInnerPowersOfOmega<BigPrimeField>(BigPrimeField omega,int basecase);
template GeneralizedFermatPrimeField* precomputeInnerPowersOfOmega<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField omega,int basecase);

// precompute the powers of omega for the run from N to K
// long int** precomputeOuterPowersOfOmega(int K, int e, long int omega, long int prime, long long int R){
// 	// setup and array allocation
//   long int** precomputed_omegas = new long int*[e];
//   long long int pP = getPp(prime,R);
// 	// set stride
// 	int N=pow(K,e);
// 	int stride=N;
// 	// for each level from e-2 -> 0 allocate space
// 	for (int l=0;l<=e-2;l++){
// 		precomputed_omegas[l]= new long int[stride];
// 		stride/=K; // update stride
// 		// Compute D_{K}{K^e-1}=D_{K}{stride}
// 		for(int i=0;i<K;i++){
// 			for(int j=0;j<stride;j++){
// 				precomputed_omegas[l][i*stride+j]=host_pow(omega,(long int)(i*j),prime,R,pP);
// 			}
// 		}
// 	}
// 	return precomputed_omegas;
// }

// precompute the powers of omega for the run from N to K
template<class FiniteField>
FiniteField** precomputeOuterPowersOfOmega(int K, int e, FiniteField omega){
	// setup and array allocation
	FiniteField** precomputed_omegas = new FiniteField*[e];
	// set stride
	int N=pow(K,e);
	int stride=N;
	// for each level from e-2 -> 0 allocate space
	for (int l=0;l<=e-2;l++){
		precomputed_omegas[l]= new FiniteField[stride];
		stride/=K; // update stride
		// Compute D_{K}{K^e-1}=D_{K}{stride}
		for(int i=0;i<K;i++){
			for(int j=0;j<stride;j++){
				precomputed_omegas[l][i*stride+j]=omega^(i*j);
			}
		}
	}
	return precomputed_omegas;
}
template SmallPrimeField** precomputeOuterPowersOfOmega<SmallPrimeField>(int K,int e,SmallPrimeField omega);
template BigPrimeField** precomputeOuterPowersOfOmega<BigPrimeField>(int K,int e,BigPrimeField omega);
template GeneralizedFermatPrimeField** precomputeOuterPowersOfOmega<GeneralizedFermatPrimeField>(int K,int e,GeneralizedFermatPrimeField omega);




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


