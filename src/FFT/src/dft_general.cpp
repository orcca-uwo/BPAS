/*
  Authored by Colin Costello to satisfy the requirements of CS4470Y
  The following is a generalized Cooley-Tukey six-step FFT.
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "../../../include/FFT/src/dft_general.h"
#include "../../../include/FFT/src/dft8.h"
#include "../../../include/FFT/src/dft16.h"
#include "../../../include/FFT/src/dft32.h"
#include "../../../include/FFT/src/dft64.h"
#include "../../../include/FFT/src/dft_utils.h"
#include <bpas.h>
#include <cilk/cilk.h>

//! Stride Permutation Function L_{m}^{mn} written as L(m,n) (which is an n X m)
// /*!
//  \param a coefficient vector
//  \param m
//  \param n
// */
// template <class FiniteField>
// void stride_permutation(FiniteField* A, int m,int n){
//   FiniteField* B = new FiniteField[m*n];
//   for (int i=0;i<n;i++)
//      for (int j=0;j<m;j++)
//         B[j*n+i]=A[i*m+j];
//   for (int i=0;i<m*n;i++)
//     A[i]=B[i];
// }

//! Twiddle Factors Function represents T_{m}^{mn} (equivalent to D_{n},{m})

/*!
 \param a coefficient vector
 \param m
 \param n
 \param omega  (omega at that level)
*/
template <class FiniteField>
void twiddle(FiniteField* vector, int m, int n,FiniteField omega_w){
  for (int j=0;j<n;j++)
    for(int i=0;i<m;i++){
      FiniteField t(omega_w);
      t=(omega_w)^(i*j);
      //point-wise multiplication
      vector[j*m+i]=vector[j*m+i]*(t);//*omega_w^(i*j);<-- compiler tried XOR instead of pow()
    }
}

//! 2 Point FFT

// /*!
//  \param a coefficient vector
//  \param omega  (omega at that level)
// */
// template <class FiniteField>
// inline void DFT_2 (FiniteField* vector,FiniteField omega_w){
//   FiniteField a0,a1;
//   a0=vector[0]+vector[1];
//   a1=vector[0]-vector[1];
//   vector[0]=a0;
//   vector[1]=a1;
// }
//! Generalized FFT

/*!
 \param a coefficient vector
 \param K basecase (8,16,32)
 \param e where K^e equals N (the size of coefficient vector)
 \param omega_w (Omega)^N=1
*/
template <class FiniteField>
void DFT_general(FiniteField* vector, int K, int e, FiniteField omega_w){
  /*       Step I       */
  for (int i=0;i<=e-2;i++){
    int loop_max_step1 = pow(K,i);
    cilk_for (int j=0;j<loop_max_step1;j++){
       stride_permutation(&vector[j*(int)pow(K,e-i)],K,(int)pow(K,e-i-1));
    }
  }
  /*       Step II       */
  int loop_max_step2 = pow(K,e-1);
  cilk_for (int l=0;l<loop_max_step2;l++){
    int idx = l*K;
    //omega_at_that_level = w^K^{e-1}
    FiniteField omega_at_that_level(omega_w);
    //omega_at_that_level = omega_at_that_level^(pow(K,e-1));
    omega_at_that_level = (omega_w)^(loop_max_step2);
    //cout << " j: " << j << " step II omega: " << omega_at_that_level << " pow: " << pow(K,e-1) << endl;
    if (K==2){
      DFT_2(&vector[idx],omega_at_that_level);
    }else if (K==8){
      DFT_8(&vector[idx],omega_at_that_level);
    }else if (K==16){
      DFT_16(&vector[idx],omega_at_that_level);
    }else if (K==32){
      DFT_32(&vector[idx],omega_at_that_level);
    }
  }

  for (int i=e-2;i>=0;i--){
    /*       Step III       */
    int stride = pow(K,e-i);
    int loop_max_step3 = pow(K,i);
    cilk_for (int m=0;m<loop_max_step3;m++){
      int index = m*stride;
      //omega_at_that_level =w^K^{i}
      FiniteField omega_at_that_level(omega_w);
      //omega_at_that_level = omega_at_that_level^(pow(K,i));
      omega_at_that_level = (omega_w)^(loop_max_step3);

      //cout << " i: " << i << " j: " << j << " step III omega: " << omega_at_that_level << " pow :" << pow(K,i) << endl;
      twiddle(&vector[index],stride/K,K,omega_at_that_level);
      stride_permutation(&vector[index],stride/K,K);
    }
    /*       Step IV       */
    int loop_max_step4 = pow(K,e-1);
    cilk_for (int o=0;o<loop_max_step4;o++){
      int idx = o*K;
      //omega_at_that_level = w^K^{e-1}
      FiniteField omega_at_that_level(omega_w);
      //omega_at_that_level = omega_at_that_level^(pow(K,e-1));
      omega_at_that_level = (omega_w)^(loop_max_step4);
      //cout << "i: " << i << " j: " << j << " step IIII omega: " << omega_at_that_level << " pow: " << pow(K,e-1) << endl;
      if (K==2){
        DFT_2(&vector[idx],omega_at_that_level);
      }else if (K==8){
        DFT_8(&vector[idx],omega_at_that_level);
      }else if (K==16){
        DFT_16(&vector[idx],omega_at_that_level);
      }else if (K==32){
        DFT_32(&vector[idx],omega_at_that_level);
      }
    }
    /*       Step V       */
    int loop_max_step5 = pow(K,i);
    cilk_for (int p=0;p<loop_max_step5;p++){
      int index = p*stride;
      stride_permutation(&vector[index],K,stride/K);
    }
  }// end of step III IV & V for_loop
}// end DFT_general

// ========================= serial_version ===============================
// void DFT_general(FiniteField* vector, int K, int e,FiniteField omega_w){
//   /*       Step I       */
//   for (int i=0;i<=e-2;i++){
//     for (int j=0;j<pow(K,i);j++){ //can cilk-for
//        stride_permutation(&vector[j*(int)pow(K,e-i)],K,(int)pow(K,e-i-1));
//     }
//   }
//   /*       Step II       */
//   for (int j=0;j<pow(K,e-1);j++){ //cilk-for
//     int idx = j*K;
//      //omega_at_that_level = w^K^{e-1}
//     FiniteField omega_at_that_level (omega_w);
//     omega_at_that_level = (omega_w)^(pow(K,e-1));
//     if (K==2){
//       DFT_2(&vector[idx],omega_at_that_level);
//     }else if (K==8){
//       DFT_8(&vector[idx],omega_at_that_level);
//     }else if (K==16){
//       DFT_16(&vector[idx],omega_at_that_level);
//     }else if (K==32){
//       DFT_32(&vector[idx],omega_at_that_level);
//     }
//   }
//
//   for (int i=e-2;i>=0;i--){
//     /*       Step III       */
//     int stride = pow(K,e-i);
//     for (int j=0;j<pow(K,i);j++){ //cilk for
//        int index = j*stride;
//        //omega_at_that_level =w^K^{i}
//        FiniteField omega_at_that_level(omega_w);
//        omega_at_that_level= (omega_w)^(pow(K,i));
//        twiddle(&vector[index],stride/K,K,omega_at_that_level);
//        stride_permutation(&vector[index],stride/K,K);
//     }
//     /*       Step IV       */
//     for (int j=0;j<pow(K,e-1);j++){ //cilk-for
//       int idx = j*K;
//       //omega_at_that_level = w^K^{e-1}
//       FiniteField omega_at_that_level(omega_w);
//       omega_at_that_level = (omega_w)^(pow(K,e-1));
//       if (K==2){
//         DFT_2(&vector[idx],omega_at_that_level);
//       }else if (K==8){
//         DFT_8(&vector[idx],omega_at_that_level);
//       }else if (K==16){
//         DFT_16(&vector[idx],omega_at_that_level);
//       }else if (K==32){
//         DFT_32(&vector[idx],omega_at_that_level);
//       }
//     }
//     /*       Step V       */
//     for (int j=0;j<pow(K,i);j++){ //cilk-for
//       int index = j*stride;
//       stride_permutation(&vector[index],K,stride/K);
//     }
//   }// end of step III IV & V for_loop
// }// end DFT_general
// ================ end serial_version ==============================

//! Inverse Generalized FFT
/*!
 \param a coefficient vector
 \param K basecase (8,16,32)
 \param e where K^e equals N (the size of coefficient vector)
 \param omega_w (Omega)^N=1
*/
template <class FiniteField>
void inverse_DFT(FiniteField* vector, int K,int e,FiniteField omega_in){
  FiniteField omega_inv = omega_in.inverse();
  DFT_general(vector,K,e,omega_inv);
  int n = pow(K,e);
  for (int i=0;i<n;i++){
    vector[i]=(vector[i]/n);
  }
}

//! Fourier Prime lookup function ouputs a fourier prime (max 32bit)
/*!
 \param e (where e is the exponent of the first operand of the prime factorization 2^{e}+x^{#}.... )
*/
long int FFTPrime(int e){
  long int prime_vector[96] = {957349889,940572673,938475521,925892609,919601153,913309697,907018241,883949569,862978049,850395137,833617921,818937857,
  802160641,800063489,770703361,745537537,718274561,655360001,605028353,581959681,531628033,493879297,468713473,447741953,409993217,399507457,361758721,
  359661569,347078657,330301441,311427073,305135617,290455553,246415361,221249537,204472321,185597953,158334977,147849217,141557761,120586241,101711873,
  70254593,28311553,962592769,950009857,924844033,899678209,824180737,799014913,786432001,740294657,715128833,710934529,648019969,639631361,635437057,
  597688321,576716801,463470593,459276289,387973121,383778817,274726913,270532609,257949697,249561089,211812353,199229441,186646529,169869313,136314881,
  132120577,111149057,81788929,69206017,943718401,935329793,918552577,683671553,666894337,415236097,230686721,163577857,155189249,138412033,113246209,
  104857601,897581057,880803841,645922817,595591169,377487361,754974721,167772161,469762049};
  int index;
  if (0<e<=20){
    index = rand() % 96;
  }else if (e==21){
    index = rand() % 52 + 44;
  }else if (e==22){
    index = rand() % 20 + 76;
  }else if (e==23){
    index = rand() % 3 + 93;
  }else if (e==24){
    index = rand() % 2 + 94;
  }else if (e==25){
    index = rand() % 1 + 95;
  }else if (e==26){
    index = 96;
  }else{
    cout << "exponent should be an integer in the range 1-26." << endl;
    return -1;
  }
  return prime_vector[index];
}

template void DFT_general<SmallPrimeField>(SmallPrimeField* A,int K,int e,SmallPrimeField omega);
template void DFT_general<BigPrimeField>(BigPrimeField* A,int K,int e,BigPrimeField omega);
template void inverse_DFT<SmallPrimeField>(SmallPrimeField* A,int K,int e,SmallPrimeField omega);
template void inverse_DFT<BigPrimeField>(BigPrimeField* A,int K,int e,BigPrimeField omega);
//template void DFT_general<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* A,int K,int e,GeneralizedFermatPrimeField omega);

#define BLOCKSIZE 16

template<class FiniteField>
void stride_permutation(FiniteField* A,int m, int n){
  int blocksize=m^((m^n)&(-(m>n)));
  blocksize=BLOCKSIZE^((BLOCKSIZE^blocksize)&(-(BLOCKSIZE>blocksize)));
  FiniteField* B = new FiniteField[m*n];
  for (int i = 0; i < n; i += blocksize) {
    for (int j = 0; j < m; j += blocksize) {
        // transpose the block beginning at [i,j]
        for (int k = i; k < i + blocksize; ++k) {
          for (int l = j; l < j + blocksize; ++l) {
              B[k+l*n] = A[l+k*m];
          }
        }
    }
  }
  for (long int i=0;i<m*n;i++)
    A[i]=B[i];

  delete[] B;
}

template void stride_permutation<SmallPrimeField>(SmallPrimeField* A,int m,int n);
template void stride_permutation<BigPrimeField>(BigPrimeField* A,int m,int n);
template void stride_permutation<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* A,int m,int n);

template<class FiniteField>
void precomputed_twiddle(FiniteField* vector,int m, int n, const FiniteField* omegas){
  FiniteField t;
  for (int j=0;j<n;j++)
    for(int i=0;i<m;i++){
			t = omegas[j*m+i];
      vector[j*m+i]=vector[j*m+i]*t;
    }
}

template void precomputed_twiddle<SmallPrimeField>(SmallPrimeField* vector,int m,int n,const SmallPrimeField* omegas);
template void precomputed_twiddle<BigPrimeField>(BigPrimeField* vector,int m,int n,const BigPrimeField* omegas);
template void precomputed_twiddle<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* vector,int m,int n,const GeneralizedFermatPrimeField* omegas);

template<class FiniteField>
void precomputed_DFT_general(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner){
  /*       Step I       */
  for (long int i=0;i<=e-2;i++){
    for (long int j=0;j<pow(K,i);j++){
       stride_permutation(&vector[j*(long int)pow(K,e-i)],K,(long int)pow(K,e-i-1));
    }
  }
  /*       Step II       */
  FiniteField omega_at_that_level_b;
  omega_at_that_level_b = powersofomega_inner[e-1];//POW(omega_w,(long int)(pow(K,e-1)),prime,R,pP);
  // omega_at_that_level_b = powersofomega_outer[e-1][0];
  FiniteField* omegas;
  omegas = precomputeInnerPowersOfOmega(omega_at_that_level_b,K);
  // omegas = powersofomega_outer[e-1];
  for (long int j=0;j<pow(K,e-1);j++){
    int idx = j*K;
    //omega_at_that_level = w^K^{e-1}
    if (K==2){
      DFT_2(&vector[idx],omega_at_that_level_b);
    }else if (K==8){
      DFT_8(&vector[idx],omegas);
    }else if (K==16){
      DFT_16(&vector[idx],omegas);
    }else if (K==32){
      DFT_32(&vector[idx],omegas);
    }else if (K==64){
      DFT_64(&vector[idx],omegas);
    }
  }

  for (int i=e-2;i>=0;i--){
    /*       Step III       */
    int stride = pow(K,e-i);
    for (int j=0;j<pow(K,i);j++){
      int index = j*stride;
      precomputed_twiddle(&vector[index],stride/K,K,powersofomega_outer[i]);
      stride_permutation(&vector[index],stride/K,K);
    }
    /*       Step IV       */
    for (int j=0;j<pow(K,e-1);j++){
      int idx = j*K;
      // omega_at_that_level = w^K^{e-1}
      if (K==2){
        DFT_2(&vector[idx],omega_at_that_level_b);
      }else if (K==8){
        DFT_8(&vector[idx],omegas);
      }else if (K==16){
        DFT_16(&vector[idx],omegas);
      }else if (K==32){
        DFT_32(&vector[idx],omegas);
      }else if (K==64){
        DFT_64(&vector[idx],omegas);
      }
    }
    /*       Step V       */
    for (int j=0;j<pow(K,i);j++){
      int index = j*stride;
      stride_permutation(&vector[index],K,stride/K);
    }
  }// end of step III IV & V for_loop
		delete[] omegas;
}// end DFT_general

template void precomputed_DFT_general<SmallPrimeField>(SmallPrimeField* vector, int K, int e, SmallPrimeField omega_w, SmallPrimeField** powersofomega_outer, SmallPrimeField* powersofomega_inner);
template void precomputed_DFT_general<BigPrimeField>(BigPrimeField* vector, int K, int e, BigPrimeField omega_w, BigPrimeField** powersofomega_outer, BigPrimeField* powersofomega_inner);
template void precomputed_DFT_general<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* vector, int K, int e, GeneralizedFermatPrimeField omega_w, GeneralizedFermatPrimeField** powersofomega_outer, GeneralizedFermatPrimeField* powersofomega_inner);

// // c functions
// void stride_permutation(long int* A,int m, int n){
//   int blocksize=m^((m^n)&(-(m>n)));
//   blocksize=BLOCKSIZE^((BLOCKSIZE^blocksize)&(-(BLOCKSIZE>blocksize)));
//   long int* B = new long int[m*n];
//   for (int i = 0; i < n; i += blocksize) {
//     for (int j = 0; j < m; j += blocksize) {
//         // transpose the block beginning at [i,j]
//         for (int k = i; k < i + blocksize; ++k) {
//           for (int l = j; l < j + blocksize; ++l) {
//               B[k+l*n] = A[l+k*m];
//           }
//         }
//     }
//   }
//   for (long int i=0;i<m*n;i++)
//     A[i]=B[i];
//
//   delete[] B;
// }
//
// void precomputed_twiddle(long int* vector,int m, int n, const long int* omegas,
// 		long int prime, long long int R, long long int pP){
//   long int t;
//   for (int j=0;j<n;j++)
//     for(int i=0;i<m;i++){
// 			t = omegas[j*m+i];
//       vector[j*m+i]=multi(vector[j*m+i],t,prime,R,pP);
//     }
// }
//
// void precomputed_DFT_general(long int* vector, int K, int e, long int omega_w, long int prime, long long int R, long int** powersofomega_outer, long int* powersofomega_inner){
//
//   long long int pP = getPp(prime,R);
//   /*       Step I       */
//   for (long int i=0;i<=e-2;i++){
//     for (long int j=0;j<pow(K,i);j++){
//        stride_permutation(&vector[j*(long int)pow(K,e-i)],K,(long int)pow(K,e-i-1));
//     }
//   }
//   /*       Step II       */
//   long int omega_at_that_level_b;
//   omega_at_that_level_b = powersofomega_inner[e-1];//POW(omega_w,(long int)(pow(K,e-1)),prime,R,pP);
//   long int* omegas;
//   omegas = precomputeInnerPowersOfOmega(omega_at_that_level_b,K,prime,R,pP);
//   for (long int j=0;j<pow(K,e-1);j++){
//     int idx = j*K;
//     //omega_at_that_level = w^K^{e-1}
//     if (K==2){
//       DFT_2(&vector[idx],omega_at_that_level_b,prime,R,pP);
//     }else if (K==8){
//       DFT_8(&vector[idx],omegas,prime,R,pP);
//     }else if (K==16){
//       DFT_16(&vector[idx],omegas,prime,R,pP);
//     }else if (K==32){
//       DFT_32(&vector[idx],omegas,prime,R,pP);
//     }else if (K==64){
//       DFT_64(&vector[idx],omegas,prime,R,pP);
//     }
//   }
//
//   for (int i=e-2;i>=0;i--){
//     /*       Step III       */
//     int stride = pow(K,e-i);
//     for (int j=0;j<pow(K,i);j++){
//       int index = j*stride;
//       precomputed_twiddle(&vector[index],stride/K,K,powersofomega_outer[i],prime,R,pP);
//       stride_permutation(&vector[index],stride/K,K);
//     }
//     /*       Step IV       */
//     for (int j=0;j<pow(K,e-1);j++){
//       int idx = j*K;
//       // omega_at_that_level = w^K^{e-1}
//       if (K==2){
//         DFT_2(&vector[idx],omega_at_that_level_b,prime,R,pP);
//       }else if (K==8){
//         DFT_8(&vector[idx],omegas,prime,R,pP);
//       }else if (K==16){
//         DFT_16(&vector[idx],omegas,prime,R,pP);
//       }else if (K==32){
//         DFT_32(&vector[idx],omegas,prime,R,pP);
//       }else if (K==64){
//         DFT_64(&vector[idx],omegas,prime,R,pP);
//       }
//     }
//     /*       Step V       */
//     for (int j=0;j<pow(K,i);j++){
//       int index = j*stride;
//       stride_permutation(&vector[index],K,stride/K);
//     }
//   }// end of step III IV & V for_loop
// 	delete[] omegas;
// }// end DFT_general
