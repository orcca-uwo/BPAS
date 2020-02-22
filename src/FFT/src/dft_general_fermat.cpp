#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <bpas.h>
#include "../../../include/FFT/src/dft_fermat_inverse_16.h"
#include "../../../include/FFT/src/dft_fermat_fwd_16.h"
#include "../../../include/FFT/src/dft_general_fermat.h"

using namespace std;
//! Helper that prints a Generalized Fermat Prime Field vector

/*!
 \param A coefficient vector
 \param n size of the vector
*/
void printArray(GeneralizedFermatPrimeField* A, int n){
  std::cout << "Array : ";
  std::cout << "[";
  for (int i=0;i<n;i++){
    if (i==(n-1)) {
      std::cout << A[i];
      std::cout << "]\n";
    }else{
      std::cout << A[i];
      std::cout << ", ";
    }
  }
}

//! Stride Permutation Function L_{m}^{mn} written as L(m,n) (which is an n X m)
/*!
 \param a coefficient vector
 \param m
 \param n
*/
void stride_permutation(GeneralizedFermatPrimeField* A, int m,int n){
  GeneralizedFermatPrimeField* B = new GeneralizedFermatPrimeField[m*n];
  for (long int i=0;i<n;i++)
     for (int j=0;j<m;j++)
        B[j*n+i]=A[i*m+j];
  for (long int i=0;i<m*n;i++)
    A[i]=B[i];
}

//! Twiddle Factors Function represents T_{m}^{mn} (equivalent to D_{n},{m}) 

/*!
 \param a coefficient vector
 \param m 
 \param n
 \param omega  (omega at that level) 
*/
void twiddle(GeneralizedFermatPrimeField* vector, int m, int n,GeneralizedFermatPrimeField omega_w){
  for (int j=0;j<n;j++)
    for(int i=0;i<m;i++){
      //point-wise multiplication
      GeneralizedFermatPrimeField t(omega_w);
      t=(omega_w)^(i*j);
      vector[j*m+i]=vector[j*m+i]*(t);
    }
}

// dft_2
inline void DFT_2 (GeneralizedFermatPrimeField* vector,GeneralizedFermatPrimeField omega_w){
  GeneralizedFermatPrimeField a0,a1;
  a0=vector[0]+vector[1];
  a1=vector[0]-vector[1];

  vector[0]=a0;
  vector[1]=a1;
}
//! Generalized FFT Generalized Fermat Prime Field

/*!
 \param vector coefficient vector
 \param K basecase (8,16,32)
 \param e where K^e equals N (the size of coefficient vector)
 \param omega_w (Omega)^N=1
*/
void dft_general_fermat(GeneralizedFermatPrimeField* vector, int K, int e, GeneralizedFermatPrimeField omega_w){
  /*       Step I       */
  for (long int i=0;i<=e-2;i++){
    for (long int j=0;j<pow(K,i);j++){
       stride_permutation(&vector[j*(long int)pow(K,e-i)],K,(long int)pow(K,e-i-1));
    }
  }
  /*       Step II       */
  GeneralizedFermatPrimeField omega_at_that_level_b(omega_w);
  omega_at_that_level_b = (omega_w)^(pow(K,e-1));
  for (long int j=0;j<pow(K,e-1);j++){
    int idx = j*K;
    if (K==16){
      dft_fermat_fwd_16(&vector[idx]);
    }
  }
  for (int i=e-2;i>=0;i--){
    /*       Step III       */
    int stride = pow(K,e-i);
    GeneralizedFermatPrimeField omega_at_that_level(omega_w);
    omega_at_that_level = (omega_w)^(pow(K,i));
    for (int j=0;j<pow(K,i);j++){
      int index = j*stride;
      twiddle(&vector[index],stride/K,K,omega_at_that_level);
      stride_permutation(&vector[index],stride/K,K);
    }
    /*       Step IV       */
    for (int j=0;j<pow(K,e-1);j++){
      int idx = j*K;
      if (K==16){
        dft_fermat_fwd_16(&vector[idx]);
      }
    }
    /*       Step V       */
    for (int j=0;j<pow(K,i);j++){
      int index = j*stride;
      stride_permutation(&vector[index],K,stride/K);
    }
  }// end of step III IV & V for_loop
}// end DFT_general

//! inverse Generalized FFT Generalized Fermat Prime Field

/*!
 \param vector coefficient vector
 \param K basecase (8,16,32)
 \param e where K^e equals N (the size of coefficient vector)
 \param omega_w (Omega)^N=1
*/
void dft_general_fermat_inv(GeneralizedFermatPrimeField* vector, int K, int e, GeneralizedFermatPrimeField omega_w){
  /*       Step I       */
  for (long int i=0;i<=e-2;i++){
    for (long int j=0;j<pow(K,i);j++){
       stride_permutation(&vector[j*(long int)pow(K,e-i)],K,(long int)pow(K,e-i-1));
    }
  }
  /*       Step II       */
  GeneralizedFermatPrimeField omega_at_that_level_b(omega_w);
  omega_at_that_level_b = (omega_w)^(pow(K,e-1));
  for (long int j=0;j<pow(K,e-1);j++){
    int idx = j*K;
    if (K==16){
      dft_fermat_inverse_16(&vector[idx]);
    }
  }
  for (int i=e-2;i>=0;i--){
    /*       Step III       */
    int stride = pow(K,e-i);
    GeneralizedFermatPrimeField omega_at_that_level(omega_w);
    omega_at_that_level = (omega_w)^(pow(K,i));
    for (int j=0;j<pow(K,i);j++){
      int index = j*stride;
      twiddle(&vector[index],stride/K,K,omega_at_that_level);
      stride_permutation(&vector[index],stride/K,K);
    }
    /*       Step IV       */
    for (int j=0;j<pow(K,e-1);j++){
      int idx = j*K;
      if (K==16){
        dft_fermat_inverse_16(&vector[idx]);
      }
    }
    /*       Step V       */
    for (int j=0;j<pow(K,i);j++){
      int index = j*stride;
      stride_permutation(&vector[index],K,stride/K);
    }
  }// end of step III IV & V for_loop
}// end DFT_general

//! Inverse FFT Generalized Fermat Prime Field

/*!
 \param vector coefficient vector
 \param K basecase (8,16,32)
 \param e where K^e equals N (the size of coefficient vector)
 \param omega_in Omega^N=1
*/
void inverse_fermat_DFT(GeneralizedFermatPrimeField* vector, int K,int e,GeneralizedFermatPrimeField omega_in){
  GeneralizedFermatPrimeField omega_inv = omega_in.inverse();
  dft_general_fermat_inv(vector,K,e,omega_inv);
  long int n = pow(K,e);
  for (long int i=0;i<n;i++){
    vector[i]=(vector[i]/n);
  }
}
