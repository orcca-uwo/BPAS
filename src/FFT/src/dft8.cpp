/*
  Authored by Colin Costello to satisfy the requirements of CS4470Y
  The following is an 8-point loop-unrolled Cooley-Tukey six-step FFT.
*/
#include <bpas.h>

#include "../../../include/FFT/src/dft8.h"
#include "../../../include/FFT/src/dft_utils.h"

//! 8 Point FFT

/*!
 \param a coefficient vector
 \param omega  (omega at that level)
*/
template <class FiniteField>
FiniteField* DFT_8(FiniteField* a, FiniteField omega){
  FiniteField b[8]; // temp array
  // step 1 :  B = L_{2}^{8} t* A   note : t* stands for tensor product
  b[0]=a[0];
  b[1]=a[2];
  b[2]=a[4];
  b[3]=a[6];
  b[4]=a[1];
  b[5]=a[3];
  b[6]=a[5];
  b[7]=a[7];
  // step 2 :  A = I_{2} t* L_{2}^{4} B
  a[0]=b[0];
  a[1]=b[2];
  a[2]=b[1];
  a[3]=b[3];
  a[4]=b[4];
  a[5]=b[6];
  a[6]=b[5];
  a[7]=b[7];
  // step 3 :  B = I_{4} t* DFT_2 A
  b[0]=a[0]+a[1];
  b[1]=a[0]-a[1];

  b[2]=a[2]+a[3];
  b[3]=a[2]-a[3];

  b[4]=a[4]+a[5];
  b[5]=a[4]-a[5];

  b[6]=a[6]+a[7];
  b[7]=a[6]-a[7];
  // step 4 :  A = I_{2} t* T_{2}^{4} B
  a[0]=b[0];
  a[1]=b[1];
  a[2]=b[2];
  a[3]=b[3]*omega*omega;
  a[4]=b[4];
  a[5]=b[5];
  a[6]=b[6];
  a[7]=b[7]*omega*omega;
  // step 5 :  B = I_{2} t* L_{2}^{4} A
  b[0]=a[0];
  b[1]=a[2];
  b[2]=a[1];
  b[3]=a[3];
  b[4]=a[4];
  b[5]=a[6];
  b[6]=a[5];
  b[7]=a[7];
  // step 6 :  A = I_{4} t* DFT_2 B
  a[0]=b[0]+b[1];
  a[1]=b[0]-b[1];
  a[2]=b[2]+b[3];
  a[3]=b[2]-b[3];
  a[4]=b[4]+b[5];
  a[5]=b[4]-b[5];
  a[6]=b[6]+b[7];
  a[7]=b[6]-b[7];
  // step 7 :  B = I_{2} t* L_{2}^{4} A
  b[0]=a[0];
  b[1]=a[2];
  b[2]=a[1];
  b[3]=a[3];
  b[4]=a[4];
  b[5]=a[6];
  b[6]=a[5];
  b[7]=a[7];
  // step 8 :  A = T_{4}^{8} t* B
  a[0]=b[0];
  a[1]=b[1];
  a[2]=b[2];
  a[3]=b[3];
  a[4]=b[4];
  a[5]=b[5]*omega;
  a[6]=b[6]*omega*omega;
  a[7]=b[7]*omega*omega*omega;
  // step 9 : B = L_{4}^{8} t* A
  b[0]=a[0];
  b[1]=a[4];
  b[2]=a[1];
  b[3]=a[5];
  b[4]=a[2];
  b[5]=a[6];
  b[6]=a[3];
  b[7]=a[7];
  // step 10 : A = I_{4} t* DFT_2
  a[0]=b[0]+b[1];
  a[1]=b[0]-b[1];
  a[2]=b[2]+b[3];
  a[3]=b[2]-b[3];
  a[4]=b[4]+b[5];
  a[5]=b[4]-b[5];
  a[6]=b[6]+b[7];
  a[7]=b[6]-b[7];
  // step 11 : B = L_{2}^{8} t* A
  b[0]=a[0];
  b[1]=a[2];
  b[2]=a[4];
  b[3]=a[6];
  b[4]=a[1];
  b[5]=a[3];
  b[6]=a[5];
  b[7]=a[7];

  // put result in a
  for (int i=0;i<8;i++){
    a[i]=b[i];
  }
  // return a.
  return a;
}

template SmallPrimeField* DFT_8<SmallPrimeField>(SmallPrimeField* A,SmallPrimeField omega);
template BigPrimeField* DFT_8<BigPrimeField>(BigPrimeField* A,BigPrimeField omega);
template GeneralizedFermatPrimeField* DFT_8<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* A,GeneralizedFermatPrimeField omega);

GeneralizedFermatPrimeField* DFT_8(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas){

  DFT2(&a[0],&a[4]); // dft on permutated indexes
  DFT2(&a[2],&a[6]);
  DFT2(&a[1],&a[5]);
  DFT2(&a[3],&a[7]);

  a[6]=a[6].MulPowR(2);//omegas[2];   //twiddle
  a[7]=a[7].MulPowR(2);//*omegas[2];

  DFT2(&a[0],&a[2]); // dft on permutated indexes
  DFT2(&a[4],&a[6]);
  DFT2(&a[1],&a[3]);
  DFT2(&a[5],&a[7]);

  a[5]=a[5].MulPowR(1);//*omegas[1]); // twiddle
  a[3]=a[3].MulPowR(2);//*omegas[2]);
  a[7]=a[7].MulPowR(3);//*omegas[3]);

  DFT2(&a[0],&a[1]);   // dft on permutated indexes
  DFT2(&a[4],&a[5]);
  DFT2(&a[2],&a[3]);
  DFT2(&a[6],&a[7]);

  swap(&a[1],&a[4]); // final permutation
  swap(&a[3],&a[6]);

  return a; // return result
}

// long int* DFT_8(long int* a,long int* omegas, long int prime,long long int R,long long int pP){
//
//   DFT2(&a[0],&a[4],prime); // dft on permutated indexes
//   DFT2(&a[1],&a[5],prime);
//   DFT2(&a[2],&a[6],prime);
//   DFT2(&a[3],&a[7],prime);
//
//   a[6]=multi(a[6],omegas[2],prime,R,pP);   // twiddle
//   a[7]=multi(a[7],omegas[2],prime,R,pP);
//
//   DFT2(&a[0],&a[2],prime); // dft on permutated indexes
//   DFT2(&a[1],&a[3],prime);
//   DFT2(&a[4],&a[6],prime);
//   DFT2(&a[5],&a[7],prime);
//
//   a[3]=multi(a[3],omegas[2],prime,R,pP);   // twiddle
//   a[5]=multi(a[5],omegas[1],prime,R,pP);
//   a[7]=multi(a[7],omegas[3],prime,R,pP);
//
//   DFT2(&a[0],&a[1],prime);   // dft on permutated indexes
//   DFT2(&a[2],&a[3],prime);
//   DFT2(&a[4],&a[5],prime);
//   DFT2(&a[6],&a[7],prime);
//
//   swap(&a[1],&a[4]); // final permutation
//   swap(&a[3],&a[6]);
//
//   return a; // return result
// }

/*!
 \param a coefficient vector
 \param omegas  (precomputed powers of omega)
*/
template<class FiniteField>
FiniteField* DFT_8(FiniteField* a,FiniteField* omegas){
  DFT2(&a[0],&a[4]); // dft on permutated indexes
  DFT2(&a[2],&a[6]);
  DFT2(&a[1],&a[5]);
  DFT2(&a[3],&a[7]);

  a[6]=a[6]*omegas[2];   //twiddle
  a[7]=a[7]*omegas[2];

  DFT2(&a[0],&a[2]); // dft on permutated indexes
  DFT2(&a[4],&a[6]);
  DFT2(&a[1],&a[3]);
  DFT2(&a[5],&a[7]);

  a[5]=a[5]*omegas[1]; // twiddle
  a[3]=a[3]*omegas[2];
  a[7]=a[7]*omegas[3];

  DFT2(&a[0],&a[1]);   // dft on permutated indexes
  DFT2(&a[4],&a[5]);
  DFT2(&a[2],&a[3]);
  DFT2(&a[6],&a[7]);

  swap(&a[1],&a[4]); // final permutation
  swap(&a[3],&a[6]);

  return a; // return result
}
template SmallPrimeField* DFT_8<SmallPrimeField>(SmallPrimeField* A,SmallPrimeField* omegas);
template BigPrimeField* DFT_8<BigPrimeField>(BigPrimeField* A,BigPrimeField* omegas);
