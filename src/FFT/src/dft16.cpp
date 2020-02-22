/*
  Authored by Colin Costello to satisfy the requirements of CS4470Y
  The following is an 16-point loop-unrolled Cooley-Tukey six-step FFT.
*/
#include <bpas.h>
#include "../../../include/FFT/src/dft16.h"
#include "../../../include/FFT/src/dft_utils.h"

///! 16-Point Loop unrolled FFT
/*!
\param a coefficient vector
\param omega (**omega at that level).
*/
template <class FiniteField>
FiniteField* DFT_16(FiniteField* a,FiniteField omega){
  // temp array
  FiniteField b[16];
  // precompute powers of omega
  FiniteField omega_pow[8];
  omega_pow[0]=1;
  for(int j=1;j<8;j++)
    omega_pow[j]=omega_pow[j-1]*omega;

  // Step 1 : B = L_{2}^{16} t* A
  b[0]=a[0];
  b[1]=a[2];
  b[2]=a[4];
  b[3]=a[6];
  b[4]=a[8];
  b[5]=a[10];
  b[6]=a[12];
  b[7]=a[14];
  b[8]=a[1];
  b[9]=a[3];
  b[10]=a[5];
  b[11]=a[7];
  b[12]=a[9];
  b[13]=a[11];
  b[14]=a[13];
  b[15]=a[15];

  // Step 2 : A = I_{2} t* L_{2}^{8} B
  a[0]=b[0];
  a[1]=b[2];
  a[2]=b[4];
  a[3]=b[6];
  a[4]=b[1];
  a[5]=b[3];
  a[6]=b[5];
  a[7]=b[7];

  a[8]=b[8];
  a[9]=b[10];
  a[10]=b[12];
  a[11]=b[14];
  a[12]=b[9];
  a[13]=b[11];
  a[14]=b[13];
  a[15]=b[15];

  // Step 3 : B = I_{4} t* L_{2}^{4} A
  b[0]=a[0];
  b[1]=a[2];
  b[2]=a[1];
  b[3]=a[3];

  b[4]=a[4];
  b[5]=a[6];
  b[6]=a[5];
  b[7]=a[7];

  b[8]=a[8];
  b[9]=a[10];
  b[10]=a[9];
  b[11]=a[11];

  b[12]=a[12];
  b[13]=a[14];
  b[14]=a[13];
  b[15]=a[15];

  // Step 4 : A = I_{8} t* DFT_{2} B
  a[0]=b[0]+b[1];
  a[1]=b[0]-b[1];
  a[2]=b[2]+b[3];
  a[3]=b[2]-b[3];
  a[4]=b[4]+b[5];
  a[5]=b[4]-b[5];
  a[6]=b[6]+b[7];
  a[7]=b[6]-b[7];
  a[8]=b[8]+b[9];
  a[9]=b[8]-b[9];
  a[10]=b[10]+b[11];
  a[11]=b[10]-b[11];
  a[12]=b[12]+b[13];
  a[13]=b[12]-b[13];
  a[14]=b[14]+b[15];
  a[15]=b[14]-b[15];
  // Step 5 : B = T_{2}^{8} t* A
  b[0]=a[0];
  b[1]=a[1];
  b[2]=a[2];
  b[3]=a[3]*omega_pow[4];//*omega*omega*omega*omega;
  b[4]=a[4];
  b[5]=a[5];
  b[6]=a[6];
  b[7]=a[7]*omega_pow[4];//*omega*omega*omega*omega;
  b[8]=a[8];
  b[9]=a[9];
  b[10]=a[10];
  b[11]=a[11]*omega_pow[4];//*omega*omega*omega*omega;
  b[12]=a[12];
  b[13]=a[13];
  b[14]=a[14];
  b[15]=a[15]*omega_pow[4];//*omega*omega*omega*omega;

  // Step 6 : A= I_{4} t* L_{2}^{4} B
  a[0]=b[0];
  a[1]=b[2];
  a[2]=b[1];
  a[3]=b[3];
  a[4]=b[4];
  a[5]=b[6];
  a[6]=b[5];
  a[7]=b[7];
  a[8]=b[8];
  a[9]=b[10];
  a[10]=b[9];
  a[11]=b[11];
  a[12]=b[12];
  a[13]=b[14];
  a[14]=b[13];
  a[15]=b[15];
  // Step 7 : B = I_{8} t* DTF_2 A
  b[0]=a[0]+a[1];
  b[1]=a[0]-a[1];
  b[2]=a[2]+a[3];
  b[3]=a[2]-a[3];
  b[4]=a[4]+a[5];
  b[5]=a[4]-a[5];
  b[6]=a[6]+a[7];
  b[7]=a[6]-a[7];
  b[8]=a[8]+a[9];
  b[9]=a[8]-a[9];
  b[10]=a[10]+a[11];
  b[11]=a[10]-a[11];
  b[12]=a[12]+a[13];
  b[13]=a[12]-a[13];
  b[14]=a[14]+a[15];
  b[15]=a[14]-a[15];
  // Step 8 : A = I_{4} t* L_{2}^{4} B
  a[0]=b[0];
  a[1]=b[2];
  a[2]=b[1];
  a[3]=b[3];
  a[4]=b[4];
  a[5]=b[6];
  a[6]=b[5];
  a[7]=b[7];
  a[8]=b[8];
  a[9]=b[10];
  a[10]=b[9];
  a[11]=b[11];
  a[12]=b[12];
  a[13]=b[14];
  a[14]=b[13];
  a[15]=b[15];
  // Step 9 : B = I_{2} t* T_{4}^{8} A
  b[0]=a[0];
  b[1]=a[1];
  b[2]=a[2];
  b[3]=a[3];
  b[4]=a[4];
  b[5]=a[5]*omega_pow[2];//*omega*omega;
  b[6]=a[6]*omega_pow[4];//*omega*omega*omega*omega;
  b[7]=a[7]*omega_pow[6];//*omega*omega*omega*omega*omega*omega;
  b[8]=a[8];
  b[9]=a[9];
  b[10]=a[10];
  b[11]=a[11];
  b[12]=a[12];
  b[13]=a[13]*omega_pow[2];//*omega*omega;
  b[14]=a[14]*omega_pow[4];//*omega*omega*omega*omega;
  b[15]=a[15]*omega_pow[6];//*omega*omega*omega*omega*omega*omega;
  // Step 10 : A = I_{2} t* L_{4}^{8} B
  a[0]=b[0];
  a[1]=b[4];
  a[2]=b[1];
  a[3]=b[5];
  a[4]=b[2];
  a[5]=b[6];
  a[6]=b[3];
  a[7]=b[7];
  a[8]=b[8];
  a[9]=b[12];
  a[10]=b[9];
  a[11]=b[13];
  a[12]=b[10];
  a[13]=b[14];
  a[14]=b[11];
  a[15]=b[15];
  // Step 11 : B = I_{8} t* DFT_2 A
  b[0]=a[0]+a[1];
  b[1]=a[0]-a[1];
  b[2]=a[2]+a[3];
  b[3]=a[2]-a[3];
  b[4]=a[4]+a[5];
  b[5]=a[4]-a[5];
  b[6]=a[6]+a[7];
  b[7]=a[6]-a[7];
  b[8]=a[8]+a[9];
  b[9]=a[8]-a[9];
  b[10]=a[10]+a[11];
  b[11]=a[10]-a[11];
  b[12]=a[12]+a[13];
  b[13]=a[12]-a[13];
  b[14]=a[14]+a[15];
  b[15]=a[14]-a[15];
  // Step 12 : A = I_{2} t* L_{2}^{8} B
  a[0]=b[0];
  a[1]=b[2];
  a[2]=b[4];
  a[3]=b[6];
  a[4]=b[1];
  a[5]=b[3];
  a[6]=b[5];
  a[7]=b[7];

  a[8]=b[8];
  a[9]=b[10];
  a[10]=b[12];
  a[11]=b[14];
  a[12]=b[9];
  a[13]=b[11];
  a[14]=b[13];
  a[15]=b[15];
  // Step 13 : B = T_{8}^{16} A
  b[0]=a[0];
  b[1]=a[1];
  b[2]=a[2];
  b[3]=a[3];
  b[4]=a[4];
  b[5]=a[5];
  b[6]=a[6];
  b[7]=a[7];
  b[8]=a[8];
  b[9]=a[9]*omega;
  b[10]=a[10]*omega_pow[2];//*omega*omega;
  b[11]=a[11]*omega_pow[3];//*omega*omega*omega;
  b[12]=a[12]*omega_pow[4];//*omega*omega*omega*omega;
  b[13]=a[13]*omega_pow[5];//*omega*omega*omega*omega*omega;
  b[14]=a[14]*omega_pow[6];//*omega*omega*omega*omega*omega*omega;
  b[15]=a[15]*omega_pow[7];//*omega*omega*omega*omega*omega*omega*omega;
  // Step 14 : A = L_{8}^{16} B
  a[0]=b[0];
  a[1]=b[8];
  a[2]=b[1];
  a[3]=b[9];
  a[4]=b[2];
  a[5]=b[10];
  a[6]=b[3];
  a[7]=b[11];
  a[8]=b[4];
  a[9]=b[12];
  a[10]=b[5];
  a[11]=b[13];
  a[12]=b[6];
  a[13]=b[14];
  a[14]=b[7];
  a[15]=b[15];
  // Step 15 : B = I_{8} t* DFT_2 A
  b[0]=a[0]+a[1];
  b[1]=a[0]-a[1];
  b[2]=a[2]+a[3];
  b[3]=a[2]-a[3];
  b[4]=a[4]+a[5];
  b[5]=a[4]-a[5];
  b[6]=a[6]+a[7];
  b[7]=a[6]-a[7];
  b[8]=a[8]+a[9];
  b[9]=a[8]-a[9];
  b[10]=a[10]+a[11];
  b[11]=a[10]-a[11];
  b[12]=a[12]+a[13];
  b[13]=a[12]-a[13];
  b[14]=a[14]+a[15];
  b[15]=a[14]-a[15];
  // Step 16 : A = L_{2}^{16} B
  a[0]=b[0];
  a[1]=b[2];
  a[2]=b[4];
  a[3]=b[6];
  a[4]=b[8];
  a[5]=b[10];
  a[6]=b[12];
  a[7]=b[14];
  a[8]=b[1];
  a[9]=b[3];
  a[10]=b[5];
  a[11]=b[7];
  a[12]=b[9];
  a[13]=b[11];
  a[14]=b[13];
  a[15]=b[15];

  return a;
}
template SmallPrimeField* DFT_16<SmallPrimeField>(SmallPrimeField* A,SmallPrimeField omega);
template BigPrimeField* DFT_16<BigPrimeField>(BigPrimeField* A,BigPrimeField omega);
template GeneralizedFermatPrimeField* DFT_16<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* A,GeneralizedFermatPrimeField omega);

// GeneralizedFermatPrimeField* DFT_16(GeneralizedFermatPrimeField* a, GeneralizedFermatPrimeField* omegas){
//
//   DFT2(&a[0],&a[8]);
//   DFT2(&a[1],&a[9]);
//   DFT2(&a[2],&a[10]);
//   DFT2(&a[3],&a[11]);
//   DFT2(&a[4],&a[12]);
//   DFT2(&a[5],&a[13]);
//   DFT2(&a[6],&a[14]);
//   DFT2(&a[7],&a[15]);
//
//   a[12]=a[12].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//   a[13]=a[13].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//   a[14]=a[14].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//   a[15]=a[15].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//
//   DFT2(&a[0],&a[4]);
//   DFT2(&a[1],&a[5]);
//   DFT2(&a[2],&a[6]);
//   DFT2(&a[3],&a[7]);
//   DFT2(&a[8],&a[12]);
//   DFT2(&a[9],&a[13]);
//   DFT2(&a[10],&a[14]);
//   DFT2(&a[11],&a[15]);
//
//   a[6]=a[6].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//   a[7]=a[7].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//   a[10]=a[10].MulPowR(2);//*omegas[2];//*omega*omega;
//   a[11]=a[11].MulPowR(2);//*omegas[2];//*omega*omega;
//   a[14]=a[14].MulPowR(6);//*omegas[6];//*omega*omega*omega*omega*omega*omega;
//   a[15]=a[15].MulPowR(6);//*omegas[6];//*omega*omega*omega*omega*omega*omega;
//
//   DFT2(&a[0],&a[2]);
//   DFT2(&a[1],&a[3]);
//   DFT2(&a[4],&a[6]);
//   DFT2(&a[5],&a[7]);
//   DFT2(&a[8],&a[10]);
//   DFT2(&a[9],&a[11]);
//   DFT2(&a[12],&a[14]);
//   DFT2(&a[13],&a[15]);
//
//   a[3]=a[3].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//   a[5]=a[5].MulPowR(2);//*omegas[2];//*omega*omega;
//   a[7]=a[7].MulPowR(6);//*omegas[6];//*omega*omega*omega*omega*omega*omega;
//   a[9]=a[9].MulPowR(1);//*omegas[1];
//   a[11]=a[11].MulPowR(5);//*omegas[5];//*omega*omega*omega*omega*omega;
//   a[13]=a[13].MulPowR(3);//*omegas[3];//*omega*omega*omega;
//   a[15]=a[15].MulPowR(7);//*omegas[7];//*omega*omega*omega*omega*omega*omega*omega;
//
//   DFT2(&a[0],&a[1]);
//   DFT2(&a[2],&a[3]);
//   DFT2(&a[4],&a[5]);
//   DFT2(&a[6],&a[7]);
//   DFT2(&a[8],&a[9]);
//   DFT2(&a[10],&a[11]);
//   DFT2(&a[12],&a[13]);
//   DFT2(&a[14],&a[15]);
//
//   swap(&a[1],&a[8]);
//   swap(&a[2],&a[4]);
//   swap(&a[3],&a[12]);
//   swap(&a[5],&a[10]);
//   swap(&a[7],&a[14]);
//   swap(&a[11],&a[13]);
//
//   return a;
// }

// long int* DFT_16(long int* a, long int* omega_pow,long int prime,long long int R,long long int pP){
//
//   DFT2(&a[0],&a[8],prime);
//   DFT2(&a[1],&a[9],prime);
//   DFT2(&a[2],&a[10],prime);
//   DFT2(&a[3],&a[11],prime);
//   DFT2(&a[4],&a[12],prime);
//   DFT2(&a[5],&a[13],prime);
//   DFT2(&a[6],&a[14],prime);
//   DFT2(&a[7],&a[15],prime);
//
//   a[12]=multi(a[12],omega_pow[4],prime,R,pP);//*omega*omega*omega*omega;
//   a[13]=multi(a[13],omega_pow[4],prime,R,pP);//*omega*omega*omega*omega;
//   a[14]=multi(a[14],omega_pow[4],prime,R,pP);//*omega*omega*omega*omega;
//   a[15]=multi(a[15],omega_pow[4],prime,R,pP);//*omega*omega*omega*omega;
//
//   DFT2(&a[0],&a[4],prime);
//   DFT2(&a[1],&a[5],prime);
//   DFT2(&a[2],&a[6],prime);
//   DFT2(&a[3],&a[7],prime);
//   DFT2(&a[8],&a[12],prime);
//   DFT2(&a[9],&a[13],prime);
//   DFT2(&a[10],&a[14],prime);
//   DFT2(&a[11],&a[15],prime);
//
//   a[6]=multi(a[6],omega_pow[4],prime,R,pP);//*omega*omega*omega*omega;
//   a[7]=multi(a[7],omega_pow[4],prime,R,pP);//*omega*omega*omega*omega;
//   a[10]=multi(a[10],omega_pow[2],prime,R,pP);//*omega*omega;
//   a[11]=multi(a[11],omega_pow[2],prime,R,pP);//*omega*omega;
//   a[14]=multi(a[14],omega_pow[6],prime,R,pP);//*omega*omega*omega*omega*omega*omega;
//   a[15]=multi(a[15],omega_pow[6],prime,R,pP);//*omega*omega*omega*omega*omega*omega;
//
//   DFT2(&a[0],&a[2],prime);
//   DFT2(&a[1],&a[3],prime);
//   DFT2(&a[4],&a[6],prime);
//   DFT2(&a[5],&a[7],prime);
//   DFT2(&a[8],&a[10],prime);
//   DFT2(&a[9],&a[11],prime);
//   DFT2(&a[12],&a[14],prime);
//   DFT2(&a[13],&a[15],prime);
//
//   a[3]=multi(a[3],omega_pow[4],prime,R,pP);//*omega*omega*omega*omega;
//   a[5]=multi(a[5],omega_pow[2],prime,R,pP);//*omega*omega;
//   a[7]=multi(a[7],omega_pow[6],prime,R,pP);//*omega*omega*omega*omega*omega*omega;
//   a[9]=multi(a[9],omega_pow[1],prime,R,pP);
//   a[13]=multi(a[13],omega_pow[3],prime,R,pP);//*omega*omega*omega;
//   a[11]=multi(a[11],omega_pow[5],prime,R,pP);//*omega*omega*omega*omega*omega;
//   a[15]=multi(a[15],omega_pow[7],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega;
//
//   DFT2(&a[0],&a[1],prime);
//   DFT2(&a[2],&a[3],prime);
//   DFT2(&a[4],&a[5],prime);
//   DFT2(&a[6],&a[7],prime);
//   DFT2(&a[8],&a[9],prime);
//   DFT2(&a[10],&a[11],prime);
//   DFT2(&a[12],&a[13],prime);
//   DFT2(&a[14],&a[15],prime);
//
//   swap(&a[1],&a[8]);
//   swap(&a[2],&a[4]);
//   swap(&a[3],&a[12]);
//   swap(&a[5],&a[10]);
//   swap(&a[7],&a[14]);
//   swap(&a[11],&a[13]);
//
//   return a;
// }

template<class FiniteField>
FiniteField* DFT_16(FiniteField* a, FiniteField* omegas){

  DFT2(&a[0],&a[8]);
  DFT2(&a[1],&a[9]);
  DFT2(&a[2],&a[10]);
  DFT2(&a[3],&a[11]);
  DFT2(&a[4],&a[12]);
  DFT2(&a[5],&a[13]);
  DFT2(&a[6],&a[14]);
  DFT2(&a[7],&a[15]);

  a[12]=a[12]*omegas[4];//*omega*omega*omega*omega;
  a[13]=a[13]*omegas[4];//*omega*omega*omega*omega;
  a[14]=a[14]*omegas[4];//*omega*omega*omega*omega;
  a[15]=a[15]*omegas[4];//*omega*omega*omega*omega;

  DFT2(&a[0],&a[4]);
  DFT2(&a[1],&a[5]);
  DFT2(&a[2],&a[6]);
  DFT2(&a[3],&a[7]);
  DFT2(&a[8],&a[12]);
  DFT2(&a[9],&a[13]);
  DFT2(&a[10],&a[14]);
  DFT2(&a[11],&a[15]);

  a[6]=a[6]*omegas[4];//*omega*omega*omega*omega;
  a[7]=a[7]*omegas[4];//*omega*omega*omega*omega;
  a[10]=a[10]*omegas[2];//*omega*omega;
  a[11]=a[11]*omegas[2];//*omega*omega;
  a[14]=a[14]*omegas[6];//*omega*omega*omega*omega*omega*omega;
  a[15]=a[15]*omegas[6];//*omega*omega*omega*omega*omega*omega;

  DFT2(&a[0],&a[2]);
  DFT2(&a[1],&a[3]);
  DFT2(&a[4],&a[6]);
  DFT2(&a[5],&a[7]);
  DFT2(&a[8],&a[10]);
  DFT2(&a[9],&a[11]);
  DFT2(&a[12],&a[14]);
  DFT2(&a[13],&a[15]);

  a[3]=a[3]*omegas[4];//*omega*omega*omega*omega;
  a[5]=a[5]*omegas[2];//*omega*omega;
  a[7]=a[7]*omegas[6];//*omega*omega*omega*omega*omega*omega;
  a[9]=a[9]*omegas[1];
  a[11]=a[11]*omegas[5];//*omega*omega*omega*omega*omega;
  a[13]=a[13]*omegas[3];//*omega*omega*omega;
  a[15]=a[15]*omegas[7];//*omega*omega*omega*omega*omega*omega*omega;

  DFT2(&a[0],&a[1]);
  DFT2(&a[2],&a[3]);
  DFT2(&a[4],&a[5]);
  DFT2(&a[6],&a[7]);
  DFT2(&a[8],&a[9]);
  DFT2(&a[10],&a[11]);
  DFT2(&a[12],&a[13]);
  DFT2(&a[14],&a[15]);

  swap(&a[1],&a[8]);
  swap(&a[2],&a[4]);
  swap(&a[3],&a[12]);
  swap(&a[5],&a[10]);
  swap(&a[7],&a[14]);
  swap(&a[11],&a[13]);

  return a;
}

template SmallPrimeField* DFT_16<SmallPrimeField>(SmallPrimeField* A,SmallPrimeField* omegas);
template BigPrimeField* DFT_16<BigPrimeField>(BigPrimeField* A,BigPrimeField* omegas);
