/*
  Authored by Colin Costello to satisfy the requirements of CS4470Y
  The following is an 32-point loop-unrolled Cooley-Tukey six-step FFT.
*/
#include <bpas.h>
#include "../../../include/FFT/src/dft32.h"
#include "../../../include/FFT/src/dft_utils.h"

//! Split into two to avoid compiler notification

/*!
 \param a coefficient vector
 \param b coefficient vector temp
 \param omega_pow current power of omega
*/
template <class FiniteField>
void firstHalf(FiniteField* a,FiniteField* b,FiniteField* omega_pow){
    b[0]=a[0];
    b[1]=a[2];
    b[2]=a[4];
    b[3]=a[6];
    b[4]=a[8];
    b[5]=a[10];
    b[6]=a[12];
    b[7]=a[14];
    b[8]=a[16];
    b[9]=a[18];
    b[10]=a[20];
    b[11]=a[22];
    b[12]=a[24];
    b[13]=a[26];
    b[14]=a[28];
    b[15]=a[30];
    b[16]=a[1];
    b[17]=a[3];
    b[18]=a[5];
    b[19]=a[7];
    b[20]=a[9];
    b[21]=a[11];
    b[22]=a[13];
    b[23]=a[15];
    b[24]=a[17];
    b[25]=a[19];
    b[26]=a[21];
    b[27]=a[23];
    b[28]=a[25];
    b[29]=a[27];
    b[30]=a[29];
    b[31]=a[31];

    // Step 2, A = I_{2}  T* L_{2}^{16} B
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

    a[16]=b[16];
    a[17]=b[18];
    a[18]=b[20];
    a[19]=b[22];
    a[20]=b[24];
    a[21]=b[26];
    a[22]=b[28];
    a[23]=b[30];
    a[24]=b[17];
    a[25]=b[19];
    a[26]=b[21];
    a[27]=b[23];
    a[28]=b[25];
    a[29]=b[27];
    a[30]=b[29];
    a[31]=b[31];

    // Step 3, B = I_{4}  T* L_{2}^{8} A
    b[0]=a[0];
    b[1]=a[2];
    b[2]=a[4];
    b[3]=a[6];
    b[4]=a[1];
    b[5]=a[3];
    b[6]=a[5];
    b[7]=a[7];

    b[8]=a[8];
    b[9]=a[10];
    b[10]=a[12];
    b[11]=a[14];
    b[12]=a[9];
    b[13]=a[11];
    b[14]=a[13];
    b[15]=a[15];

    b[16]=a[16];
    b[17]=a[18];
    b[18]=a[20];
    b[19]=a[22];
    b[20]=a[17];
    b[21]=a[19];
    b[22]=a[21];
    b[23]=a[23];

    b[24]=a[24];
    b[25]=a[26];
    b[26]=a[28];
    b[27]=a[30];
    b[28]=a[25];
    b[29]=a[27];
    b[30]=a[29];
    b[31]=a[31];

    // Step 4, A = I_{8}  T* L_{2}^{4} B
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

    a[16]=b[16];
    a[17]=b[18];
    a[18]=b[17];
    a[19]=b[19];

    a[20]=b[20];
    a[21]=b[22];
    a[22]=b[21];
    a[23]=b[23];

    a[24]=b[24];
    a[25]=b[26];
    a[26]=b[25];
    a[27]=b[27];

    a[28]=b[28];
    a[29]=b[30];
    a[30]=b[29];
    a[31]=b[31];

    // Step 5, B = I_{16} T* DFT_2 A
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
    b[16]=a[16]+a[17];
    b[17]=a[16]-a[17];
    b[18]=a[18]+a[19];
    b[19]=a[18]-a[19];
    b[20]=a[20]+a[21];
    b[21]=a[20]-a[21];
    b[22]=a[22]+a[23];
    b[23]=a[22]-a[23];
    b[24]=a[24]+a[25];
    b[25]=a[24]-a[25];
    b[26]=a[26]+a[27];
    b[27]=a[26]-a[27];
    b[28]=a[28]+a[29];
    b[29]=a[28]-a[29];
    b[30]=a[30]+a[31];
    b[31]=a[30]-a[31];

    // Step 6, A = I_{8}  T* T_{2}^{4} B
    a[0]=b[0];
    a[1]=b[1];
    a[2]=b[2];
    a[3]=b[3]*omega_pow[8]; //omega*omega*omega*omega;

    a[4]=b[4];
    a[5]=b[5];
    a[6]=b[6];
    a[7]=b[7]*omega_pow[8];//*omega*omega*omega*omega;

    a[8]=b[8];
    a[9]=b[9];
    a[10]=b[10];
    a[11]=b[11]*omega_pow[8];//*omega*omega*omega*omega;

    a[12]=b[12];
    a[13]=b[13];
    a[14]=b[14];
    a[15]=b[15]*omega_pow[8];//*omega*omega*omega*omega;

    a[16]=b[16];
    a[17]=b[17];
    a[18]=b[18];
    a[19]=b[19]*omega_pow[8];//*omega*omega*omega*omega;

    a[20]=b[20];
    a[21]=b[21];
    a[22]=b[22];
    a[23]=b[23]*omega_pow[8];//*omega*omega*omega*omega;

    a[24]=b[24];
    a[25]=b[25];
    a[26]=b[26];
    a[27]=b[27]*omega_pow[8];//*omega*omega*omega*omega;

    a[28]=b[28];
    a[29]=b[29];
    a[30]=b[30];
    a[31]=b[31]*omega_pow[8];//*omega*omega*omega*omega;

    // Step 7, B = I_{8}  T* L_{2}^{4} A
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

    b[16]=a[16];
    b[17]=a[18];
    b[18]=a[17];
    b[19]=a[19];

    b[20]=a[20];
    b[21]=a[22];
    b[22]=a[21];
    b[23]=a[23];

    b[24]=a[24];
    b[25]=a[26];
    b[26]=a[25];
    b[27]=a[27];

    b[28]=a[28];
    b[29]=a[30];
    b[30]=a[29];
    b[31]=a[31];

    // Step 8, A = I_{16} T* DFT_2 B
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
    a[16]=b[16]+b[17];
    a[17]=b[16]-b[17];
    a[18]=b[18]+b[19];
    a[19]=b[18]-b[19];
    a[20]=b[20]+b[21];
    a[21]=b[20]-b[21];
    a[22]=b[22]+b[23];
    a[23]=b[22]-b[23];
    a[24]=b[24]+b[25];
    a[25]=b[24]-b[25];
    a[26]=b[26]+b[27];
    a[27]=b[26]-b[27];
    a[28]=b[28]+b[29];
    a[29]=b[28]-b[29];
    a[30]=b[30]+b[31];
    a[31]=b[30]-b[31];
}
//! Split into two to avoid compiler notification

/*!
 \param a coefficient vector
 \param b coefficient vector temp
 \param omega_pow current power of omega
*/
template <class FiniteField>
void secondHalf(FiniteField* a,FiniteField* b,FiniteField* omega_pow){
    // Step 9, B = I_{8}  T* L_{2}^{4} A
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

    b[16]=a[16];
    b[17]=a[18];
    b[18]=a[17];
    b[19]=a[19];

    b[20]=a[20];
    b[21]=a[22];
    b[22]=a[21];
    b[23]=a[23];

    b[24]=a[24];
    b[25]=a[26];
    b[26]=a[25];
    b[27]=a[27];

    b[28]=a[28];
    b[29]=a[30];
    b[30]=a[29];
    b[31]=a[31];

    // Step 10, A = I_{4} T* T_{4}^{8} B
    a[0]=b[0];
    a[1]=b[1];
    a[2]=b[2];
    a[3]=b[3];
    a[4]=b[4];
    a[5]=b[5]*omega_pow[4];//*omega*omega;
    a[6]=b[6]*omega_pow[8];//*omega*omega*omega*omega;
    a[7]=b[7]*omega_pow[12];//omega*omega*omega*omega*omega*omega;

    a[8]=b[8];
    a[9]=b[9];
    a[10]=b[10];
    a[11]=b[11];
    a[12]=b[12];
    a[13]=b[13]*omega_pow[4];//*omega*omega;
    a[14]=b[14]*omega_pow[8];//*omega*omega*omega*omega;
    a[15]=b[15]*omega_pow[12];//*omega*omega*omega*omega*omega*omega;

    a[16]=b[16];
    a[17]=b[17];
    a[18]=b[18];
    a[19]=b[19];
    a[20]=b[20];
    a[21]=b[21]*omega_pow[4];//*omega*omega;
    a[22]=b[22]*omega_pow[8];//*omega*omega*omega*omega;
    a[23]=b[23]*omega_pow[12];//*omega*omega*omega*omega*omega*omega;

    a[24]=b[24];
    a[25]=b[25];
    a[26]=b[26];
    a[27]=b[27];
    a[28]=b[28];
    a[29]=b[29]*omega_pow[4];//*omega*omega;
    a[30]=b[30]*omega_pow[8];//*omega*omega*omega*omega;
    a[31]=b[31]*omega_pow[12];//*omega*omega*omega*omega*omega*omega;

    // Step 11, B = I_{4} T* L_{4}^{8} A
    b[0]=a[0];
    b[1]=a[4];
    b[2]=a[1];
    b[3]=a[5];
    b[4]=a[2];
    b[5]=a[6];
    b[6]=a[3];
    b[7]=a[7];

    b[8]=a[8];
    b[9]=a[12];
    b[10]=a[9];
    b[11]=a[13];
    b[12]=a[10];
    b[13]=a[14];
    b[14]=a[11];
    b[15]=a[15];

    b[16]=a[16];
    b[17]=a[20];
    b[18]=a[17];
    b[19]=a[21];
    b[20]=a[18];
    b[21]=a[22];
    b[22]=a[19];
    b[23]=a[23];

    b[24]=a[24];
    b[25]=a[28];
    b[26]=a[25];
    b[27]=a[29];
    b[28]=a[26];
    b[29]=a[30];
    b[30]=a[27];
    b[31]=a[31];

    // Step 12, A = I_{16} T* DFT_2 B
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
    a[16]=b[16]+b[17];
    a[17]=b[16]-b[17];
    a[18]=b[18]+b[19];
    a[19]=b[18]-b[19];
    a[20]=b[20]+b[21];
    a[21]=b[20]-b[21];
    a[22]=b[22]+b[23];
    a[23]=b[22]-b[23];
    a[24]=b[24]+b[25];
    a[25]=b[24]-b[25];
    a[26]=b[26]+b[27];
    a[27]=b[26]-b[27];
    a[28]=b[28]+b[29];
    a[29]=b[28]-b[29];
    a[30]=b[30]+b[31];
    a[31]=b[30]-b[31];

    // Step 13, B = I_{4} T* L_{2}^{8} A
    b[0]=a[0];
    b[1]=a[2];
    b[2]=a[4];
    b[3]=a[6];
    b[4]=a[1];
    b[5]=a[3];
    b[6]=a[5];
    b[7]=a[7];

    b[8]=a[8];
    b[9]=a[10];
    b[10]=a[12];
    b[11]=a[14];
    b[12]=a[9];
    b[13]=a[11];
    b[14]=a[13];
    b[15]=a[15];

    b[16]=a[16];
    b[17]=a[18];
    b[18]=a[20];
    b[19]=a[22];
    b[20]=a[17];
    b[21]=a[19];
    b[22]=a[21];
    b[23]=a[23];

    b[24]=a[24];
    b[25]=a[26];
    b[26]=a[28];
    b[27]=a[30];
    b[28]=a[25];
    b[29]=a[27];
    b[30]=a[29];
    b[31]=a[31];

    // Step 14, A = I_{2} T* T_{8}^{16} B
    a[0]=b[0];
    a[1]=b[1];
    a[2]=b[2];
    a[3]=b[3];
    a[4]=b[4];
    a[5]=b[5];
    a[6]=b[6];
    a[7]=b[7];
    a[8]=b[8];
    a[9]=b[9]*omega_pow[2];//*omega*omega;
    a[10]=b[10]*omega_pow[4];//*omega*omega*omega*omega;
    a[11]=b[11]*omega_pow[6];//*omega*omega*omega*omega*omega*omega;
    a[12]=b[12]*omega_pow[8];//*omega*omega*omega*omega*omega*omega*omega*omega;
    a[13]=b[13]*omega_pow[10];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[14]=b[14]*omega_pow[12];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[15]=b[15]*omega_pow[14];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;

    a[16]=b[16];
    a[17]=b[17];
    a[18]=b[18];
    a[19]=b[19];
    a[20]=b[20];
    a[21]=b[21];
    a[22]=b[22];
    a[23]=b[23];
    a[24]=b[24];
    a[25]=b[25]*omega_pow[2];//*omega*omega;
    a[26]=b[26]*omega_pow[4];//*omega*omega*omega*omega;
    a[27]=b[27]*omega_pow[6];//*omega*omega*omega*omega*omega*omega;
    a[28]=b[28]*omega_pow[8];//*omega*omega*omega*omega*omega*omega*omega*omega;
    a[29]=b[29]*omega_pow[10];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[30]=b[30]*omega_pow[12];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[31]=b[31]*omega_pow[14];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;

    // Step 15, B = I_{2} T* L_{8}^{16} A
    b[0]=a[0];
    b[1]=a[8];
    b[2]=a[1];
    b[3]=a[9];
    b[4]=a[2];
    b[5]=a[10];
    b[6]=a[3];
    b[7]=a[11];
    b[8]=a[4];
    b[9]=a[12];
    b[10]=a[5];
    b[11]=a[13];
    b[12]=a[6];
    b[13]=a[14];
    b[14]=a[7];
    b[15]=a[15];

    b[16]=a[16];
    b[17]=a[24];
    b[18]=a[17];
    b[19]=a[25];
    b[20]=a[18];
    b[21]=a[26];
    b[22]=a[19];
    b[23]=a[27];
    b[24]=a[20];
    b[25]=a[28];
    b[26]=a[21];
    b[27]=a[29];
    b[28]=a[22];
    b[29]=a[30];
    b[30]=a[23];
    b[31]=a[31];

    // Step 16, A = I_{16} T* DFT_2 B
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
    a[16]=b[16]+b[17];
    a[17]=b[16]-b[17];
    a[18]=b[18]+b[19];
    a[19]=b[18]-b[19];
    a[20]=b[20]+b[21];
    a[21]=b[20]-b[21];
    a[22]=b[22]+b[23];
    a[23]=b[22]-b[23];
    a[24]=b[24]+b[25];
    a[25]=b[24]-b[25];
    a[26]=b[26]+b[27];
    a[27]=b[26]-b[27];
    a[28]=b[28]+b[29];
    a[29]=b[28]-b[29];
    a[30]=b[30]+b[31];
    a[31]=b[30]-b[31];

    // Step 17, B = I_{2} T* L_{2}^{16} A
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

    b[16]=a[16];
    b[17]=a[18];
    b[18]=a[20];
    b[19]=a[22];
    b[20]=a[24];
    b[21]=a[26];
    b[22]=a[28];
    b[23]=a[30];
    b[24]=a[17];
    b[25]=a[19];
    b[26]=a[21];
    b[27]=a[23];
    b[28]=a[25];
    b[29]=a[27];
    b[30]=a[29];
    b[31]=a[31];

    // Step 18, A = T_{16}^{32} T* B
    a[0]=b[0];
    a[1]=b[1];
    a[2]=b[2];
    a[3]=b[3];
    a[4]=b[4];
    a[5]=b[5];
    a[6]=b[6];
    a[7]=b[7];
    a[8]=b[8];
    a[9]=b[9];
    a[10]=b[10];
    a[11]=b[11];
    a[12]=b[12];
    a[13]=b[13];
    a[14]=b[14];
    a[15]=b[15];
    //--------------
    a[16]=b[16];
    a[17]=b[17]*omega_pow[1];//*omega
    a[18]=b[18]*omega_pow[2];//*omega*omega;
    a[19]=b[19]*omega_pow[3];//*omega*omega*omega;
    a[20]=b[20]*omega_pow[4];//*omega*omega*omega*omega;
    a[21]=b[21]*omega_pow[5];//*omega*omega*omega*omega*omega;
    a[22]=b[22]*omega_pow[6];//*omega*omega*omega*omega*omega*omega;
    a[23]=b[23]*omega_pow[7];//*omega*omega*omega*omega*omega*omega*omega;
    a[24]=b[24]*omega_pow[8];//*omega*omega*omega*omega*omega*omega*omega*omega;
    a[25]=b[25]*omega_pow[9];//*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[26]=b[26]*omega_pow[10];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[27]=b[27]*omega_pow[11];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[28]=b[28]*omega_pow[12];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[29]=b[29]*omega_pow[13];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[30]=b[30]*omega_pow[14];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    a[31]=b[31]*omega_pow[15];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;

    // Step 19, B = L_{16}^{32} T* A
    b[0]=a[0];
    b[1]=a[16];
    b[2]=a[1];
    b[3]=a[17];
    b[4]=a[2];
    b[5]=a[18];
    b[6]=a[3];
    b[7]=a[19];
    b[8]=a[4];
    b[9]=a[20];
    b[10]=a[5];
    b[11]=a[21];
    b[12]=a[6];
    b[13]=a[22];
    b[14]=a[7];
    b[15]=a[23];
    b[16]=a[8];
    b[17]=a[24];
    b[18]=a[9];
    b[19]=a[25];
    b[20]=a[10];
    b[21]=a[26];
    b[22]=a[11];
    b[23]=a[27];
    b[24]=a[12];
    b[25]=a[28];
    b[26]=a[13];
    b[27]=a[29];
    b[28]=a[14];
    b[29]=a[30];
    b[30]=a[15];
    b[31]=a[31];

    // Step 20, A = I_{16} T* DFT_2 B
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
    a[16]=b[16]+b[17];
    a[17]=b[16]-b[17];
    a[18]=b[18]+b[19];
    a[19]=b[18]-b[19];
    a[20]=b[20]+b[21];
    a[21]=b[20]-b[21];
    a[22]=b[22]+b[23];
    a[23]=b[22]-b[23];
    a[24]=b[24]+b[25];
    a[25]=b[24]-b[25];
    a[26]=b[26]+b[27];
    a[27]=b[26]-b[27];
    a[28]=b[28]+b[29];
    a[29]=b[28]-b[29];
    a[30]=b[30]+b[31];
    a[31]=b[30]-b[31];

    // Step 21, B = L_{2}^{32} T* A
    b[0]=a[0];
    b[1]=a[2];
    b[2]=a[4];
    b[3]=a[6];
    b[4]=a[8];
    b[5]=a[10];
    b[6]=a[12];
    b[7]=a[14];
    b[8]=a[16];
    b[9]=a[18];
    b[10]=a[20];
    b[11]=a[22];
    b[12]=a[24];
    b[13]=a[26];
    b[14]=a[28];
    b[15]=a[30];
    b[16]=a[1];
    b[17]=a[3];
    b[18]=a[5];
    b[19]=a[7];
    b[20]=a[9];
    b[21]=a[11];
    b[22]=a[13];
    b[23]=a[15];
    b[24]=a[17];
    b[25]=a[19];
    b[26]=a[21];
    b[27]=a[23];
    b[28]=a[25];
    b[29]=a[27];
    b[30]=a[29];
    b[31]=a[31];

    for (int i=0;i<32;i++){
      a[i]=b[i];
    }
}
//! 32-point FFT (Split into two to avoid compiler notification)

/*!
 \param a coefficient vector
 \param omega (omega at that level)
*/
template <class FiniteField>
FiniteField* DFT_32(FiniteField* a,FiniteField omega){
    // temp array
    FiniteField b[32];
    // precompute powers of omega
    FiniteField omega_pow[16];
    omega_pow[0]=1;
    //omeag_pow[1]=(omega)%prime_here;
    for(int j=1;j<16;j++)
      omega_pow[j]=omega_pow[j-1]*omega; // PrimeField should multiply mod prime
    // Step 1, B = L_{2}^{32} A
    firstHalf(a,b,omega_pow);
    secondHalf(a,b,omega_pow);
    return a;
}


template SmallPrimeField* DFT_32<SmallPrimeField>(SmallPrimeField* A,SmallPrimeField omega);
template BigPrimeField* DFT_32<BigPrimeField>(BigPrimeField* A,BigPrimeField omega);
template GeneralizedFermatPrimeField* DFT_32<GeneralizedFermatPrimeField>(GeneralizedFermatPrimeField* A,GeneralizedFermatPrimeField omega);

// GeneralizedFermatPrimeField* DFT_32(GeneralizedFermatPrimeField* A, GeneralizedFermatPrimeField* omegas){
//
//     DFT2(&A[0],&A[16]);
//     DFT2(&A[1],&A[17]);
//     DFT2(&A[2],&A[18]);
//     DFT2(&A[3],&A[19]);
//     DFT2(&A[4],&A[20]);
//     DFT2(&A[5],&A[21]);
//     DFT2(&A[6],&A[22]);
//     DFT2(&A[7],&A[23]);
//     DFT2(&A[8],&A[24]);
//     DFT2(&A[9],&A[25]);
//     DFT2(&A[10],&A[26]);
//     DFT2(&A[11],&A[27]);
//     DFT2(&A[12],&A[28]);
//     DFT2(&A[13],&A[29]);
//     DFT2(&A[14],&A[30]);
//     DFT2(&A[15],&A[31]);
//
//     A[24]=A[24].MulPowR(8);//*omegas[8];
//     A[25]=A[25].MulPowR(8);//*omegas[8];
//     A[26]=A[26].MulPowR(8);//*omegas[8];
//     A[27]=A[27].MulPowR(8);//*omegas[8];
//     A[28]=A[28].MulPowR(8);//*omegas[8];
//     A[29]=A[29].MulPowR(8);//*omegas[8];
//     A[30]=A[30].MulPowR(8);//*omegas[8];
//     A[31]=A[31].MulPowR(8);//*omegas[8];
//
//     DFT2(&A[0],&A[8]);
//     DFT2(&A[1],&A[9]);
//     DFT2(&A[2],&A[10]);
//     DFT2(&A[3],&A[11]);
//     DFT2(&A[4],&A[12]);
//     DFT2(&A[5],&A[13]);
//     DFT2(&A[6],&A[14]);
//     DFT2(&A[7],&A[15]);
//     DFT2(&A[16],&A[24]);
//     DFT2(&A[17],&A[25]);
//     DFT2(&A[18],&A[26]);
//     DFT2(&A[19],&A[27]);
//     DFT2(&A[20],&A[28]);
//     DFT2(&A[21],&A[29]);
//     DFT2(&A[22],&A[30]);
//     DFT2(&A[23],&A[31]);
//
//     A[12]=A[12].MulPowR(8);//*omegas[8];
//     A[13]=A[13].MulPowR(8);//*omegas[8];
//     A[14]=A[14].MulPowR(8);//*omegas[8];
//     A[15]=A[15].MulPowR(8);//*omegas[8];
//     A[20]=A[20].MulPowR(4);//*omegas[4];
//     A[21]=A[21].MulPowR(4);//*omegas[4];
//     A[22]=A[22].MulPowR(4);//*omegas[4];
//     A[23]=A[23].MulPowR(4);//*omegas[4];
//     A[28]=A[28].MulPowR(12);//*omegas[12];
//     A[29]=A[29].MulPowR(12);//*omegas[12];
//     A[30]=A[30].MulPowR(12);//*omegas[12];
//     A[31]=A[31].MulPowR(12);//*omegas[12];
//
//     DFT2(&A[0],&A[4]);
//     DFT2(&A[1],&A[5]);
//     DFT2(&A[2],&A[6]);
//     DFT2(&A[3],&A[7]);
//     DFT2(&A[8],&A[12]);
//     DFT2(&A[9],&A[13]);
//     DFT2(&A[10],&A[14]);
//     DFT2(&A[11],&A[15]);
//     DFT2(&A[16],&A[20]);
//     DFT2(&A[17],&A[21]);
//     DFT2(&A[18],&A[22]);
//     DFT2(&A[19],&A[23]);
//     DFT2(&A[24],&A[28]);
//     DFT2(&A[25],&A[29]);
//     DFT2(&A[26],&A[30]);
//     DFT2(&A[27],&A[31]);
//
//     A[6]=A[6].MulPowR(8);//*omegas[8];
//     A[7]=A[7].MulPowR(8);//*omegas[8];
//     A[10]=A[10].MulPowR(4);//*omegas[4];
//     A[11]=A[11].MulPowR(4);//*omegas[4];
//     A[14]=A[14].MulPowR(12);//*omegas[12];
//     A[15]=A[15].MulPowR(12);//*omegas[12];
//     A[18]=A[18].MulPowR(2);//*omegas[2];
//     A[19]=A[19].MulPowR(2);//*omegas[2];
//     A[22]=A[22].MulPowR(10);//*omegas[10];
//     A[23]=A[23].MulPowR(10);//*omegas[10];
//     A[26]=A[26].MulPowR(6);//*omegas[6];
//     A[27]=A[27].MulPowR(6);//*omegas[6];
//     A[30]=A[30].MulPowR(14);//*omegas[14];
//     A[31]=A[31].MulPowR(14);//*omegas[14];
//
//     DFT2(&A[0],&A[2]);
//     DFT2(&A[1],&A[3]);
//     DFT2(&A[4],&A[6]);
//     DFT2(&A[5],&A[7]);
//     DFT2(&A[8],&A[10]);
//     DFT2(&A[9],&A[11]);
//     DFT2(&A[12],&A[14]);
//     DFT2(&A[13],&A[15]);
//     DFT2(&A[16],&A[18]);
//     DFT2(&A[17],&A[19]);
//     DFT2(&A[20],&A[22]);
//     DFT2(&A[21],&A[23]);
//     DFT2(&A[24],&A[26]);
//     DFT2(&A[25],&A[27]);
//     DFT2(&A[28],&A[30]);
//     DFT2(&A[29],&A[31]);
//
//     A[3]=A[3].MulPowR(8);//*omegas[8];//*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[5]=A[5].MulPowR(4);//*omegas[4];//*omega*omega*omega*omega;
//     A[7]=A[7].MulPowR(12);//*omegas[12];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[9]=A[9].MulPowR(2);//*omegas[2];//*omega*omega;
//     A[11]=A[11].MulPowR(10);//*omegas[10];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[13]=A[13].MulPowR(6);//*omegas[6];//*omega*omega*omega*omega*omega*omega;
//     A[15]=A[15].MulPowR(14);//*omegas[14];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[17]=A[17].MulPowR(1);//*omegas[1];//*omega
//     A[19]=A[19].MulPowR(9);//*omegas[9];//*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[21]=A[21].MulPowR(5);//*omegas[5];//*omega*omega*omega*omega*omega;
//     A[23]=A[23].MulPowR(13);//*omegas[13];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[25]=A[25].MulPowR(3);//*omegas[3];//*omega*omega*omega;
//     A[27]=A[27].MulPowR(11)//;*omegas[11];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[29]=A[29].MulPowR(7);//*omegas[7];//*omega*omega*omega*omega*omega*omega*omega;
//     A[31]=A[31].MulPowR(15);//*omegas[15];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//
//     DFT2(&A[0],&A[1]);
//     DFT2(&A[2],&A[3]);
//     DFT2(&A[4],&A[5]);
//     DFT2(&A[6],&A[7]);
//     DFT2(&A[8],&A[9]);
//     DFT2(&A[10],&A[11]);
//     DFT2(&A[12],&A[13]);
//     DFT2(&A[14],&A[15]);
//     DFT2(&A[16],&A[17]);
//     DFT2(&A[18],&A[19]);
//     DFT2(&A[20],&A[21]);
//     DFT2(&A[22],&A[23]);
//     DFT2(&A[24],&A[25]);
//     DFT2(&A[26],&A[27]);
//     DFT2(&A[28],&A[29]);
//     DFT2(&A[30],&A[31]);
//
//     swap(&A[1],&A[16]);
//     swap(&A[2],&A[8]);
//     swap(&A[3],&A[24]);
//     swap(&A[5],&A[20]);
//     swap(&A[6],&A[12]);
//     swap(&A[7],&A[28]);
//     swap(&A[9],&A[18]);
//     swap(&A[11],&A[26]);
//     swap(&A[13],&A[22]);
//     swap(&A[15],&A[30]);
//     swap(&A[19],&A[25]);
//     swap(&A[23],&A[29]);
//
//     return A;
// }

// long int* DFT_32(long int* A, long int* omegas, long int prime,long long int R,long long int pP){
//   /*// precompute powers of omega
//   // long int omegas[16];
//   // omegas[0]=1;
//   // omeag_pow[1]=omega;
//   // for(int j=2;j<16;j++)
//   //   omegas[j]=multi(omegas[j-1],omega,prime,R,pP);
//   */
//     DFT2(&A[0],&A[16],prime);
//     DFT2(&A[1],&A[17],prime);
//     DFT2(&A[2],&A[18],prime);
//     DFT2(&A[3],&A[19],prime);
//     DFT2(&A[4],&A[20],prime);
//     DFT2(&A[5],&A[21],prime);
//     DFT2(&A[6],&A[22],prime);
//     DFT2(&A[7],&A[23],prime);
//     DFT2(&A[8],&A[24],prime);
//     DFT2(&A[9],&A[25],prime);
//     DFT2(&A[10],&A[26],prime);
//     DFT2(&A[11],&A[27],prime);
//     DFT2(&A[12],&A[28],prime);
//     DFT2(&A[13],&A[29],prime);
//     DFT2(&A[14],&A[30],prime);
//     DFT2(&A[15],&A[31],prime);
//
//     A[24]=multi(A[24],omegas[8],prime,R,pP);
//     A[25]=multi(A[25],omegas[8],prime,R,pP);
//     A[26]=multi(A[26],omegas[8],prime,R,pP);
//     A[27]=multi(A[27],omegas[8],prime,R,pP);
//     A[28]=multi(A[28],omegas[8],prime,R,pP);
//     A[29]=multi(A[29],omegas[8],prime,R,pP);
//     A[30]=multi(A[30],omegas[8],prime,R,pP);
//     A[31]=multi(A[31],omegas[8],prime,R,pP);
//
//     DFT2(&A[0],&A[8],prime);
//     DFT2(&A[1],&A[9],prime);
//     DFT2(&A[2],&A[10],prime);
//     DFT2(&A[3],&A[11],prime);
//     DFT2(&A[4],&A[12],prime);
//     DFT2(&A[5],&A[13],prime);
//     DFT2(&A[6],&A[14],prime);
//     DFT2(&A[7],&A[15],prime);
//     DFT2(&A[16],&A[24],prime);
//     DFT2(&A[17],&A[25],prime);
//     DFT2(&A[18],&A[26],prime);
//     DFT2(&A[19],&A[27],prime);
//     DFT2(&A[20],&A[28],prime);
//     DFT2(&A[21],&A[29],prime);
//     DFT2(&A[22],&A[30],prime);
//     DFT2(&A[23],&A[31],prime);
//
//     A[12]=multi(A[12],omegas[8],prime,R,pP);
//     A[13]=multi(A[13],omegas[8],prime,R,pP);
//     A[14]=multi(A[14],omegas[8],prime,R,pP);
//     A[15]=multi(A[15],omegas[8],prime,R,pP);
//     A[20]=multi(A[20],omegas[4],prime,R,pP);
//     A[21]=multi(A[21],omegas[4],prime,R,pP);
//     A[22]=multi(A[22],omegas[4],prime,R,pP);
//     A[23]=multi(A[23],omegas[4],prime,R,pP);
//     A[28]=multi(A[28],omegas[12],prime,R,pP);
//     A[29]=multi(A[29],omegas[12],prime,R,pP);
//     A[30]=multi(A[30],omegas[12],prime,R,pP);
//     A[31]=multi(A[31],omegas[12],prime,R,pP);
//
//     DFT2(&A[0],&A[4],prime);
//     DFT2(&A[1],&A[5],prime);
//     DFT2(&A[2],&A[6],prime);
//     DFT2(&A[3],&A[7],prime);
//     DFT2(&A[8],&A[12],prime);
//     DFT2(&A[9],&A[13],prime);
//     DFT2(&A[10],&A[14],prime);
//     DFT2(&A[11],&A[15],prime);
//     DFT2(&A[16],&A[20],prime);
//     DFT2(&A[17],&A[21],prime);
//     DFT2(&A[18],&A[22],prime);
//     DFT2(&A[19],&A[23],prime);
//     DFT2(&A[24],&A[28],prime);
//     DFT2(&A[25],&A[29],prime);
//     DFT2(&A[26],&A[30],prime);
//     DFT2(&A[27],&A[31],prime);
//
//     A[6]=multi(A[6],omegas[8],prime,R,pP);
//     A[7]=multi(A[7],omegas[8],prime,R,pP);
//     A[10]=multi(A[10],omegas[4],prime,R,pP);
//     A[11]=multi(A[11],omegas[4],prime,R,pP);
//     A[14]=multi(A[14],omegas[12],prime,R,pP);
//     A[15]=multi(A[15],omegas[12],prime,R,pP);
//     A[18]=multi(A[18],omegas[2],prime,R,pP);
//     A[19]=multi(A[19],omegas[2],prime,R,pP);
//     A[22]=multi(A[22],omegas[10],prime,R,pP);
//     A[23]=multi(A[23],omegas[10],prime,R,pP);
//     A[26]=multi(A[26],omegas[6],prime,R,pP);
//     A[27]=multi(A[27],omegas[6],prime,R,pP);
//     A[30]=multi(A[30],omegas[14],prime,R,pP);
//     A[31]=multi(A[31],omegas[14],prime,R,pP);
//
//     DFT2(&A[0],&A[2],prime);
//     DFT2(&A[1],&A[3],prime);
//     DFT2(&A[4],&A[6],prime);
//     DFT2(&A[5],&A[7],prime);
//     DFT2(&A[8],&A[10],prime);
//     DFT2(&A[9],&A[11],prime);
//     DFT2(&A[12],&A[14],prime);
//     DFT2(&A[13],&A[15],prime);
//     DFT2(&A[16],&A[18],prime);
//     DFT2(&A[17],&A[19],prime);
//     DFT2(&A[20],&A[22],prime);
//     DFT2(&A[21],&A[23],prime);
//     DFT2(&A[24],&A[26],prime);
//     DFT2(&A[25],&A[27],prime);
//     DFT2(&A[28],&A[30],prime);
//     DFT2(&A[29],&A[31],prime);
//
//     A[3]=multi(A[3],omegas[8],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[5]=multi(A[5],omegas[4],prime,R,pP);//*omega*omega*omega*omega;
//     A[7]=multi(A[7],omegas[12],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[9]=multi(A[9],omegas[2],prime,R,pP);//*omega*omega;
//     A[11]=multi(A[11],omegas[10],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[13]=multi(A[13],omegas[6],prime,R,pP);//*omega*omega*omega*omega*omega*omega;
//     A[15]=multi(A[15],omegas[14],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[17]=multi(A[17],omegas[1],prime,R,pP);//*omega
//     A[19]=multi(A[19],omegas[9],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[21]=multi(A[21],omegas[5],prime,R,pP);//*omega*omega*omega*omega*omega;
//     A[23]=multi(A[23],omegas[13],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[25]=multi(A[25],omegas[3],prime,R,pP);//*omega*omega*omega;
//     A[27]=multi(A[27],omegas[11],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//     A[29]=multi(A[29],omegas[7],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega;
//     A[31]=multi(A[31],omegas[15],prime,R,pP);//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
//
//     DFT2(&A[0],&A[1],prime);
//     DFT2(&A[2],&A[3],prime);
//     DFT2(&A[4],&A[5],prime);
//     DFT2(&A[6],&A[7],prime);
//     DFT2(&A[8],&A[9],prime);
//     DFT2(&A[10],&A[11],prime);
//     DFT2(&A[12],&A[13],prime);
//     DFT2(&A[14],&A[15],prime);
//     DFT2(&A[16],&A[17],prime);
//     DFT2(&A[18],&A[19],prime);
//     DFT2(&A[20],&A[21],prime);
//     DFT2(&A[22],&A[23],prime);
//     DFT2(&A[24],&A[25],prime);
//     DFT2(&A[26],&A[27],prime);
//     DFT2(&A[28],&A[29],prime);
//     DFT2(&A[30],&A[31],prime);
//
//     swap(&A[1],&A[16]);
//     swap(&A[2],&A[8]);
//     swap(&A[3],&A[24]);
//     swap(&A[5],&A[20]);
//     swap(&A[6],&A[12]);
//     swap(&A[7],&A[28]);
//     swap(&A[9],&A[18]);
//     swap(&A[11],&A[26]);
//     swap(&A[13],&A[22]);
//     swap(&A[15],&A[30]);
//     swap(&A[19],&A[25]);
//     swap(&A[23],&A[29]);
//
//     return A;
// }

template<class FiniteField>
FiniteField* DFT_32(FiniteField* A, FiniteField* omegas){

    DFT2(&A[0],&A[16]);
    DFT2(&A[1],&A[17]);
    DFT2(&A[2],&A[18]);
    DFT2(&A[3],&A[19]);
    DFT2(&A[4],&A[20]);
    DFT2(&A[5],&A[21]);
    DFT2(&A[6],&A[22]);
    DFT2(&A[7],&A[23]);
    DFT2(&A[8],&A[24]);
    DFT2(&A[9],&A[25]);
    DFT2(&A[10],&A[26]);
    DFT2(&A[11],&A[27]);
    DFT2(&A[12],&A[28]);
    DFT2(&A[13],&A[29]);
    DFT2(&A[14],&A[30]);
    DFT2(&A[15],&A[31]);

    A[24]=A[24]*omegas[8];
    A[25]=A[25]*omegas[8];
    A[26]=A[26]*omegas[8];
    A[27]=A[27]*omegas[8];
    A[28]=A[28]*omegas[8];
    A[29]=A[29]*omegas[8];
    A[30]=A[30]*omegas[8];
    A[31]=A[31]*omegas[8];

    DFT2(&A[0],&A[8]);
    DFT2(&A[1],&A[9]);
    DFT2(&A[2],&A[10]);
    DFT2(&A[3],&A[11]);
    DFT2(&A[4],&A[12]);
    DFT2(&A[5],&A[13]);
    DFT2(&A[6],&A[14]);
    DFT2(&A[7],&A[15]);
    DFT2(&A[16],&A[24]);
    DFT2(&A[17],&A[25]);
    DFT2(&A[18],&A[26]);
    DFT2(&A[19],&A[27]);
    DFT2(&A[20],&A[28]);
    DFT2(&A[21],&A[29]);
    DFT2(&A[22],&A[30]);
    DFT2(&A[23],&A[31]);

    A[12]=A[12]*omegas[8];
    A[13]=A[13]*omegas[8];
    A[14]=A[14]*omegas[8];
    A[15]=A[15]*omegas[8];
    A[20]=A[20]*omegas[4];
    A[21]=A[21]*omegas[4];
    A[22]=A[22]*omegas[4];
    A[23]=A[23]*omegas[4];
    A[28]=A[28]*omegas[12];
    A[29]=A[29]*omegas[12];
    A[30]=A[30]*omegas[12];
    A[31]=A[31]*omegas[12];

    DFT2(&A[0],&A[4]);
    DFT2(&A[1],&A[5]);
    DFT2(&A[2],&A[6]);
    DFT2(&A[3],&A[7]);
    DFT2(&A[8],&A[12]);
    DFT2(&A[9],&A[13]);
    DFT2(&A[10],&A[14]);
    DFT2(&A[11],&A[15]);
    DFT2(&A[16],&A[20]);
    DFT2(&A[17],&A[21]);
    DFT2(&A[18],&A[22]);
    DFT2(&A[19],&A[23]);
    DFT2(&A[24],&A[28]);
    DFT2(&A[25],&A[29]);
    DFT2(&A[26],&A[30]);
    DFT2(&A[27],&A[31]);

    A[6]=A[6]*omegas[8];
    A[7]=A[7]*omegas[8];
    A[10]=A[10]*omegas[4];
    A[11]=A[11]*omegas[4];
    A[14]=A[14]*omegas[12];
    A[15]=A[15]*omegas[12];
    A[18]=A[18]*omegas[2];
    A[19]=A[19]*omegas[2];
    A[22]=A[22]*omegas[10];
    A[23]=A[23]*omegas[10];
    A[26]=A[26]*omegas[6];
    A[27]=A[27]*omegas[6];
    A[30]=A[30]*omegas[14];
    A[31]=A[31]*omegas[14];

    DFT2(&A[0],&A[2]);
    DFT2(&A[1],&A[3]);
    DFT2(&A[4],&A[6]);
    DFT2(&A[5],&A[7]);
    DFT2(&A[8],&A[10]);
    DFT2(&A[9],&A[11]);
    DFT2(&A[12],&A[14]);
    DFT2(&A[13],&A[15]);
    DFT2(&A[16],&A[18]);
    DFT2(&A[17],&A[19]);
    DFT2(&A[20],&A[22]);
    DFT2(&A[21],&A[23]);
    DFT2(&A[24],&A[26]);
    DFT2(&A[25],&A[27]);
    DFT2(&A[28],&A[30]);
    DFT2(&A[29],&A[31]);

    A[3]=A[3]*omegas[8];//*omega*omega*omega*omega*omega*omega*omega*omega;
    A[5]=A[5]*omegas[4];//*omega*omega*omega*omega;
    A[7]=A[7]*omegas[12];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    A[9]=A[9]*omegas[2];//*omega*omega;
    A[11]=A[11]*omegas[10];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    A[13]=A[13]*omegas[6];//*omega*omega*omega*omega*omega*omega;
    A[15]=A[15]*omegas[14];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    A[17]=A[17]*omegas[1];//*omega
    A[19]=A[19]*omegas[9];//*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    A[21]=A[21]*omegas[5];//*omega*omega*omega*omega*omega;
    A[23]=A[23]*omegas[13];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    A[25]=A[25]*omegas[3];//*omega*omega*omega;
    A[27]=A[27]*omegas[11];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
    A[29]=A[29]*omegas[7];//*omega*omega*omega*omega*omega*omega*omega;
    A[31]=A[31]*omegas[15];//*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;

    DFT2(&A[0],&A[1]);
    DFT2(&A[2],&A[3]);
    DFT2(&A[4],&A[5]);
    DFT2(&A[6],&A[7]);
    DFT2(&A[8],&A[9]);
    DFT2(&A[10],&A[11]);
    DFT2(&A[12],&A[13]);
    DFT2(&A[14],&A[15]);
    DFT2(&A[16],&A[17]);
    DFT2(&A[18],&A[19]);
    DFT2(&A[20],&A[21]);
    DFT2(&A[22],&A[23]);
    DFT2(&A[24],&A[25]);
    DFT2(&A[26],&A[27]);
    DFT2(&A[28],&A[29]);
    DFT2(&A[30],&A[31]);

    swap(&A[1],&A[16]);
    swap(&A[2],&A[8]);
    swap(&A[3],&A[24]);
    swap(&A[5],&A[20]);
    swap(&A[6],&A[12]);
    swap(&A[7],&A[28]);
    swap(&A[9],&A[18]);
    swap(&A[11],&A[26]);
    swap(&A[13],&A[22]);
    swap(&A[15],&A[30]);
    swap(&A[19],&A[25]);
    swap(&A[23],&A[29]);

    return A;
}

template SmallPrimeField* DFT_32<SmallPrimeField>(SmallPrimeField* A,SmallPrimeField* omegas);
template BigPrimeField* DFT_32<BigPrimeField>(BigPrimeField* A,BigPrimeField* omegas);
