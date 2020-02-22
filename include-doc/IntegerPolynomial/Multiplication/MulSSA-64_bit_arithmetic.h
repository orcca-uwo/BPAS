#pragma once
/**
 	Implementation of polynomial multiplication using Schonhage - Strassen Algorithm.
 
 	@author Yuzhen Xie
 */
#include "gmp.h"
#include "math.h"

#include "Mul.h"
#include "../BivariatePoly.h"

#include "../../FFT/src/general_routine.h"
#include "../../FFT/src/basic_routine.h"
#if FURER
#include "../../../include/FFT/src/fft_furer1.h"
#else
#include "../../../include/FFT/src/fft_iter1.h"
#endif


class MulSSA: public Mul{
 public:
  int size;
  int* goodN;
  int* goodM;
  int* goodK;	

  MONTP_OPT2_AS_GENE *MontP1;
  MONTP_OPT2_AS_GENE *MontP2;
  MONTP_OPT2_AS_GENE *MontP;

  int* RevBidMap; //bit reversal for the base case iterative DFT

  static const int DFTBASESIZE = 1024; //size of base case iterative DFT

  static const sfixn P = 4179340454199820289;
//  static const sfixn P = 3799912185593857;
  static const sfixn HALF_P = P>>1; //P/2

  static const sfixn P1 = 4179340454199820289;
  static const sfixn P2 = 2485986994308513793;
  static const sfixn U1 = -740506764262110493;
  static const sfixn U2 = 1244909922527606046;  
  static const sfixn U2_R1_sft = 1422754197174406896;


  //static const sfixn P1 = 3799912185593857;
  //static const sfixn P2 = 3553621580972033;
  //static const sfixn U1 = -507660225853162;

  __int128 P1_P2, HALF_P1_P2, N_HALF_P1_P2;

  int N;
  int K;
  int M;
  unsigned long LIMB_BITS; //tight bound for u_i and v_i
  
  MulSSA() {
	P1_P2 = (__int128) P1 * P2;
	HALF_P1_P2 = P1_P2 >> 1;
	N_HALF_P1_P2 = -HALF_P1_P2;

	//init1();
	init2();
  }

  ~MulSSA() {
	delete [] goodN;
	delete [] goodM;
	delete [] goodK;

	//delete [] RevBidMap;

	//my_free(MontP1);
	//my_free(MontP2);
	//my_free(MontP);
  }
  /**
   * Multiplying two univariate polynomial
   *
   * @param a first polynomial
   * @param b second polynomial
   * @return result polynomial
   */
  virtual UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
	//int csize = a->getSize() + b->getSize() ;
	//int es2 = logceiling(csize);
	//int dims2 = 1<<es2;
	//if(dims2 != csize)
		UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
	//else
	//	UnivariateIntegerPolynomial c(RESULT_REAL_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
 
    mul2C(a, b, &c);
    
    return c;
  };
  
  /**
   * Multiplying two univariate polynomial
   *
   * @param a first polynomial
   * @param b second polynomial
   * @param c result polynomial
   */
  virtual void multiply(UnivariateIntegerPolynomial* a, UnivariateIntegerPolynomial* b, UnivariateIntegerPolynomial* c){
	mul2C(a, b, c);
  };
  
 private:

  
  void inline init2() 
  {//---------------------
    //std::cout<<"init2\n";
    //MontP1, MontP2
    MontP1 = (MONTP_OPT2_AS_GENE *)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
    MontP2 = (MONTP_OPT2_AS_GENE *)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
    MontP = (MONTP_OPT2_AS_GENE *)my_malloc(sizeof(MONTP_OPT2_AS_GENE)); 
   
    EX_MontP_Init_OPT2_AS_GENE(MontP1, P1);
    EX_MontP_Init_OPT2_AS_GENE(MontP2, P2);
    EX_MontP_Init_OPT2_AS_GENE(MontP, P); 

    //if (MontP1==NULL) 
    //  std::cout<<"MontP1: "<<MontP1<<std::endl;
    //else
    //  std::cout<<"MontP1 is not null: "<<MontP1->P<<std::endl;

    //std::cout<<"after init Mont\n";

    //RevBidMap
    RevBidMap = new int[DFTBASESIZE];
    for(int i=0; i<DFTBASESIZE; ++i)   
      RevBidMap[i] = i;
    PBPAS::RevBitInd(DFTBASESIZE, RevBidMap);
    //std::cout<<"after init RevBidMap\n";

    //N,K,M
	size = 432;
	goodN = new int[size];
	goodM = new int[size];
	goodK = new int[size];

	goodN[0] = 2; goodK[0] = 1; goodM[0] = 1;
	goodN[1] = 4; goodK[1] = 1; goodM[1] = 2;
	goodN[2] = 6; goodK[2] = 1; goodM[2] = 3;
	goodN[3] = 8; goodK[3] = 1; goodM[3] = 4;
	goodN[4] = 10; goodK[4] = 1; goodM[4] = 5;
	goodN[5] = 12; goodK[5] = 1; goodM[5] = 6;
	goodN[6] = 14; goodK[6] = 1; goodM[6] = 7;
	goodN[7] = 16; goodK[7] = 1; goodM[7] = 8;
	goodN[8] = 18; goodK[8] = 1; goodM[8] = 9;
	goodN[9] = 20; goodK[9] = 1; goodM[9] = 10;
	goodN[10] = 22; goodK[10] = 1; goodM[10] = 11;
	goodN[11] = 24; goodK[11] = 1; goodM[11] = 12;
	goodN[12] = 26; goodK[12] = 1; goodM[12] = 13;
	goodN[13] = 28; goodK[13] = 1; goodM[13] = 14;
	goodN[14] = 30; goodK[14] = 1; goodM[14] = 15;
	goodN[15] = 32; goodK[15] = 1; goodM[15] = 16;
	goodN[16] = 34; goodK[16] = 1; goodM[16] = 17;
	goodN[17] = 36; goodK[17] = 1; goodM[17] = 18;
	goodN[18] = 38; goodK[18] = 1; goodM[18] = 19;
	goodN[19] = 40; goodK[19] = 1; goodM[19] = 20;
	goodN[20] = 42; goodK[20] = 1; goodM[20] = 21;
	goodN[21] = 44; goodK[21] = 1; goodM[21] = 22;
	goodN[22] = 46; goodK[22] = 1; goodM[22] = 23;
	goodN[23] = 48; goodK[23] = 1; goodM[23] = 24;
	goodN[24] = 50; goodK[24] = 1; goodM[24] = 25;
	goodN[25] = 52; goodK[25] = 1; goodM[25] = 26;
	goodN[26] = 54; goodK[26] = 1; goodM[26] = 27;
	goodN[27] = 56; goodK[27] = 1; goodM[27] = 28;
	goodN[28] = 58; goodK[28] = 1; goodM[28] = 29;
	goodN[29] = 60; goodK[29] = 1; goodM[29] = 30;
	goodN[30] = 62; goodK[30] = 1; goodM[30] = 31;
	goodN[31] = 64; goodK[31] = 1; goodM[31] = 32;
	goodN[32] = 66; goodK[32] = 1; goodM[32] = 33;
	goodN[33] = 68; goodK[33] = 1; goodM[33] = 34;
	goodN[34] = 70; goodK[34] = 1; goodM[34] = 35;
	goodN[35] = 72; goodK[35] = 1; goodM[35] = 36;
	goodN[36] = 74; goodK[36] = 1; goodM[36] = 37;
	goodN[37] = 76; goodK[37] = 1; goodM[37] = 38;
	goodN[38] = 78; goodK[38] = 1; goodM[38] = 39;
	goodN[39] = 80; goodK[39] = 1; goodM[39] = 40;
	goodN[40] = 82; goodK[40] = 1; goodM[40] = 41;
	goodN[41] = 84; goodK[41] = 1; goodM[41] = 42;
	goodN[42] = 86; goodK[42] = 1; goodM[42] = 43;
	goodN[43] = 88; goodK[43] = 1; goodM[43] = 44;
	goodN[44] = 90; goodK[44] = 1; goodM[44] = 45;
	goodN[45] = 92; goodK[45] = 1; goodM[45] = 46;
	goodN[46] = 94; goodK[46] = 1; goodM[46] = 47;
	goodN[47] = 96; goodK[47] = 1; goodM[47] = 48;
	goodN[48] = 98; goodK[48] = 1; goodM[48] = 49;
	goodN[49] = 100; goodK[49] = 1; goodM[49] = 50;
	goodN[50] = 102; goodK[50] = 1; goodM[50] = 51;
	goodN[51] = 104; goodK[51] = 1; goodM[51] = 52;
	goodN[52] = 106; goodK[52] = 1; goodM[52] = 53;
	goodN[53] = 108; goodK[53] = 1; goodM[53] = 54;
	goodN[54] = 110; goodK[54] = 1; goodM[54] = 55;
	goodN[55] = 112; goodK[55] = 1; goodM[55] = 56;
	goodN[56] = 114; goodK[56] = 1; goodM[56] = 57;
	goodN[57] = 116; goodK[57] = 2; goodM[57] = 29;
	goodN[58] = 120; goodK[58] = 2; goodM[58] = 30;
	goodN[59] = 124; goodK[59] = 2; goodM[59] = 31;
	goodN[60] = 128; goodK[60] = 2; goodM[60] = 32;
	goodN[61] = 132; goodK[61] = 2; goodM[61] = 33;
	goodN[62] = 136; goodK[62] = 2; goodM[62] = 34;
	goodN[63] = 140; goodK[63] = 2; goodM[63] = 35;
	goodN[64] = 144; goodK[64] = 2; goodM[64] = 36;
	goodN[65] = 148; goodK[65] = 2; goodM[65] = 37;
	goodN[66] = 152; goodK[66] = 2; goodM[66] = 38;
	goodN[67] = 156; goodK[67] = 2; goodM[67] = 39;
	goodN[68] = 160; goodK[68] = 2; goodM[68] = 40;
	goodN[69] = 164; goodK[69] = 2; goodM[69] = 41;
	goodN[70] = 168; goodK[70] = 2; goodM[70] = 42;
	goodN[71] = 172; goodK[71] = 2; goodM[71] = 43;
	goodN[72] = 176; goodK[72] = 2; goodM[72] = 44;
	goodN[73] = 180; goodK[73] = 2; goodM[73] = 45;
	goodN[74] = 184; goodK[74] = 2; goodM[74] = 46;
	goodN[75] = 188; goodK[75] = 2; goodM[75] = 47;
	goodN[76] = 192; goodK[76] = 2; goodM[76] = 48;
	goodN[77] = 196; goodK[77] = 2; goodM[77] = 49;
	goodN[78] = 200; goodK[78] = 2; goodM[78] = 50;
	goodN[79] = 204; goodK[79] = 2; goodM[79] = 51;
	goodN[80] = 208; goodK[80] = 2; goodM[80] = 52;
	goodN[81] = 212; goodK[81] = 2; goodM[81] = 53;
	goodN[82] = 216; goodK[82] = 2; goodM[82] = 54;
	goodN[83] = 220; goodK[83] = 2; goodM[83] = 55;
	goodN[84] = 224; goodK[84] = 2; goodM[84] = 56;
	goodN[85] = 232; goodK[85] = 3; goodM[85] = 29;
	goodN[86] = 240; goodK[86] = 3; goodM[86] = 30;
	goodN[87] = 248; goodK[87] = 3; goodM[87] = 31;
	goodN[88] = 256; goodK[88] = 3; goodM[88] = 32;
	goodN[89] = 264; goodK[89] = 3; goodM[89] = 33;
	goodN[90] = 272; goodK[90] = 3; goodM[90] = 34;
	goodN[91] = 280; goodK[91] = 3; goodM[91] = 35;
	goodN[92] = 288; goodK[92] = 3; goodM[92] = 36;
	goodN[93] = 296; goodK[93] = 3; goodM[93] = 37;
	goodN[94] = 304; goodK[94] = 3; goodM[94] = 38;
	goodN[95] = 312; goodK[95] = 3; goodM[95] = 39;
	goodN[96] = 320; goodK[96] = 3; goodM[96] = 40;
	goodN[97] = 328; goodK[97] = 3; goodM[97] = 41;
	goodN[98] = 336; goodK[98] = 3; goodM[98] = 42;
	goodN[99] = 344; goodK[99] = 3; goodM[99] = 43;
	goodN[100] = 352; goodK[100] = 3; goodM[100] = 44;
	goodN[101] = 360; goodK[101] = 3; goodM[101] = 45;
	goodN[102] = 368; goodK[102] = 3; goodM[102] = 46;
	goodN[103] = 376; goodK[103] = 3; goodM[103] = 47;
	goodN[104] = 384; goodK[104] = 3; goodM[104] = 48;
	goodN[105] = 392; goodK[105] = 3; goodM[105] = 49;
	goodN[106] = 400; goodK[106] = 3; goodM[106] = 50;
	goodN[107] = 408; goodK[107] = 3; goodM[107] = 51;
	goodN[108] = 416; goodK[108] = 3; goodM[108] = 52;
	goodN[109] = 424; goodK[109] = 3; goodM[109] = 53;
	goodN[110] = 432; goodK[110] = 3; goodM[110] = 54;
	goodN[111] = 440; goodK[111] = 3; goodM[111] = 55;
	goodN[112] = 448; goodK[112] = 4; goodM[112] = 28;
	goodN[113] = 464; goodK[113] = 4; goodM[113] = 29;
	goodN[114] = 480; goodK[114] = 4; goodM[114] = 30;
	goodN[115] = 496; goodK[115] = 4; goodM[115] = 31;
	goodN[116] = 512; goodK[116] = 4; goodM[116] = 32;
	goodN[117] = 528; goodK[117] = 4; goodM[117] = 33;
	goodN[118] = 544; goodK[118] = 4; goodM[118] = 34;
	goodN[119] = 560; goodK[119] = 4; goodM[119] = 35;
	goodN[120] = 576; goodK[120] = 4; goodM[120] = 36;
	goodN[121] = 592; goodK[121] = 4; goodM[121] = 37;
	goodN[122] = 608; goodK[122] = 4; goodM[122] = 38;
	goodN[123] = 624; goodK[123] = 4; goodM[123] = 39;
	goodN[124] = 640; goodK[124] = 4; goodM[124] = 40;
	goodN[125] = 656; goodK[125] = 4; goodM[125] = 41;
	goodN[126] = 672; goodK[126] = 4; goodM[126] = 42;
	goodN[127] = 688; goodK[127] = 4; goodM[127] = 43;
	goodN[128] = 704; goodK[128] = 4; goodM[128] = 44;
	goodN[129] = 720; goodK[129] = 4; goodM[129] = 45;
	goodN[130] = 736; goodK[130] = 4; goodM[130] = 46;
	goodN[131] = 752; goodK[131] = 4; goodM[131] = 47;
	goodN[132] = 768; goodK[132] = 4; goodM[132] = 48;
	goodN[133] = 784; goodK[133] = 4; goodM[133] = 49;
	goodN[134] = 800; goodK[134] = 4; goodM[134] = 50;
	goodN[135] = 816; goodK[135] = 4; goodM[135] = 51;
	goodN[136] = 832; goodK[136] = 4; goodM[136] = 52;
	goodN[137] = 848; goodK[137] = 4; goodM[137] = 53;
	goodN[138] = 864; goodK[138] = 4; goodM[138] = 54;
	goodN[139] = 896; goodK[139] = 5; goodM[139] = 28;
	goodN[140] = 928; goodK[140] = 5; goodM[140] = 29;
	goodN[141] = 960; goodK[141] = 5; goodM[141] = 30;
	goodN[142] = 992; goodK[142] = 5; goodM[142] = 31;
	goodN[143] = 1024; goodK[143] = 5; goodM[143] = 32;
	goodN[144] = 1056; goodK[144] = 5; goodM[144] = 33;
	goodN[145] = 1088; goodK[145] = 5; goodM[145] = 34;
	goodN[146] = 1120; goodK[146] = 5; goodM[146] = 35;
	goodN[147] = 1152; goodK[147] = 5; goodM[147] = 36;
	goodN[148] = 1184; goodK[148] = 5; goodM[148] = 37;
	goodN[149] = 1216; goodK[149] = 5; goodM[149] = 38;
	goodN[150] = 1248; goodK[150] = 5; goodM[150] = 39;
	goodN[151] = 1280; goodK[151] = 5; goodM[151] = 40;
	goodN[152] = 1312; goodK[152] = 5; goodM[152] = 41;
	goodN[153] = 1344; goodK[153] = 5; goodM[153] = 42;
	goodN[154] = 1376; goodK[154] = 5; goodM[154] = 43;
	goodN[155] = 1408; goodK[155] = 5; goodM[155] = 44;
	goodN[156] = 1440; goodK[156] = 5; goodM[156] = 45;
	goodN[157] = 1472; goodK[157] = 5; goodM[157] = 46;
	goodN[158] = 1504; goodK[158] = 5; goodM[158] = 47;
	goodN[159] = 1536; goodK[159] = 5; goodM[159] = 48;
	goodN[160] = 1568; goodK[160] = 5; goodM[160] = 49;
	goodN[161] = 1600; goodK[161] = 5; goodM[161] = 50;
	goodN[162] = 1632; goodK[162] = 5; goodM[162] = 51;
	goodN[163] = 1664; goodK[163] = 5; goodM[163] = 52;
	goodN[164] = 1696; goodK[164] = 5; goodM[164] = 53;
	goodN[165] = 1728; goodK[165] = 6; goodM[165] = 27;
	goodN[166] = 1792; goodK[166] = 6; goodM[166] = 28;
	goodN[167] = 1856; goodK[167] = 6; goodM[167] = 29;
	goodN[168] = 1920; goodK[168] = 6; goodM[168] = 30;
	goodN[169] = 1984; goodK[169] = 6; goodM[169] = 31;
	goodN[170] = 2048; goodK[170] = 6; goodM[170] = 32;
	goodN[171] = 2112; goodK[171] = 6; goodM[171] = 33;
	goodN[172] = 2176; goodK[172] = 6; goodM[172] = 34;
	goodN[173] = 2240; goodK[173] = 6; goodM[173] = 35;
	goodN[174] = 2304; goodK[174] = 6; goodM[174] = 36;
	goodN[175] = 2368; goodK[175] = 6; goodM[175] = 37;
	goodN[176] = 2432; goodK[176] = 6; goodM[176] = 38;
	goodN[177] = 2496; goodK[177] = 6; goodM[177] = 39;
	goodN[178] = 2560; goodK[178] = 6; goodM[178] = 40;
	goodN[179] = 2624; goodK[179] = 6; goodM[179] = 41;
	goodN[180] = 2688; goodK[180] = 6; goodM[180] = 42;
	goodN[181] = 2752; goodK[181] = 6; goodM[181] = 43;
	goodN[182] = 2816; goodK[182] = 6; goodM[182] = 44;
	goodN[183] = 2880; goodK[183] = 6; goodM[183] = 45;
	goodN[184] = 2944; goodK[184] = 6; goodM[184] = 46;
	goodN[185] = 3008; goodK[185] = 6; goodM[185] = 47;
	goodN[186] = 3072; goodK[186] = 6; goodM[186] = 48;
	goodN[187] = 3136; goodK[187] = 6; goodM[187] = 49;
	goodN[188] = 3200; goodK[188] = 6; goodM[188] = 50;
	goodN[189] = 3264; goodK[189] = 6; goodM[189] = 51;
	goodN[190] = 3328; goodK[190] = 6; goodM[190] = 52;
	goodN[191] = 3456; goodK[191] = 7; goodM[191] = 27;
	goodN[192] = 3584; goodK[192] = 7; goodM[192] = 28;
	goodN[193] = 3712; goodK[193] = 7; goodM[193] = 29;
	goodN[194] = 3840; goodK[194] = 7; goodM[194] = 30;
	goodN[195] = 3968; goodK[195] = 7; goodM[195] = 31;
	goodN[196] = 4096; goodK[196] = 7; goodM[196] = 32;
	goodN[197] = 4224; goodK[197] = 7; goodM[197] = 33;
	goodN[198] = 4352; goodK[198] = 7; goodM[198] = 34;
	goodN[199] = 4480; goodK[199] = 7; goodM[199] = 35;
	goodN[200] = 4608; goodK[200] = 7; goodM[200] = 36;
	goodN[201] = 4736; goodK[201] = 7; goodM[201] = 37;
	goodN[202] = 4864; goodK[202] = 7; goodM[202] = 38;
	goodN[203] = 4992; goodK[203] = 7; goodM[203] = 39;
	goodN[204] = 5120; goodK[204] = 7; goodM[204] = 40;
	goodN[205] = 5248; goodK[205] = 7; goodM[205] = 41;
	goodN[206] = 5376; goodK[206] = 7; goodM[206] = 42;
	goodN[207] = 5504; goodK[207] = 7; goodM[207] = 43;
	goodN[208] = 5632; goodK[208] = 7; goodM[208] = 44;
	goodN[209] = 5760; goodK[209] = 7; goodM[209] = 45;
	goodN[210] = 5888; goodK[210] = 7; goodM[210] = 46;
	goodN[211] = 6016; goodK[211] = 7; goodM[211] = 47;
	goodN[212] = 6144; goodK[212] = 7; goodM[212] = 48;
	goodN[213] = 6272; goodK[213] = 7; goodM[213] = 49;
	goodN[214] = 6400; goodK[214] = 7; goodM[214] = 50;
	goodN[215] = 6528; goodK[215] = 7; goodM[215] = 51;
	goodN[216] = 6656; goodK[216] = 8; goodM[216] = 26;
	goodN[217] = 6912; goodK[217] = 8; goodM[217] = 27;
	goodN[218] = 7168; goodK[218] = 8; goodM[218] = 28;
	goodN[219] = 7424; goodK[219] = 8; goodM[219] = 29;
	goodN[220] = 7680; goodK[220] = 8; goodM[220] = 30;
	goodN[221] = 7936; goodK[221] = 8; goodM[221] = 31;
	goodN[222] = 8192; goodK[222] = 8; goodM[222] = 32;
	goodN[223] = 8448; goodK[223] = 8; goodM[223] = 33;
	goodN[224] = 8704; goodK[224] = 8; goodM[224] = 34;
	goodN[225] = 8960; goodK[225] = 8; goodM[225] = 35;
	goodN[226] = 9216; goodK[226] = 8; goodM[226] = 36;
	goodN[227] = 9472; goodK[227] = 8; goodM[227] = 37;
	goodN[228] = 9728; goodK[228] = 8; goodM[228] = 38;
	goodN[229] = 9984; goodK[229] = 8; goodM[229] = 39;
	goodN[230] = 10240; goodK[230] = 8; goodM[230] = 40;
	goodN[231] = 10496; goodK[231] = 8; goodM[231] = 41;
	goodN[232] = 10752; goodK[232] = 8; goodM[232] = 42;
	goodN[233] = 11008; goodK[233] = 8; goodM[233] = 43;
	goodN[234] = 11264; goodK[234] = 8; goodM[234] = 44;
	goodN[235] = 11520; goodK[235] = 8; goodM[235] = 45;
	goodN[236] = 11776; goodK[236] = 8; goodM[236] = 46;
	goodN[237] = 12032; goodK[237] = 8; goodM[237] = 47;
	goodN[238] = 12288; goodK[238] = 8; goodM[238] = 48;
	goodN[239] = 12544; goodK[239] = 8; goodM[239] = 49;
	goodN[240] = 12800; goodK[240] = 8; goodM[240] = 50;
	goodN[241] = 13312; goodK[241] = 9; goodM[241] = 26;
	goodN[242] = 13824; goodK[242] = 9; goodM[242] = 27;
	goodN[243] = 14336; goodK[243] = 9; goodM[243] = 28;
	goodN[244] = 14848; goodK[244] = 9; goodM[244] = 29;
	goodN[245] = 15360; goodK[245] = 9; goodM[245] = 30;
	goodN[246] = 15872; goodK[246] = 9; goodM[246] = 31;
	goodN[247] = 16384; goodK[247] = 9; goodM[247] = 32;
	goodN[248] = 16896; goodK[248] = 9; goodM[248] = 33;
	goodN[249] = 17408; goodK[249] = 9; goodM[249] = 34;
	goodN[250] = 17920; goodK[250] = 9; goodM[250] = 35;
	goodN[251] = 18432; goodK[251] = 9; goodM[251] = 36;
	goodN[252] = 18944; goodK[252] = 9; goodM[252] = 37;
	goodN[253] = 19456; goodK[253] = 9; goodM[253] = 38;
	goodN[254] = 19968; goodK[254] = 9; goodM[254] = 39;
	goodN[255] = 20480; goodK[255] = 9; goodM[255] = 40;
	goodN[256] = 20992; goodK[256] = 9; goodM[256] = 41;
	goodN[257] = 21504; goodK[257] = 9; goodM[257] = 42;
	goodN[258] = 22016; goodK[258] = 9; goodM[258] = 43;
	goodN[259] = 22528; goodK[259] = 9; goodM[259] = 44;
	goodN[260] = 23040; goodK[260] = 9; goodM[260] = 45;
	goodN[261] = 23552; goodK[261] = 9; goodM[261] = 46;
	goodN[262] = 24064; goodK[262] = 9; goodM[262] = 47;
	goodN[263] = 24576; goodK[263] = 9; goodM[263] = 48;
	goodN[264] = 25088; goodK[264] = 9; goodM[264] = 49;
	goodN[265] = 25600; goodK[265] = 10; goodM[265] = 25;
	goodN[266] = 26624; goodK[266] = 10; goodM[266] = 26;
	goodN[267] = 27648; goodK[267] = 10; goodM[267] = 27;
	goodN[268] = 28672; goodK[268] = 10; goodM[268] = 28;
	goodN[269] = 29696; goodK[269] = 10; goodM[269] = 29;
	goodN[270] = 30720; goodK[270] = 10; goodM[270] = 30;
	goodN[271] = 31744; goodK[271] = 10; goodM[271] = 31;
	goodN[272] = 32768; goodK[272] = 10; goodM[272] = 32;
	goodN[273] = 33792; goodK[273] = 10; goodM[273] = 33;
	goodN[274] = 34816; goodK[274] = 10; goodM[274] = 34;
	goodN[275] = 35840; goodK[275] = 10; goodM[275] = 35;
	goodN[276] = 36864; goodK[276] = 10; goodM[276] = 36;
	goodN[277] = 37888; goodK[277] = 10; goodM[277] = 37;
	goodN[278] = 38912; goodK[278] = 10; goodM[278] = 38;
	goodN[279] = 39936; goodK[279] = 10; goodM[279] = 39;
	goodN[280] = 40960; goodK[280] = 10; goodM[280] = 40;
	goodN[281] = 41984; goodK[281] = 10; goodM[281] = 41;
	goodN[282] = 43008; goodK[282] = 10; goodM[282] = 42;
	goodN[283] = 44032; goodK[283] = 10; goodM[283] = 43;
	goodN[284] = 45056; goodK[284] = 10; goodM[284] = 44;
	goodN[285] = 46080; goodK[285] = 10; goodM[285] = 45;
	goodN[286] = 47104; goodK[286] = 10; goodM[286] = 46;
	goodN[287] = 48128; goodK[287] = 10; goodM[287] = 47;
	goodN[288] = 49152; goodK[288] = 10; goodM[288] = 48;
	goodN[289] = 51200; goodK[289] = 11; goodM[289] = 25;
	goodN[290] = 53248; goodK[290] = 11; goodM[290] = 26;
	goodN[291] = 55296; goodK[291] = 11; goodM[291] = 27;
	goodN[292] = 57344; goodK[292] = 11; goodM[292] = 28;
	goodN[293] = 59392; goodK[293] = 11; goodM[293] = 29;
	goodN[294] = 61440; goodK[294] = 11; goodM[294] = 30;
	goodN[295] = 63488; goodK[295] = 11; goodM[295] = 31;
	goodN[296] = 65536; goodK[296] = 11; goodM[296] = 32;
	goodN[297] = 67584; goodK[297] = 11; goodM[297] = 33;
	goodN[298] = 69632; goodK[298] = 11; goodM[298] = 34;
	goodN[299] = 71680; goodK[299] = 11; goodM[299] = 35;
	goodN[300] = 73728; goodK[300] = 11; goodM[300] = 36;
	goodN[301] = 75776; goodK[301] = 11; goodM[301] = 37;
	goodN[302] = 77824; goodK[302] = 11; goodM[302] = 38;
	goodN[303] = 79872; goodK[303] = 11; goodM[303] = 39;
	goodN[304] = 81920; goodK[304] = 11; goodM[304] = 40;
	goodN[305] = 83968; goodK[305] = 11; goodM[305] = 41;
	goodN[306] = 86016; goodK[306] = 11; goodM[306] = 42;
	goodN[307] = 88064; goodK[307] = 11; goodM[307] = 43;
	goodN[308] = 90112; goodK[308] = 11; goodM[308] = 44;
	goodN[309] = 92160; goodK[309] = 11; goodM[309] = 45;
	goodN[310] = 94208; goodK[310] = 11; goodM[310] = 46;
	goodN[311] = 96256; goodK[311] = 11; goodM[311] = 47;
	goodN[312] = 98304; goodK[312] = 12; goodM[312] = 24;
	goodN[313] = 102400; goodK[313] = 12; goodM[313] = 25;
	goodN[314] = 106496; goodK[314] = 12; goodM[314] = 26;
	goodN[315] = 110592; goodK[315] = 12; goodM[315] = 27;
	goodN[316] = 114688; goodK[316] = 12; goodM[316] = 28;
	goodN[317] = 118784; goodK[317] = 12; goodM[317] = 29;
	goodN[318] = 122880; goodK[318] = 12; goodM[318] = 30;
	goodN[319] = 126976; goodK[319] = 12; goodM[319] = 31;
	goodN[320] = 131072; goodK[320] = 12; goodM[320] = 32;
	goodN[321] = 135168; goodK[321] = 12; goodM[321] = 33;
	goodN[322] = 139264; goodK[322] = 12; goodM[322] = 34;
	goodN[323] = 143360; goodK[323] = 12; goodM[323] = 35;
	goodN[324] = 147456; goodK[324] = 12; goodM[324] = 36;
	goodN[325] = 151552; goodK[325] = 12; goodM[325] = 37;
	goodN[326] = 155648; goodK[326] = 12; goodM[326] = 38;
	goodN[327] = 159744; goodK[327] = 12; goodM[327] = 39;
	goodN[328] = 163840; goodK[328] = 12; goodM[328] = 40;
	goodN[329] = 167936; goodK[329] = 12; goodM[329] = 41;
	goodN[330] = 172032; goodK[330] = 12; goodM[330] = 42;
	goodN[331] = 176128; goodK[331] = 12; goodM[331] = 43;
	goodN[332] = 180224; goodK[332] = 12; goodM[332] = 44;
	goodN[333] = 184320; goodK[333] = 12; goodM[333] = 45;
	goodN[334] = 188416; goodK[334] = 12; goodM[334] = 46;
	goodN[335] = 196608; goodK[335] = 13; goodM[335] = 24;
	goodN[336] = 204800; goodK[336] = 13; goodM[336] = 25;
	goodN[337] = 212992; goodK[337] = 13; goodM[337] = 26;
	goodN[338] = 221184; goodK[338] = 13; goodM[338] = 27;
	goodN[339] = 229376; goodK[339] = 13; goodM[339] = 28;
	goodN[340] = 237568; goodK[340] = 13; goodM[340] = 29;
	goodN[341] = 245760; goodK[341] = 13; goodM[341] = 30;
	goodN[342] = 253952; goodK[342] = 13; goodM[342] = 31;
	goodN[343] = 262144; goodK[343] = 13; goodM[343] = 32;
	goodN[344] = 270336; goodK[344] = 13; goodM[344] = 33;
	goodN[345] = 278528; goodK[345] = 13; goodM[345] = 34;
	goodN[346] = 286720; goodK[346] = 13; goodM[346] = 35;
	goodN[347] = 294912; goodK[347] = 13; goodM[347] = 36;
	goodN[348] = 303104; goodK[348] = 13; goodM[348] = 37;
	goodN[349] = 311296; goodK[349] = 13; goodM[349] = 38;
	goodN[350] = 319488; goodK[350] = 13; goodM[350] = 39;
	goodN[351] = 327680; goodK[351] = 13; goodM[351] = 40;
	goodN[352] = 335872; goodK[352] = 13; goodM[352] = 41;
	goodN[353] = 344064; goodK[353] = 13; goodM[353] = 42;
	goodN[354] = 352256; goodK[354] = 13; goodM[354] = 43;
	goodN[355] = 360448; goodK[355] = 13; goodM[355] = 44;
	goodN[356] = 368640; goodK[356] = 13; goodM[356] = 45;
	goodN[357] = 376832; goodK[357] = 14; goodM[357] = 23;
	goodN[358] = 393216; goodK[358] = 14; goodM[358] = 24;
	goodN[359] = 409600; goodK[359] = 14; goodM[359] = 25;
	goodN[360] = 425984; goodK[360] = 14; goodM[360] = 26;
	goodN[361] = 442368; goodK[361] = 14; goodM[361] = 27;
	goodN[362] = 458752; goodK[362] = 14; goodM[362] = 28;
	goodN[363] = 475136; goodK[363] = 14; goodM[363] = 29;
	goodN[364] = 491520; goodK[364] = 14; goodM[364] = 30;
	goodN[365] = 507904; goodK[365] = 14; goodM[365] = 31;
	goodN[366] = 524288; goodK[366] = 14; goodM[366] = 32;
	goodN[367] = 540672; goodK[367] = 14; goodM[367] = 33;
	goodN[368] = 557056; goodK[368] = 14; goodM[368] = 34;
	goodN[369] = 573440; goodK[369] = 14; goodM[369] = 35;
	goodN[370] = 589824; goodK[370] = 14; goodM[370] = 36;
	goodN[371] = 606208; goodK[371] = 14; goodM[371] = 37;
	goodN[372] = 622592; goodK[372] = 14; goodM[372] = 38;
	goodN[373] = 638976; goodK[373] = 14; goodM[373] = 39;
	goodN[374] = 655360; goodK[374] = 14; goodM[374] = 40;
	goodN[375] = 671744; goodK[375] = 14; goodM[375] = 41;
	goodN[376] = 688128; goodK[376] = 14; goodM[376] = 42;
	goodN[377] = 704512; goodK[377] = 14; goodM[377] = 43;
	goodN[378] = 720896; goodK[378] = 14; goodM[378] = 44;
	goodN[379] = 753664; goodK[379] = 15; goodM[379] = 23;
	goodN[380] = 786432; goodK[380] = 15; goodM[380] = 24;
	goodN[381] = 819200; goodK[381] = 15; goodM[381] = 25;
	goodN[382] = 851968; goodK[382] = 15; goodM[382] = 26;
	goodN[383] = 884736; goodK[383] = 15; goodM[383] = 27;
	goodN[384] = 917504; goodK[384] = 15; goodM[384] = 28;
	goodN[385] = 950272; goodK[385] = 15; goodM[385] = 29;
	goodN[386] = 983040; goodK[386] = 15; goodM[386] = 30;
	goodN[387] = 1015808; goodK[387] = 15; goodM[387] = 31;
	goodN[388] = 1048576; goodK[388] = 15; goodM[388] = 32;
	goodN[389] = 1081344; goodK[389] = 15; goodM[389] = 33;
	goodN[390] = 1114112; goodK[390] = 15; goodM[390] = 34;
	goodN[391] = 1146880; goodK[391] = 15; goodM[391] = 35;
	goodN[392] = 1179648; goodK[392] = 15; goodM[392] = 36;
	goodN[393] = 1212416; goodK[393] = 15; goodM[393] = 37;
	goodN[394] = 1245184; goodK[394] = 15; goodM[394] = 38;
	goodN[395] = 1277952; goodK[395] = 15; goodM[395] = 39;
	goodN[396] = 1310720; goodK[396] = 15; goodM[396] = 40;
	goodN[397] = 1343488; goodK[397] = 15; goodM[397] = 41;
	goodN[398] = 1376256; goodK[398] = 15; goodM[398] = 42;
	goodN[399] = 1409024; goodK[399] = 15; goodM[399] = 43;
	goodN[400] = 1441792; goodK[400] = 16; goodM[400] = 22;
	goodN[401] = 1507328; goodK[401] = 16; goodM[401] = 23;
	goodN[402] = 1572864; goodK[402] = 16; goodM[402] = 24;
	goodN[403] = 1638400; goodK[403] = 16; goodM[403] = 25;
	goodN[404] = 1703936; goodK[404] = 16; goodM[404] = 26;
	goodN[405] = 1769472; goodK[405] = 16; goodM[405] = 27;
	goodN[406] = 1835008; goodK[406] = 16; goodM[406] = 28;
	goodN[407] = 1900544; goodK[407] = 16; goodM[407] = 29;
	goodN[408] = 1966080; goodK[408] = 16; goodM[408] = 30;
	goodN[409] = 2031616; goodK[409] = 16; goodM[409] = 31;
	goodN[410] = 2097152; goodK[410] = 16; goodM[410] = 32;
	goodN[411] = 2162688; goodK[411] = 16; goodM[411] = 33;
	goodN[412] = 2228224; goodK[412] = 16; goodM[412] = 34;
	goodN[413] = 2293760; goodK[413] = 16; goodM[413] = 35;
	goodN[414] = 2359296; goodK[414] = 16; goodM[414] = 36;
	goodN[415] = 2424832; goodK[415] = 16; goodM[415] = 37;
	goodN[416] = 2490368; goodK[416] = 16; goodM[416] = 38;
	goodN[417] = 2555904; goodK[417] = 16; goodM[417] = 39;
	goodN[418] = 2621440; goodK[418] = 16; goodM[418] = 40;
	goodN[419] = 2686976; goodK[419] = 16; goodM[419] = 41;
	goodN[420] = 2752512; goodK[420] = 16; goodM[420] = 42;
	goodN[421] = 2883584; goodK[421] = 17; goodM[421] = 22;
	goodN[422] = 3014656; goodK[422] = 17; goodM[422] = 23;
	goodN[423] = 3145728; goodK[423] = 17; goodM[423] = 24;
	goodN[424] = 3276800; goodK[424] = 17; goodM[424] = 25;
	goodN[425] = 3407872; goodK[425] = 17; goodM[425] = 26;
	goodN[426] = 3538944; goodK[426] = 17; goodM[426] = 27;
	goodN[427] = 3670016; goodK[427] = 17; goodM[427] = 28;
	goodN[428] = 3801088; goodK[428] = 17; goodM[428] = 29;
	goodN[429] = 3932160; goodK[429] = 17; goodM[429] = 30;
	goodN[430] = 4063232; goodK[430] = 17; goodM[430] = 31;
	goodN[431] = 4194304; goodK[431] = 17; goodM[431] = 32;

  }


  inline void init1() {
	size = 173;
	goodN = new int[size];
	goodM = new int[size];
	goodK = new int[size];

	goodN[0] = 2; goodK[0] = 1; goodM[0] = 1;
	goodN[1] = 4; goodK[1] = 1; goodM[1] = 2;
	goodN[2] = 6; goodK[2] = 1; goodM[2] = 3;
	goodN[3] = 8; goodK[3] = 1; goodM[3] = 4;
	goodN[4] = 10; goodK[4] = 1; goodM[4] = 5;
	goodN[5] = 12; goodK[5] = 1; goodM[5] = 6;
	goodN[6] = 14; goodK[6] = 1; goodM[6] = 7;
	goodN[7] = 16; goodK[7] = 1; goodM[7] = 8;
	goodN[8] = 18; goodK[8] = 1; goodM[8] = 9;
	goodN[9] = 20; goodK[9] = 1; goodM[9] = 10;
	goodN[10] = 22; goodK[10] = 1; goodM[10] = 11;
	goodN[11] = 24; goodK[11] = 1; goodM[11] = 12;
	goodN[12] = 26; goodK[12] = 1; goodM[12] = 13;
	goodN[13] = 28; goodK[13] = 1; goodM[13] = 14;
	goodN[14] = 30; goodK[14] = 1; goodM[14] = 15;
	goodN[15] = 32; goodK[15] = 1; goodM[15] = 16;
	goodN[16] = 34; goodK[16] = 1; goodM[16] = 17;
	goodN[17] = 36; goodK[17] = 1; goodM[17] = 18;
	goodN[18] = 38; goodK[18] = 1; goodM[18] = 19;
	goodN[19] = 40; goodK[19] = 1; goodM[19] = 20;
	goodN[20] = 42; goodK[20] = 1; goodM[20] = 21;
	goodN[21] = 44; goodK[21] = 1; goodM[21] = 22;
	goodN[22] = 46; goodK[22] = 1; goodM[22] = 23;
	goodN[23] = 48; goodK[23] = 1; goodM[23] = 24;
	goodN[24] = 50; goodK[24] = 1; goodM[24] = 25;
	goodN[25] = 52; goodK[25] = 1; goodM[25] = 26;
	goodN[26] = 56; goodK[26] = 2; goodM[26] = 14;
	goodN[27] = 60; goodK[27] = 2; goodM[27] = 15;
	goodN[28] = 64; goodK[28] = 2; goodM[28] = 16;
	goodN[29] = 68; goodK[29] = 2; goodM[29] = 17;
	goodN[30] = 72; goodK[30] = 2; goodM[30] = 18;
	goodN[31] = 76; goodK[31] = 2; goodM[31] = 19;
	goodN[32] = 80; goodK[32] = 2; goodM[32] = 20;
	goodN[33] = 84; goodK[33] = 2; goodM[33] = 21;
	goodN[34] = 88; goodK[34] = 2; goodM[34] = 22;
	goodN[35] = 92; goodK[35] = 2; goodM[35] = 23;
	goodN[36] = 96; goodK[36] = 2; goodM[36] = 24;
	goodN[37] = 100; goodK[37] = 2; goodM[37] = 25;
	goodN[38] = 104; goodK[38] = 3; goodM[38] = 13;
	goodN[39] = 112; goodK[39] = 3; goodM[39] = 14;
	goodN[40] = 120; goodK[40] = 3; goodM[40] = 15;
	goodN[41] = 128; goodK[41] = 3; goodM[41] = 16;
	goodN[42] = 136; goodK[42] = 3; goodM[42] = 17;
	goodN[43] = 144; goodK[43] = 3; goodM[43] = 18;
	goodN[44] = 152; goodK[44] = 3; goodM[44] = 19;
	goodN[45] = 160; goodK[45] = 3; goodM[45] = 20;
	goodN[46] = 168; goodK[46] = 3; goodM[46] = 21;
	goodN[47] = 176; goodK[47] = 3; goodM[47] = 22;
	goodN[48] = 184; goodK[48] = 3; goodM[48] = 23;
	goodN[49] = 192; goodK[49] = 3; goodM[49] = 24;
	goodN[50] = 208; goodK[50] = 4; goodM[50] = 13;
	goodN[51] = 224; goodK[51] = 4; goodM[51] = 14;
	goodN[52] = 240; goodK[52] = 4; goodM[52] = 15;
	goodN[53] = 256; goodK[53] = 4; goodM[53] = 16;
	goodN[54] = 272; goodK[54] = 4; goodM[54] = 17;
	goodN[55] = 288; goodK[55] = 4; goodM[55] = 18;
	goodN[56] = 304; goodK[56] = 4; goodM[56] = 19;
	goodN[57] = 320; goodK[57] = 4; goodM[57] = 20;
	goodN[58] = 336; goodK[58] = 4; goodM[58] = 21;
	goodN[59] = 352; goodK[59] = 4; goodM[59] = 22;
	goodN[60] = 368; goodK[60] = 4; goodM[60] = 23;
	goodN[61] = 384; goodK[61] = 5; goodM[61] = 12;
	goodN[62] = 416; goodK[62] = 5; goodM[62] = 13;
	goodN[63] = 448; goodK[63] = 5; goodM[63] = 14;
	goodN[64] = 480; goodK[64] = 5; goodM[64] = 15;
	goodN[65] = 512; goodK[65] = 5; goodM[65] = 16;
	goodN[66] = 544; goodK[66] = 5; goodM[66] = 17;
	goodN[67] = 576; goodK[67] = 5; goodM[67] = 18;
	goodN[68] = 608; goodK[68] = 5; goodM[68] = 19;
	goodN[69] = 640; goodK[69] = 5; goodM[69] = 20;
	goodN[70] = 672; goodK[70] = 5; goodM[70] = 21;
	goodN[71] = 704; goodK[71] = 5; goodM[71] = 22;
	goodN[72] = 768; goodK[72] = 6; goodM[72] = 12;
	goodN[73] = 832; goodK[73] = 6; goodM[73] = 13;
	goodN[74] = 896; goodK[74] = 6; goodM[74] = 14;
	goodN[75] = 960; goodK[75] = 6; goodM[75] = 15;
	goodN[76] = 1024; goodK[76] = 6; goodM[76] = 16;
	goodN[77] = 1088; goodK[77] = 6; goodM[77] = 17;
	goodN[78] = 1152; goodK[78] = 6; goodM[78] = 18;
	goodN[79] = 1216; goodK[79] = 6; goodM[79] = 19;
	goodN[80] = 1280; goodK[80] = 6; goodM[80] = 20;
	goodN[81] = 1344; goodK[81] = 6; goodM[81] = 21;
	goodN[82] = 1408; goodK[82] = 7; goodM[82] = 11;
	goodN[83] = 1536; goodK[83] = 7; goodM[83] = 12;
	goodN[84] = 1664; goodK[84] = 7; goodM[84] = 13;
	goodN[85] = 1792; goodK[85] = 7; goodM[85] = 14;
	goodN[86] = 1920; goodK[86] = 7; goodM[86] = 15;
	goodN[87] = 2048; goodK[87] = 7; goodM[87] = 16;
	goodN[88] = 2176; goodK[88] = 7; goodM[88] = 17;
	goodN[89] = 2304; goodK[89] = 7; goodM[89] = 18;
	goodN[90] = 2432; goodK[90] = 7; goodM[90] = 19;
	goodN[91] = 2560; goodK[91] = 7; goodM[91] = 20;
	goodN[92] = 2816; goodK[92] = 8; goodM[92] = 11;
	goodN[93] = 3072; goodK[93] = 8; goodM[93] = 12;
	goodN[94] = 3328; goodK[94] = 8; goodM[94] = 13;
	goodN[95] = 3584; goodK[95] = 8; goodM[95] = 14;
	goodN[96] = 3840; goodK[96] = 8; goodM[96] = 15;
	goodN[97] = 4096; goodK[97] = 8; goodM[97] = 16;
	goodN[98] = 4352; goodK[98] = 8; goodM[98] = 17;
	goodN[99] = 4608; goodK[99] = 8; goodM[99] = 18;
	goodN[100] = 4864; goodK[100] = 8; goodM[100] = 19;
	goodN[101] = 5120; goodK[101] = 9; goodM[101] = 10;
	goodN[102] = 5632; goodK[102] = 9; goodM[102] = 11;
	goodN[103] = 6144; goodK[103] = 9; goodM[103] = 12;
	goodN[104] = 6656; goodK[104] = 9; goodM[104] = 13;
	goodN[105] = 7168; goodK[105] = 9; goodM[105] = 14;
	goodN[106] = 7680; goodK[106] = 9; goodM[106] = 15;
	goodN[107] = 8192; goodK[107] = 9; goodM[107] = 16;
	goodN[108] = 8704; goodK[108] = 9; goodM[108] = 17;
	goodN[109] = 9216; goodK[109] = 9; goodM[109] = 18;
	goodN[110] = 10240; goodK[110] = 10; goodM[110] = 10;
	goodN[111] = 11264; goodK[111] = 10; goodM[111] = 11;
	goodN[112] = 12288; goodK[112] = 10; goodM[112] = 12;
	goodN[113] = 13312; goodK[113] = 10; goodM[113] = 13;
	goodN[114] = 14336; goodK[114] = 10; goodM[114] = 14;
	goodN[115] = 15360; goodK[115] = 10; goodM[115] = 15;
	goodN[116] = 16384; goodK[116] = 10; goodM[116] = 16;
	goodN[117] = 17408; goodK[117] = 10; goodM[117] = 17;
	goodN[118] = 18432; goodK[118] = 11; goodM[118] = 9;
	goodN[119] = 20480; goodK[119] = 11; goodM[119] = 10;
	goodN[120] = 22528; goodK[120] = 11; goodM[120] = 11;
	goodN[121] = 24576; goodK[121] = 11; goodM[121] = 12;
	goodN[122] = 26624; goodK[122] = 11; goodM[122] = 13;
	goodN[123] = 28672; goodK[123] = 11; goodM[123] = 14;
	goodN[124] = 30720; goodK[124] = 11; goodM[124] = 15;
	goodN[125] = 32768; goodK[125] = 11; goodM[125] = 16;
	goodN[126] = 36864; goodK[126] = 12; goodM[126] = 9;
	goodN[127] = 40960; goodK[127] = 12; goodM[127] = 10;
	goodN[128] = 45056; goodK[128] = 12; goodM[128] = 11;
	goodN[129] = 49152; goodK[129] = 12; goodM[129] = 12;
	goodN[130] = 53248; goodK[130] = 12; goodM[130] = 13;
	goodN[131] = 57344; goodK[131] = 12; goodM[131] = 14;
	goodN[132] = 61440; goodK[132] = 12; goodM[132] = 15;
	goodN[133] = 65536; goodK[133] = 12; goodM[133] = 16;
	goodN[134] = 73728; goodK[134] = 13; goodM[134] = 9;
	goodN[135] = 81920; goodK[135] = 13; goodM[135] = 10;
	goodN[136] = 90112; goodK[136] = 13; goodM[136] = 11;
	goodN[137] = 98304; goodK[137] = 13; goodM[137] = 12;
	goodN[138] = 106496; goodK[138] = 13; goodM[138] = 13;
	goodN[139] = 114688; goodK[139] = 13; goodM[139] = 14;
	goodN[140] = 122880; goodK[140] = 13; goodM[140] = 15;
	goodN[141] = 131072; goodK[141] = 14; goodM[141] = 8;
	goodN[142] = 147456; goodK[142] = 14; goodM[142] = 9;
	goodN[143] = 163840; goodK[143] = 14; goodM[143] = 10;
	goodN[144] = 180224; goodK[144] = 14; goodM[144] = 11;
	goodN[145] = 196608; goodK[145] = 14; goodM[145] = 12;
	goodN[146] = 212992; goodK[146] = 14; goodM[146] = 13;
	goodN[147] = 229376; goodK[147] = 14; goodM[147] = 14;
	goodN[148] = 262144; goodK[148] = 15; goodM[148] = 8;
	goodN[149] = 294912; goodK[149] = 15; goodM[149] = 9;
	goodN[150] = 327680; goodK[150] = 15; goodM[150] = 10;
	goodN[151] = 360448; goodK[151] = 15; goodM[151] = 11;
	goodN[152] = 393216; goodK[152] = 15; goodM[152] = 12;
	goodN[153] = 425984; goodK[153] = 15; goodM[153] = 13;
	goodN[154] = 458752; goodK[154] = 16; goodM[154] = 7;
	goodN[155] = 524288; goodK[155] = 16; goodM[155] = 8;
	goodN[156] = 589824; goodK[156] = 16; goodM[156] = 9;
	goodN[157] = 655360; goodK[157] = 16; goodM[157] = 10;
	goodN[158] = 720896; goodK[158] = 16; goodM[158] = 11;
	goodN[159] = 786432; goodK[159] = 16; goodM[159] = 12;
	goodN[160] = 917504; goodK[160] = 17; goodM[160] = 7;
	goodN[161] = 1048576; goodK[161] = 17; goodM[161] = 8;
	goodN[162] = 1179648; goodK[162] = 17; goodM[162] = 9;
	goodN[163] = 1310720; goodK[163] = 17; goodM[163] = 10;
	goodN[164] = 1441792; goodK[164] = 17; goodM[164] = 11;
	goodN[165] = 1572864; goodK[165] = 18; goodM[165] = 6;
	goodN[166] = 1835008; goodK[166] = 18; goodM[166] = 7;
	goodN[167] = 2097152; goodK[167] = 18; goodM[167] = 8;
	goodN[168] = 2359296; goodK[168] = 18; goodM[168] = 9;
	goodN[169] = 2621440; goodK[169] = 18; goodM[169] = 10;
	goodN[170] = 3145728; goodK[170] = 19; goodM[170] = 6;
	goodN[171] = 3670016; goodK[171] = 19; goodM[171] = 7;
	goodN[172] = 4194304; goodK[172] = 19; goodM[172] = 8;

  }

  /*
  **
   * Find the next good N which is in the form of N = 2*k * M.
   * The global variables of M, K and N will be set based on 
   * the index of the good N.
   *
   * @param n the integer that might not be good (not in the form
   * 		of N = 2*k * M)
   */
  void determineGoodN(int n, sfixn d);
  
  /**
   * Search for the nearest good N to the input integer
   * using binary search on the goodN array.
   */
  int binarySearch(int n, int start, int end);
  
  /**
   * Multiply two univariate polynomials using the 
   * Schonhage-Strassen Algorithm.
   */
  void mul2C(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c);
  void mul2C1(UnivariateIntegerPolynomial*, UnivariateIntegerPolynomial*, UnivariateIntegerPolynomial*);
  void mul2C2(UnivariateIntegerPolynomial*, UnivariateIntegerPolynomial*, UnivariateIntegerPolynomial*);
  
  /**
   * for both C^+ (ncc) and C^- (cc), for each coefficient of y^i,
   * which is a large integer encoded as a polynomial in x, i.e.
   * ncc_i(x), cc_i(x), 
   * (1) convert ncc_i and cc_i to GMP, get u_i and v_i
   * (2) compute c_i = (u_i+v_i)/2+(-u_i+v_i)/2*2^N
   */
  void RecoverProduct(sfixn*, int, int, mpz_class*);
  void RecoverProduct(sfixn*, sfixn*, int, int, mpz_class*);  

  void ToMPZ(mpz_t, sfixn*);
  void CRTtoMPZ(mpz_t, sfixn*, sfixn*, int);
  
  /**
   * Convert a univariate polynomial with large coefficients to 
   * a bivariate polynomial representation with relatively 
   * small coefficients.
   *
   * univariate --> a1    * y ^ 0 + a2    * y ^ 1 + ... + ad    * y ^ (d-1)
   * bivariate  --> A1(x) * y ^ 0 + A2(x) * y ^ 1 + ... + Ad(x) * y ^ (d-1)
   *
   * Ai(x) = b1 * x ^ 0 + b2 * x ^ 1 + ... + bK * x ^ (K-1)
   */
  BivariatePolynomial * ToBivarMod(mpz_class*, int, sfixn);
  BivariatePolynomial * ToBivarMod(mpz_class*, int, sfixn, sfixn);

  void mpzToPolyMod(sfixn*, mpz_t, int, sfixn);
  void mpzToPolyMod(sfixn*, sfixn*, mpz_t, int, sfixn, sfixn);
};

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


