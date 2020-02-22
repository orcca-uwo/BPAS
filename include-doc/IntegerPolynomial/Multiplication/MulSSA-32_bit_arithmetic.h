#pragma once
/**
   Implementation of polynomial multiplication using Schonhage - Strassen Algorithm.
 
   @author Farnam Mansouri and Yuzhen Xie
 */

#include "Mul.h"
#include "../BivariatePoly.h"

#include "../../FFT/src/general_routine.h"
#include "../../FFT/src/basic_routine.h"

class MulSSA: public Mul{
 public:
  int size;
  int* goodN;
  int* goodM;
  int* goodK;		
  
  static const sfixn prime1 = 962592769; // The first prime for CRT
  static const sfixn prime2 = 950009857; // The second prime for CRT
  static const sfixn prime3 = 940572673; // The third prime for CRT
  static const sfixn U1 = 475004853;
  static const sfixn U2 = -481296308;
  
  long int P1_P2;        //prime1 * prime2;
  long int HALF_P1_P2;   //P1_P2 / 2;
  long int N_HALF_P1_P2; //HALF_P1_P2;
  long int P1_U1;        //prime1*U1;
  long int P2_U2;        //prime2*U2;
  unsigned long int LIMB_BITS; 	//bound for u_i and v_i

  long int P2_P3, P3_P1;
  int U23, U31, U12;
  __int128 P1_P2_P3;	// prime1 * prime2 * prime3
  __int128 HALF_P1_P2_P3;
  __int128 N_HALF_P1_P2_P3;

  //static const sfixn P = 180143985094819841;
  
  int N;
  int K;
  int M;
  
  void init () {
	size = 179;
	goodN = new int[size];
	goodM = new int[size];
	goodK = new int[size];

	// 2 primes
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
	goodN[24] = 52; goodK[24] = 2; goodM[24] = 13;
	goodN[25] = 56; goodK[25] = 2; goodM[25] = 14;
	goodN[26] = 60; goodK[26] = 2; goodM[26] = 15;
	goodN[27] = 64; goodK[27] = 2; goodM[27] = 16;
	goodN[28] = 68; goodK[28] = 2; goodM[28] = 17;
	goodN[29] = 72; goodK[29] = 2; goodM[29] = 18;
	goodN[30] = 76; goodK[30] = 2; goodM[30] = 19;
	goodN[31] = 80; goodK[31] = 2; goodM[31] = 20;
	goodN[32] = 84; goodK[32] = 2; goodM[32] = 21;
	goodN[33] = 88; goodK[33] = 2; goodM[33] = 22;
	goodN[34] = 92; goodK[34] = 2; goodM[34] = 23;
	goodN[35] = 96; goodK[35] = 3; goodM[35] = 12;
	goodN[36] = 104; goodK[36] = 3; goodM[36] = 13;
	goodN[37] = 112; goodK[37] = 3; goodM[37] = 14;
	goodN[38] = 120; goodK[38] = 3; goodM[38] = 15;
	goodN[39] = 128; goodK[39] = 3; goodM[39] = 16;
	goodN[40] = 136; goodK[40] = 3; goodM[40] = 17;
	goodN[41] = 144; goodK[41] = 3; goodM[41] = 18;
	goodN[42] = 152; goodK[42] = 3; goodM[42] = 19;
	goodN[43] = 160; goodK[43] = 3; goodM[43] = 20;
	goodN[44] = 168; goodK[44] = 3; goodM[44] = 21;
	goodN[45] = 176; goodK[45] = 3; goodM[45] = 22;
	goodN[46] = 192; goodK[46] = 4; goodM[46] = 12;
	goodN[47] = 208; goodK[47] = 4; goodM[47] = 13;
	goodN[48] = 224; goodK[48] = 4; goodM[48] = 14;
	goodN[49] = 240; goodK[49] = 4; goodM[49] = 15;
	goodN[50] = 256; goodK[50] = 4; goodM[50] = 16;
	goodN[51] = 272; goodK[51] = 4; goodM[51] = 17;
	goodN[52] = 288; goodK[52] = 4; goodM[52] = 18;
	goodN[53] = 304; goodK[53] = 4; goodM[53] = 19;
	goodN[54] = 320; goodK[54] = 4; goodM[54] = 20;
	goodN[55] = 336; goodK[55] = 4; goodM[55] = 21;
	goodN[56] = 352; goodK[56] = 5; goodM[56] = 11;
	goodN[57] = 384; goodK[57] = 5; goodM[57] = 12;
	goodN[58] = 416; goodK[58] = 5; goodM[58] = 13;
	goodN[59] = 448; goodK[59] = 5; goodM[59] = 14;
	goodN[60] = 480; goodK[60] = 5; goodM[60] = 15;
	goodN[61] = 512; goodK[61] = 5; goodM[61] = 16;
	goodN[62] = 544; goodK[62] = 5; goodM[62] = 17;
	goodN[63] = 576; goodK[63] = 5; goodM[63] = 18;
	goodN[64] = 608; goodK[64] = 5; goodM[64] = 19;
	goodN[65] = 640; goodK[65] = 5; goodM[65] = 20;
	goodN[66] = 704; goodK[66] = 6; goodM[66] = 11;
	goodN[67] = 768; goodK[67] = 6; goodM[67] = 12;
	goodN[68] = 832; goodK[68] = 6; goodM[68] = 13;
	goodN[69] = 896; goodK[69] = 6; goodM[69] = 14;
	goodN[70] = 960; goodK[70] = 6; goodM[70] = 15;
	goodN[71] = 1024; goodK[71] = 6; goodM[71] = 16;
	goodN[72] = 1088; goodK[72] = 6; goodM[72] = 17;
	goodN[73] = 1152; goodK[73] = 6; goodM[73] = 18;
	goodN[74] = 1216; goodK[74] = 6; goodM[74] = 19;
	goodN[75] = 1280; goodK[75] = 7; goodM[75] = 10;
	goodN[76] = 1408; goodK[76] = 7; goodM[76] = 11;
	goodN[77] = 1536; goodK[77] = 7; goodM[77] = 12;
	goodN[78] = 1664; goodK[78] = 7; goodM[78] = 13;
	goodN[79] = 1792; goodK[79] = 7; goodM[79] = 14;
	goodN[80] = 1920; goodK[80] = 7; goodM[80] = 15;
	goodN[81] = 2048; goodK[81] = 7; goodM[81] = 16;
	goodN[82] = 2176; goodK[82] = 7; goodM[82] = 17;
	goodN[83] = 2304; goodK[83] = 7; goodM[83] = 18;
	goodN[84] = 2560; goodK[84] = 8; goodM[84] = 10;
	goodN[85] = 2816; goodK[85] = 8; goodM[85] = 11;
	goodN[86] = 3072; goodK[86] = 8; goodM[86] = 12;
	goodN[87] = 3328; goodK[87] = 8; goodM[87] = 13;
	goodN[88] = 3584; goodK[88] = 8; goodM[88] = 14;
	goodN[89] = 3840; goodK[89] = 8; goodM[89] = 15;
	goodN[90] = 4096; goodK[90] = 8; goodM[90] = 16;
	// 3 primes
	goodN[91] = 4352; goodK[91] = 8; goodM[91] = 17;
	goodN[92] = 4608; goodK[92] = 8; goodM[92] = 18;
	goodN[93] = 4864; goodK[93] = 8; goodM[93] = 19;
	goodN[94] = 5120; goodK[94] = 8; goodM[94] = 20;
	goodN[95] = 5376; goodK[95] = 8; goodM[95] = 21;
	goodN[96] = 5632; goodK[96] = 8; goodM[96] = 22;
	goodN[97] = 5888; goodK[97] = 8; goodM[97] = 23;
	goodN[98] = 6144; goodK[98] = 8; goodM[98] = 24;
	goodN[99] = 6400; goodK[99] = 8; goodM[99] = 25;
	goodN[100] = 6656; goodK[100] = 8; goodM[100] = 26;
	goodN[101] = 6912; goodK[101] = 8; goodM[101] = 27;
	goodN[102] = 7168; goodK[102] = 8; goodM[102] = 28;
	goodN[103] = 7424; goodK[103] = 8; goodM[103] = 29;
	goodN[104] = 7680; goodK[104] = 9; goodM[104] = 15;
	goodN[105] = 8192; goodK[105] = 9; goodM[105] = 16;
	goodN[106] = 8704; goodK[106] = 9; goodM[106] = 17;
	goodN[107] = 9216; goodK[107] = 9; goodM[107] = 18;
	goodN[108] = 9728; goodK[108] = 9; goodM[108] = 19;
	goodN[109] = 10240; goodK[109] = 9; goodM[109] = 20;
	goodN[110] = 10752; goodK[110] = 9; goodM[110] = 21;
	goodN[111] = 11264; goodK[111] = 9; goodM[111] = 22;
	goodN[112] = 11776; goodK[112] = 9; goodM[112] = 23;
	goodN[113] = 12288; goodK[113] = 9; goodM[113] = 24;
	goodN[114] = 12800; goodK[114] = 9; goodM[114] = 25;
	goodN[115] = 13312; goodK[115] = 9; goodM[115] = 26;
	goodN[116] = 13824; goodK[116] = 9; goodM[116] = 27;
	goodN[117] = 14336; goodK[117] = 9; goodM[117] = 28;
	goodN[118] = 14848; goodK[118] = 9; goodM[118] = 29;
	goodN[119] = 15360; goodK[119] = 10; goodM[119] = 15;
	goodN[120] = 16384; goodK[120] = 10; goodM[120] = 16;
	goodN[121] = 17408; goodK[121] = 10; goodM[121] = 17;
	goodN[122] = 18432; goodK[122] = 10; goodM[122] = 18;
	goodN[123] = 19456; goodK[123] = 10; goodM[123] = 19;
	goodN[124] = 20480; goodK[124] = 10; goodM[124] = 20;
	goodN[125] = 21504; goodK[125] = 10; goodM[125] = 21;
	goodN[126] = 22528; goodK[126] = 10; goodM[126] = 22;
	goodN[127] = 23552; goodK[127] = 10; goodM[127] = 23;
	goodN[128] = 24576; goodK[128] = 10; goodM[128] = 24;
	goodN[129] = 25600; goodK[129] = 10; goodM[129] = 25;
	goodN[130] = 26624; goodK[130] = 10; goodM[130] = 26;
	goodN[131] = 27648; goodK[131] = 10; goodM[131] = 27;
	goodN[132] = 28672; goodK[132] = 10; goodM[132] = 28;
	goodN[133] = 29696; goodK[133] = 10; goodM[133] = 29;
	goodN[134] = 30720; goodK[134] = 11; goodM[134] = 15;
	goodN[135] = 32768; goodK[135] = 11; goodM[135] = 16;
	goodN[136] = 34816; goodK[136] = 11; goodM[136] = 17;
	goodN[137] = 36864; goodK[137] = 11; goodM[137] = 18;
	goodN[138] = 38912; goodK[138] = 11; goodM[138] = 19;
	goodN[139] = 40960; goodK[139] = 11; goodM[139] = 20;
	goodN[140] = 43008; goodK[140] = 11; goodM[140] = 21;
	goodN[141] = 45056; goodK[141] = 11; goodM[141] = 22;
	goodN[142] = 47104; goodK[142] = 11; goodM[142] = 23;
	goodN[143] = 49152; goodK[143] = 11; goodM[143] = 24;
	goodN[144] = 51200; goodK[144] = 11; goodM[144] = 25;
	goodN[145] = 53248; goodK[145] = 11; goodM[145] = 26;
	goodN[146] = 55296; goodK[146] = 11; goodM[146] = 27;
	goodN[147] = 57344; goodK[147] = 11; goodM[147] = 28;
	goodN[148] = 59392; goodK[148] = 11; goodM[148] = 29;
	goodN[149] = 61440; goodK[149] = 12; goodM[149] = 15;
	goodN[150] = 65536; goodK[150] = 12; goodM[150] = 16;
	goodN[151] = 69632; goodK[151] = 12; goodM[151] = 17;
	goodN[152] = 73728; goodK[152] = 12; goodM[152] = 18;
	goodN[153] = 77824; goodK[153] = 12; goodM[153] = 19;
	goodN[154] = 81920; goodK[154] = 12; goodM[154] = 20;
	goodN[155] = 86016; goodK[155] = 12; goodM[155] = 21;
	goodN[156] = 90112; goodK[156] = 12; goodM[156] = 22;
	goodN[157] = 94208; goodK[157] = 12; goodM[157] = 23;
	goodN[158] = 98304; goodK[158] = 12; goodM[158] = 24;
	goodN[159] = 102400; goodK[159] = 12; goodM[159] = 25;
	goodN[160] = 106496; goodK[160] = 12; goodM[160] = 26;
	goodN[161] = 110592; goodK[161] = 12; goodM[161] = 27;
	goodN[162] = 114688; goodK[162] = 12; goodM[162] = 28;
	goodN[163] = 122880; goodK[163] = 13; goodM[163] = 15;
	goodN[164] = 131072; goodK[164] = 13; goodM[164] = 16;
	goodN[165] = 139264; goodK[165] = 13; goodM[165] = 17;
	goodN[166] = 147456; goodK[166] = 13; goodM[166] = 18;
	goodN[167] = 155648; goodK[167] = 13; goodM[167] = 19;
	goodN[168] = 163840; goodK[168] = 13; goodM[168] = 20;
	goodN[169] = 172032; goodK[169] = 13; goodM[169] = 21;
	goodN[170] = 180224; goodK[170] = 13; goodM[170] = 22;
	goodN[171] = 188416; goodK[171] = 13; goodM[171] = 23;
	goodN[172] = 196608; goodK[172] = 13; goodM[172] = 24;
	goodN[173] = 204800; goodK[173] = 13; goodM[173] = 25;
	goodN[174] = 212992; goodK[174] = 13; goodM[174] = 26;
	goodN[175] = 221184; goodK[175] = 13; goodM[175] = 27;
	goodN[176] = 229376; goodK[176] = 14; goodM[176] = 14;
	goodN[177] = 245760; goodK[177] = 14; goodM[177] = 15;
	goodN[178] = 262144; goodK[178] = 14; goodM[178] = 16;
  }

  void init2 () {
	size = 136;
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
	goodN[24] = 52; goodK[24] = 2; goodM[24] = 13;
	goodN[25] = 56; goodK[25] = 2; goodM[25] = 14;
	goodN[26] = 60; goodK[26] = 2; goodM[26] = 15;
	goodN[27] = 64; goodK[27] = 2; goodM[27] = 16;
	goodN[28] = 68; goodK[28] = 2; goodM[28] = 17;
	goodN[29] = 72; goodK[29] = 2; goodM[29] = 18;
	goodN[30] = 76; goodK[30] = 2; goodM[30] = 19;
	goodN[31] = 80; goodK[31] = 2; goodM[31] = 20;
	goodN[32] = 84; goodK[32] = 2; goodM[32] = 21;
	goodN[33] = 88; goodK[33] = 2; goodM[33] = 22;
	goodN[34] = 92; goodK[34] = 2; goodM[34] = 23;
	goodN[35] = 96; goodK[35] = 3; goodM[35] = 12;
	goodN[36] = 104; goodK[36] = 3; goodM[36] = 13;
	goodN[37] = 112; goodK[37] = 3; goodM[37] = 14;
	goodN[38] = 120; goodK[38] = 3; goodM[38] = 15;
	goodN[39] = 128; goodK[39] = 3; goodM[39] = 16;
	goodN[40] = 136; goodK[40] = 3; goodM[40] = 17;
	goodN[41] = 144; goodK[41] = 3; goodM[41] = 18;
	goodN[42] = 152; goodK[42] = 3; goodM[42] = 19;
	goodN[43] = 160; goodK[43] = 3; goodM[43] = 20;
	goodN[44] = 168; goodK[44] = 3; goodM[44] = 21;
	goodN[45] = 176; goodK[45] = 3; goodM[45] = 22;
	goodN[46] = 192; goodK[46] = 4; goodM[46] = 12;
	goodN[47] = 208; goodK[47] = 4; goodM[47] = 13;
	goodN[48] = 224; goodK[48] = 4; goodM[48] = 14;
	goodN[49] = 240; goodK[49] = 4; goodM[49] = 15;
	goodN[50] = 256; goodK[50] = 4; goodM[50] = 16;
	goodN[51] = 272; goodK[51] = 4; goodM[51] = 17;
	goodN[52] = 288; goodK[52] = 4; goodM[52] = 18;
	goodN[53] = 304; goodK[53] = 4; goodM[53] = 19;
	goodN[54] = 320; goodK[54] = 4; goodM[54] = 20;
	goodN[55] = 336; goodK[55] = 4; goodM[55] = 21;
	goodN[56] = 352; goodK[56] = 5; goodM[56] = 11;
	goodN[57] = 384; goodK[57] = 5; goodM[57] = 12;
	goodN[58] = 416; goodK[58] = 5; goodM[58] = 13;
	goodN[59] = 448; goodK[59] = 5; goodM[59] = 14;
	goodN[60] = 480; goodK[60] = 5; goodM[60] = 15;
	goodN[61] = 512; goodK[61] = 5; goodM[61] = 16;
	goodN[62] = 544; goodK[62] = 5; goodM[62] = 17;
	goodN[63] = 576; goodK[63] = 5; goodM[63] = 18;
	goodN[64] = 608; goodK[64] = 5; goodM[64] = 19;
	goodN[65] = 640; goodK[65] = 5; goodM[65] = 20;
	goodN[66] = 704; goodK[66] = 6; goodM[66] = 11;
	goodN[67] = 768; goodK[67] = 6; goodM[67] = 12;
	goodN[68] = 832; goodK[68] = 6; goodM[68] = 13;
	goodN[69] = 896; goodK[69] = 6; goodM[69] = 14;
	goodN[70] = 960; goodK[70] = 6; goodM[70] = 15;
	goodN[71] = 1024; goodK[71] = 6; goodM[71] = 16;
	goodN[72] = 1088; goodK[72] = 6; goodM[72] = 17;
	goodN[73] = 1152; goodK[73] = 6; goodM[73] = 18;
	goodN[74] = 1216; goodK[74] = 6; goodM[74] = 19;
	goodN[75] = 1280; goodK[75] = 7; goodM[75] = 10;
	goodN[76] = 1408; goodK[76] = 7; goodM[76] = 11;
	goodN[77] = 1536; goodK[77] = 7; goodM[77] = 12;
	goodN[78] = 1664; goodK[78] = 7; goodM[78] = 13;
	goodN[79] = 1792; goodK[79] = 7; goodM[79] = 14;
	goodN[80] = 1920; goodK[80] = 7; goodM[80] = 15;
	goodN[81] = 2048; goodK[81] = 7; goodM[81] = 16;
	goodN[82] = 2176; goodK[82] = 7; goodM[82] = 17;
	goodN[83] = 2304; goodK[83] = 7; goodM[83] = 18;
	goodN[84] = 2560; goodK[84] = 8; goodM[84] = 10;
	goodN[85] = 2816; goodK[85] = 8; goodM[85] = 11;
	goodN[86] = 3072; goodK[86] = 8; goodM[86] = 12;
	goodN[87] = 3328; goodK[87] = 8; goodM[87] = 13;
	goodN[88] = 3584; goodK[88] = 8; goodM[88] = 14;
	goodN[89] = 3840; goodK[89] = 8; goodM[89] = 15;
	goodN[90] = 4096; goodK[90] = 8; goodM[90] = 16;
	goodN[91] = 4352; goodK[91] = 8; goodM[91] = 17;
	goodN[92] = 4608; goodK[92] = 9; goodM[92] = 9;
	goodN[93] = 5120; goodK[93] = 9; goodM[93] = 10;
	goodN[94] = 5632; goodK[94] = 9; goodM[94] = 11;
	goodN[95] = 6144; goodK[95] = 9; goodM[95] = 12;
	goodN[96] = 6656; goodK[96] = 9; goodM[96] = 13;
	goodN[97] = 7168; goodK[97] = 9; goodM[97] = 14;
	goodN[98] = 7680; goodK[98] = 9; goodM[98] = 15;
	goodN[99] = 8192; goodK[99] = 9; goodM[99] = 16;
	goodN[100] = 9216; goodK[100] = 10; goodM[100] = 9;
	goodN[101] = 10240; goodK[101] = 10; goodM[101] = 10;
	goodN[102] = 11264; goodK[102] = 10; goodM[102] = 11;
	goodN[103] = 12288; goodK[103] = 10; goodM[103] = 12;
	goodN[104] = 13312; goodK[104] = 10; goodM[104] = 13;
	goodN[105] = 14336; goodK[105] = 10; goodM[105] = 14;
	goodN[106] = 15360; goodK[106] = 10; goodM[106] = 15;
	goodN[107] = 16384; goodK[107] = 10; goodM[107] = 16;
	goodN[108] = 18432; goodK[108] = 11; goodM[108] = 9;
	goodN[109] = 20480; goodK[109] = 11; goodM[109] = 10;
	goodN[110] = 22528; goodK[110] = 11; goodM[110] = 11;
	goodN[111] = 24576; goodK[111] = 11; goodM[111] = 12;
	goodN[112] = 26624; goodK[112] = 11; goodM[112] = 13;
	goodN[113] = 28672; goodK[113] = 11; goodM[113] = 14;
	goodN[114] = 30720; goodK[114] = 11; goodM[114] = 15;
	goodN[115] = 32768; goodK[115] = 12; goodM[115] = 8;
	goodN[116] = 36864; goodK[116] = 12; goodM[116] = 9;
	goodN[117] = 40960; goodK[117] = 12; goodM[117] = 10;
	goodN[118] = 45056; goodK[118] = 12; goodM[118] = 11;
	goodN[119] = 49152; goodK[119] = 12; goodM[119] = 12;
	goodN[120] = 53248; goodK[120] = 12; goodM[120] = 13;
	goodN[121] = 57344; goodK[121] = 12; goodM[121] = 14;
	goodN[122] = 65536; goodK[122] = 13; goodM[122] = 8;
	goodN[123] = 73728; goodK[123] = 13; goodM[123] = 9;
	goodN[124] = 81920; goodK[124] = 13; goodM[124] = 10;
	goodN[125] = 90112; goodK[125] = 13; goodM[125] = 11;
	goodN[126] = 98304; goodK[126] = 13; goodM[126] = 12;
	goodN[127] = 106496; goodK[127] = 13; goodM[127] = 13;
	goodN[128] = 114688; goodK[128] = 14; goodM[128] = 7;
	goodN[129] = 131072; goodK[129] = 14; goodM[129] = 8;
	goodN[130] = 147456; goodK[130] = 14; goodM[130] = 9;
	goodN[131] = 163840; goodK[131] = 14; goodM[131] = 10;
	goodN[132] = 180224; goodK[132] = 14; goodM[132] = 11;
	goodN[133] = 196608; goodK[133] = 14; goodM[133] = 12;
	goodN[134] = 229376; goodK[134] = 15; goodM[134] = 7;
	goodN[135] = 262144; goodK[135] = 15; goodM[135] = 8;
  }

  void init3 () {
	size = 209;
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
	goodN[29] = 60; goodK[29] = 2; goodM[29] = 15;
	goodN[30] = 64; goodK[30] = 2; goodM[30] = 16;
	goodN[31] = 68; goodK[31] = 2; goodM[31] = 17;
	goodN[32] = 72; goodK[32] = 2; goodM[32] = 18;
	goodN[33] = 76; goodK[33] = 2; goodM[33] = 19;
	goodN[34] = 80; goodK[34] = 2; goodM[34] = 20;
	goodN[35] = 84; goodK[35] = 2; goodM[35] = 21;
	goodN[36] = 88; goodK[36] = 2; goodM[36] = 22;
	goodN[37] = 92; goodK[37] = 2; goodM[37] = 23;
	goodN[38] = 96; goodK[38] = 2; goodM[38] = 24;
	goodN[39] = 100; goodK[39] = 2; goodM[39] = 25;
	goodN[40] = 104; goodK[40] = 2; goodM[40] = 26;
	goodN[41] = 108; goodK[41] = 2; goodM[41] = 27;
	goodN[42] = 112; goodK[42] = 2; goodM[42] = 28;
	goodN[43] = 116; goodK[43] = 2; goodM[43] = 29;
	goodN[44] = 120; goodK[44] = 3; goodM[44] = 15;
	goodN[45] = 128; goodK[45] = 3; goodM[45] = 16;
	goodN[46] = 136; goodK[46] = 3; goodM[46] = 17;
	goodN[47] = 144; goodK[47] = 3; goodM[47] = 18;
	goodN[48] = 152; goodK[48] = 3; goodM[48] = 19;
	goodN[49] = 160; goodK[49] = 3; goodM[49] = 20;
	goodN[50] = 168; goodK[50] = 3; goodM[50] = 21;
	goodN[51] = 176; goodK[51] = 3; goodM[51] = 22;
	goodN[52] = 184; goodK[52] = 3; goodM[52] = 23;
	goodN[53] = 192; goodK[53] = 3; goodM[53] = 24;
	goodN[54] = 200; goodK[54] = 3; goodM[54] = 25;
	goodN[55] = 208; goodK[55] = 3; goodM[55] = 26;
	goodN[56] = 216; goodK[56] = 3; goodM[56] = 27;
	goodN[57] = 224; goodK[57] = 3; goodM[57] = 28;
	goodN[58] = 232; goodK[58] = 3; goodM[58] = 29;
	goodN[59] = 240; goodK[59] = 4; goodM[59] = 15;
	goodN[60] = 256; goodK[60] = 4; goodM[60] = 16;
	goodN[61] = 272; goodK[61] = 4; goodM[61] = 17;
	goodN[62] = 288; goodK[62] = 4; goodM[62] = 18;
	goodN[63] = 304; goodK[63] = 4; goodM[63] = 19;
	goodN[64] = 320; goodK[64] = 4; goodM[64] = 20;
	goodN[65] = 336; goodK[65] = 4; goodM[65] = 21;
	goodN[66] = 352; goodK[66] = 4; goodM[66] = 22;
	goodN[67] = 368; goodK[67] = 4; goodM[67] = 23;
	goodN[68] = 384; goodK[68] = 4; goodM[68] = 24;
	goodN[69] = 400; goodK[69] = 4; goodM[69] = 25;
	goodN[70] = 416; goodK[70] = 4; goodM[70] = 26;
	goodN[71] = 432; goodK[71] = 4; goodM[71] = 27;
	goodN[72] = 448; goodK[72] = 4; goodM[72] = 28;
	goodN[73] = 464; goodK[73] = 4; goodM[73] = 29;
	goodN[74] = 480; goodK[74] = 5; goodM[74] = 15;
	goodN[75] = 512; goodK[75] = 5; goodM[75] = 16;
	goodN[76] = 544; goodK[76] = 5; goodM[76] = 17;
	goodN[77] = 576; goodK[77] = 5; goodM[77] = 18;
	goodN[78] = 608; goodK[78] = 5; goodM[78] = 19;
	goodN[79] = 640; goodK[79] = 5; goodM[79] = 20;
	goodN[80] = 672; goodK[80] = 5; goodM[80] = 21;
	goodN[81] = 704; goodK[81] = 5; goodM[81] = 22;
	goodN[82] = 736; goodK[82] = 5; goodM[82] = 23;
	goodN[83] = 768; goodK[83] = 5; goodM[83] = 24;
	goodN[84] = 800; goodK[84] = 5; goodM[84] = 25;
	goodN[85] = 832; goodK[85] = 5; goodM[85] = 26;
	goodN[86] = 864; goodK[86] = 5; goodM[86] = 27;
	goodN[87] = 896; goodK[87] = 5; goodM[87] = 28;
	goodN[88] = 928; goodK[88] = 5; goodM[88] = 29;
	goodN[89] = 960; goodK[89] = 6; goodM[89] = 15;
	goodN[90] = 1024; goodK[90] = 6; goodM[90] = 16;
	goodN[91] = 1088; goodK[91] = 6; goodM[91] = 17;
	goodN[92] = 1152; goodK[92] = 6; goodM[92] = 18;
	goodN[93] = 1216; goodK[93] = 6; goodM[93] = 19;
	goodN[94] = 1280; goodK[94] = 6; goodM[94] = 20;
	goodN[95] = 1344; goodK[95] = 6; goodM[95] = 21;
	goodN[96] = 1408; goodK[96] = 6; goodM[96] = 22;
	goodN[97] = 1472; goodK[97] = 6; goodM[97] = 23;
	goodN[98] = 1536; goodK[98] = 6; goodM[98] = 24;
	goodN[99] = 1600; goodK[99] = 6; goodM[99] = 25;
	goodN[100] = 1664; goodK[100] = 6; goodM[100] = 26;
	goodN[101] = 1728; goodK[101] = 6; goodM[101] = 27;
	goodN[102] = 1792; goodK[102] = 6; goodM[102] = 28;
	goodN[103] = 1856; goodK[103] = 6; goodM[103] = 29;
	goodN[104] = 1920; goodK[104] = 7; goodM[104] = 15;
	goodN[105] = 2048; goodK[105] = 7; goodM[105] = 16;
	goodN[106] = 2176; goodK[106] = 7; goodM[106] = 17;
	goodN[107] = 2304; goodK[107] = 7; goodM[107] = 18;
	goodN[108] = 2432; goodK[108] = 7; goodM[108] = 19;
	goodN[109] = 2560; goodK[109] = 7; goodM[109] = 20;
	goodN[110] = 2688; goodK[110] = 7; goodM[110] = 21;
	goodN[111] = 2816; goodK[111] = 7; goodM[111] = 22;
	goodN[112] = 2944; goodK[112] = 7; goodM[112] = 23;
	goodN[113] = 3072; goodK[113] = 7; goodM[113] = 24;
	goodN[114] = 3200; goodK[114] = 7; goodM[114] = 25;
	goodN[115] = 3328; goodK[115] = 7; goodM[115] = 26;
	goodN[116] = 3456; goodK[116] = 7; goodM[116] = 27;
	goodN[117] = 3584; goodK[117] = 7; goodM[117] = 28;
	goodN[118] = 3712; goodK[118] = 7; goodM[118] = 29;
	goodN[119] = 3840; goodK[119] = 8; goodM[119] = 15;
	goodN[120] = 4096; goodK[120] = 8; goodM[120] = 16;
	goodN[121] = 4352; goodK[121] = 8; goodM[121] = 17;
	goodN[122] = 4608; goodK[122] = 8; goodM[122] = 18;
	goodN[123] = 4864; goodK[123] = 8; goodM[123] = 19;
	goodN[124] = 5120; goodK[124] = 8; goodM[124] = 20;
	goodN[125] = 5376; goodK[125] = 8; goodM[125] = 21;
	goodN[126] = 5632; goodK[126] = 8; goodM[126] = 22;
	goodN[127] = 5888; goodK[127] = 8; goodM[127] = 23;
	goodN[128] = 6144; goodK[128] = 8; goodM[128] = 24;
	goodN[129] = 6400; goodK[129] = 8; goodM[129] = 25;
	goodN[130] = 6656; goodK[130] = 8; goodM[130] = 26;
	goodN[131] = 6912; goodK[131] = 8; goodM[131] = 27;
	goodN[132] = 7168; goodK[132] = 8; goodM[132] = 28;
	goodN[133] = 7424; goodK[133] = 8; goodM[133] = 29;
	goodN[134] = 7680; goodK[134] = 9; goodM[134] = 15;
	goodN[135] = 8192; goodK[135] = 9; goodM[135] = 16;
	goodN[136] = 8704; goodK[136] = 9; goodM[136] = 17;
	goodN[137] = 9216; goodK[137] = 9; goodM[137] = 18;
	goodN[138] = 9728; goodK[138] = 9; goodM[138] = 19;
	goodN[139] = 10240; goodK[139] = 9; goodM[139] = 20;
	goodN[140] = 10752; goodK[140] = 9; goodM[140] = 21;
	goodN[141] = 11264; goodK[141] = 9; goodM[141] = 22;
	goodN[142] = 11776; goodK[142] = 9; goodM[142] = 23;
	goodN[143] = 12288; goodK[143] = 9; goodM[143] = 24;
	goodN[144] = 12800; goodK[144] = 9; goodM[144] = 25;
	goodN[145] = 13312; goodK[145] = 9; goodM[145] = 26;
	goodN[146] = 13824; goodK[146] = 9; goodM[146] = 27;
	goodN[147] = 14336; goodK[147] = 9; goodM[147] = 28;
	goodN[148] = 14848; goodK[148] = 9; goodM[148] = 29;
	goodN[149] = 15360; goodK[149] = 10; goodM[149] = 15;
	goodN[150] = 16384; goodK[150] = 10; goodM[150] = 16;
	goodN[151] = 17408; goodK[151] = 10; goodM[151] = 17;
	goodN[152] = 18432; goodK[152] = 10; goodM[152] = 18;
	goodN[153] = 19456; goodK[153] = 10; goodM[153] = 19;
	goodN[154] = 20480; goodK[154] = 10; goodM[154] = 20;
	goodN[155] = 21504; goodK[155] = 10; goodM[155] = 21;
	goodN[156] = 22528; goodK[156] = 10; goodM[156] = 22;
	goodN[157] = 23552; goodK[157] = 10; goodM[157] = 23;
	goodN[158] = 24576; goodK[158] = 10; goodM[158] = 24;
	goodN[159] = 25600; goodK[159] = 10; goodM[159] = 25;
	goodN[160] = 26624; goodK[160] = 10; goodM[160] = 26;
	goodN[161] = 27648; goodK[161] = 10; goodM[161] = 27;
	goodN[162] = 28672; goodK[162] = 10; goodM[162] = 28;
	goodN[163] = 29696; goodK[163] = 10; goodM[163] = 29;
	goodN[164] = 30720; goodK[164] = 11; goodM[164] = 15;
	goodN[165] = 32768; goodK[165] = 11; goodM[165] = 16;
	goodN[166] = 34816; goodK[166] = 11; goodM[166] = 17;
	goodN[167] = 36864; goodK[167] = 11; goodM[167] = 18;
	goodN[168] = 38912; goodK[168] = 11; goodM[168] = 19;
	goodN[169] = 40960; goodK[169] = 11; goodM[169] = 20;
	goodN[170] = 43008; goodK[170] = 11; goodM[170] = 21;
	goodN[171] = 45056; goodK[171] = 11; goodM[171] = 22;
	goodN[172] = 47104; goodK[172] = 11; goodM[172] = 23;
	goodN[173] = 49152; goodK[173] = 11; goodM[173] = 24;
	goodN[174] = 51200; goodK[174] = 11; goodM[174] = 25;
	goodN[175] = 53248; goodK[175] = 11; goodM[175] = 26;
	goodN[176] = 55296; goodK[176] = 11; goodM[176] = 27;
	goodN[177] = 57344; goodK[177] = 11; goodM[177] = 28;
	goodN[178] = 59392; goodK[178] = 11; goodM[178] = 29;
	goodN[179] = 61440; goodK[179] = 12; goodM[179] = 15;
	goodN[180] = 65536; goodK[180] = 12; goodM[180] = 16;
	goodN[181] = 69632; goodK[181] = 12; goodM[181] = 17;
	goodN[182] = 73728; goodK[182] = 12; goodM[182] = 18;
	goodN[183] = 77824; goodK[183] = 12; goodM[183] = 19;
	goodN[184] = 81920; goodK[184] = 12; goodM[184] = 20;
	goodN[185] = 86016; goodK[185] = 12; goodM[185] = 21;
	goodN[186] = 90112; goodK[186] = 12; goodM[186] = 22;
	goodN[187] = 94208; goodK[187] = 12; goodM[187] = 23;
	goodN[188] = 98304; goodK[188] = 12; goodM[188] = 24;
	goodN[189] = 102400; goodK[189] = 12; goodM[189] = 25;
	goodN[190] = 106496; goodK[190] = 12; goodM[190] = 26;
	goodN[191] = 110592; goodK[191] = 12; goodM[191] = 27;
	goodN[192] = 114688; goodK[192] = 12; goodM[192] = 28;
	goodN[193] = 122880; goodK[193] = 13; goodM[193] = 15;
	goodN[194] = 131072; goodK[194] = 13; goodM[194] = 16;
	goodN[195] = 139264; goodK[195] = 13; goodM[195] = 17;
	goodN[196] = 147456; goodK[196] = 13; goodM[196] = 18;
	goodN[197] = 155648; goodK[197] = 13; goodM[197] = 19;
	goodN[198] = 163840; goodK[198] = 13; goodM[198] = 20;
	goodN[199] = 172032; goodK[199] = 13; goodM[199] = 21;
	goodN[200] = 180224; goodK[200] = 13; goodM[200] = 22;
	goodN[201] = 188416; goodK[201] = 13; goodM[201] = 23;
	goodN[202] = 196608; goodK[202] = 13; goodM[202] = 24;
	goodN[203] = 204800; goodK[203] = 13; goodM[203] = 25;
	goodN[204] = 212992; goodK[204] = 13; goodM[204] = 26;
	goodN[205] = 221184; goodK[205] = 13; goodM[205] = 27;
	goodN[206] = 229376; goodK[206] = 14; goodM[206] = 14;
	goodN[207] = 245760; goodK[207] = 14; goodM[207] = 15;
	goodN[208] = 262144; goodK[208] = 14; goodM[208] = 16;
  }

  MulSSA(){
    //precomputed data for CRT
    P1_P2 = 914472618826924033; 
    HALF_P1_P2 = 457236309413462016; 
    N_HALF_P1_P2 = -457236309413462016;
    P1_U1 = 457236236737707957; 
    P2_U2 = -457236236737707956; 

     P2_P3 = 893553310574837761;
     P3_P1 = 905388453748801537;
     U23 = -412536414;
     U31 = 316662352;
     U12 = 89582607;

     P1_P2_P3 = (__int128) P1_P2 * prime3;
     HALF_P1_P2_P3 = P1_P2_P3 >> 1;
     N_HALF_P1_P2_P3 = -HALF_P1_P2_P3;

     // Hard Code the Good N, K, M    
     init();
     //init2(); 
     //init3();
   }
   ~MulSSA() {
	delete [] goodN;
	delete [] goodM;
	delete [] goodK;
   }
  
  /**
   * Multiplying two univariate polynomial
   *
   * @param a first polynomial
   * @param b second polynomial
   * @return result polynomial
   */
  virtual UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
    UnivariateIntegerPolynomial c(RESULT_DEGREE_P(a, b), MAX_DIGITS_P(a, b));
    
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
  virtual void multiply(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){
    mul2C(a, b, c);
  };
  
 private:
  
  /**
   * Multiply two univariate polynomials using the 
   * Schonhage-Strassen Algorithm.
   *
   * @param a the first polynomial
   * @param b the second polynomial
   * @param c the result polynomial
   */
  void mul2C(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c);
  void mul2C2(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c);
  void mul2C3(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c);

  /**
   * for both C^+ (ncc) and C^- (cc), for each coefficient of y^i,
   * which is a large integer encoded as a polynomial in x, i.e.
   * ncc1_i(x), ncc2_i(x), cc1_i(x), cc2_i(x)
   * (1) for each coefficient of x, apply CRT for the two prime numbers, 
   *     get ncc_i and cc_i
   * (2) convert ncc_i and cc_i to GMP, get u_i and v_i
   * (3) compute c_i = (u_i+v_i)/2+(-u_i+v_i)/2*2^N
   * @param ncc1 result from TwoConvolutionMod
   * @param d1 size of y in the first int poly
   * @param d2 size of y in the second int poly
   * @param c coefficient vector of the product
   */
  void CRTtoMPZ(mpz_t z, sfixn* a, sfixn* b, int size);
  void CRTtoMPZ(mpz_t z, sfixn* a, sfixn* b, sfixn* c, int size);
  void recoverMPZ(mpz_t res, mpz_t u, mpz_t v);
  void CRT_ToGMP_Recovering2(sfixn *ncc1, int d1, int d2, mpz_class *c);
  void CRT_ToGMP_Recovering(sfixn *ncc1, int d1, int d2, mpz_class *c);
		
  /**
   * Convert a univariate polynomial with large coefficients to 
   * a bivariate polynomial representation with relatively 
   * small coefficients.
   *
   * univariate --> a1    * y ^ 0 + a2    * y ^ 1 + ... + ad    * y ^ (d-1)
   * bivariate  --> A1(x) * y ^ 0 + A2(x) * y ^ 1 + ... + Ad(x) * y ^ (d-1)
   *
   * Ai(x) = b1 * x ^ 0 + b2 * x ^ 1 + ... + bK * x ^ (K-1)
   *
   * @param coeff the coefficients of large integers
   * @param d partial degree of the univariate polynomial plus one
   */
  BivariatePolynomial * ToBivarTwoMod(mpz_class *coeff, int d, int p1, int p2);
  BivariatePolynomial* toBivariateMod3(mpz_class *coeff, int d, sfixn p1, sfixn p2, sfixn p3);
	
  void mpzToPolyTwoMod(mpz_t z, int M, 
		       sfixn prime1, sfixn prime2,
		       sfixn *X1, sfixn *X2); 
  void mpzToPolyMod3(sfixn* a1, sfixn* a2, sfixn* a3, mpz_t coef, int M, sfixn p1, sfixn p2, sfixn p3);

  BivariatePolynomial * ToBivarTwoMod0(mpz_class *coeff, int d, int p1, int p2);

  /**
   * Find the next good N which is in the form of N = 2*k * M.
   * The global variables of M, K and N will be set based on 
   * the index of the good N.
   *
   * @param n the integer that might not be good (not in the form
   * 		of N = 2*k * M)
   * @param d, the maximum degree among a and b
   */
  void determineGoodN(int n, sfixn d);
  
  /**
   * Search for the nearest good N to the input integer
   * using binary search on the goodN array.
   *
   * @param n the integer
   * @param start index of the array
   * @param end end index of the array
   * @return the index
   */
  int binarySearch(int n, int start, int end);
  	
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


