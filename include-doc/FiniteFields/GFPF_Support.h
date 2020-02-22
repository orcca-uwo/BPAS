#ifndef _GFPF_SUPPORT_H_
#define _GFPF_SUPPORT_H_

// #ifdef __cplusplus
// extern "C" { 
// #endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <unistd.h>

/**************************************/
#ifndef VERIFICATION_ENABLED
#define VERIFICATION_ENABLED 0
#endif
/**************************************/

#ifndef LOOP_UNROLLING_ENABLED
#define LOOP_UNROLLING_ENABLED
#endif
/**************************************/

#define CONVOLUTION_CACHE_SIZE 4

typedef long long int sfixn64;
typedef unsigned long long int usfixn64;
//typedef __int128 longfixnum;

/**************************************/


typedef struct
{
	usfixn64 m1;
	usfixn64 m2;
	usfixn64 p1;
	usfixn64 p2;

	usfixn64 p1_inv_m;
	usfixn64 p1_inv_q;
	usfixn64 p2_inv_m;
	usfixn64 p2_inv_q;
//	__int128 p1p2;
	usfixn64 p1p2_q;
	usfixn64 p1p2_m;

	usfixn64 qb;
	usfixn64 mb;

	usfixn64 radix;
	usfixn64 r_inv_0;
	usfixn64 r_inv_1;

} crt_u192_data;

//#ifndef t_crt_data_global
extern crt_u192_data t_crt_data_global;
//#endif

/**************************************/

typedef struct
{
	sfixn64 radix;
	int k;
} srgfn_prime;

/**************************************/
/**************************************/

/**
 * the following constants are used for montogmery
 * multiplication, addition, sub.
 *
 * R=(2**64)
 * prime = p
 * prime_inv = (-1/p) mod R
 * prime * prime_inv == -1 mod R
 * R1 = (2**128) mod p
 */

#define U64_MASK 18446744073709551615L

#define BASE_63 63

#define CONVOLUTION_PRIME_1 4179340454199820289 //29*2^57 + 1
#define INV_CONVOLUTION_PRIME_1 4179340454199820287   //-p^-1 mod R;
//#define RP1 1729382256910270460  //R = 2^64; RP = R ^ -1 mod p
#define R_11 1878466934230121386 //(2^64)^2 mod MY_PRIME1

#define CONVOLUTION_PRIME_2 2485986994308513793
#define INV_CONVOLUTION_PRIME_2 2485986994308513791
//#define RP2 1044835113549955065
#define R_22 1974795801822054070 // (2^64)^2 mod MY_PRIME2

#define R8  ((1L<<59)+(1L<<16))
#define R16 ((1L<<58)+(1L<<10))
#define R32 ((1L<<56)+(1L<<21)) //prime = R^K + 1;
#define R64 ((1L<<47)+(1L<<32))

#define r_inv_0_8  18446744073642442752L
#define r_inv_1_8  31L
#define r_inv_0_16 18446744073705357312L
#define r_inv_1_16 63L
#define r_inv_0_32 18446743936270598147L
#define r_inv_1_32 255L
#define r_inv_0_64 2251731096305600L
#define r_inv_1_64 131068L

/**************************************/
#define MAX_TMP_ARRAY_SIZE 64

extern sfixn64 global_tmp_a[MAX_TMP_ARRAY_SIZE];


extern sfixn64 global_tmp_b[MAX_TMP_ARRAY_SIZE];


extern sfixn64 global_tmp_c[MAX_TMP_ARRAY_SIZE];


extern sfixn64 global_tmp_dft2[MAX_TMP_ARRAY_SIZE];




// ////conv8_step1
extern const sfixn64 w_conv8_step1;
//3324705732702508476; //8th primitive root of unity WRT prime1
extern const sfixn64 winv_conv8_step1;
//1324460247237108937; //w^-1 mod prime1
//sfixn64 theta = 2840515125011340363; //16th primitive root of unity WRT prime1
//theta^(i) mod p, i=0,...,7

extern const sfixn64 theta_list[8];


//theta^(-i) mod p, i=0,...,7
extern const sfixn64 thetainv_list[8];

extern const sfixn64 ninv_conv8_step1;
//3656922897424842753; // 1/32 mod prime

////conv8_step2
/**************************************/
extern const sfixn64 w_conv8_step2;
//433351031471881650; //8th primitive root of unity WRT prime2
extern const sfixn64 winv_conv8_step2;
//906946876039850912; //w^-1 mod prime2
//sfixn64 theta = 408596721056289947; //16th primitive root of unity WRT prime2
//theta^(i) mod p, i=0,...,7
extern const sfixn64 theta_list_conv8_step1[8];

//theta^(-i) mod p, i=0,...,7
extern const sfixn64 thetainv_list_conv8_step1[8];
extern const sfixn64 ninv_conv8_step2;

//2175238620019949569; // 1/32 mod prime

/**************************************/
////conv16_step1
/**************************************/
extern const sfixn64 w_conv16_step1; 
//1300873456667227704; //8th primitive root of unity WRT prime2
extern const sfixn64 winv_conv16_step1;
//2840515125011340363; //w^-1 mod prime2
//sfixn64 theta = 2691566718071997393; //16th primitive root of unity WRT prime2
//theta^(i) mod p, i=0,...,7
extern const sfixn64 theta_list_conv16_step1[16];

//theta^(-i) mod p, i=0,...,7
extern const sfixn64 thetainv_list_conv16_step1[16];

extern const sfixn64 ninv_conv16_step1;
//3918131675812331521; // 1/32 mod prime

/**************************************/
////conv16_step2
/**************************************/
extern const sfixn64 w_conv16_step2;
//2425116930690660708; //8th primitive root of unity WRT prime2
extern const sfixn64 winv_conv16_step2;
//2077390273252223846; //w^-1 mod prime2
//sfixn64 theta = 1028685717932621537; //16th primitive root of unity WRT prime2
//theta^(i) mod p, i=0,...,7
extern const sfixn64 theta_list_conv16_step2[16];

//theta^(-i) mod p, i=0,...,7
extern const sfixn64 thetainv_list_conv16_step2[16];

extern const sfixn64 ninv_conv16_step2;
//2330612807164231681; // 1/32 mod prime

/**************************************/
////conv32_step1
/**************************************/
extern sfixn64 w_conv32_step1;
//3652694528069969887; //32th primitive root of unity WRT prime1
extern sfixn64 winv_conv32_step1;
//486130157111821104; //w^-1 mod prime1
//sfixn64 theta = 2421024014398278637; //64th primitive root of unity WRT prime1
//theta^(i) mod p, i=0,...,31
extern sfixn64 theta_list_conv32_step1[32];
//theta^(-i) mod p, i=0,...,31

extern sfixn64 thetainv_list_conv32_step1[32];

extern sfixn64 ninv_conv32_step1;
//4048736065006075905; // 1/32 mod prime

/**************************************/
////conv32_step2
/**************************************/
extern sfixn64 w_conv32_step2;
//1028685717932621537; //32th primitive root of unity WRT prime2
extern sfixn64 winv_conv32_step2;
//788787063340104355; //w^-1 mod prime2
//sfixn64 theta = 806891853415459562; //64th primitive root of unity WRT prime1
//theta^(i) mod p, i=0,...,31
extern sfixn64 theta_list_conv32_step2[32];
//theta^(-i) mod p, i=0,...,31
//#ifndef thetainv_list_conv32_step2
extern sfixn64 thetainv_list_conv32_step2[32];
//#endif
extern sfixn64 ninv_conv32_step2;
//2408299900736372737; // 1/32 mod prime
/**************************************/
////conv64_step1
/**************************************/
extern sfixn64 w_conv64_step1;
//163657867345391920; //32th primitive root of unity WRT prime1
extern sfixn64 winv_conv64_step1;
//752526341387853577; //w^-1 mod prime1
//sfixn64 theta = 2535834845514761989; //64th primitive root of unity WRT prime1
//theta^(i) mod p, i=0,...,31
extern sfixn64 theta_list_conv64_step1[64];

//theta^(-i) mod p, i=0,...,31
extern sfixn64 thetainv_list_conv64_step1[64];

extern sfixn64 ninv_conv64_step1;
//4114038259602948097; // 1/32 mod prime

/**************************************/
////conv64_step2
/**************************************/
extern sfixn64 w_conv64_step2;
//2428251324049660672; //32th primitive root of unity WRT prime1
extern sfixn64 winv_conv64_step2;
//537092545712996909; //w^-1 mod prime1
//sfixn64 theta = 106594563439452883; //64th primitive root of unity WRT prime1
//theta^(i) mod p, i=0,...,31
extern sfixn64 theta_list_conv64_step2[64];

//theta^(-i) mod p, i=0,...,31
extern sfixn64 thetainv_list_conv64_step2[64];

extern sfixn64 ninv_conv64_step2;

//pointers for precoputed powers of omega in various steps of 
//convolution (fft-based mult)
extern sfixn64 * vec_pow_omega_conv_step1;

extern sfixn64 * vec_pow_omega_inv_conv_step1;

extern sfixn64 * vec_pow_omega_conv_step2;

extern sfixn64 * vec_pow_omega_inv_conv_step2;

#define MAX_CONVOLUTION_SIZE 64
extern sfixn64 global_x1[MAX_CONVOLUTION_SIZE];

extern sfixn64 global_x2[MAX_CONVOLUTION_SIZE];
extern sfixn64 global_y1[MAX_CONVOLUTION_SIZE];
extern sfixn64 global_y2[MAX_CONVOLUTION_SIZE];
extern usfixn64 global_s0[MAX_CONVOLUTION_SIZE];
extern usfixn64 global_s1[MAX_CONVOLUTION_SIZE];
extern sfixn64 global_l_vec[MAX_CONVOLUTION_SIZE];
extern sfixn64 global_h_vec[MAX_CONVOLUTION_SIZE];
extern sfixn64 global_c_vec[MAX_CONVOLUTION_SIZE];


extern char global_post[MAX_CONVOLUTION_SIZE];

#define LHC_NEGATIVE_SIGN 1

//bigint -> vector in radix-based representation
void mpz_to_radix_based_s64(sfixn64 *vector, const mpz_t bigint, sfixn64 radix,
		int max_input_vector_size);


//convert radix-based representation -> vector
void radix_based_to_mpz_s64(mpz_t & bigint, sfixn64 *vector, sfixn64 radix,
		int input_vector_size);

void  GFPFMultiplication(sfixn64 *x, sfixn64 *y, int k,
		const crt_u192_data & t_crt_data);

void compute_srgfn_p_gmp(mpz_t &p, usfixn64 radix, int coefficient_size);

void init_gfpf_mult_data(crt_u192_data& t_crt_data, usfixn64 r_inv_0,
		usfixn64 r_inv_1, usfixn64 r);

//x=x+y
//addtion for GFPF
void addition_big_elements(sfixn64 * x, sfixn64 *y, const int k,
		const sfixn64 r);

//x=x-y
// subtraction for GFPF
void subtraction_big_elements(sfixn64 *x, sfixn64 *y, const int k,
		const sfixn64 r);

//x=x*(r^s)
void __inline__ mult_pow_R(sfixn64 *x, int s, const int k, const sfixn64 r);

void plain_mult_gmp_big_elements(sfixn64 *x, sfixn64 *y, const int k,
		const sfixn64 r, const mpz_t prime);

void precompute_pow_omega(sfixn64* pow_omega, const sfixn64 omega, int n,
		sfixn64 prime, sfixn64 pP);


// #ifdef __cplusplus
// }
// #endif

#endif/* This file is part of the BPAS library http://www.bpaslib.org

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


