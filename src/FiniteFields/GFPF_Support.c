#include "FiniteFields/GFPF_Support.h"
#include "FiniteFields/gfpf_gmp_fft.h"
#include "FiniteFields/gfpf_gmp_tools.h"
#include "FiniteFields/gfpf_arithmetic.h"
#include "FiniteFields/gfpf_six_step_fft.h"
#include "FiniteFields/small_prime_field_fft.h"
// /*******************************************************************************
// This is the C functions for Generalized Fermat Prime Fields. 
// For now, only the following prime numbers can be used(in the 
// form of p = r^k + 1)
// 	p = (2^59 + 2^16)^8 + 1;
// 	p = (2^58 + 2^10)^16 + 1;
// 	p = (2^56 + 2^21)^32 + 1;
// 	p = (2^47 + 2^32)^64 + 1;

// Before using the function, the following steps need to be done.
// Please specify the k value as needed.


//   int k; 
//   mpz_t p0;
//   mpz_init(p0);
//   sfixn64 R;

//   srgfn_prime p8_k16_59_16;
//   srgfn_prime p16_k32_58_10;
//   srgfn_prime p32_k64_56_21;
//   srgfn_prime p64_k128_47_32;

//   p8_k16_59_16.k = 8;
//   p8_k16_59_16.radix = ((1L << 59) + (1L << 16));

//   p16_k32_58_10.k = 16;
//   p16_k32_58_10.radix = ((1L << 58) + (1L << 10));

//   p32_k64_56_21.k = 32;
//   p32_k64_56_21.radix = ((1L << 56) + (1L << 21));

//   p64_k128_47_32.k = 64;
//   p64_k128_47_32.radix = ((1L << 47) + (1L << 32));

//   switch (k)
//   {
//   case 8:
//     compute_srgfn_p_gmp(p0, p8_k16_59_16.radix, p8_k16_59_16.k);
//     R = p8_k16_59_16.radix;
//     break;

//   case 16:
//     compute_srgfn_p_gmp(p0, p16_k32_58_10.radix, p16_k32_58_10.k);
//     R = p16_k32_58_10.radix;
//     break;

//   case 32:
//     compute_srgfn_p_gmp(p0, p32_k64_56_21.radix, p32_k64_56_21.k);
//     R = p32_k64_56_21.radix;
//     break;

//   case 64:
//     compute_srgfn_p_gmp(p0, p64_k128_47_32.radix, p64_k128_47_32.k);
//     R = p64_k128_47_32.radix;
//     break;

//   }


//   usfixn64 r_inv_0, r_inv_1;

//   switch (k)
//   {
//   case 8:
//     r_inv_0 = r_inv_0_8;
//     r_inv_1 = r_inv_1_8;
//     break;
//   case 16:
//     r_inv_0 = r_inv_0_16;
//     r_inv_1 = r_inv_1_16;
//     break;
//   case 32:
//     r_inv_0 = r_inv_0_32;
//     r_inv_1 = r_inv_1_32;
//     break;

//   case 64:
//     r_inv_0 = r_inv_0_64;
//     r_inv_1 = r_inv_1_64;
//     break;
//   }
//   init_gfpf_mult_data(t_crt_data_global, r_inv_0, r_inv_1, R);

//   vec_pow_omega_conv_step1 = (sfixn64*) malloc(
//       k * sizeof(sfixn64));
//   vec_pow_omega_inv_conv_step1 = (sfixn64*) malloc(
//       k * sizeof(sfixn64));
//   vec_pow_omega_conv_step2 = (sfixn64*) malloc(
//       k * sizeof(sfixn64));
//   vec_pow_omega_inv_conv_step2 = (sfixn64*) malloc(
//       k * sizeof(sfixn64));

//   sfixn64 tmp_w_conv_step1, tmp_winv_conv_step1;
//   sfixn64 tmp_w_conv_step2, tmp_winv_conv_step2;

//   switch (k)
//   {
//   case 8:
//     tmp_w_conv_step1 = w_conv8_step1;
//     tmp_winv_conv_step1 = winv_conv8_step1;
//     tmp_w_conv_step2 = w_conv8_step2;
//     tmp_winv_conv_step2 = winv_conv8_step2;
//     break;

//   case 16:
//     tmp_w_conv_step1 = w_conv16_step1;
//     tmp_winv_conv_step1 = winv_conv16_step1;
//     tmp_w_conv_step2 = w_conv16_step2;
//     tmp_winv_conv_step2 = winv_conv16_step2;
//     break;

//   case 32:
//     tmp_w_conv_step1 = w_conv32_step1;
//     tmp_winv_conv_step1 = winv_conv32_step1;
//     tmp_w_conv_step2 = w_conv32_step2;
//     tmp_winv_conv_step2 = winv_conv32_step2;
//     break;

//   case 64:
//     tmp_w_conv_step1 = w_conv64_step1;
//     tmp_winv_conv_step1 = winv_conv64_step1;
//     tmp_w_conv_step2 = w_conv64_step2;
//     tmp_winv_conv_step2 = winv_conv64_step2;
//     break;

//   }

//   precompute_pow_omega(vec_pow_omega_conv_step1, tmp_w_conv_step1,
//       k, CONVOLUTION_PRIME_1, INV_CONVOLUTION_PRIME_1);

//   precompute_pow_omega(vec_pow_omega_inv_conv_step1, tmp_winv_conv_step1,
//       k, CONVOLUTION_PRIME_1, INV_CONVOLUTION_PRIME_1);

//   precompute_pow_omega(vec_pow_omega_conv_step2, tmp_w_conv_step2,
//       k, CONVOLUTION_PRIME_2, INV_CONVOLUTION_PRIME_2);

//   precompute_pow_omega(vec_pow_omega_inv_conv_step2, tmp_winv_conv_step2,
//       k, CONVOLUTION_PRIME_2, INV_CONVOLUTION_PRIME_2);

// **************************************************************************/

// #include "GFPF_Support.h" 

// //conv8_step1
// const sfixn64 w_conv8_step1 = 457531513967587773;
// //3324705732702508476; //8th primitive root of unity WRT prime1

// const sfixn64 winv_conv8_step1 = 2080441839125523171;
// //1324460247237108937; //w^-1 mod prime1
// //sfixn64 theta = 2840515125011340363; //16th primitive root of unity WRT prime1
// //theta^(i) mod p, i=0,...,7

// const sfixn64 theta_list[8] =
// { 1729382256910270460, 1585842362418304202, 3721808940232232516,
// 		3511859793675373719, 2751416685589087298, 484245702337590844,
// 		2080441839125523171, 409524332690061911 };

// //theta^(-i) mod p, i=0,...,7
// const sfixn64 thetainv_list[8] =
// { 1729382256910270460, 3769816121509758378, 2098898615074297118,
// 		3695094751862229445, 1427923768610732991, 667480660524446570,
// 		457531513967587773, 2593498091781516087 };

// const sfixn64 ninv_conv8_step1 = 2305843009213693952;
// //3656922897424842753; // 1/32 mod prime

// ////conv8_step2
// /**************************************/

// const sfixn64 w_conv8_step2 = 1540210943987252404;
// //433351031471881650; //8th primitive root of unity WRT prime2

// const sfixn64 winv_conv8_step2 = 1088440249747011973;
// //906946876039850912; //w^-1 mod prime2
// //sfixn64 theta = 408596721056289947; //16th primitive root of unity WRT prime2
// //theta^(i) mod p, i=0,...,7
// const sfixn64 theta_list_conv8_step1[8] =
// { 1044835113549955065, 1219491944167886590, 945776050321261389,
// 		975757521232241402, 1371793355634421488, 2355896704721207921,
// 		1088440249747011973, 1244307049958378014 };

// //theta^(-i) mod p, i=0,...,7
// const sfixn64 thetainv_list_conv8_step1[8] =
// { 1044835113549955065, 1241679944350135779, 1397546744561501820,
// 		130090289587305872, 1114193638674092305, 1510229473076272391,
// 		1540210943987252404, 1266495050140627203 };

// const sfixn64 ninv_conv8_step2 = 2305843009213693952;
// //2175238620019949569; // 1/32 mod prime

// /**************************************/
// ////conv16_step1
// /**************************************/
// const sfixn64 w_conv16_step1 = 3769816121509758378;
// //1300873456667227704; //8th primitive root of unity WRT prime2

// const sfixn64 winv_conv16_step1 = 1585842362418304202;
// //2840515125011340363; //w^-1 mod prime2
// //sfixn64 theta = 2691566718071997393; //16th primitive root of unity WRT prime2
// //theta^(i) mod p, i=0,...,7
// const sfixn64 theta_list_conv16_step1[16] =
// { 1729382256910270460, 2819730566040207006, 3511859793675373719,
// 		2362357316524691841, 2080441839125523171, 4021483989832543751,
// 		2593498091781516087, 3151870990171459294, 1427923768610732991,
// 		1457201767058013997, 3769816121509758378, 1792593325179242211,
// 		3721808940232232516, 3528745245532525794, 484245702337590844,
// 		4118802380293107367 };

// //theta^(-i) mod p, i=0,...,7
// const sfixn64 thetainv_list_conv16_step1[16] =
// { 1729382256910270460, 60538073906712922, 3695094751862229445,
// 		650595208667294495, 457531513967587773, 2386747129020578078,
// 		409524332690061911, 2722138687141806292, 2751416685589087298,
// 		1027469464028360995, 1585842362418304202, 157856464367276538,
// 		2098898615074297118, 1816983137675128448, 667480660524446570,
// 		1359609888159613283 };

// const sfixn64 ninv_conv16_step1 = 1152921504606846976;
// //3918131675812331521; // 1/32 mod prime

// /**************************************/
// ////conv16_step2
// /**************************************/

// const sfixn64 w_conv16_step2 = 1244307049958378014;
// //2425116930690660708; //8th primitive root of unity WRT prime2

// const sfixn64 winv_conv16_step2 = 1266495050140627203;
// //2077390273252223846; //w^-1 mod prime2
// //sfixn64 theta = 1028685717932621537; //16th primitive root of unity WRT prime2
// //theta^(i) mod p, i=0,...,7
// const sfixn64 theta_list_conv16_step2[16] =
// { 1044835113549955065, 2022571404248977820, 2355896704721207921,
// 		1532822032730435630, 1540210943987252404, 904657620539112591,
// 		1241679944350135779, 172337847095355835, 1371793355634421488,
// 		2121652943084085323, 1266495050140627203, 1356725222271650252,
// 		1397546744561501820, 401594649585253768, 975757521232241402,
// 		1133256234430697433 };

// //theta^(-i) mod p, i=0,...,7
// const sfixn64 thetainv_list_conv16_step2[16] =
// { 1044835113549955065, 1352730759877816360, 1510229473076272391,
// 		2084392344723260025, 1088440249747011973, 1129261772036863541,
// 		1219491944167886590, 364334051224428470, 1114193638674092305,
// 		2313649147213157958, 1244307049958378014, 1581329373769401202,
// 		945776050321261389, 953164961578078163, 130090289587305872,
// 		463415590059535973 };

// const sfixn64 ninv_conv16_step2 = 1152921504606846976;
// //2330612807164231681; // 1/32 mod prime

// /**************************************/
// ////conv32_step1
// /**************************************/

// sfixn64 w_conv32_step1 = 1027469464028360995;
// //3652694528069969887; //32th primitive root of unity WRT prime1

// sfixn64 winv_conv32_step1 = 1457201767058013997;
// //486130157111821104; //w^-1 mod prime1
// //sfixn64 theta = 2421024014398278637; //64th primitive root of unity WRT prime1
// //theta^(i) mod p, i=0,...,31
// sfixn64 theta_list_conv32_step1[32] =
// { 1729382256910270460, 2716813393020814529, 157856464367276538,
// 		4089334770727709972, 3769816121509758378, 3163687887632643281,
// 		60538073906712922, 3383072288340994664, 2098898615074297118,
// 		1563935707605160488, 1457201767058013997, 1255185265220861680,
// 		3695094751862229445, 2875458431980415028, 1816983137675128448,
// 		487605760032951608, 1427923768610732991, 3708624070457302916,
// 		650595208667294495, 1117245083902298516, 667480660524446570,
// 		3783947300525059822, 3151870990171459294, 1873774784217575639,
// 		457531513967587773, 3941932662964242749, 1359609888159613283,
// 		1154697530441711263, 2593498091781516087, 1023233648403010707,
// 		2386747129020578078, 3143170229016239819 };

// //theta^(-i) mod p, i=0,...,31
// sfixn64 thetainv_list_conv32_step1[32] =
// { 1729382256910270460, 1036170225183580470, 1792593325179242211,
// 		3156106805796809582, 1585842362418304202, 3024642923758109026,
// 		2819730566040207006, 237407791235577540, 3721808940232232516,
// 		2305565669982244650, 1027469464028360995, 395393153674760467,
// 		3511859793675373719, 3062095370297521773, 3528745245532525794,
// 		470716383742517373, 2751416685589087298, 3691734694166868681,
// 		2362357316524691841, 1303882022219405261, 484245702337590844,
// 		2924155188978958609, 2722138687141806292, 2615404746594659801,
// 		2080441839125523171, 796268165858825625, 4118802380293107367,
// 		1015652566567177008, 409524332690061911, 90005683472110317,
// 		4021483989832543751, 1462527061179005760 };

// sfixn64 ninv_conv32_step1 = 576460752303423488;
// //4048736065006075905; // 1/32 mod prime

// /**************************************/
// ////conv32_step2
// /**************************************/
// sfixn64 w_conv32_step2 = 2022571404248977820;
// //1028685717932621537; //32th primitive root of unity WRT prime2
// sfixn64 winv_conv32_step2 = 1352730759877816360;
// //788787063340104355; //w^-1 mod prime2
// //sfixn64 theta = 806891853415459562; //64th primitive root of unity WRT prime1
// //theta^(i) mod p, i=0,...,31
// sfixn64 theta_list_conv32_step2[32] =
// { 1044835113549955065, 2047050301189612755, 1356725222271650252,
// 		644588190411360590, 1244307049958378014, 298217946831434736,
// 		2022571404248977820, 1310169535590195091, 1397546744561501820,
// 		1685516285716994721, 2313649147213157958, 425239377029723869,
// 		2355896704721207921, 1569584657395041006, 401594649585253768,
// 		319767642474155278, 1114193638674092305, 750036283956181994,
// 		1532822032730435630, 558366276389926692, 975757521232241402,
// 		1897145380361432392, 364334051224428470, 759232621735111511,
// 		1540210943987252404, 1267062247725635992, 1133256234430697433,
// 		2309309006264702774, 1219491944167886590, 1031021012634149078,
// 		904657620539112591, 562437447091789940 };
// //theta^(-i) mod p, i=0,...,31
// sfixn64 thetainv_list_conv32_step2[32] =
// { 1044835113549955065, 1923549547216723853, 1581329373769401202,
// 		1454965981674364715, 1266495050140627203, 176677988043811019,
// 		1352730759877816360, 1218924746582877801, 945776050321261389,
// 		1726754372573402282, 2121652943084085323, 588841613947081401,
// 		1510229473076272391, 1927620717918587101, 953164961578078163,
// 		1735950710352331799, 1371793355634421488, 2166219351834358515,
// 		2084392344723260025, 916402336913472787, 130090289587305872,
// 		2060747617278789924, 172337847095355835, 800470708591519072,
// 		1088440249747011973, 1175817458718318702, 463415590059535973,
// 		2187769047477079057, 1241679944350135779, 1841398803897153203,
// 		1129261772036863541, 438936693118901038 };

// sfixn64 ninv_conv32_step2 = 576460752303423488;
// //2408299900736372737; // 1/32 mod prime
// //    for (int i = 0; i < 32; ++i)
// /**************************************/
// ////conv64_step1
// /**************************************/

// sfixn64 w_conv64_step1 = 3024642923758109026;
// //163657867345391920; //32th primitive root of unity WRT prime1

// sfixn64 winv_conv64_step1 = 3163687887632643281;
// //752526341387853577; //w^-1 mod prime1
// //sfixn64 theta = 2535834845514761989; //64th primitive root of unity WRT prime1
// //theta^(i) mod p, i=0,...,31
// sfixn64 theta_list_conv64_step1[64] =
// { 1729382256910270460, 2065946950361859724, 1154697530441711263,

// 2476108590791941094, 1027469464028360995, 2059359081148152507,

// 3708624070457302916, 538099369645060078, 484245702337590844,

// 917771333061254578, 3383072288340994664, 304891985036234466,

// 4021483989832543751, 3598980956079120329, 3156106805796809582,

// 374626057821032220, 457531513967587773, 1795117943575646548,

// 3062095370297521773, 979325666219533743, 1816983137675128448,

// 3894101812361035661, 2615404746594659801, 3915434929025088050,

// 3769816121509758378, 733637758923106799, 3143170229016239819,

// 442673444517998762, 2819730566040207006, 1735758398812268926,

// 3783947300525059822, 1004802775640331475, 2751416685589087298,

// 2841785773942353161, 1255185265220861680, 253743403371786327,

// 4118802380293107367, 1359278460547709117, 2716813393020814529,

// 2703879931380536893, 2593498091781516087, 3866216525844242598,

// 2305565669982244650, 2102926903966095368, 650595208667294495,

// 1795897946482823557, 1303882022219405261, 266908249491553749,

// 2098898615074297118, 1674470114139062335, 90005683472110317,

// 2471180623826992337, 1792593325179242211, 612119759564674832,

// 3941932662964242749, 4122878476577186715, 3511859793675373719,

// 2933410441910542137, 487605760032951608, 3550631267640149211,

// 2722138687141806292, 1797000015247851304, 3163687887632643281,

// 754332124317356728 };

// //theta^(-i) mod p, i=0,...,31
// sfixn64 thetainv_list_conv64_step1[64] =
// { 1729382256910270460, 3425008329882463561, 1015652566567177008,

// 2382340438951968985, 1457201767058013997, 628709186559671078,

// 3691734694166868681, 1245930012289278152, 667480660524446570,

// 56461977622633574, 237407791235577540, 3567220694635145457,

// 2386747129020578078, 1708159830372827952, 4089334770727709972,

// 2504870340060757954, 2080441839125523171, 3912432204708266540,

// 2875458431980415028, 2383442507716996732, 3528745245532525794,

// 2076413550233724921, 1873774784217575639, 313123928355577691,

// 1585842362418304202, 1475460522819283396, 1462527061179005760,

// 2820061993652111172, 60538073906712922, 3925597050828033962,

// 2924155188978958609, 1337554680257467128, 1427923768610732991,

// 3174537678559488814, 395393153674760467, 2443582055387551363,

// 1359609888159613283, 3736667009681821527, 1036170225183580470,

// 3445702695276713490, 409524332690061911, 263905525174732239,

// 1563935707605160488, 285238641838784628, 2362357316524691841,

// 3200014787980286546, 1117245083902298516, 2384222510624173741,

// 3721808940232232516, 3804714396378788069, 1023233648403010707,

// 580359498120699960, 157856464367276538, 3874448469163585823,

// 796268165858825625, 3261569121138565711, 3695094751862229445,

// 3641241084554760211, 470716383742517373, 2119981373051667782,

// 3151870990171459294, 1703231863407879195, 3024642923758109026,

// 2113393503837960565 };

// sfixn64 ninv_conv64_step1 = 288230376151711744;
// //4114038259602948097; // 1/32 mod prime

// /**************************************/
// ////conv64_step2
// /**************************************/

// sfixn64 w_conv64_step2 = 644588190411360590;
// //2428251324049660672; //32th primitive root of unity WRT prime1

// sfixn64 winv_conv64_step2 = 1454965981674364715;
// //537092545712996909; //w^-1 mod prime1
// //sfixn64 theta = 106594563439452883; //64th primitive root of unity WRT prime1
// //theta^(i) mod p, i=0,...,31
// sfixn64 theta_list_conv64_step2[64] =
// { 1044835113549955065, 2199427594878358967, 2309309006264702774,

// 1265010417846975610, 2121652943084085323, 131538619409147699,

// 750036283956181994, 1879907367050655700, 130090289587305872,

// 1950222215232092632, 1310169535590195091, 1725250145434944198,

// 1129261772036863541, 2158995144080862303, 1454965981674364715,

// 2022703779642475294, 1540210943987252404, 1185256421254406077,

// 1927620717918587101, 2214140757189243146, 401594649585253768,

// 783342093435844151, 800470708591519072, 1162758175855657432,

// 1244307049958378014, 907813655099824896, 562437447091789940,

// 1572759043687496971, 1352730759877816360, 1920694390471811445,

// 1897145380361432392, 2118295947336690355, 1371793355634421488,

// 2164083013020691108, 425239377029723869, 954709305002152387,

// 463415590059535973, 609158136575339544, 2047050301189612755,

// 1421123125906192443, 1219491944167886590, 1613135571073605938,

// 1726754372573402282, 450813148394426698, 1532822032730435630,

// 624848380044809920, 916402336913472787, 836886091165719047,

// 1397546744561501820, 1729036124443350336, 1841398803897153203,

// 714389610640452211, 1581329373769401202, 266557363513473354,

// 1267062247725635992, 1807805997174689797, 1510229473076272391,

// 2436809563503831582, 319767642474155278, 2166131625677576363,

// 172337847095355835, 575908515216647102, 298217946831434736,

// 2309739332721813385 };

// //theta^(-i) mod p, i=0,...,31
// sfixn64 thetainv_list_conv64_step2[64] =
// { 1044835113549955065, 176247661586700408, 2187769047477079057,

// 1910078479091866691, 2313649147213157958, 319855368630937430,

// 2166219351834358515, 49177430804682211, 975757521232241402,

// 678180997133823996, 1218924746582877801, 2219429630795040439,

// 904657620539112591, 1771597383668061582, 644588190411360590,

// 756950869865163457, 1088440249747011973, 1649100903142794746,

// 1569584657395041006, 1861138614263703873, 953164961578078163,

// 2035173845914087095, 759232621735111511, 872851423234907855,

// 1266495050140627203, 1064863868402321350, 438936693118901038,

// 1876828857733174249, 2022571404248977820, 1531277689306361406,

// 2060747617278789924, 321903981287822685, 1114193638674092305,

// 367691046971823438, 588841613947081401, 565292603836702348,

// 1133256234430697433, 913227950621016822, 1923549547216723853,

// 1578173339208688897, 1241679944350135779, 1323228818452856361,

// 1685516285716994721, 1702644900872669642, 2084392344723260025,

// 271846237119270647, 558366276389926692, 1300730573054107716,

// 945776050321261389, 463283214666038499, 1031021012634149078,

// 326991850227651490, 1356725222271650252, 760736848873569595,

// 1175817458718318702, 535764779076421161, 2355896704721207921,

// 606079627257858093, 1735950710352331799, 2354448374899366094,

// 364334051224428470, 1220976576461538183, 176677988043811019,

// 286559399430154826 };

// sfixn64 ninv_conv64_step2 = 288230376151711744;

// //pointers for precoputed powers of omega in various steps of 
// //convolution (fft-based mult)
// sfixn64 * vec_pow_omega_conv_step1, *vec_pow_omega_inv_conv_step1,
// 		*vec_pow_omega_conv_step2, *vec_pow_omega_inv_conv_step2;

// sfixn64 global_x1[MAX_CONVOLUTION_SIZE], global_x2[MAX_CONVOLUTION_SIZE],
// 		global_y1[MAX_CONVOLUTION_SIZE], global_y2[MAX_CONVOLUTION_SIZE];
// usfixn64 global_s0[MAX_CONVOLUTION_SIZE], global_s1[MAX_CONVOLUTION_SIZE];
// sfixn64 global_l_vec[MAX_CONVOLUTION_SIZE], global_h_vec[MAX_CONVOLUTION_SIZE],
// 		global_c_vec[MAX_CONVOLUTION_SIZE];

// char global_post[MAX_CONVOLUTION_SIZE];

// crt_u192_data t_crt_data_global;

// inline sfixn64 AddModMont(sfixn64 a, sfixn64 b, const sfixn64 MY_PRIME)
// {
// 	sfixn64 r = a + b;
// 	r -= MY_PRIME;
// 	r += (r >> BASE_63) & MY_PRIME;
// 	return r;
// }

// /**************************************/

// inline sfixn64 SubModMont(const sfixn64 a, const sfixn64 b, const sfixn64 MY_PRIME)
// {
// 	sfixn64 r = a - b;
// 	r += (r >> BASE_63) & MY_PRIME;
// 	return r;
// }

// /**************************************/

// sfixn64 MulModMont(sfixn64 a, sfixn64 b,
// 		sfixn64 MY_PRIME, sfixn64 INV_PRIME)
// {
// 	__asm__ volatile (
// 			"mulq %2\n\t"
// 			"movq %%rax,%%rsi\n\t"
// 			"movq %%rdx,%%rdi\n\t"
// 			"imulq %3,%%rax\n\t"
// 			"mulq %4\n\t"
// 			"add %%rsi,%%rax\n\t"
// 			"adc %%rdi,%%rdx\n\t"
// 			"subq %4,%%rdx\n\t"
// 			"mov %%rdx,%%rax\n\t"
// 			"sar $63,%%rax\n\t"
// 			"andq %4,%%rax\n\t"
// 			"addq %%rax,%%rdx\n\t"
// 			: "=&d" (a)
// 			: "a"(a),"rm"(b),"b"((sfixn64) INV_PRIME),"c"((sfixn64) MY_PRIME)
// 			:"rsi","rdi");
// 	return a;
// }

// /**************************************/


// //a*b mod n;
// inline sfixn64  convertIn(const sfixn64 a, const sfixn64 r, const sfixn64 MY_PRIME,
// 		const sfixn64 INV_PRIME)
// {
// 	return MulModMont(a, r, MY_PRIME, INV_PRIME);
// }

// /**************************************/

// inline sfixn64  convertOut(const sfixn64 r, const sfixn64 MY_PRIME,
// 		const sfixn64 INV_PRIME)
// {
// 	return MulModMont(r, 1, MY_PRIME, INV_PRIME);
// }

// /**************************************/

// inline void  GFPF_DFT_2(sfixn64* a0, sfixn64* a1, const sfixn64 prime)
// {
// 	sfixn64 sum;
// 	sum = AddModMont(*a0, *a1, prime);
// 	*a1 = SubModMont(*a0, *a1, prime);
// 	*a0 = sum;
// }

// /**************************************/

// inline void swap(sfixn64* a, sfixn64* b)
// {
// 	sfixn64 tmp;
// 	tmp = *a;
// 	*a = *b;
// 	*b = tmp;
// }


// /**************************************/

// //pow_omega is already allocated and is of size n*sizeof(sfixn64);
// void precompute_pow_omega(sfixn64* pow_omega, const sfixn64 omega, int n,
// 		sfixn64 prime, sfixn64 pP)
// {

// 	memset(pow_omega, 0x00, n * sizeof(sfixn64));
// 	pow_omega[0] = 1;
// 	pow_omega[1] = omega;
// 	//ToDO: should it include i=1?
// 	//maybe for arithmetic considerations.
// //	for (int i = 1; i < n; i++)
// 	for (int i = 2; i < n; i++)
// 	{
// 		pow_omega[i] = MulModMont(pow_omega[i - 1],
// 				omega, prime, pP);
// //		pow_omega[i] = convertIn(prime, )
// //		R2, CONVOLUTION_PRIME_2,
// //						INV_CONVOLUTION_PRIME_2
// 	}
// }

// /**************************************/

// sfixn64* GFPF_DFT_8(sfixn64* a, const sfixn64 *omegas, const sfixn64 prime,
// 		const sfixn64 pP)
// {


// 	GFPF_DFT_2(&a[0], &a[4], prime); // dft on permutated indexes
// 	GFPF_DFT_2(&a[2], &a[6], prime);
// 	GFPF_DFT_2(&a[1], &a[5], prime);
// 	GFPF_DFT_2(&a[3], &a[7], prime);

// 	a[6] = MulModMont(a[6], omegas[2], prime, pP); //twiddle
// 	a[7] = MulModMont(a[7], omegas[2], prime, pP);

// 	GFPF_DFT_2(&a[0], &a[2], prime); // dft on permutated indexes
// 	GFPF_DFT_2(&a[1], &a[3], prime);
// 	GFPF_DFT_2(&a[4], &a[6], prime);
// 	GFPF_DFT_2(&a[5], &a[7], prime);

// 	a[5] = MulModMont(a[5], omegas[1], prime, pP); // twiddle
// 	a[3] = MulModMont(a[3], omegas[2], prime, pP);
// 	a[7] = MulModMont(a[7], omegas[3], prime, pP);

// 	GFPF_DFT_2(&a[0], &a[1], prime);   // dft on permutated indexes
// 	GFPF_DFT_2(&a[2], &a[3], prime);
// 	GFPF_DFT_2(&a[4], &a[5], prime);
// 	GFPF_DFT_2(&a[6], &a[7], prime);

// 	//ToDO: use swap(,) here.
// 	sfixn64 tmp;   // final permutation
// 	tmp = a[1];
// 	a[1] = a[4];
// 	a[4] = tmp;
// 	tmp = a[3];
// 	a[3] = a[6];
// 	a[6] = tmp;
// 	// return result
// 	return a;
// }

// /**************************************/

// void GFPF_DFT_16(sfixn64* a, const sfixn64 * omega_pow, const sfixn64 prime,
// 		const sfixn64 pP)
// {


// 	GFPF_DFT_2(&a[0], &a[8], prime);
// 	GFPF_DFT_2(&a[1], &a[9], prime);
// 	GFPF_DFT_2(&a[2], &a[10], prime);
// 	GFPF_DFT_2(&a[3], &a[11], prime);
// 	GFPF_DFT_2(&a[4], &a[12], prime);
// 	GFPF_DFT_2(&a[5], &a[13], prime);
// 	GFPF_DFT_2(&a[6], &a[14], prime);
// 	GFPF_DFT_2(&a[7], &a[15], prime);

// 	a[12] = MulModMont(a[12], omega_pow[4], prime, pP); //*omega*omega*omega*omega;
// 	a[14] = MulModMont(a[14], omega_pow[4], prime, pP); //*omega*omega*omega*omega;
// 	a[13] = MulModMont(a[13], omega_pow[4], prime, pP); //*omega*omega*omega*omega;
// 	a[15] = MulModMont(a[15], omega_pow[4], prime, pP); //*omega*omega*omega*omega;

// 	GFPF_DFT_2(&a[0], &a[4], prime);
// 	GFPF_DFT_2(&a[1], &a[5], prime);
// 	GFPF_DFT_2(&a[2], &a[6], prime);
// 	GFPF_DFT_2(&a[3], &a[7], prime);
// 	GFPF_DFT_2(&a[8], &a[12], prime);
// 	GFPF_DFT_2(&a[9], &a[13], prime);
// 	GFPF_DFT_2(&a[10], &a[14], prime);
// 	GFPF_DFT_2(&a[11], &a[15], prime);

// 	a[6] = MulModMont(a[6], omega_pow[4], prime, pP); //*omega*omega*omega*omega;
// 	a[7] = MulModMont(a[7], omega_pow[4], prime, pP); //*omega*omega*omega*omega;

// 	a[10] = MulModMont(a[10], omega_pow[2], prime, pP); //*omega*omega;
// 	a[11] = MulModMont(a[11], omega_pow[2], prime, pP); //*omega*omega;


// 	a[14] = MulModMont(a[14], omega_pow[6], prime, pP); //*omega*omega*omega*omega*omega*omega;
// 	a[15] = MulModMont(a[15], omega_pow[6], prime, pP); //*omega*omega*omega*omega*omega*omega;

// 	GFPF_DFT_2(&a[0], &a[2], prime);
// 	GFPF_DFT_2(&a[1], &a[3], prime);
// 	GFPF_DFT_2(&a[4], &a[6], prime);
// 	GFPF_DFT_2(&a[5], &a[7], prime);
// 	GFPF_DFT_2(&a[8], &a[10], prime);
// 	GFPF_DFT_2(&a[9], &a[11], prime);
// 	GFPF_DFT_2(&a[12], &a[14], prime);
// 	GFPF_DFT_2(&a[13], &a[15], prime);

// 	a[3] = MulModMont(a[3], omega_pow[4], prime, pP); //*omega*omega*omega*omega;
// 	a[5] = MulModMont(a[5], omega_pow[2], prime, pP); //*omega*omega;
// 	a[7] = MulModMont(a[7], omega_pow[6], prime, pP); //*omega*omega*omega*omega*omega*omega;
// 	a[9] = MulModMont(a[9], omega_pow[1], prime, pP);
// 	a[11] = MulModMont(a[11], omega_pow[5], prime, pP); //*omega*omega*omega*omega*omega;
// 	a[13] = MulModMont(a[13], omega_pow[3], prime, pP); //*omega*omega*omega;
// 	a[15] = MulModMont(a[15], omega_pow[7], prime, pP); //*omega*omega*omega*omega*omega*omega*omega;

// 	GFPF_DFT_2(&a[0], &a[1], prime);
// 	GFPF_DFT_2(&a[2], &a[3], prime);
// 	GFPF_DFT_2(&a[4], &a[5], prime);
// 	GFPF_DFT_2(&a[6], &a[7], prime);
// 	GFPF_DFT_2(&a[8], &a[9], prime);
// 	GFPF_DFT_2(&a[10], &a[11], prime);
// 	GFPF_DFT_2(&a[12], &a[13], prime);
// 	GFPF_DFT_2(&a[14], &a[15], prime);

// //	sfixn64 tmp;
// //	tmp = a[1];
// //	a[1] = a[8];
// //	a[8] = tmp;

// 	swap(&a[1],&a[8]);
// //	tmp = a[2];
// //	a[2] = a[4];
// //	a[4] = tmp;
// 	swap(&a[2],&a[4]);

// //	tmp = a[3];
// //	a[3] = a[12];
// //	a[12] = tmp;

// 	swap(&a[3],&a[12]);

// //	tmp = a[5];
// //	a[5] = a[10];
// //	a[10] = tmp;
// 	swap(&a[5],&a[10]);

// //	tmp = a[7];
// //	a[7] = a[14];
// //	a[14] = tmp;
// 	swap(&a[7],&a[14]);
// 	swap(&a[11],&a[13]);

// //	return a;
// }

// /**************************************/

// sfixn64* GFPF_DFT_32(sfixn64* A, const sfixn64 *omega_pow, const sfixn64 prime,
// 		const sfixn64 pP)
// {

// 	GFPF_DFT_2(&A[0], &A[16], prime);
// 	GFPF_DFT_2(&A[1], &A[17], prime);
// 	GFPF_DFT_2(&A[2], &A[18], prime);
// 	GFPF_DFT_2(&A[3], &A[19], prime);
// 	GFPF_DFT_2(&A[4], &A[20], prime);
// 	GFPF_DFT_2(&A[5], &A[21], prime);
// 	GFPF_DFT_2(&A[6], &A[22], prime);
// 	GFPF_DFT_2(&A[7], &A[23], prime);
// 	GFPF_DFT_2(&A[8], &A[24], prime);
// 	GFPF_DFT_2(&A[9], &A[25], prime);
// 	GFPF_DFT_2(&A[10], &A[26], prime);
// 	GFPF_DFT_2(&A[11], &A[27], prime);
// 	GFPF_DFT_2(&A[12], &A[28], prime);
// 	GFPF_DFT_2(&A[13], &A[29], prime);
// 	GFPF_DFT_2(&A[14], &A[30], prime);
// 	GFPF_DFT_2(&A[15], &A[31], prime);

// 	A[24] = MulModMont(A[24], omega_pow[8], prime, pP);
// 	A[28] = MulModMont(A[28], omega_pow[8], prime, pP);
// 	A[26] = MulModMont(A[26], omega_pow[8], prime, pP);
// 	A[30] = MulModMont(A[30], omega_pow[8], prime, pP);
// 	A[25] = MulModMont(A[25], omega_pow[8], prime, pP);
// 	A[29] = MulModMont(A[29], omega_pow[8], prime, pP);
// 	A[27] = MulModMont(A[27], omega_pow[8], prime, pP);
// 	A[31] = MulModMont(A[31], omega_pow[8], prime, pP);

// 	GFPF_DFT_2(&A[0], &A[8], prime);
// 	GFPF_DFT_2(&A[1], &A[9], prime);
// 	GFPF_DFT_2(&A[2], &A[10], prime);
// 	GFPF_DFT_2(&A[3], &A[11], prime);
// 	GFPF_DFT_2(&A[4], &A[12], prime);
// 	GFPF_DFT_2(&A[5], &A[13], prime);
// 	GFPF_DFT_2(&A[6], &A[14], prime);
// 	GFPF_DFT_2(&A[7], &A[15], prime);
// 	GFPF_DFT_2(&A[16], &A[24], prime);
// 	GFPF_DFT_2(&A[17], &A[25], prime);
// 	GFPF_DFT_2(&A[18], &A[26], prime);
// 	GFPF_DFT_2(&A[19], &A[27], prime);
// 	GFPF_DFT_2(&A[20], &A[28], prime);
// 	GFPF_DFT_2(&A[21], &A[29], prime);
// 	GFPF_DFT_2(&A[22], &A[30], prime);
// 	GFPF_DFT_2(&A[23], &A[31], prime);

// 	A[20] = MulModMont(A[20], omega_pow[4], prime, pP);
// 	A[12] = MulModMont(A[12], omega_pow[8], prime, pP);
// 	A[28] = MulModMont(A[28], omega_pow[12], prime, pP);

// 	A[22] = MulModMont(A[22], omega_pow[4], prime, pP);
// 	A[14] = MulModMont(A[14], omega_pow[8], prime, pP);
// 	A[30] = MulModMont(A[30], omega_pow[12], prime, pP);

// 	A[21] = MulModMont(A[21], omega_pow[4], prime, pP);
// 	A[13] = MulModMont(A[13], omega_pow[8], prime, pP);
// 	A[29] = MulModMont(A[29], omega_pow[12], prime, pP);

// 	A[23] = MulModMont(A[23], omega_pow[4], prime, pP);
// 	A[15] = MulModMont(A[15], omega_pow[8], prime, pP);
// 	A[31] = MulModMont(A[31], omega_pow[12], prime, pP);

// 	GFPF_DFT_2(&A[0], &A[4], prime);
// 	GFPF_DFT_2(&A[16], &A[20], prime);
// 	GFPF_DFT_2(&A[8], &A[12], prime);
// 	GFPF_DFT_2(&A[24], &A[28], prime);
// 	GFPF_DFT_2(&A[2], &A[6], prime);
// 	GFPF_DFT_2(&A[18], &A[22], prime);
// 	GFPF_DFT_2(&A[10], &A[14], prime);
// 	GFPF_DFT_2(&A[26], &A[30], prime);
// 	GFPF_DFT_2(&A[1], &A[5], prime);
// 	GFPF_DFT_2(&A[17], &A[21], prime);
// 	GFPF_DFT_2(&A[9], &A[13], prime);
// 	GFPF_DFT_2(&A[25], &A[29], prime);
// 	GFPF_DFT_2(&A[3], &A[7], prime);
// 	GFPF_DFT_2(&A[19], &A[23], prime);
// 	GFPF_DFT_2(&A[11], &A[15], prime);
// 	GFPF_DFT_2(&A[27], &A[31], prime);

// 	A[18] = MulModMont(A[18], omega_pow[2], prime, pP);
// 	A[10] = MulModMont(A[10], omega_pow[4], prime, pP);
// 	A[26] = MulModMont(A[26], omega_pow[6], prime, pP);
// 	A[6] = MulModMont(A[6], omega_pow[8], prime, pP);
// 	A[22] = MulModMont(A[22], omega_pow[10], prime, pP);
// 	A[14] = MulModMont(A[14], omega_pow[12], prime, pP);
// 	A[30] = MulModMont(A[30], omega_pow[14], prime, pP);

// 	A[19] = MulModMont(A[19], omega_pow[2], prime, pP);
// 	A[11] = MulModMont(A[11], omega_pow[4], prime, pP);
// 	A[27] = MulModMont(A[27], omega_pow[6], prime, pP);
// 	A[7] = MulModMont(A[7], omega_pow[8], prime, pP);
// 	A[23] = MulModMont(A[23], omega_pow[10], prime, pP);
// 	A[15] = MulModMont(A[15], omega_pow[12], prime, pP);
// 	A[31] = MulModMont(A[31], omega_pow[14], prime, pP);

// 	GFPF_DFT_2(&A[0], &A[2], prime);
// 	GFPF_DFT_2(&A[1], &A[3], prime);
// 	GFPF_DFT_2(&A[4], &A[6], prime);
// 	GFPF_DFT_2(&A[5], &A[7], prime);
// 	GFPF_DFT_2(&A[8], &A[10], prime);
// 	GFPF_DFT_2(&A[9], &A[11], prime);
// 	GFPF_DFT_2(&A[12], &A[14], prime);
// 	GFPF_DFT_2(&A[13], &A[15], prime);
// 	GFPF_DFT_2(&A[16], &A[18], prime);
// 	GFPF_DFT_2(&A[17], &A[19], prime);
// 	GFPF_DFT_2(&A[20], &A[22], prime);
// 	GFPF_DFT_2(&A[21], &A[23], prime);
// 	GFPF_DFT_2(&A[24], &A[26], prime);
// 	GFPF_DFT_2(&A[25], &A[27], prime);
// 	GFPF_DFT_2(&A[28], &A[30], prime);
// 	GFPF_DFT_2(&A[29], &A[31], prime);

// 	A[17] = MulModMont(A[17], omega_pow[1], prime, pP); //*omega
// 	A[9] = MulModMont(A[9], omega_pow[2], prime, pP); //*omega*omega;
// 	A[25] = MulModMont(A[25], omega_pow[3], prime, pP); //*omega*omega*omega;
// 	A[5] = MulModMont(A[5], omega_pow[4], prime, pP); //*omega*omega*omega*omega;
// 	A[21] = MulModMont(A[21], omega_pow[5], prime, pP); //*omega*omega*omega*omega*omega;
// 	A[13] = MulModMont(A[13], omega_pow[6], prime, pP); //*omega*omega*omega*omega*omega*omega;
// 	A[29] = MulModMont(A[29], omega_pow[7], prime, pP); //*omega*omega*omega*omega*omega*omega*omega;
// 	A[3] = MulModMont(A[3], omega_pow[8], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega;
// 	A[19] = MulModMont(A[19], omega_pow[9], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega*omega;
// 	A[11] = MulModMont(A[11], omega_pow[10], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
// 	A[27] = MulModMont(A[27], omega_pow[11], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
// 	A[7] = MulModMont(A[7], omega_pow[12], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
// 	A[23] = MulModMont(A[23], omega_pow[13], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
// 	A[15] = MulModMont(A[15], omega_pow[14], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;
// 	A[31] = MulModMont(A[31], omega_pow[15], prime, pP); //*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega*omega;

// 	GFPF_DFT_2(&A[0], &A[1], prime);
// 	GFPF_DFT_2(&A[2], &A[3], prime);
// 	GFPF_DFT_2(&A[4], &A[5], prime);
// 	GFPF_DFT_2(&A[6], &A[7], prime);
// 	GFPF_DFT_2(&A[8], &A[9], prime);
// 	GFPF_DFT_2(&A[10], &A[11], prime);
// 	GFPF_DFT_2(&A[12], &A[13], prime);
// 	GFPF_DFT_2(&A[14], &A[15], prime);
// 	GFPF_DFT_2(&A[16], &A[17], prime);
// 	GFPF_DFT_2(&A[18], &A[19], prime);
// 	GFPF_DFT_2(&A[20], &A[21], prime);
// 	GFPF_DFT_2(&A[22], &A[23], prime);
// 	GFPF_DFT_2(&A[24], &A[25], prime);
// 	GFPF_DFT_2(&A[26], &A[27], prime);
// 	GFPF_DFT_2(&A[28], &A[29], prime);
// 	GFPF_DFT_2(&A[30], &A[31], prime);

// 	sfixn64 tmp;
// 	tmp = A[1];
// 	A[1] = A[16];
// 	A[16] = tmp;
// 	tmp = A[2];
// 	A[2] = A[8];
// 	A[8] = tmp;
// 	tmp = A[3];
// 	A[3] = A[24];
// 	A[24] = tmp;
// 	tmp = A[5];
// 	A[5] = A[20];
// 	A[20] = tmp;
// 	tmp = A[6];
// 	A[6] = A[12];
// 	A[12] = tmp;
// 	tmp = A[7];
// 	A[7] = A[28];
// 	A[28] = tmp;
// 	tmp = A[9];
// 	A[9] = A[18];
// 	A[18] = tmp;
// 	tmp = A[11];
// 	A[11] = A[26];
// 	A[26] = tmp;
// 	tmp = A[13];
// 	A[13] = A[22];
// 	A[22] = tmp;
// 	tmp = A[15];
// 	A[15] = A[30];
// 	A[30] = tmp;
// 	tmp = A[19];
// 	A[19] = A[25];
// 	A[25] = tmp;
// 	tmp = A[23];
// 	A[23] = A[29];
// 	A[29] = tmp;

// 	return A;
// }

// /**************************************/

// sfixn64* GFPF_DFT_64(sfixn64* A, const sfixn64 * omegas, sfixn64 prime, sfixn64 pP)
// {

// 	GFPF_DFT_2(&A[0], &A[32], prime);
// 	GFPF_DFT_2(&A[16], &A[48], prime);
// 	GFPF_DFT_2(&A[8], &A[40], prime);
// 	GFPF_DFT_2(&A[24], &A[56], prime);
// 	GFPF_DFT_2(&A[4], &A[36], prime);
// 	GFPF_DFT_2(&A[20], &A[52], prime);
// 	GFPF_DFT_2(&A[12], &A[44], prime);
// 	GFPF_DFT_2(&A[28], &A[60], prime);
// 	GFPF_DFT_2(&A[2], &A[34], prime);
// 	GFPF_DFT_2(&A[18], &A[50], prime);
// 	GFPF_DFT_2(&A[10], &A[42], prime);
// 	GFPF_DFT_2(&A[26], &A[58], prime);
// 	GFPF_DFT_2(&A[6], &A[38], prime);
// 	GFPF_DFT_2(&A[22], &A[54], prime);
// 	GFPF_DFT_2(&A[14], &A[46], prime);
// 	GFPF_DFT_2(&A[30], &A[62], prime);

// 	GFPF_DFT_2(&A[1], &A[33], prime);
// 	GFPF_DFT_2(&A[17], &A[49], prime);
// 	GFPF_DFT_2(&A[9], &A[41], prime);
// 	GFPF_DFT_2(&A[25], &A[57], prime);
// 	GFPF_DFT_2(&A[5], &A[37], prime);
// 	GFPF_DFT_2(&A[21], &A[53], prime);
// 	GFPF_DFT_2(&A[13], &A[45], prime);
// 	GFPF_DFT_2(&A[29], &A[61], prime);
// 	GFPF_DFT_2(&A[3], &A[35], prime);
// 	GFPF_DFT_2(&A[19], &A[51], prime);
// 	GFPF_DFT_2(&A[11], &A[43], prime);
// 	GFPF_DFT_2(&A[27], &A[59], prime);
// 	GFPF_DFT_2(&A[7], &A[39], prime);
// 	GFPF_DFT_2(&A[23], &A[55], prime);
// 	GFPF_DFT_2(&A[15], &A[47], prime);
// 	GFPF_DFT_2(&A[31], &A[63], prime);

// 	//T_{2}^{4}
// 	A[48] = MulModMont(A[48], omegas[16], prime, pP);
// 	A[56] = MulModMont(A[56], omegas[16], prime, pP);
// 	A[52] = MulModMont(A[52], omegas[16], prime, pP);
// 	A[60] = MulModMont(A[60], omegas[16], prime, pP);
// 	A[50] = MulModMont(A[50], omegas[16], prime, pP);
// 	A[58] = MulModMont(A[58], omegas[16], prime, pP);
// 	A[54] = MulModMont(A[54], omegas[16], prime, pP);
// 	A[62] = MulModMont(A[62], omegas[16], prime, pP);
// 	A[49] = MulModMont(A[49], omegas[16], prime, pP);
// 	A[57] = MulModMont(A[57], omegas[16], prime, pP);
// 	A[53] = MulModMont(A[53], omegas[16], prime, pP);
// 	A[61] = MulModMont(A[61], omegas[16], prime, pP);
// 	A[51] = MulModMont(A[51], omegas[16], prime, pP);
// 	A[59] = MulModMont(A[59], omegas[16], prime, pP);
// 	A[55] = MulModMont(A[55], omegas[16], prime, pP);
// 	A[63] = MulModMont(A[63], omegas[16], prime, pP);

// 	//I_{32} \tx GFPF_DFT_2
// 	GFPF_DFT_2(&A[0], &A[16], prime);
// 	GFPF_DFT_2(&A[32], &A[48], prime);
// 	GFPF_DFT_2(&A[8], &A[24], prime);
// 	GFPF_DFT_2(&A[40], &A[56], prime);
// 	GFPF_DFT_2(&A[4], &A[20], prime);
// 	GFPF_DFT_2(&A[36], &A[52], prime);
// 	GFPF_DFT_2(&A[12], &A[28], prime);
// 	GFPF_DFT_2(&A[44], &A[60], prime);
// 	GFPF_DFT_2(&A[2], &A[18], prime);
// 	GFPF_DFT_2(&A[34], &A[50], prime);
// 	GFPF_DFT_2(&A[10], &A[26], prime);
// 	GFPF_DFT_2(&A[42], &A[58], prime);
// 	GFPF_DFT_2(&A[6], &A[22], prime);
// 	GFPF_DFT_2(&A[38], &A[54], prime);
// 	GFPF_DFT_2(&A[14], &A[30], prime);
// 	GFPF_DFT_2(&A[46], &A[62], prime);
// 	GFPF_DFT_2(&A[1], &A[17], prime);
// 	GFPF_DFT_2(&A[33], &A[49], prime);
// 	GFPF_DFT_2(&A[9], &A[25], prime);
// 	GFPF_DFT_2(&A[41], &A[57], prime);
// 	GFPF_DFT_2(&A[5], &A[21], prime);
// 	GFPF_DFT_2(&A[37], &A[53], prime);
// 	GFPF_DFT_2(&A[13], &A[29], prime);
// 	GFPF_DFT_2(&A[45], &A[61], prime);
// 	GFPF_DFT_2(&A[3], &A[19], prime);
// 	GFPF_DFT_2(&A[35], &A[51], prime);
// 	GFPF_DFT_2(&A[11], &A[27], prime);
// 	GFPF_DFT_2(&A[43], &A[59], prime);
// 	GFPF_DFT_2(&A[7], &A[23], prime);
// 	GFPF_DFT_2(&A[39], &A[55], prime);
// 	GFPF_DFT_2(&A[15], &A[31], prime);
// 	GFPF_DFT_2(&A[47], &A[63], prime);

// 	// I_{8} \tx T_{4}^{8}
// 	// A[0]
// 	// A[32]
// 	// A[16]
// 	// A[48]
// 	// A[8]
// 	A[40] = MulModMont(A[40], omegas[8], prime, pP);
// 	A[24] = MulModMont(A[24], omegas[16], prime, pP);
// 	A[56] = MulModMont(A[56], omegas[24], prime, pP);

// 	// A[4]
// 	// A[36]
// 	// A[20]
// 	// A[52]
// 	// A[12]
// 	A[44] = MulModMont(A[44], omegas[8], prime, pP);
// 	A[28] = MulModMont(A[28], omegas[16], prime, pP);
// 	A[60] = MulModMont(A[60], omegas[24], prime, pP);

// 	// A[2]
// 	// A[34]
// 	// A[18]
// 	// A[50]
// 	// A[10]
// 	A[42] = MulModMont(A[42], omegas[8], prime, pP);
// 	A[26] = MulModMont(A[26], omegas[16], prime, pP);
// 	A[58] = MulModMont(A[58], omegas[24], prime, pP);

// 	// A[6]
// 	// A[38]
// 	// A[22]
// 	// A[54]
// 	// A[14]
// 	A[46] = MulModMont(A[46], omegas[8], prime, pP);
// 	A[30] = MulModMont(A[30], omegas[16], prime, pP);
// 	A[62] = MulModMont(A[62], omegas[24], prime, pP);

// 	// A[1]
// 	// A[33]
// 	// A[17]
// 	// A[49]
// 	// A[9]
// 	A[41] = MulModMont(A[41], omegas[8], prime, pP);
// 	A[25] = MulModMont(A[25], omegas[16], prime, pP);
// 	A[57] = MulModMont(A[57], omegas[24], prime, pP);

// 	// A[5]
// 	// A[37]
// 	// A[21]
// 	// A[53]
// 	// A[13]
// 	A[45] = MulModMont(A[45], omegas[8], prime, pP);
// 	A[29] = MulModMont(A[29], omegas[16], prime, pP);
// 	A[61] = MulModMont(A[61], omegas[24], prime, pP);

// 	// A[3]
// 	// A[35]
// 	// A[19]
// 	// A[51]
// 	// A[11]
// 	A[43] = MulModMont(A[43], omegas[8], prime, pP);
// 	A[27] = MulModMont(A[27], omegas[16], prime, pP);
// 	A[59] = MulModMont(A[59], omegas[24], prime, pP);

// 	// A[7]
// 	// A[39]
// 	// A[23]
// 	// A[55]
// 	// A[15]
// 	A[47] = MulModMont(A[47], omegas[8], prime, pP);
// 	A[31] = MulModMont(A[31], omegas[16], prime, pP);
// 	A[63] = MulModMont(A[63], omegas[24], prime, pP);

// 	// I_{32} \tx GFPF_DFT_2
// 	GFPF_DFT_2(&A[0], &A[8], prime);
// 	GFPF_DFT_2(&A[32], &A[40], prime);
// 	GFPF_DFT_2(&A[16], &A[24], prime);
// 	GFPF_DFT_2(&A[48], &A[56], prime);
// 	GFPF_DFT_2(&A[4], &A[12], prime);
// 	GFPF_DFT_2(&A[36], &A[44], prime);
// 	GFPF_DFT_2(&A[20], &A[28], prime);
// 	GFPF_DFT_2(&A[52], &A[60], prime);

// 	GFPF_DFT_2(&A[2], &A[10], prime);
// 	GFPF_DFT_2(&A[34], &A[42], prime);
// 	GFPF_DFT_2(&A[18], &A[26], prime);
// 	GFPF_DFT_2(&A[50], &A[58], prime);
// 	GFPF_DFT_2(&A[6], &A[14], prime);
// 	GFPF_DFT_2(&A[38], &A[46], prime);
// 	GFPF_DFT_2(&A[22], &A[30], prime);
// 	GFPF_DFT_2(&A[54], &A[62], prime);

// 	GFPF_DFT_2(&A[1], &A[9], prime);
// 	GFPF_DFT_2(&A[33], &A[41], prime);
// 	GFPF_DFT_2(&A[17], &A[25], prime);
// 	GFPF_DFT_2(&A[49], &A[57], prime);
// 	GFPF_DFT_2(&A[5], &A[13], prime);
// 	GFPF_DFT_2(&A[37], &A[45], prime);
// 	GFPF_DFT_2(&A[21], &A[29], prime);
// 	GFPF_DFT_2(&A[53], &A[61], prime);

// 	GFPF_DFT_2(&A[3], &A[11], prime);
// 	GFPF_DFT_2(&A[35], &A[43], prime);
// 	GFPF_DFT_2(&A[19], &A[27], prime);
// 	GFPF_DFT_2(&A[51], &A[59], prime);
// 	GFPF_DFT_2(&A[7], &A[15], prime);
// 	GFPF_DFT_2(&A[39], &A[47], prime);
// 	GFPF_DFT_2(&A[23], &A[31], prime);
// 	GFPF_DFT_2(&A[55], &A[63], prime);

// 	// T_{8}^{16}
// 	// A[0]
// 	// A[32]
// 	// A[16]
// 	// A[48]
// 	// A[8]
// 	// A[40]
// 	// A[24]
// 	// A[56]
// 	//
// 	// A[4]
// 	A[36] = MulModMont(A[36], omegas[4], prime, pP);
// 	A[20] = MulModMont(A[20], omegas[8], prime, pP);
// 	A[52] = MulModMont(A[52], omegas[12], prime, pP);
// 	A[12] = MulModMont(A[12], omegas[16], prime, pP);
// 	A[44] = MulModMont(A[44], omegas[20], prime, pP);
// 	A[28] = MulModMont(A[28], omegas[24], prime, pP);
// 	A[60] = MulModMont(A[60], omegas[28], prime, pP);

// 	A[38] = MulModMont(A[38], omegas[4], prime, pP);
// 	A[22] = MulModMont(A[22], omegas[8], prime, pP);
// 	A[54] = MulModMont(A[54], omegas[12], prime, pP);
// 	A[14] = MulModMont(A[14], omegas[16], prime, pP);
// 	A[46] = MulModMont(A[46], omegas[20], prime, pP);
// 	A[30] = MulModMont(A[30], omegas[24], prime, pP);
// 	A[62] = MulModMont(A[62], omegas[28], prime, pP);

// 	A[37] = MulModMont(A[37], omegas[4], prime, pP);
// 	A[21] = MulModMont(A[21], omegas[8], prime, pP);
// 	A[53] = MulModMont(A[53], omegas[12], prime, pP);
// 	A[13] = MulModMont(A[13], omegas[16], prime, pP);
// 	A[45] = MulModMont(A[45], omegas[20], prime, pP);
// 	A[29] = MulModMont(A[29], omegas[24], prime, pP);
// 	A[61] = MulModMont(A[61], omegas[28], prime, pP);

// 	A[39] = MulModMont(A[39], omegas[4], prime, pP);
// 	A[23] = MulModMont(A[23], omegas[8], prime, pP);
// 	A[55] = MulModMont(A[55], omegas[12], prime, pP);
// 	A[15] = MulModMont(A[15], omegas[16], prime, pP);
// 	A[47] = MulModMont(A[47], omegas[20], prime, pP);
// 	A[31] = MulModMont(A[31], omegas[24], prime, pP);
// 	A[63] = MulModMont(A[63], omegas[28], prime, pP);

// 	// I_{32} \ts GFPF_DFT_2
// 	GFPF_DFT_2(&A[0], &A[4], prime);
// 	GFPF_DFT_2(&A[32], &A[36], prime);
// 	GFPF_DFT_2(&A[16], &A[20], prime);
// 	GFPF_DFT_2(&A[48], &A[52], prime);
// 	GFPF_DFT_2(&A[8], &A[12], prime);
// 	GFPF_DFT_2(&A[40], &A[44], prime);
// 	GFPF_DFT_2(&A[24], &A[28], prime);
// 	GFPF_DFT_2(&A[56], &A[60], prime);

// 	GFPF_DFT_2(&A[2], &A[6], prime);
// 	GFPF_DFT_2(&A[34], &A[38], prime);
// 	GFPF_DFT_2(&A[18], &A[22], prime);
// 	GFPF_DFT_2(&A[50], &A[54], prime);
// 	GFPF_DFT_2(&A[10], &A[14], prime);
// 	GFPF_DFT_2(&A[42], &A[46], prime);
// 	GFPF_DFT_2(&A[26], &A[30], prime);
// 	GFPF_DFT_2(&A[58], &A[62], prime);

// 	GFPF_DFT_2(&A[1], &A[5], prime);
// 	GFPF_DFT_2(&A[33], &A[37], prime);
// 	GFPF_DFT_2(&A[17], &A[21], prime);
// 	GFPF_DFT_2(&A[49], &A[53], prime);
// 	GFPF_DFT_2(&A[9], &A[13], prime);
// 	GFPF_DFT_2(&A[41], &A[45], prime);
// 	GFPF_DFT_2(&A[25], &A[29], prime);
// 	GFPF_DFT_2(&A[57], &A[61], prime);

// 	GFPF_DFT_2(&A[3], &A[7], prime);
// 	GFPF_DFT_2(&A[35], &A[39], prime);
// 	GFPF_DFT_2(&A[19], &A[23], prime);
// 	GFPF_DFT_2(&A[51], &A[55], prime);
// 	GFPF_DFT_2(&A[11], &A[15], prime);
// 	GFPF_DFT_2(&A[43], &A[47], prime);
// 	GFPF_DFT_2(&A[27], &A[31], prime);
// 	GFPF_DFT_2(&A[59], &A[63], prime);

// 	A[34] = MulModMont(A[34], omegas[2], prime, pP);
// 	A[18] = MulModMont(A[18], omegas[4], prime, pP);
// 	A[50] = MulModMont(A[50], omegas[6], prime, pP);
// 	A[10] = MulModMont(A[10], omegas[8], prime, pP);
// 	A[42] = MulModMont(A[42], omegas[10], prime, pP);
// 	A[26] = MulModMont(A[26], omegas[12], prime, pP);
// 	A[58] = MulModMont(A[58], omegas[14], prime, pP);
// 	A[6] = MulModMont(A[6], omegas[16], prime, pP);
// 	A[38] = MulModMont(A[38], omegas[18], prime, pP);
// 	A[22] = MulModMont(A[22], omegas[20], prime, pP);
// 	A[54] = MulModMont(A[54], omegas[22], prime, pP);
// 	A[14] = MulModMont(A[14], omegas[24], prime, pP);
// 	A[46] = MulModMont(A[46], omegas[26], prime, pP);
// 	A[30] = MulModMont(A[30], omegas[28], prime, pP);
// 	A[62] = MulModMont(A[62], omegas[30], prime, pP);

// 	A[35] = MulModMont(A[35], omegas[2], prime, pP);
// 	A[19] = MulModMont(A[19], omegas[4], prime, pP);
// 	A[51] = MulModMont(A[51], omegas[6], prime, pP);
// 	A[11] = MulModMont(A[11], omegas[8], prime, pP);
// 	A[43] = MulModMont(A[43], omegas[10], prime, pP);
// 	A[27] = MulModMont(A[27], omegas[12], prime, pP);
// 	A[59] = MulModMont(A[59], omegas[14], prime, pP);
// 	A[7] = MulModMont(A[7], omegas[16], prime, pP);
// 	A[39] = MulModMont(A[39], omegas[18], prime, pP);
// 	A[23] = MulModMont(A[23], omegas[20], prime, pP);
// 	A[55] = MulModMont(A[55], omegas[22], prime, pP);
// 	A[15] = MulModMont(A[15], omegas[24], prime, pP);
// 	A[47] = MulModMont(A[47], omegas[26], prime, pP);
// 	A[31] = MulModMont(A[31], omegas[28], prime, pP);
// 	A[63] = MulModMont(A[63], omegas[30], prime, pP);

// 	// I_{2} I_{16} \tx GFPF_DFT_2
// 	GFPF_DFT_2(&A[0], &A[2], prime);
// 	GFPF_DFT_2(&A[32], &A[34], prime);
// 	GFPF_DFT_2(&A[16], &A[18], prime);
// 	GFPF_DFT_2(&A[48], &A[50], prime);
// 	GFPF_DFT_2(&A[8], &A[10], prime);
// 	GFPF_DFT_2(&A[40], &A[42], prime);
// 	GFPF_DFT_2(&A[24], &A[26], prime);
// 	GFPF_DFT_2(&A[56], &A[58], prime);

// 	GFPF_DFT_2(&A[4], &A[6], prime);
// 	GFPF_DFT_2(&A[36], &A[38], prime);
// 	GFPF_DFT_2(&A[20], &A[22], prime);
// 	GFPF_DFT_2(&A[52], &A[54], prime);
// 	GFPF_DFT_2(&A[12], &A[14], prime);
// 	GFPF_DFT_2(&A[44], &A[46], prime);
// 	GFPF_DFT_2(&A[28], &A[30], prime);
// 	GFPF_DFT_2(&A[60], &A[62], prime);

// 	GFPF_DFT_2(&A[1], &A[3], prime);
// 	GFPF_DFT_2(&A[33], &A[35], prime);
// 	GFPF_DFT_2(&A[17], &A[19], prime);
// 	GFPF_DFT_2(&A[49], &A[51], prime);
// 	GFPF_DFT_2(&A[9], &A[11], prime);
// 	GFPF_DFT_2(&A[41], &A[43], prime);
// 	GFPF_DFT_2(&A[25], &A[27], prime);
// 	GFPF_DFT_2(&A[57], &A[59], prime);

// 	GFPF_DFT_2(&A[5], &A[7], prime);
// 	GFPF_DFT_2(&A[37], &A[39], prime);
// 	GFPF_DFT_2(&A[21], &A[23], prime);
// 	GFPF_DFT_2(&A[53], &A[55], prime);
// 	GFPF_DFT_2(&A[13], &A[15], prime);
// 	GFPF_DFT_2(&A[45], &A[47], prime);
// 	GFPF_DFT_2(&A[29], &A[31], prime);
// 	GFPF_DFT_2(&A[61], &A[63], prime);

// 	A[33] = MulModMont(A[33], omegas[1], prime, pP);
// 	A[17] = MulModMont(A[17], omegas[2], prime, pP);
// 	A[49] = MulModMont(A[49], omegas[3], prime, pP);
// 	A[9] = MulModMont(A[9], omegas[4], prime, pP);
// 	A[41] = MulModMont(A[41], omegas[5], prime, pP);
// 	A[25] = MulModMont(A[25], omegas[6], prime, pP);
// 	A[57] = MulModMont(A[57], omegas[7], prime, pP);
// 	A[5] = MulModMont(A[5], omegas[8], prime, pP);
// 	A[37] = MulModMont(A[37], omegas[9], prime, pP);
// 	A[21] = MulModMont(A[21], omegas[10], prime, pP);
// 	A[53] = MulModMont(A[53], omegas[11], prime, pP);
// 	A[13] = MulModMont(A[13], omegas[12], prime, pP);
// 	A[45] = MulModMont(A[45], omegas[13], prime, pP);
// 	A[29] = MulModMont(A[29], omegas[14], prime, pP);
// 	A[61] = MulModMont(A[61], omegas[15], prime, pP);
// 	A[3] = MulModMont(A[3], omegas[16], prime, pP);
// 	A[35] = MulModMont(A[35], omegas[17], prime, pP);
// 	A[19] = MulModMont(A[19], omegas[18], prime, pP);
// 	A[51] = MulModMont(A[51], omegas[19], prime, pP);
// 	A[11] = MulModMont(A[11], omegas[20], prime, pP);
// 	A[43] = MulModMont(A[43], omegas[21], prime, pP);
// 	A[27] = MulModMont(A[27], omegas[22], prime, pP);
// 	A[59] = MulModMont(A[59], omegas[23], prime, pP);
// 	A[7] = MulModMont(A[7], omegas[24], prime, pP);
// 	A[39] = MulModMont(A[39], omegas[25], prime, pP);
// 	A[23] = MulModMont(A[23], omegas[26], prime, pP);
// 	A[55] = MulModMont(A[55], omegas[27], prime, pP);
// 	A[15] = MulModMont(A[15], omegas[28], prime, pP);
// 	A[47] = MulModMont(A[47], omegas[29], prime, pP);
// 	A[31] = MulModMont(A[31], omegas[30], prime, pP);
// 	A[63] = MulModMont(A[63], omegas[31], prime, pP);

// 	// I_{32} \tx GFPF_DFT_2
// 	GFPF_DFT_2(&A[0], &A[1], prime);
// 	GFPF_DFT_2(&A[32], &A[33], prime);
// 	GFPF_DFT_2(&A[16], &A[17], prime);
// 	GFPF_DFT_2(&A[48], &A[49], prime);
// 	GFPF_DFT_2(&A[8], &A[9], prime);
// 	GFPF_DFT_2(&A[40], &A[41], prime);
// 	GFPF_DFT_2(&A[24], &A[25], prime);
// 	GFPF_DFT_2(&A[56], &A[57], prime);

// 	GFPF_DFT_2(&A[4], &A[5], prime);
// 	GFPF_DFT_2(&A[36], &A[37], prime);
// 	GFPF_DFT_2(&A[20], &A[21], prime);
// 	GFPF_DFT_2(&A[52], &A[53], prime);
// 	GFPF_DFT_2(&A[12], &A[13], prime);
// 	GFPF_DFT_2(&A[44], &A[45], prime);
// 	GFPF_DFT_2(&A[28], &A[29], prime);
// 	GFPF_DFT_2(&A[60], &A[61], prime);

// 	GFPF_DFT_2(&A[2], &A[3], prime);
// 	GFPF_DFT_2(&A[34], &A[35], prime);
// 	GFPF_DFT_2(&A[18], &A[19], prime);
// 	GFPF_DFT_2(&A[50], &A[51], prime);
// 	GFPF_DFT_2(&A[10], &A[11], prime);
// 	GFPF_DFT_2(&A[42], &A[43], prime);
// 	GFPF_DFT_2(&A[26], &A[27], prime);
// 	GFPF_DFT_2(&A[58], &A[59], prime);

// 	GFPF_DFT_2(&A[6], &A[7], prime);
// 	GFPF_DFT_2(&A[38], &A[39], prime);
// 	GFPF_DFT_2(&A[22], &A[23], prime);
// 	GFPF_DFT_2(&A[54], &A[55], prime);
// 	GFPF_DFT_2(&A[14], &A[15], prime);
// 	GFPF_DFT_2(&A[46], &A[47], prime);
// 	GFPF_DFT_2(&A[30], &A[31], prime);
// 	GFPF_DFT_2(&A[62], &A[63], prime);

// 	// Final permutation
// 	swap(&A[1], &A[32]);
// 	swap(&A[2], &A[16]);
// 	swap(&A[3], &A[48]);
// 	swap(&A[4], &A[8]);
// 	swap(&A[5], &A[40]);
// 	swap(&A[6], &A[24]);
// 	swap(&A[7], &A[56]);
// 	swap(&A[9], &A[36]);
// 	swap(&A[10], &A[20]);
// 	swap(&A[11], &A[52]);
// 	swap(&A[13], &A[44]);
// 	swap(&A[14], &A[28]);
// 	swap(&A[15], &A[60]);
// 	swap(&A[17], &A[34]);
// 	swap(&A[19], &A[50]);
// 	swap(&A[21], &A[42]);
// 	swap(&A[22], &A[26]);
// 	swap(&A[23], &A[58]);
// 	swap(&A[25], &A[38]);
// 	swap(&A[27], &A[54]);
// 	swap(&A[29], &A[46]);
// 	swap(&A[31], &A[62]);
// 	swap(&A[35], &A[49]);
// 	swap(&A[37], &A[41]);
// 	swap(&A[39], &A[57]);
// 	swap(&A[43], &A[53]);
// 	swap(&A[47], &A[61]);
// 	swap(&A[55], &A[59]);

// 	return A;
// }


// //bigint -> vector in radix-based representation
// void mpz_to_radix_based_s64(sfixn64 *vector, const mpz_t bigint, sfixn64 radix,
// 		int max_input_vector_size)
// {

// 	mpz_t q, r;
// 	mpz_init_set_ui(q, 0);
// 	mpz_init_set_ui(r, 0);

// 	mpz_set(q, bigint);
// 	for (int i = 0; i < max_input_vector_size; i++)
// 	{
// 		//      mpz_fdiv_qr (q, r, bigint, base);
// 		mpz_fdiv_qr_ui(q, r, q, radix);
// 		vector[i] = mpz_get_ui(r);
// 		//      mpz_set (bigint, q);
// 	}

// 	mpz_clear(q);
// 	mpz_clear(r);
// }

// /**************************************/


// //convert radix-based representation -> vector
// void radix_based_to_mpz_s64(mpz_t bigint, sfixn64 *vector, sfixn64 radix,
// 		int input_vector_size)
// {
// //	char cmd_str[1024];
// 	mpz_init_set_ui(bigint, 0);
// //Horner method for evaluation
// 	mpz_init_set_ui(bigint, vector[input_vector_size - 1]);
// 	for (int i = input_vector_size - 1; i > 0; i--)
// 	{
// 		mpz_mul_ui(bigint, bigint, radix);
// 		mpz_add_ui(bigint, bigint, vector[i - 1]);
// 	}
// }

// /**************************************/

// void compute_srgfn_p_gmp(mpz_t p, usfixn64 radix, int coefficient_size)
// {
// 	mpz_t m;
// 	mpz_t r;

// 	mpz_init_set_ui(m, 0);
// 	mpz_init_set_ui(r, radix);

// 	mpz_pow_ui(m, r, coefficient_size);
// 	mpz_add_ui(p, m, 1);
// 	mpz_clear(m);
// }


// /**************************************/

// int verify_mult_u64_u64_gmp(const usfixn64 a, const usfixn64 b,
// 		const usfixn64 s0, const usfixn64 s1)
// {
// //	printf("a=%lu\n", a);
// //	printf("b=%lu\n", b);
// //	printf("s0in=%lu\n", s0);
// //	printf("s1in=%lu\n", s1);
// 	mpz_t x_zz;
// 	mpz_init_set_ui(x_zz, a);
// 	mpz_mul_ui(x_zz, x_zz, b);

// 	mpz_t u64_zz;
// 	mpz_init_set_ui(u64_zz, U64_MASK);
// 	mpz_add_ui(u64_zz, u64_zz, 1);

// 	mpz_t s0_zz, s1_zz;
// 	mpz_init(s0_zz);
// 	mpz_init(s1_zz);

// 	mpz_tdiv_qr(s1_zz, s0_zz, x_zz, u64_zz);

// 	if (mpz_get_ui(s0_zz) != s0)
// 	{
// 		printf("mismatch in s0!\n");
// 		printf("s0     =%lu\n", s0);
// 		printf("s0[gmp]=%lu\n", mpz_get_ui(s0_zz));
// 		return -1;
// 	}
// 	if (mpz_get_ui(s1_zz) != s1)
// 	{
// 		printf("mismatch in s1!\n");
// 		printf("s1=%lu\n", s1);
// 		printf("s1[gmp]=%lu\n", mpz_get_ui(s1_zz));
// 		return -1;
// 	}
// 	return 0;
// }

// /**************************************/

// int verify_mult_u64_u128_gmp(const usfixn64 a, const usfixn64 b0,
// 		const usfixn64 b1, const usfixn64 s0, const usfixn64 s1)
// {
// //	printf("a=%lu\n", a);
// //	printf("b=%lu\n", b);
// //	printf("s0in=%lu\n", s0);
// //	printf("s1in=%lu\n", s1);
// 	mpz_t x_zz;
// 	mpz_init_set_ui(x_zz, a);

// 	//u64
// 	mpz_t u64_zz;
// 	mpz_init_set_ui(u64_zz, U64_MASK);
// 	mpz_add_ui(u64_zz, u64_zz, 1);

// 	//b=b0+b1.u64
// 	mpz_t b_zz;
// 	mpz_init_set_ui(b_zz, b1);
// 	mpz_mul(b_zz, b_zz, u64_zz);
// 	mpz_add_ui(b_zz, b_zz, b0);

// 	//x=a.b
// 	mpz_mul(x_zz, x_zz, b_zz);

// 	//x -> (s0,s1)
// 	mpz_t s0_zz, s1_zz;
// 	mpz_init(s0_zz);
// 	mpz_init(s1_zz);

// 	mpz_tdiv_qr(s1_zz, s0_zz, x_zz, u64_zz);

// 	if (mpz_get_ui(s0_zz) != s0)
// 	{
// 		printf("mismatch in s0!\n");
// 		printf("s0     =%lu\n", s0);
// 		printf("s0[gmp]=%lu\n", mpz_get_ui(s0_zz));
// 		return -1;
// 	}
// 	if (mpz_get_ui(s1_zz) != s1)
// 	{
// 		printf("mismatch in s1!\n");
// 		printf("s1=%lu\n", s1);
// 		printf("s1[gmp]=%lu\n", mpz_get_ui(s1_zz));
// 		return -1;
// 	}
// 	return 0;
// }

// /**************************************/

// int verify_sub_u128_u128_gmp(const usfixn64 a0, const usfixn64 a1,
// 		const usfixn64 b0, const usfixn64 b1, const usfixn64 s0,
// 		const usfixn64 s1)
// {
// //	printf("a=%lu\n", a);
// //	printf("b=%lu\n", b);
// //	printf("s0in=%lu\n", s0);
// //	printf("s1in=%lu\n", s1);

// //u64
// 	mpz_t u64_zz;
// 	mpz_init_set_ui(u64_zz, U64_MASK);
// 	mpz_add_ui(u64_zz, u64_zz, 1);

// 	//b=b0+b1.u64
// 	mpz_t b_zz;
// 	mpz_init_set_ui(b_zz, b1);
// 	mpz_mul(b_zz, b_zz, u64_zz);
// 	mpz_add_ui(b_zz, b_zz, b0);

// 	mpz_t x_zz;
// 	mpz_init_set_ui(x_zz, a1);
// 	mpz_mul(x_zz, x_zz, u64_zz);
// 	mpz_add_ui(x_zz, x_zz, a0);

// 	//x=a-b
// 	mpz_sub(x_zz, x_zz, b_zz);

// 	//x -> (s0,s1)
// 	mpz_t s0_zz, s1_zz;
// 	mpz_init(s0_zz);
// 	mpz_init(s1_zz);

// 	mpz_tdiv_qr(s1_zz, s0_zz, x_zz, u64_zz);

// 	if (mpz_get_ui(s0_zz) != s0)
// 	{
// 		printf("mismatch in s0!\n");
// 		printf("s0     =%lu\n", s0);
// 		printf("s0[gmp]=%lu\n", mpz_get_ui(s0_zz));
// 		return -1;
// 	}
// 	if (mpz_get_ui(s1_zz) != s1)
// 	{
// 		printf("mismatch in s1!\n");
// 		printf("s1=%lu\n", s1);
// 		printf("s1[gmp]=%lu\n", mpz_get_ui(s1_zz));
// 		return -1;
// 	}
// 	return 0;
// }

// /**************************************/

// int verify_add_u128_u128_gmp(const usfixn64 a0, const usfixn64 a1,
// 		const usfixn64 b0, const usfixn64 b1, const usfixn64 s0,
// 		const usfixn64 s1, const usfixn64 s2)
// {
// //	printf("a0=%lu\n", a0);
// //	printf("a1=%lu\n", a1);
// //	printf("b0=%lu\n", b0);
// //	printf("b1=%lu\n", b1);
// //	printf("s0in=%lu\n", s0);
// //	printf("s1in=%lu\n", s1);
// //	printf("s2in=%lu\n", s2);

// //u64
// 	mpz_t u64_zz;
// 	mpz_init_set_ui(u64_zz, U64_MASK);
// 	mpz_add_ui(u64_zz, u64_zz, 1);

// 	//b=b0+b1.u64
// 	mpz_t b_zz;
// 	mpz_init_set_ui(b_zz, b1);
// 	mpz_mul(b_zz, b_zz, u64_zz);
// 	mpz_add_ui(b_zz, b_zz, b0);

// 	mpz_t x_zz;
// 	mpz_init_set_ui(x_zz, a1);
// 	mpz_mul(x_zz, x_zz, u64_zz);
// 	mpz_add_ui(x_zz, x_zz, a0);

// 	//x=a-b
// 	mpz_add(x_zz, x_zz, b_zz);

// 	//x -> (s0,s1)
// 	mpz_t s0_zz, s1_zz, s2_zz;
// 	mpz_init(s0_zz);
// 	mpz_init(s1_zz);
// 	mpz_init(s2_zz);

// 	mpz_tdiv_qr(s2_zz, s0_zz, x_zz, u64_zz);
// 	mpz_tdiv_qr(s2_zz, s1_zz, s2_zz, u64_zz);

// 	if (mpz_get_ui(s0_zz) != s0)
// 	{
// 		printf("mismatch in s0!\n");
// 		printf("s0     =%lu\n", s0);
// 		printf("s0[gmp]=%lu\n", mpz_get_ui(s0_zz));
// 		return -1;
// 	}
// 	if (mpz_get_ui(s1_zz) != s1)
// 	{
// 		printf("mismatch in s1!\n");
// 		printf("s1=%lu\n", s1);
// 		printf("s1[gmp]=%lu\n", mpz_get_ui(s1_zz));
// 		return -1;
// 	}

// 	if (mpz_get_ui(s2_zz) != s2)
// 	{
// 		printf("mismatch in s2!\n");
// 		printf("s2=%lu\n", s2);
// 		printf("s2[gmp]=%lu\n", mpz_get_ui(s2_zz));
// 		return -1;
// 	}
// 	return 0;
// }

// /**************************************/

// /**************************************/
// // faster version of mult_u64_u64
// // relies on __int128 support of gcc for x64 processors
// // a and be must be references; any other types will
// // severely affect the correctness of assembly.
// inline void mult_u64_u64(const usfixn64 a, const usfixn64 b, usfixn64* s0,
// 		usfixn64 *s1)
// {
// ////	equivalent version using gcc128
// //	__int128 mult = (__int128) a * (__int128) b;
// //	s0 = mult & (U64_MASK);
// //	s1 = mult >> 64;

// //	printf("a=%lu\n", a);
// //	printf("b=%lu\n", b);
// //	__asm__ (
// //			"movq  %2, %%rax;\n\t"          // rax = a
// //			"mulq  %3;\n\t"// rdx:rax = a * b
// //			"movq  %%rax, %0;\n\t"// s0 = rax
// //			"movq  %%rdx, %1;\n\t"// s1 = rdx
// //			: "=&rm" (s0),"=&rm"(s1)
// //			: "rm"(a), "rm"(b)
// //			: "%rax", "%rdx");

// 	asm volatile(
// 			"movq  %2, %%rax;\n\t"          // rax = a
// 		//	"movq  %3, %%rdx;\n\t"// rax = a
// 		//	"mulq  %%rdx;\n\t"// rdx:rax = a * b
// 			"mulq %3;\n\t"
//       "movq  %%rax, %0;\n\t"// s0 = rax
// 			"movq  %%rdx, %1;\n\t"// s1 = rdx
// 			: "=&rm" (*s0),"=&rm"(*s1)
// 			: "rm"(a), "rm"(b)
// 			: "%rax", "%rdx");//, "memory", "cc");

// #if VERIFICATION_ENABLED == 1
// 	printf("IN MULT!\n");

// 	if (verify_mult_u64_u64_gmp(a, b, *s0, *s1) != 0)
// 	{
// 		printf("FAILED! @mult_u64_u64\n");
// 		exit(0);
// 		return;
// 	}
// 	else
// 	{
// 		printf("VERIFIED! @mult_u64_u64\n");
// 	}
// 	printf("=========\n");
// #endif
// }

// /**************************************/

// //compute a*(b0+b1.u64);
// //results are correct only if the whole product
// //is less than u128.
// inline void mult_u64_u128(const usfixn64 a, const usfixn64 b0, const usfixn64 b1,
// 		usfixn64* s0, usfixn64* s1)
// {
// ////	equivalent version using gcc128
// //	__int128 mult = (__int128) a * (__int128) b;
// //	s0 = mult & (U64_MASK);
// //	s1 = mult >> 64;

// //	__asm__ (
// //			"movq  %2, %%rax;\n\t"          // rax = a
// //			"mulq  %3;\n\t"// rdx:rax = a * b
// //			"movq  %%rax, %0;\n\t"// s0 = rax
// //			"movq  %%rdx, %1;\n\t"// s1 = rdx
// //			: "=&rm" (s0),"=&rm"(s1)
// //			: "rm"(a), "rm"(b)
// //			: "%rax", "%rdx");

// 	//ab0+ab1u64

// //	printf("a=%lu\n", a);
// //	printf("b0=%lu\n", b0);
// //	printf("b1=%lu\n", b1);

// 	//(s0,s1)=a.b0
// 	asm volatile(
// 			"movq  %2, %%rax;\n\t"          // rax = a
// 			//"movq  %3, %%rdx;\n\t"// rax = a
// 			//"mulq  %%rdx;\n\t"// rdx:rax = a * b
//       "mulq %3;\n\t"
// 			"movq  %%rax, %0;\n\t"// s0 = rax
// 			"movq  %%rdx, %1;\n\t"// s1 = rdx

// 			"movq  %2, %%rax;\n\t"// rax = a
// 			//"movq  %4, %%rdx;\n\t"// rax = a
//       "xor %%rdx,%%rdx;\n\t"
//       "mulq %4;\n\t"
// 			"mulq  %%rdx; \n\t"
// 			"addq  %%rax, %1;\n\t"
// 			: "+&rm" (*s0),"+&rm"(*s1)
// 			: "rm"(a), "rm"(b0), "rm"(b1)
// 			: "%rax", "%rdx");//, "memory", "cc");

// #if VERIFICATION_ENABLED == 1
// 	//(s0,s1)=(s0,s1)+a.b1.u64
// 	printf("@mult_u64_u128!\n");

// 	if (verify_mult_u64_u128_gmp(a, b0, b1, *s0, *s1) != 0)
// 	{
// 		printf("FAILED! @mult_u64_u128\n");
// 		exit(0);
// 		return;
// 	}
// 	else
// 	{
// 		printf("VERIFIED! @mult_u64_u128\n");
// 	}
// 	printf("=========\n");
// #endif
// }

// /**************************************/

// //(x0,x1)=(x0,x1)-(y0,y1);
// inline void sub_u128_u128(const usfixn64 y0, const usfixn64 y1, usfixn64* x0,
// 		usfixn64* x1)
// {
// ////	equivalent version using gcc128
// //	__int128 mult = (__int128) a * (__int128) b;
// //	s0 = mult & (U64_MASK);
// //	s1 = mult >> 64;

// //	__asm__ (
// //			"movq  %2, %%rax;\n\t"          // rax = a
// //			"mulq  %3;\n\t"// rdx:rax = a * b
// //			"movq  %%rax, %0;\n\t"// s0 = rax
// //			"movq  %%rdx, %1;\n\t"// s1 = rdx
// //			: "=&rm" (s0),"=&rm"(s1)
// //			: "rm"(a), "rm"(b)
// //			: "%rax", "%rdx");

// 	//ab0+ab1u64

// 	//(s0,s1)=a.b0

// 	usfixn64 a0, a1;
// 	a0 = *x0;
// 	a1 = *x1;

// 	usfixn64 m0, m1;
// 	m0 = 0 - 1 - y0;
// 	m1 = 0 - 1 - y1;
// 	asm volatile(
// 			"addq  %2, %0;\n\t"          // rax = a
// 			"adcq  %3, %1;\n\t"// rax = a
// 			"addq  $1, %0;\n\t"// rax = a
// 			"adcq  $0, %1;\n\t"// rax = a
// 			: "+&rm" (*x0),"+&rm"(*x1)
// 			: "rm"(m0), "rm"(m1) );
// //			: "memory", "cc");

// #if VERIFICATION_ENABLED == 1
// 	//(s0,s1)=(s0,s1)+a.b1.u64
// 	printf("@sub_u128_u128!\n");

// 	if (verify_sub_u128_u128_gmp(a0, a1, y0, y1, *x0, *x1) != 0)
// 	{
// 		printf("FAILED! @sub_u128_u128\n");
// 		exit(0);
// 		return;
// 	}
// 	else
// 	{
// 		printf("VERIFIED! @sub_u128_u128\n");
// 	}
// 	printf("=========\n");
// #endif
// }

// /**************************************/

// inline void add_u128_u128(const usfixn64 x0, const usfixn64 x1, const usfixn64 y0,
// 		const usfixn64 y1, usfixn64* s0, usfixn64* s1, usfixn64* s2)
// {
// //	__int128 mult = (__int128) a * (__int128) b;
// //	s0 = mult & (U64_MASK);
// //	s1 = mult >> 64;

// //	s0 = 0;
// //	s1 = 0;
// //	s2 = 0;

// //	usfixn64 a0, a1, a2;
// //	__asm__ (
// //			"movq  %3, %0;\n\t" //r0=x0
// //			"addq  %5, %0;\n\t"//r0+=y0
// //			"adcq  %4, %1;\n\t"
// //			"adcq  $0, %2;\n\t"
// //			"addq  %5, %1;\n\t"
// //			"adcq  $0, %2;\n\t"
// //			: "=rm" (s0),"=rm"(s1),"=rm"(s2)
// //			: "rm"(x0), "rm"(x1),"rm"(y0),"rm"(y1)
// //			: "%rax", "%rdx", "memory");
// //	__asm__ volatile (
// //			"movq  %3, %0;\n\t" //r0=x0
// //			"addq  %5, %0;\n\t"//r0+=y0
// //			"adcq  %4, %1;\n\t"
// //			"adcq  $0, %2;\n\t"
// //			"addq  %5, %1;\n\t"
// //			"adcq  $0, %2;\n\t"
// //			: "+&q" (a0),"+&q"(a1),"+&q"(a2)
// //			: "q"(x0), "q"(x1),"q"(y0),"q"(y1)
// //			: "memory");

// 	// a0 = 0;
// 	// a1 = 0;
// 	// a2 = 0;
//   *s0 = x0;
//   *s1 = x1;
// 	__asm__ volatile (
// 			//"movq  %3, %0;\n\t" //r0=x0
// 			//"movq  %4, %1;\n\t"//r0=x0
// 			"addq  %5, %0;\n\t"//r0+=y0
// 			"adcq  %6, %1;\n\t"
// 			"adcq  $0, %2;\n\t"
// 			: "+&rm" (*s0),"+&rm"(*s1),"+&rm"(*s2)
// 			: "rm"(x0), "rm"(x1),"rm"(y0),"rm"(y1));
// //			: "memory");
// 	// *s0 = a0;
// 	// *s1 = a1;
// 	// *s2 = a2;

// #if VERIFICATION_ENABLED == 1
// 	if (verify_add_u128_u128_gmp(x0, x1, y0, y1, *s0, *s1, *s2) != 0)
// 	{
// 		printf("FAILED! @add_u128_u128!\n");
// 		exit(0);
// 	}
// 	else
// 	{

// 		printf("VERIFIED! @add_u128_u128!\n");
// 	}
// #endif
// //	printf("x0=%lu\n",x0);
// //	printf("x1=%lu\n",x1);
// //	printf("x0=%lu\n",y0);
// //	printf("x1=%lu\n",y1);
// //	printf("s0=%lu\n",s0);
// //	printf("s1=%lu\n",s1);
// //	printf("s2=%lu\n",s2);
// }

// /**************************************/
// /**************************************/

// // computes (a_u64)*(b0_u64, b1_u64)
// // a=a_u64;
// // b=b0_u64 + (b1_u64)*u64
// /*************************************************/

// inline void  u64_mod_u64(usfixn64* a, const usfixn64* n)
// {

// //	a = a % n;
// //	return;
// 	double ninv = 1 / (double) *n;
// 	usfixn64 q = (usfixn64) ((((double) *a)) * ninv);
// 	usfixn64 res;
// 	res = *a - q * (*n);
// 	*a = res & (U64_MASK);
// }

// /**************************************/

// //mult_u128_u128_hi128: returns ((x0+x1.u64)*(y0+y1.u64))>>128
// inline void  mult_u128_u128_hi128(const usfixn64* x0, const usfixn64* x1,
// 		const usfixn64* y0, const usfixn64* y1, usfixn64* q)
// {
// 	usfixn64 s0, s1, s2;
// 	usfixn64 c1;

// 	*q = 0;

// 	//	s0 = (__int128) x0 * (__int128) y0;
// 	// 	s0>>=64;
// 	__asm__ volatile(
// 			"movq  %1, %%rax;\n\t"          // rax = a
// 			"mulq  %2;\n\t"// rdx:rax = a * b
// 			"movq  %%rdx, %0;\n\t"// s1 = rdx
// 			: "=rm" (s0)
// 			: "rm"(*x0), "rm"(*y0)
// 			: "%rax", "%rdx");//, "memory");

// //	s1 = (__int128) x1 * (__int128) y0;
// //	c1 = (s1 >> 64);
// //	s1 = s1 & (U64_MASK);

// 	__asm__ volatile(
// 			"movq  %2, %%rax;\n\t"          // rax = a
// 			"mulq  %3;\n\t"// rdx:rax = a * b
// 			"movq  %%rax, %0;\n\t"// s1 = rdx
// 			"movq  %%rdx, %1;\n\t"// s1 = rdx
// 			: "=rm" (s1),"=rm"(c1)
// 			: "rm"(*x1), "rm"(*y0)
// 			: "%rax", "%rdx");//, "memory");

// //	s2 = (__int128) x0 * (__int128) y1;
// //	c2 = (s2 >> 64);
// //	s2 = s2 & (U64_MASK);

// 	__asm__ volatile(
// 			"movq  %2, %%rax;\n\t"          // rax = a
// 			"mulq  %3;\n\t"// rdx:rax = a * b
// 			"movq  %%rax, %0;\n\t"// s1 = rdx
// //			"movq  %%rdx, %1;\n\t"// s1 = rdx
// 			"addq  %%rdx, %1;\n\t"// s1 = rdx
// 			: "=rm" (s2),"+rm"(c1)
// 			: "rm"(*x0), "rm"(*y1)
// 			: "%rax", "%rdx");//,"memory");

// //	c1+=c2;
// 	*q += c1;
// //	s3 = (__int128) x1 * (__int128) y1;
// 	*q += (*x1) * (*y1);

// 	__asm__ volatile(
// 			"movq  %1, %%rax;\n\t"          // rax = a
// 			"addq  %2, %%rax;\n\t"// rdx:rax = a * b
// 			"adcq  $0x0, %0;\n\t"
// 			"addq  %3, %%rax;\n\t"// rdx:rax = a * b
// 			"adcq  $0x0, %0;\n\t"
// //			"movq  %%rax, %0;\n\t"// s1 = rdx
// 			: "+rm"(*q)
// 			: "rm" (s0), "rm"(s1), "rm"(s2)
// 			: "%rax");//,"memory");
// }

// /**************************************/

// //inv_p_mod_u128: returns [p_inv_m, p_inv_q];
// // p_inv = p_inv_m + p_inv_q.u64  = int(u128/p);
// inline void inv_p_mod_u128(const usfixn64* p, usfixn64* p_inv_m, usfixn64* p_inv_q)
// {
// 	mpz_t u128;
// 	mpz_t u64;
// 	mpz_t q, m;

// 	mpz_init(u128);
// 	mpz_init(u64);
// 	mpz_init(q);
// 	mpz_init(m);

// 	mpz_set_ui(u64, U64_MASK);
// 	mpz_add_ui(u64, u64, 1);
// 	mpz_mul(u128, u64, u64);

// 	mpz_tdiv_q_ui(u128, u128, *p);
// 	mpz_tdiv_qr(q, m, u128, u64);
// 	*p_inv_m = mpz_get_ui(m);
// 	*p_inv_q = mpz_get_ui(q);
// }

// /**************************************/

// // - compute [a1m2p1+a2m1p2] and return [u0,u1,u2]
// int verify_crt_mult_sub_u192_with_reduction_gmp(const usfixn64 a1,
// 		const usfixn64 a2, crt_u192_data* data, const usfixn64 s0,
// 		const usfixn64 s1)
// {

// 	usfixn64 p1 = data->p1;
// 	usfixn64 p2 = data->p2;

// 	usfixn64 m1 = data->m1;
// 	usfixn64 m2 = data->m2;

// 	/*  compute the following using gmp:
// 	 * 	s=((a1*m2)%p1)*p2+((a2*m1)%p2)*p1
// 	 *  rewrite s==s0+s1->u64 as two machine words (s0,s1)
// 	 */

// 	mpz_t t0, t1, p1p2, u64;
// 	mpz_init(t0);
// 	mpz_init(t1);
// 	mpz_init(p1p2);
// 	mpz_init(u64);

// 	//u64=2^64
// 	mpz_set_ui(u64, U64_MASK);
// 	mpz_add_ui(u64, u64, 1);

// //	char * tmp = (char*) malloc(1024);

// 	//p1p2=p1*p2
// 	mpz_set_ui(p1p2, p1);
// 	//	mpz_get_str(tmp, 10, p1p2);
// 	////	printf("p1p2=%s\n", tmp);
// 	mpz_mul_ui(p1p2, p1p2, p2);

// 	//	mpz_get_str(tmp, 10, p1p2);
// 	//	printf("p1p2=%s\n", tmp);

// 	mpz_set_ui(t0, a1);
// 	mpz_mul_ui(t0, t0, m2);
// 	mpz_mul_ui(t0, t0, p2);

// 	//	mpz_get_str(tmp, 10, t0);
// 	//	printf("a1m2p2=%s\n", tmp);

// 	mpz_set_ui(t1, a2);
// 	mpz_mul_ui(t1, t1, m1);
// 	mpz_mul_ui(t1, t1, p1);

// 	//	mpz_get_str(tmp, 10, t1);
// 	//	printf("a2m1p1=%s\n", tmp);

// 	mpz_add(t0, t0, t1);
// 	mpz_mod(t0, t0, p1p2);

// 	//	mpz_get_str(tmp, 10, t0);
// 	//	printf("t0+t1=%s\n", tmp);

// 	mpz_t u0, u1;
// 	mpz_init(u0);
// 	mpz_init(u1);

// 	mpz_tdiv_qr(u1, u0, t0, u64);

// 	//	printf("s0=%lu \n", s0);
// 	//	printf("u0=%lu \n", mpz_get_ui(u0));
// 	//	printf("s1=%lu \n", s1);
// 	//	printf("u1=%lu \n", mpz_get_ui(u1));

// //	mpz_get_str(tmp, 10, t0);
// 	//	printf("t0=%s\n", tmp);

// 	if (mpz_get_ui(u0) != s0)
// 	{
// 		printf("mismatch in s0!\n");
// 		printf("s0     =%lu\n", s0);
// 		printf("s0[gmp]=%lu\n", mpz_get_ui(u0));
// 		return -1;
// 	}

// 	if (mpz_get_ui(u1) != s1)
// 	{
// 		printf("mismatch in s1!\n");
// 		printf("s1     =%lu\n", s1);
// 		printf("s1[gmp]=%lu\n", mpz_get_ui(u1));
// 		return -1;
// 	}

// //	if ((s0 == mpz_get_ui(u0)) && (s1 == mpz_get_ui(u1)))
// //	{
// //
// //		return 0;
// //
// //	}
// //	else
// //	{
// //		return -1;
// //	}

// 	return 0;

// }

// /**************************************/

// /* - compute [a1+m1p1+a2m2p2] and return [u0,u1,u2]
//  *   in base u64;
//  * - it is assumed that m2 is negative.
//  */

// inline void crt_mult_sub_u192_with_reduction(const usfixn64 a1, const usfixn64 a2,
// 		const crt_u192_data* data, usfixn64* s0, usfixn64* s1)
// {
// 	usfixn64 t[4];
// 	usfixn64 q[2];
// //	__int128 r[2];

// //	printf("in crt: mults\n");
// 	mult_u64_u64(a1, data->m2, &t[0], &t[1]);
// 	mult_u64_u64(a2, data->m1, &t[2], &t[3]);

// //	printf("u128 crt: mults\n");
// 	mult_u128_u128_hi128(&t[0], &t[1], &data->p1_inv_m, &data->p1_inv_q, &q[0]);
// 	mult_u128_u128_hi128(&t[2], &t[3], &data->p2_inv_m, &data->p2_inv_q, &q[1]);
// 	usfixn64 m0, m1;

// 	__asm__ volatile(
// 			"movq  %2, %%rax;\n\t"		// rax = a
// 			"mulq  %3;\n\t"// rdx:rax = a * b
// 			"movq  %%rax, %0;\n\t"// s1 = rdx
// 			"movq  %%rdx, %1;\n\t"// s1 = rdx
// 			: "=rm" (m0),"=rm"(m1)
// 			: "rm"(q[0]), "rm"(data->p1)
// 			: "%rax", "%rdx");//,"memory");

// 	m0 = U64_MASK - m0;
// 	m1 = U64_MASK - m1;

// 	__asm__ volatile(
// 			"addq %2, %0; \n\t"
// 			"adcq %3, %1; \n\t"
// 			"addq $0x1, %0; \n\t"
// 			"adcq $0x0, %1; \n\t"
// 			: "+rm" (t[0]),"+rm"(t[1])
// 			: "rm"(m0), "rm"(m1));
// //			: "memory");

// 	///////////////////////////////////

// 	m0 = 0;
// 	m1 = 0;

// 	__asm__ volatile(
// 			"movq  %2, %%rax;\n\t"		// rax = a
// 			"mulq  %3;\n\t"// rdx:rax = a * b
// 			"movq  %%rax, %0;\n\t"// s1 = rdx
// 			"movq  %%rdx, %1;\n\t"// s1 = rdx
// 			: "=rm" (m0),"=rm"(m1)
// 			: "rm"(q[1]), "rm"(data->p2)
// 			: "%rax", "%rdx");//,"memory");

// 	m0 = U64_MASK - m0;
// 	m1 = U64_MASK - m1;

// 	__asm__ volatile (
// 			"addq %2, %0; \n\t"
// 			"adcq %3, %1; \n\t"
// 			"addq $0x1, %0; \n\t"
// 			"adcq $0x0, %1; \n\t"
// 			: "+rm" (t[2]),"+rm"(t[3])
// 			: "rm"(m0), "rm"(m1));
// //			: "memory");

// 	///////////////////////////////////

// 	if (t[0] >= data->p1)
// 		t[0] -= data->p1;
// 	if (t[2] >= data->p2)
// 		t[2] -= data->p2;

// //	r[0] = (__int128) t[0] ;///+ ((__int128) t[1] << 64);
// //	r[1] = (__int128) t[2] ;//+ ((__int128) t[3] << 64);

// //	printf("third round of mults\n");
// 	mult_u64_u64(t[0], data->p2, &t[0], &t[1]);
// 	mult_u64_u64(t[2], data->p1, &t[2], &t[3]);
// //could be removed
// 	m0 = t[0];
// 	m1 = t[1];
// 	__asm__ volatile(
// 			"addq %2, %0; \n\t"
// 			"adcq %3, %1; \n\t"
// //					"addq $0x1, %0; \n\t"
// //					"adcq $0x0, %1; \n\t"
// 			: "+rm" (t[0]),"+rm"(t[1])
// 			: "rm"(t[2]), "rm"(t[3]));
// //			: "memory");
// 	if ((t[1] > data->p1p2_q) || ((t[1] == data->p1p2_q) && (t[0] > data->p1p2_m)))
// 	{
// 		m0 = U64_MASK - data->p1p2_m;
// 		m1 = U64_MASK - data->p1p2_q;
// 		__asm__ volatile(
// 				"addq %2, %0; \n\t"
// 				"adcq %3, %1; \n\t"
// 				"addq $0x1, %0; \n\t"
// 				"adcq $0x0, %1; \n\t"
// 				: "+rm" (t[0]),"+rm"(t[1])
// 				: "rm"(m0), "rm"(m1));
// //				: "memory");

// 	}
// 	*s0 = t[0];
// 	*s1 = t[1];

// #if VERIFICATION_ENABLED == 1
// 	if (verify_crt_mult_sub_u192_with_reduction_gmp(a1, a2, data, *s0, *s1) != 0)
// 	{
// 		printf("FAILED! @crt_mult_sub_u192_with_reduction\n");
// 		exit(0);
// 	}
// 	else
// 	{
// 		printf("VERIFIED! @crt_mult_sub_u192_with_reduction\n");
// 	}
// #endif

// }

// /**************************************/


// ///**************************************/
// int verify_div_by_const_R_gmp(const usfixn64 x0_u64, const usfixn64 x1_u64,
// 		const usfixn64 q, const usfixn64 m, sfixn64 R)
// {
// 	mpz_t x_zz;
// 	mpz_init_set_ui(x_zz, U64_MASK);
// 	mpz_add_ui(x_zz, x_zz, 1);
// 	mpz_mul_ui(x_zz, x_zz, x1_u64);
// 	mpz_add_ui(x_zz, x_zz, x0_u64);

// 	mpz_t q_zz, m_zz;
// 	mpz_init(q_zz);
// 	mpz_init(m_zz);
// 	mpz_tdiv_qr_ui(q_zz, m_zz, x_zz, R);

// 	if (mpz_get_ui(q_zz) != q)
// 	{
// 		printf("mismatch q!\n");
// 		printf("q =%lu\n", q);
// 		printf("q[gmp]=%lu\n", mpz_get_ui(q_zz));
// 		return -1;
// 	}

// 	if (mpz_get_ui(m_zz) != m)
// 	{
// 		printf("mismatch m!\n");
// 		printf("m =%lu\n", m);
// 		printf("m[gmp]=%lu\n", mpz_get_ui(m_zz));
// 		return -1;
// 	}
// 	return 0;

// }

// ///**************************************/

// inline void  div_by_const_R(const usfixn64 x0_u64, const usfixn64 x1_u64,
// 		const usfixn64 r0, const usfixn64 r1, usfixn64* q, usfixn64* m,
// 		sfixn64 R)
// {
// //	r_inv= (u128/r_in);
// //	r1,r0=[r_inv/u64, r_inv%u64];
// //	x1,x0=[x/u64, x%u64];
// //	v0=x0*r0;
// //	v1=x0*r1;
// //	v2=x1*r0;

// //	__int128 v1, v2, q0;
// //	usfixn64 x0 = x0_u64;
// //	usfixn64 x1 = x1_u64;
// //	v0 = 0;
// //	v1 = 0;
// //	v2 = 0;

// 	usfixn64 v0_lo = 0, v0_hi = 0;
// 	usfixn64 v1_lo = 0, v1_hi = 0;
// 	usfixn64 v2_lo = 0, v2_hi = 0;

// //	v0 = (__int128) x0 * (__int128) r0;
// 	mult_u64_u64(x0_u64, r0, &v0_lo, &v0_hi);
// //	v1 = (__int128) x0 * (__int128) r1;
// //	v2 = (__int128) x1 * (__int128) r0;

// 	mult_u64_u64(x0_u64, r1, &v1_lo, &v1_hi);
// 	mult_u64_u64(x1_u64, r0, &v2_lo, &v2_hi);

// //	v0 >>= 64;
// //	v1 += (v0_hi);
// //	v2 += v1;
// //	v2 >>= 64;

// 	// (v1_lo,v1_hi)+=(v0_hi) //add with carry
// //	__asm__ (
// //			"movq  %2, %%rax;\n\t"          // rax = a
// //			"addq  %%rax, %0;\n\t"// rdx:rax = a * b
// //			"adcq  $0x0, %1;\n\t"
// //			: "+&rm"(v1_lo),"+&rm"(v1_hi)
// //			: "rm" (v0_hi)
// //			: "%rax");
// 	__asm__ volatile(
// 			"addq  %2, %0;\n\t"		// rdx:rax = a * b
// 			"adcq  $0x0, %1;\n\t"
// 			: "+rm"(v1_lo),"+rm"(v1_hi)
// 			: "rm" (v0_hi));
// //			: "memory");

// 	// (v2_lo, v2_hi) += (v1_lo, v1_hi);
// //	__asm__ volatile(
// //			"movq  %2, %%rax;\n\t"          // rax = a
// //			"addq  %%rax, %0;\n\t"// rdx:rax = a * b
// //			"adcq  $0x0, %1;\n\t"
// //			"movq  %3, %%rax;\n\t"// rax = a
// //			"addq  %%rax, %1;\n\t"// rdx:rax = a * b
// //			: "+q"(v2_lo),"+q"(v2_hi)
// //			: "q" (v1_lo),"q" (v1_hi)
// //			: "%rax");

// 	__asm__ volatile(
// 			"addq  %2, %0;\n\t"		// rdx:rax = a * b
// 			"adcq  %3, %1;\n\t"
// 			: "+rm"(v2_lo),"+rm"(v2_hi)
// 			: "rm" (v1_lo),"rm" (v1_hi));
// //			:"memory");

// //	q0 = 0;
// //	q0 = (__int128) x1 * (__int128) r1;

// 	usfixn64 q0_hi = 0, q0_lo = 0;

// 	mult_u64_u64(x1_u64, r1, &q0_lo, &q0_hi);

// 	//at this point, only v2_hi is required.
// //	q0 += v2;
// //	q0+=v2_hi;

// //	__asm__ (
// //			"movq  %2, %%rax;\n\t"          // rax = a
// //			"addq  %%rax, %0;\n\t"// rdx:rax = a * b
// //			"adcq  $0x0, %1;\n\t"
// //			: "+&rm"(q0_lo),"+&rm"(q0_hi)
// //			: "rm" (v2_hi)
// //			: "%rax");
// 	__asm__ volatile (
// 			"addq  %2, %0;\n\t"		// rdx:rax = a * b
// 			"adcq  $0x0, %1;\n\t"
// 			: "+rm"(q0_lo),"+rm"(q0_hi)
// 			: "rm" (v2_hi));
// //			: "memory");

// //	printf("(q0_l0,hi)+v2= %lu , %lu\n", q0_lo, q0_hi);

// 	usfixn64 x0 = x0_u64;
// 	usfixn64 x1 = x1_u64;
// //	m0=x-(q0*r_in);
// //	__int128 m0;
// //	m0 = (__int128)(1L << 64);
// //	m0 *= (__int128) x1;
// //	m0 += (__int128) x0;
// //
// //	__int128 q0;
// //	q0 = (__int128)(1L << 64);
// //	q0 *= (__int128) q0_hi;
// //	q0 += (__int128) q0_lo;
// //	m0 = m0 - q0 * (__int128)(R);

// 	usfixn64 m0_lo, m0_hi;
// 	x0 = q0_lo;
// 	x1 = q0_hi;
// 	mult_u64_u128(R, q0_lo, q0_hi, &m0_lo, &m0_hi);

// 	x0 = x0_u64;
// 	x1 = x1_u64;
// 	sub_u128_u128(m0_lo, m0_hi, &x0, &x1);

// //	usfixn64 m0_lo, m0_hi;

// //	if (m0 >= (__int128)(R))
// //	{
// //		//	printf("carry\n");
// //		m0 -= (__int128) R;
// //		q0 += 1;
// //	}

// 	if ((x1 >= 1) || ((x1 == 0) && (x0 >= R)))
// 	{
// #if VERIFICATION_ENABLED ==1
// 		printf("carry\n");
// #endif
//     usfixn64 zero = 0;
// 		sub_u128_u128(R, zero, &x0, &x1);
// 		q0_lo++;
// 	}

// //	m = (usfixn64) (m0 & U64_MASK);
// //	q = (usfixn64) (q0 & (U64_MASK));

// 	*m = x0;
// 	*q = q0_lo;

// 	if ((q0_hi) > 0)
// 	{
// 		printf("WARNING: q >= u64!\n");
// 		exit(0);
// 	}

// #if VERIFICATION_ENABLED == 1
// 	printf("in div by const r\n"
// 			"x0 = %lu\n"
// 			"x1 = %lu\n"
// 			"R  = %lu\n"
// 			"m  = %lu\n"
// 			"q  = %lu\n", x0_u64, x1_u64, R, *m, *q);
// 	if (verify_div_by_const_R_gmp(x0_u64, x1_u64, *q, *m, R) != 0)
// 	{
// 		printf("FAILED! @div_by_const_R!\n");
// 		exit(0);
// 		return;
// 	}
// 	else
// 	{
// 		printf("VERIFIED! @div_by_const_R!\n");
// 	}
// 	printf("=========\n");
// #endif

// }

// /**************************************/

// int verify_add_lhc_gmp(const usfixn64 *lhc_0, const usfixn64 *lhc_1,
// 		const usfixn64 s0, const usfixn64 s1, const usfixn64 s2, usfixn64 r)
// {
// 	usfixn64 l0, h0, c0;
// 	usfixn64 l1, h1, c1;

// 	l0 = lhc_0[0];
// 	h0 = lhc_0[1];
// 	c0 = lhc_0[2];

// 	l1 = lhc_1[0];
// 	h1 = lhc_1[1];
// 	c1 = lhc_1[2];

// //	printf("r = %lu\n", r);
// //	printf("l0 = %lu \n", l0);
// //	printf("h0 = %lu \n", h0);
// //	printf("c0 = %lu \n", c0);
// //	printf("l1 = %lu \n", l1);
// //	printf("h1 = %lu \n", h1);
// //	printf("c1 = %lu \n", c1);
// //	printf("s0 = %lu \n", s0);
// //	printf("s1 = %lu \n", s1);
// //	printf("s2 = %lu \n", s2);

// 	mpz_t a_zz, b_zz, u64_zz;
// 	mpz_t s0_zz, s1_zz, s2_zz;

// 	mpz_init(s0_zz);
// 	mpz_init(s1_zz);
// 	mpz_init(s2_zz);

// 	mpz_init_set_ui(u64_zz, U64_MASK);
// 	mpz_add_ui(u64_zz, u64_zz, 1);

// 	mpz_init_set_ui(a_zz, c0);
// 	mpz_mul_ui(a_zz, a_zz, r);
// 	mpz_add_ui(a_zz, a_zz, h0);
// 	mpz_mul_ui(a_zz, a_zz, r);
// 	mpz_add_ui(a_zz, a_zz, l0);

// 	mpz_init_set_ui(b_zz, c1);
// 	mpz_mul_ui(b_zz, b_zz, r);
// 	mpz_add_ui(b_zz, b_zz, h1);
// 	mpz_mul_ui(b_zz, b_zz, r);
// 	mpz_add_ui(b_zz, b_zz, l1);

// 	mpz_add(a_zz, a_zz, b_zz);

// 	mpz_tdiv_qr_ui(a_zz, s0_zz, a_zz, r);
// 	mpz_tdiv_qr_ui(s2_zz, s1_zz, a_zz, r);

// 	if (mpz_get_ui(s0_zz) != s0)
// 	{
// 		printf("failed at s0!\n");
// 		printf("s0     = %lu\n", s0);
// 		printf("s0[gmp]= %lu\n", mpz_get_ui(s0_zz));
// 		return -1;
// 	}
// 	if (mpz_get_ui(s1_zz) != s1)
// 	{
// 		printf("failed at s1!\n");
// 		printf("s1     = %lu\n", s1);
// 		printf("s1[gmp]= %lu\n", mpz_get_ui(s1_zz));
// 		return -1;
// 	}
// 	if (mpz_get_ui(s2_zz) != s2)
// 	{
// 		printf("failed at s2!\n");
// 		printf("s2     = %lu\n", s2);
// 		printf("s2[gmp]= %lu\n", mpz_get_ui(s2_zz));
// 		return -1;
// 	}

// 	return 0;

// }

// /**************************************/

// /**************************************/
// // none of the input and output arrays should
// // point to the same memory address, namely,
// // *lhc_0, *lhc_1, *lhc_ans
// inline void add_lhc(const usfixn64 *lhc_0, const usfixn64* lhc_1,
// 		usfixn64 * lhc_ans, sfixn64 r)
// {
//   //  usfixn64 s3[3];
//   usfixn64 c_in = 0, c = 0, s = 0;
//  // s = lhc_ans[0];
//   for (int i = 0; i < 3; i++)
//     {
//       __asm__ __volatile__(
// //      "xorq  %0, %0;\n\t" //s=0
// //      "xorq  %1, %1;\n\t"//c=0
//     "movq  %2, %0;\n\t"//s=i0
//     "addq  %3, %0;\n\t"//s=i0+i1
//     "adcq  $0, %1;\n\t"//c+=carry;
//     "addq  %5, %0;\n\t"//s=i0+i1
//     "adcq  $0, %1;\n\t"//c+=carry;
//     "cmpq  %6, %0;\n\t"//c+=carry;
//     "jl L0_%=;\n\t"//c+=carry;
//     "subq %6, %0;\n\t"
//     "inc %1;\n\t"
//     "L0_%=:;\n\t"
//     "addq  %4, %0;\n\t"//s=i0+i1
//     "adcq  $0, %1;\n\t"//c+=carry;
//     "cmpq  %6, %0;\n\t"//c+=carry;
//     "jl L1_%=;\n\t"//c+=carry;
//     "subq %6, %0;\n\t"
//     "inc %1;\n\t"
//     "L1_%=:;\n\t"
//     : "+g" (s),"+g"(c)
//     : "g"(lhc_ans[i]), "g"(lhc_0[i]),"g"(lhc_1[i]), "g"(c_in),"g"(r)
//     : );
//       lhc_ans[i] = s;
//       c_in = c;
//       c = 0;
//     }
// //def add_lhc(lhc_0,lhc_1):

// // 	usfixn64 l0, h0, c0;
// // 	usfixn64 l1, h1, c1;
// // 	usfixn64 s0, s1, s2;
// // 	usfixn64 c = 0;
// // 	usfixn64 s = 0;
// // 	usfixn64 m0, m1, m2;
// // 	usfixn64 R = (usfixn64) r;

// // 	l0 = lhc_0[0];
// // 	h0 = lhc_0[1];
// // 	c0 = lhc_0[2];
// // //	printf("h0=%lu\n", h0);

// // 	l1 = lhc_1[0];
// // 	h1 = lhc_1[1];
// // 	c1 = lhc_1[2];

// // 	m0 = 0;
// // 	m1 = 0;
// // 	m2 = 0;

// // 	//	[s0,s1,s2]=[0,0,0]
// // 	s0 = 0;
// // 	s1 = 0;
// // 	s2 = 0;

// // 	s = (l0 + l1);
// // 	c = (s < l0) | (s < l1) | (s >= R);
// // 	s -= ((0 - (s >= R))) & R;
// // 	s0 = s;
// // 	m0 = c;

// // //	printf("for l0 + l1 \n");
// // //	printf("s=%lu\n", s);
// // //	printf("c=%lu\n", c);
// // //	##############################
// // 	s = (h0 + h1);
// // 	c = (s < h0) | (s < h1) | (s >= R);
// // //	printf("h0=%lu\n", h0);
// // //	printf("h1=%lu\n", h1);
// // //	printf("s=%lu\n", s);
// // //	printf("step 0 : S>=R = %d\n", (s >= R));
// // 	s = s - (((0 - (s >= R))) & R);
// // 	s1 = s;
// // 	m1 = c;

// // //	printf("for h0 +h1 \n");
// // //	printf("s=%lu\n", s);
// // //	printf("c=%lu\n", c);
// // //	##############################
// // 	s = (c0 + c1);
// // 	c = (s < c0) | (s < c1) | (s >= R);
// // 	s -= ((0 - (s >= R))) & R;
// // 	s2 = s;
// // 	m2 = c;

// // //	printf("for c0+c1 \n");
// // //	printf("s=%lu\n", s);
// // //	printf("c=%lu\n", c);
// // //	##############################
// // 	s = (s1 + m0);
// // 	c = (s < s1) | (s < m0) | (s >= R);
// // 	s -= ((0 - (s >= R))) & R;
// // 	s1 = s;
// // 	m1 += c;

// // //	printf("carry handling 1 \n");
// // //	printf("s=%lu\n", s);
// // //	printf("c=%lu\n", c);
// // //	##############################
// // 	s = (s2 + m1);
// // 	c = (s < s2) | (s < m1) | (s >= R);
// // 	s -= ((0 - (s >= R))) & R;
// // 	s2 = s;
// // 	m2 += c;

// // //	printf("carry handling 2 \n");
// // //	printf("s=%lu\n", s);
// // //	printf("c=%lu\n", c);

// // 	lhc_ans[0] = s0;
// // 	lhc_ans[1] = s1;
// // 	lhc_ans[2] = s2;
// // //
// // //	printf("h0=%lu\n", h0);
// // //
// // //	printf("r = %lu\n", r);
// // //	printf("l0 = %lu \n", l0);
// // //	printf("h0 = %lu \n", h0);
// // //	printf("c0 = %lu \n", c0);
// // //	printf("l1 = %lu \n", l1);
// // //	printf("h1 = %lu \n", h1);
// // //	printf("c1 = %lu \n", c1);
// // //	printf("s0 = %lu \n", s0);
// // //	printf("s1 = %lu \n", s1);
// // //	printf("s2 = %lu \n", s2);

// #if VERIFICATION_ENABLED == 1
// 	printf("@add_lhc\n");
// 	if (verify_add_lhc_gmp(lhc_0, lhc_1, s0, s1, s2, r) != 0)
// 	{
// 		printf("FAILED: @add_lhc!\n");
// 		exit(0);
// 	}
// 	else
// 	{
// 		printf("VERIFIED: @add_lhc!\n");
// 	}
// 	printf("=========\n");
// #endif
// }


// /**************************************/

// int check_lhc_gmp(const usfixn64 s0, const usfixn64 s1, const usfixn64 s2,
// 		const usfixn64 r)
// {
// 	mpz_t x_zz, r3;
// 	mpz_t u64_zz;

// 	mpz_init_set_ui(u64_zz, U64_MASK);
// 	mpz_add_ui(u64_zz, u64_zz, 1);

// 	mpz_init_set_ui(x_zz, s2);
// 	mpz_mul(x_zz, x_zz, u64_zz);
// 	mpz_add_ui(x_zz, x_zz, s1);
// 	mpz_mul(x_zz, x_zz, u64_zz);
// 	mpz_add_ui(x_zz, x_zz, s0);

// 	mpz_init_set_ui(r3, r);
// 	mpz_pow_ui(r3, r3, 3);

// 	if (mpz_cmp(x_zz, r3) >= 1)
// 	{
// 		printf("the sum is greater than r^3!\n");
// 		printf("lhc failed.\n");
// 		return -1;
// 		exit(0);
// 	}
// 	else
// 	{
// 		return 0;
// 	}

// }

// /**************************************/

// int verify_lhc_by_R_u128_gmp(const usfixn64 x0, const usfixn64 x1,
// 		const usfixn64 s0, const usfixn64 s1, const usfixn64 s2,
// 		const usfixn64 r)
// {

// //	printf("x0=%lu\n", x0);
// //	printf("x1=%lu\n", x1);
// //	printf("R =%lu\n",r);
// //	printf("s0=%lu\n",s0);
// //	printf("s1=%lu\n",s1);
// //	printf("s2=%lu\n",s2);
// 	mpz_t x_zz;
// 	mpz_t u64_zz;
// 	mpz_t rinv_zz;

// 	mpz_init_set_ui(u64_zz, U64_MASK);
// 	mpz_add_ui(u64_zz, u64_zz, 1);

// 	mpz_init_set_ui(x_zz, x1);
// 	mpz_mul(x_zz, x_zz, u64_zz);
// 	mpz_add_ui(x_zz, x_zz, x0);

// 	mpz_t s0_zz, s1_zz, s2_zz;
// 	mpz_init(s0_zz);
// 	mpz_init(s1_zz);
// 	mpz_init(s2_zz);

// 	mpz_tdiv_qr_ui(s2_zz, s0_zz, x_zz, r);
// 	mpz_tdiv_qr_ui(s2_zz, s1_zz, s2_zz, r);

// 	if (mpz_get_ui(s0_zz) != s0)
// 	{
// 		printf("lhc fail at s0!\n");
// 		printf("s0     =%lu\n", s0);
// 		printf("s0[gmp]=%lu\n", mpz_get_ui(s0_zz));
// 		return -1;
// 		exit(0);
// 	}
// 	else
// 	{
// 		printf("verified s0!\n");
// 	}

// 	if (mpz_get_ui(s1_zz) != s1)
// 	{
// 		printf("lhc fail at s1!\n");
// 		printf("s1     =%lu\n", s1);
// 		printf("s1[gmp]=%lu\n", mpz_get_ui(s1_zz));
// 		return -1;
// 		exit(0);
// 	}
// 	else
// 	{
// 		printf("verified s1!\n");
// 	}

// 	if (mpz_get_ui(s2_zz) != s2)
// 	{
// 		printf("lhc fail at s2!\n");
// 		printf("s2     =%lu\n", s2);
// 		printf("s2[gmp]=%lu\n", mpz_get_ui(s2_zz));
// 		return -1;
// 		exit(0);
// 	}
// 	else
// 	{
// 		printf("verified s2!\n");
// 	}
// 	return 0;
// }

// /**************************************/
// inline void lhc_by_R_u128(const usfixn64 x0, const usfixn64 x1,
// 		const usfixn64 r0, const usfixn64 r1, const usfixn64 u64_mod_R_q,
// 		const usfixn64 u64_mod_R_m, usfixn64* s0, usfixn64* s1, usfixn64* s2,
// 		usfixn64 R)
// {
// //def div_by_R_u128(x0,x1,x2=0,v=1):
// //
// //	usfixn64 x2 = 0;
// 	usfixn64 qb, mb;
// 	usfixn64 m0, q0;
// 	usfixn64 m1, q1;
// //	### should be precomputed
// //	[mb,qb]=div_by_const_R(u64);
// //	div_by_const_R(0, 1, r0, r1, qb, mb);

// 	qb = u64_mod_R_q;
// 	mb = u64_mod_R_m;
// //	[m0,q0]=div_by_const_R(x0);
// 	div_by_const_R(x0, 0, r0, r1, &q0, &m0, R);

// //	# q1=x1/R;s
// //	# m1=x1%R;
// //	[m1,q1]=div_by_const_R(x1);
// 	div_by_const_R(x1, 0, r0, r1, &q1, &m1, R);

// //	__int128 l0, l1;
// //	__int128 h0, h1;
// //	__int128 c0, c1;

// 	usfixn64 l0_lo;
// 	usfixn64 l1_lo, l1_hi;
// 	usfixn64 h0_lo;
// 	usfixn64 h1_lo, h1_hi;
// 	usfixn64 c0_lo, c0_hi;
// 	usfixn64 c1_lo, c1_hi;

// //	l0 = m0;
// //	l1 = (__int128) m1 * (__int128) mb;
// 	l0_lo = m0;
// 	mult_u64_u64(m1, mb, &l1_lo, &l1_hi);

// //	#l2=0;#x2*R0_u128;

// //	h0 = (__int128) q0;
// //	h1 = (__int128) q1 * (__int128) mb + (__int128) qb * (__int128) m1;

// //	h0 = (__int128) q0;
// //	h1 = (__int128) q1 * (__int128) mb + (__int128) qb * (__int128) m1;
// 	h0_lo = q0;

// 	usfixn64 t0, t1, t2, t3, t4;

// 	t0 = 0;
// 	t1 = 0;
// 	t2 = 0;
// 	t3 = 0;
// 	t4 = 0;
// 	mult_u64_u64(q1, mb, &t0, &t1);
// 	mult_u64_u64(qb, m1, &t2, &t3);
// 	add_u128_u128(t0, t1, t2, t3, &h1_lo, &h1_hi, &t4);

// #if VERIFICATION_ENABLED == 1
// 	if (t4 != 0)
// 	{
// 		printf("FAILED: @lhc_by_R_u128!\n");
// 		printf("t4!=0 in add_u128_u128(t0, t1, t2, t3, h1_lo, h1_hi, t4)\n");
// 		printf("t0=%lu \n", t0);
// 		printf("t1=%lu \n", t1);
// 		printf("t2=%lu \n", t2);
// 		printf("t3=%lu \n", t3);
// 		printf("t4=%lu \n", t4);
// 		exit(0);
// 	}
// #endif

// //	#h2=0;#x2*R1_u128;

// //	printf("t4=%lu\n",t4);
// //	exit(0);
// //	c0 = 0;
// //	c1 = (__int128) q1 * (__int128) qb;
// 	c0_lo = 0;
// //	c1 = (__int128) q1 * (__int128) qb;
// 	mult_u64_u64(q1, qb, &c1_lo, &c1_hi);
// //	#c2=0;#x2*R2_u128

// 	usfixn64 lhc_l0h0c0[3] =
// 	{ l0_lo, h0_lo, c0_lo };
// //	lhc_l0h0c0=[l0,h0,c0];

// 	usfixn64 lhc_l1c1[3] =
// 	{ 0, 0, c1_hi };
// 	usfixn64 lhc_l2c2[3] =
// 	{ 0, 0, 0 };
// 	usfixn64 lhc_h1h2[3] =
// 	{ 0, 0, 0 };
// 	usfixn64 lhc_ans[3] =
// 	{ 0, 0, 0 };

// //	printf("first div by R\n");
// //	div_by_const_R(l1 & U64_MASK, l1 >> 64, r0, r1, lhc_l1c1[1], lhc_l1c1[0],
// //			R);
// //	div_by_const_R(l1_lo, l1_hi, r0, r1, lhc_l1c1[1], lhc_l1c1[0], R);

// 	div_by_const_R(l1_lo, l1_hi, r0, r1, &t1, &t0, R);
// 	lhc_l1c1[1] = t1;
// 	lhc_l1c1[0] = t0;
// //	printf("second div by R\n");
// //	div_by_const_R(h1 & U64_MASK, h1 >> 64, r0, r1, lhc_h1h2[2], lhc_h1h2[1],
// //			R);

// //	printf("h1_lo=%lu\n",h1_lo);
// //	printf("h1_hi=%lu\n",h1_hi);
// //	div_by_const_R(h1_lo, h1_hi, r0, r1, lhc_h1h2[2], lhc_h1h2[1], R);
// 	div_by_const_R(h1_lo, h1_hi, r0, r1, &t1, &t0, R);
// 	lhc_h1h2[2] = t1;
// 	lhc_h1h2[1] = t0;

// //	lhc_ans=add_lhc(lhc_l0h0c0,lhc_l1c1)
// //	lhc_ans=add_lhc(lhc_ans,lhc_l2c2)

// 	usfixn64 lhc_ans_tmp[3] =
// 	{ 0, 0, 0 };
// 	add_lhc(lhc_l0h0c0, lhc_l1c1, lhc_ans_tmp, R);
// 	add_lhc(lhc_ans_tmp, lhc_h1h2, lhc_ans, R);

// //	printf("l=%lu, h=%lu, c=%lu\n", lhc_ans[0], lhc_ans[1], lhc_ans[2]);

// 	*s0 = lhc_ans[0];
// 	*s1 = lhc_ans[1];
// 	*s2 = lhc_ans[2];

// #if VERIFICATION_ENABLED == 1
// 	if (check_lhc_gmp(*s0, *s1, *s2, R) != 0)
// 	{
// 		printf("FAILED: @lhc_by_R_u128\n");
// 		return;
// 	}

// 	if (verify_lhc_by_R_u128_gmp(x0, x1, *s0, *s1, *s2, R) != 0)
// 	{
// 		printf("FAILED: @lhc_by_R_u128\n");
// 		exit(0);
// 		return;
// 	}
// 	else
// 	{
// 		printf("VERIFIED: @lhc_by_R_u128\n");
// 	}
// 	printf("=========\n");
// #endif

// }


// /**************************************/
// #define ARITHMETIC_CACHE_SIZE 4

// void Add(sfixn64 *x, sfixn64 *y, int k, sfixn64 r)
// {

// 	short c = 0;
// 	short post = 0;
// 	sfixn64 sum = 0;
// 	int i = 0;

// 	for (int b = 0; b < k; b += ARITHMETIC_CACHE_SIZE)
// //	for (i = 0; i < k; i++)
// 		for (i = b; i < b + ARITHMETIC_CACHE_SIZE; i++)
// 		{
// 			sum = x[i] + y[i] + c;

// 			if (sum < -r)
// 			{
// 				c = -2;
// 				x[i] = sum + 2 * r;
// 			}
// 			else if (sum < 0)
// 			{
// 				c = -1;
// 				x[i] = sum + r;
// 			}
// 			else if (sum >= r)
// 			{
// 				x[i] = sum - r;
// 				c = 1;
// 			}
// 			else
// 			{
// 				x[i] = sum;
// 				c = 0;
// 			}
// //		cout << x[i]  << " " << c << endl;
// 		}

// 	if (c == 1)
// 	{
// 		post = -1;

// 		for (i = 0; i < k; i++)
// 		{
// 			if (x[i] != 0)
// 			{
// 				post = i;
// 				break;
// 			}
// 		}

// 		if (post >= 0)
// 		{
// 			for (i = 0; i < post; i++)
// 			{
// 				x[i] = r - 1;
// 			}
// 			x[post]--;
// 		}
// 		else
// 		{
// 			x[k - 1] = r;

// //			for (i = 0; i < k - 1; i++)
// //			{
// //				x[i] = 0;
// //			}
// 			memset(x, 0x00, (k - 1) * sizeof(sfixn64));
// 		}
// 	}
// 	else if (c < 0)
// 	{
// 		x[0] -= c;
// 		if (x[0] > r - 1)
// 		{
// 			x[0] -= r;
// 			post = -1;

// 			for (i = 1; i < k; i++)
// 			{
// 				if (x[i] < r - 1)
// 				{
// 					post = i;
// 					break;
// 				}
// 			}
// 			if (post > 0)
// 			{
// //				for (i = 1; i < post; i++)
// //				{
// //					x[i] = 0;
// //				}
// 				memset(&x[1], 0x00, (post - 1) * sizeof(sfixn64));
// 				x[post]++;
// 			}
// 			else
// 			{
// 				if (x[0] = 1)
// 					x[0]--;
// 				else
// 				{
// 					x[0] = r - 1;
// 					x[1]--;
// 				}
// 			}
// 		}

// 	}
// //	return x;
// }

// /**************************************/

// inline void  oneShiftRight(sfixn64 * xs, int  k)
// {
// 	sfixn64 tmp;
// 	tmp = xs[k - 1];

// 	for (int b = k - 1; b > 0; b -= 4)
// //	for (short i = k - 1; i > 0; i--)
// 		for (int i = b; i > b - 4; i--)
// 		{
// 			xs[i] = xs[i - 1];
// 			if (i == 0)
// 				break;
// 		}
// 	xs[0] = -tmp;
// }

// /**************************************/

// inline void  twoShiftRight(sfixn64 * xs, int  k)
// {
// 	sfixn64 tmp0, tmp1;
// 	tmp1 = xs[k - 1];
// 	tmp0 = xs[k - 2];

// //	for (short i = k - 1; i > 1; i--)
// 	for (int b = k - 1; b > 0; b -= 4)
// 		//	for (short i = k - 1; i > 0; i--)
// 		for (int i = b; i > b - 4; i--)
// 		{
// 			xs[i] = xs[i - 2];
// 			if (i == 0)
// 				break;
// 		}
// 	xs[0] = -tmp0;
// 	xs[1] = -tmp1;
// }
// /**************************************/


// //x=x+y
// //addtion for GFPF
// void addition_big_elements(sfixn64 * x, sfixn64 *y, const int k,
// 		const sfixn64 r)
// {

// 	short c = 0;
// 	short post = 0;
// 	sfixn64 sum = 0;
// 	int i = 0;

// 	for (i = 0; i < k; i++)
// 	{
// 		sum = x[i] + y[i] + c;

// 		if (sum >= r)
// 		{
// 			c = 1;
// 			x[i] = sum - r;
// 		}
// 		else
// 		{
// 			x[i] = sum;
// 			c = 0;
// 		}
// 	}

// 	if (c > 0)
// 	{
// 		post = -1;

// 		for (i = 0; i < k; i++)
// 		{
// 			if (x[i] != 0)
// 			{
// 				post = i;
// 				break;
// 			}
// 		}

// 		if (post >= 0)
// 		{
// 			for (i = 0; i < post; i++)
// 			{
// 				x[i] = r - 1;
// 			}
// 			x[post]--;
// 		}
// 		else
// 		{
// 			x[k - 1] = r;

// 			for (i = 0; i < k - 1; i++)
// 			{
// 				x[i] = 0;
// 			}
// 		}
// 	}
// }

// /**************************************/

// //x=x-y
// // subtraction for GFPF
// void subtraction_big_elements(sfixn64 *x, sfixn64 *y, const int k,
// 		const sfixn64 r)
// {

// 	int c = 0;
// 	int post = 0;
// 	sfixn64 sub = 0;
// 	int i = 0;
// 	for (i = 0; i < k; i++)
// 	{
// 		sub = y[i] + c;

// 		if (x[i] < sub)
// 		{
// 			c = 1;
// 			x[i] = r - sub + x[i];
// 		}
// 		else
// 		{
// 			c = 0;
// 			x[i] = x[i] - sub;
// 		}
// 	}

// 	if (c > 0)
// 	{
// 		post = -1;
// 		for (i = 0; i < k; i++)
// 		{
// 			if (x[i] < (r - 1))
// 			{
// 				post = i;
// 				break;
// 			}
// 		}

// 		if (post >= 0)
// 		{
// 			for (i = 0; i < post; i++)
// 			{
// 				x[i] = 0;
// 			}

// 			x[post]++;
// 		}
// 		else
// 		{
// 			x[k - 1] = r;

// 			for (i = 0; i < k - 1; i++)
// 			{
// 				x[i] = 0;
// 			}
// 		}
// 	}
// }

// /**************************************/

// //x=x*(r^s)
// void  mult_pow_R(sfixn64 *x, int s, const int k, const sfixn64 r)
// {
// #if PROFILING_ENABLED == 1
// 	inc_profiling_counter(n_mult_pow_R_called);
// #endif
// 	int n_bytes = k * sizeof(sfixn64);
// 	sfixn64 *a = (sfixn64*) malloc(sizeof(sfixn64)* k);
// 	sfixn64 *b = (sfixn64*) malloc(sizeof(sfixn64)* k);
// 	sfixn64 *c = (sfixn64*) malloc(sizeof(sfixn64)* k);

// 	// sfixn64 *a = global_tmp_a;
// 	// sfixn64 *b = global_tmp_b;
// 	// sfixn64 *c = global_tmp_c;
// 	memset(a, 0x0, n_bytes);
// 	memset(b, 0x0, n_bytes);
// 	memset(c, 0x0, n_bytes);
// //	s = s % (2 * k);
// 	s = s & ((k << 1) - 1);

// //	if ((s == 0) || (s > (k << 1)))
// 	if ((s == 0))
// 	{
// 		return;
// 	}
// //	if
// //		{
// //			printf("FAILED!!");
// //					return;
// //		}
// 	else if (s == k)
// 	{
// 		subtraction_big_elements(c, x, k, r);
// 		memcpy(x, c, n_bytes);
// //		free(a);
// //		free(b);
// //		free(c);
// 		return;
// 	}
// 	else if ((s > k) && (s < (2 * k)))
// //	else if ((s > k) && (s < (k << 1)))
// 	{
// 		s = s - k;
// //		mult_pow_R(x, s, k ,r);
// 		subtraction_big_elements(c, x, k, r);
// 		memcpy(x, c, n_bytes);
// 	}

// 	for (int i = 0; i < (k - s); i++)
// 	{
// 		b[i + s] = x[i];
// 	}
// 	for (int i = k - s; i < k; i++)
// 	{
// 		a[i - (k - s)] = x[i];
// 	}
// 	if (x[k - 1] == r)
// 	{
// 		a[s - 1] -= r;
// 		a[s]++;
// 	}
// 	subtraction_big_elements(b, a, k, r);
// 	memcpy(x, b, n_bytes);
// //	free(a);
// //	free(b);
// //	free(c);
// 	return;
// }

// /**************************************/

// void plain_mult_gmp_big_elements(sfixn64 *x, sfixn64 *y, const int k,
// 		const sfixn64 r, const mpz_t prime)
// {

// //	printf("in plain mult\n");
// 	// ToDo:
// 	// verification required.
// 	// also, compare against the c++ version
// 	// important to note how the input is converted to gmp

// //	sfixn64 *res = (sfixn64*) malloc(k * sizeof(sfixn64));

// //	mpz_class X = number(x, k, r);
// //	mpz_class Y = number(y, k, r);
// 	mpz_t X, Y;
// 	mpz_init(X);
// 	mpz_init(Y);
// 	radix_based_to_mpz_s64(X, x, r, k);
// 	radix_based_to_mpz_s64(Y, y, r, k);

// //	char xstr[2048], ystr[2048];
// //	mpz_get_str(xstr, 10, X);
// //	mpz_get_str(ystr, 10, Y);
// //
// //	printf ("xstr= %s\n", xstr);
// //	printf ("ystr= %s\n", ystr);
// //	printf("==================\n");

// //	X = (X * Y) % prime;
// 	mpz_mul(X, X, Y);
// 	mpz_mod(X, X, prime);

// //	TODO: proper change of representation to be done;
// //	ConvertToGFPF(res, X, prime, k, r);
// 	mpz_to_radix_based_s64(x, X, r, k);
// //	return res;
// //	mempcpy(x, res, k*sizeof(sfixn64));
// //	free(res);
// 	mpz_clear(X);
// 	mpz_clear(Y);

// }

// /**************************************/

// void test_mult_u64_u64()
// {
// //	usfixn64 x = 10;
// //	usfixn64 y = 20;
// 	usfixn64 s0, s1;
// 	s0 = 0;
// 	s1 = 0;

// 	for (usfixn64 x = U64_MASK; x > U64_MASK - 1024; x--)
// 		for (usfixn64 y = U64_MASK; y > U64_MASK - 1024; y--)
// 			mult_u64_u64(x, y, &s0, &s1);
// }

// /**************************************/

// void test_div_by_const_R()
// {

// //	usfixn64 x, y;
// //	x = 10;
// //	y = 20;
// 	usfixn64 q, m;
// 	for (usfixn64 x = U64_MASK; x > U64_MASK - 1024; x--)
// 		for (usfixn64 y = U64_MASK; y > U64_MASK - 1024; y--)
// 			div_by_const_R(x, y, r_inv_0_16, r_inv_1_16, &q, &m, R16);
// }

// void convolution8_step1(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{
// 		x[i] = MulModMont(x[i], theta_list[i], MY_PRIME,
// 				INV_PRIME);
// 	}

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{
// 		y[i] = MulModMont(y[i], theta_list[i], MY_PRIME,
// 				INV_PRIME);
// 	}

// 	//	GFPF_DFT_8(&x[0], w_conv8_step1, MY_PRIME, INV_PRIME);
// 	//	GFPF_DFT_8(&y[0], w_conv8_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_8(&x[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_8(&y[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{
// 		x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 				INV_PRIME);
// 	}
// 	//	GFPF_DFT_8(&x[0], winv_conv8_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_8(&x[0], vec_pow_omega_inv_conv_step1, MY_PRIME, INV_PRIME);

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{
// 		x[i] = MulModMont(x[i], thetainv_list[i],
// 				MY_PRIME, INV_PRIME);
// 		x[i] = MulModMont(x[i], ninv_conv8_step1,
// 				MY_PRIME, INV_PRIME);
// 	}
// 	//	return x;
// }

// /**************************************/

// void convolution8_step2(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{
// 		x[i] = MulModMont(x[i],
// 				theta_list_conv8_step1[i], MY_PRIME, INV_PRIME);
// 	}

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{

// 		y[i] = MulModMont(y[i],
// 				theta_list_conv8_step1[i], MY_PRIME, INV_PRIME);
// 	}

// 	//	GFPF_DFT_8(&x[0], w_conv8_step2, MY_PRIME, INV_PRIME);
// 	//	GFPF_DFT_8(&y[0], w_conv8_step2, MY_PRIME, INV_PRIME);

// 	GFPF_DFT_8(&x[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_8(&y[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{
// 		x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 				INV_PRIME);
// 	}
// 	//	GFPF_DFT_8(&x[0], winv_conv8_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_8(&x[0], vec_pow_omega_inv_conv_step2, MY_PRIME, INV_PRIME);

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 8; ++i)
// 	{
// 		x[i] = MulModMont(x[i],
// 				thetainv_list_conv8_step1[i], MY_PRIME, INV_PRIME);
// 		x[i] = MulModMont(x[i], ninv_conv8_step2,
// 				MY_PRIME, INV_PRIME);
// 	}
// 	//	return x;
// }

// /**************************************/

// void convolution16_step1(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// 	//	printf ("myprime = %lu\n", MY_PRIME);
// 	//	printf ("inv_prime = %lu\n", INV_PRIME);
// 	//	printf("before first mont !\n");
// 	//	for (int i = 0; i < 16; i++)
// 	//	{
// 	//		printf("x[%d]=%lu\n", i, x[i]);
// 	//	}

// //	CALLGRIND_START_INSTRUMENTATION
// //	;
// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 16; ++i)
// //	{
// //		x[i] = MulModMont(x[i],
// //				theta_list_conv16_step1[i], MY_PRIME, INV_PRIME);
// //	}

// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 			x[i] = MulModMont(x[i],
// 					theta_list_conv16_step1[i], MY_PRIME, INV_PRIME);

// //	for (int i = 0; i < 16; ++i)
// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			y[i] = MulModMont(y[i],
// 					theta_list_conv16_step1[i], MY_PRIME, INV_PRIME);
// 		}

// 	//	GFPF_DFT_16(&x[0], w_conv16_step1, MY_PRIME, INV_PRIME);
// 	//	GFPF_DFT_16(&y[0], w_conv16_step1, MY_PRIME, INV_PRIME);

// 	//	printf("before GFPF_DFT_16!\n");
// 	//	for (int i = 0; i < 16; i++)
// 	//	{
// 	//		printf("x[%d]=%lu\n", i, x[i]);
// 	//	}

// 	GFPF_DFT_16(&x[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);

// 	GFPF_DFT_16(&y[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);
// 	//	printf("after GFPF_DFT_16!\n");
// 	//	for (int i = 0; i < 16; i++)
// 	//	{
// 	//		printf("x[%d]=%lu\n", i, x[i]);
// 	//	}

// //	for (int i = 0; i < 16; ++i)
// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 					INV_PRIME);
// 		}

// 	//	printf("after mont!\n");
// 	//	for (int i = 0; i < 16; i++)
// 	//	{
// 	//		printf("x[%d]=%lu\n", i, x[i]);
// 	//	}
// 	//	GFPF_DFT_16(&x[0], winv_conv16_step1, MY_PRIME, INV_PRIME);

// 	GFPF_DFT_16(&x[0], vec_pow_omega_inv_conv_step1, MY_PRIME, INV_PRIME);

// 	//
// 	//	printf("after GFPF_DFT_16back!\n");
// 	//	for (int i = 0; i < 16; i++)
// 	//	{
// 	//		printf("x[%d]=%lu\n", i, x[i]);
// 	//	}

// //	CALLGRIND_START_INSTRUMENTATION
// //		;

// //	for (int i = 0; i < 16; ++i)
// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i],
// 					thetainv_list_conv16_step1[i], MY_PRIME, INV_PRIME);
// 			x[i] = MulModMont(x[i], ninv_conv16_step1,
// 					MY_PRIME, INV_PRIME);
// 		}

// //	CALLGRIND_STOP_INSTRUMENTATION
// //		;
// //		CALLGRIND_DUMP_STATS
// //		;

// 	//	printf("after mont mont !\n");
// 	//	for (int i = 0; i < 16; i++)
// 	//	{
// 	//		printf("x[%d]=%lu\n", i, x[i]);
// 	//	}
// 	//	return x;
// }

// /**************************************/

// void convolution16_step2(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// 	//	printf ("s2 - myprime = %lu\n", MY_PRIME);
// 	//	printf ("s2 -inv_prime = %lu\n", INV_PRIME);
// //	for (int i = 0; i < 16; ++i)
// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED

// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i],
// 					theta_list_conv16_step2[i], MY_PRIME, INV_PRIME);
// 		}

// //	for (int i = 0; i < 16; ++i)
// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED

// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			y[i] = MulModMont(y[i],
// 					theta_list_conv16_step2[i], MY_PRIME, INV_PRIME);
// 		}

// 	//	GFPF_DFT_16(&x[0], w_conv16_step2, MY_PRIME, INV_PRIME);
// 	//	GFPF_DFT_16(&y[0], w_conv16_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_16(&x[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_16(&y[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);

// //	for (int i = 0; i < 16; ++i)
// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 					INV_PRIME);
// 		}
// 	//	GFPF_DFT_16(&x[0], winv_conv16_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_16(&x[0], vec_pow_omega_inv_conv_step2, MY_PRIME, INV_PRIME);

// //	for (int i = 0; i < 16; ++i)
// 	for (int b = 0; b < 16; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i],
// 					thetainv_list_conv16_step2[i], MY_PRIME, INV_PRIME);
// 			x[i] = MulModMont(x[i], ninv_conv16_step2,
// 					MY_PRIME, INV_PRIME);
// 		}
// 	//	return x;
// }

// /**************************************/

// //sfixn64*
// void convolution32_step1(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// 	//    for (int i = 0; i < 32; ++i)
// 	// {
// 	// 	theta_list[i] = convertToMontMulR(theta_list[i],RP,MY_PRIME);
// 	// 	thetainv_list[i] = convertToMontMulR(thetainv_list[i],RP,MY_PRIME);
// 	// //	cout << x[i] << " ";
// 	// }

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i],
// 					theta_list_conv32_step1[i], MY_PRIME, INV_PRIME);
// 		}

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			y[i] = MulModMont(y[i],
// 					theta_list_conv32_step1[i], MY_PRIME, INV_PRIME);
// 		}
// 	//w = convertToMontMulR(w,RP,MY_PRIME);

// 	//	GFPF_DFT_32(&x[0], w_conv32_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_32(&x[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout << convertOut(x[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;
// 	//	GFPF_DFT_32(&y[0], w_conv32_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_32(&y[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout <<  convertOut(y[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;
// //
// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 					INV_PRIME);
// 		}

// 	//	winv = convertToMontMulR(winv,RP,MY_PRIME);
// 	//	GFPF_DFT_32(&x[0], winv_conv32_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_32(&x[0], vec_pow_omega_inv_conv_step1, MY_PRIME, INV_PRIME);

// 	//	ninv = convertToMontMulR(ninv,RP,MY_PRIME);

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i],
// 					thetainv_list_conv32_step1[i], MY_PRIME, INV_PRIME);
// 			x[i] = MulModMont(x[i], ninv_conv32_step1,
// 					MY_PRIME, INV_PRIME);
// 		}
// 	//	return x;
// }

// /**************************************/

// //sfixn64*
// void convolution32_step2(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// 	// {
// 	// 	theta_list[i] = convertToMontMulR(theta_list[i],RP,MY_PRIME);
// 	// 	thetainv_list[i] = convertToMontMulR(thetainv_list[i],RP,MY_PRIME);
// 	// //	cout << x[i] << " ";
// 	// }

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i],
// 					theta_list_conv32_step2[i], MY_PRIME, INV_PRIME);
// 		}

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			y[i] = MulModMont(y[i],
// 					theta_list_conv32_step2[i], MY_PRIME, INV_PRIME);
// 		}
// 	//	w = convertToMontMulR(w,RP,MY_PRIME);
// 	//cout << convertToMontMulR(4) << endl;
// 	//	GFPF_DFT_32(&x[0], w_conv32_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_32(&x[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout << convertOut(x[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;
// 	//	GFPF_DFT_32(&y[0], w_conv32_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_32(&y[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout <<  convertOut(y[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 					INV_PRIME);
// 		}

// 	//	winv = convertToMontMulR(winv,RP,MY_PRIME);
// 	//	GFPF_DFT_32(&x[0], winv_conv32_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_32(&x[0], vec_pow_omega_inv_conv_step2, MY_PRIME, INV_PRIME);

// 	//	ninv = convertToMontMulR(ninv,RP,MY_PRIME);

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < 32; ++i)
// 	for (int b = 0; b < 32; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			x[i] = MulModMont(x[i],
// 					thetainv_list_conv32_step2[i], MY_PRIME, INV_PRIME);
// 			x[i] = MulModMont(x[i], ninv_conv32_step2,
// 					MY_PRIME, INV_PRIME);
// 		}
// 	//	return x;
// }

// /**************************************/

// //sfixn64*
// void convolution64_step1(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// 	//    for (int i = 0; i < 32; ++i)
// 	// {
// 	// 	theta_list[i] = convertToMontMulR(theta_list[i],RP,MY_PRIME);
// 	// 	thetainv_list[i] = convertToMontMulR(thetainv_list[i],RP,MY_PRIME);
// 	// //	cout << x[i] << " ";
// 	// }

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		x[i] = MulModMont(x[i],
// 				theta_list_conv64_step1[i], MY_PRIME, INV_PRIME);
// 	}

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		y[i] = MulModMont(y[i],
// 				theta_list_conv64_step1[i], MY_PRIME, INV_PRIME);
// 	}
// 	//w = convertToMontMulR(w,RP,MY_PRIME);

// 	//	GFPF_DFT_64(&x[0], w_conv64_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_64(&x[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout << convertOut(x[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;

// 	//	GFPF_DFT_64(&y[0], w_conv64_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_64(&y[0], vec_pow_omega_conv_step1, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout <<  convertOut(y[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 				INV_PRIME);
// 	}

// 	//	winv = convertToMontMulR(winv,RP,MY_PRIME);
// 	//	GFPF_DFT_64(&x[0], winv_conv64_step1, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_64(&x[0], vec_pow_omega_inv_conv_step1, MY_PRIME, INV_PRIME);

// 	//	ninv = convertToMontMulR(ninv,RP,MY_PRIME);

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		x[i] = MulModMont(x[i],
// 				thetainv_list_conv64_step1[i], MY_PRIME, INV_PRIME);
// 		x[i] = MulModMont(x[i], ninv_conv64_step1,
// 				MY_PRIME, INV_PRIME);
// 	}
// 	//	return x;
// }

// /**************************************/

// //sfixn64*
// void convolution64_step2(sfixn64 *x, sfixn64 *y, sfixn64 MY_PRIME,
// 		sfixn64 INV_PRIME)
// {

// 	//2447143447522443265; // 1/32 mod prime
// 	//    for (int i = 0; i < 32; ++i)
// 	// {
// 	// 	theta_list[i] = convertToMontMulR(theta_list[i],RP,MY_PRIME);
// 	// 	thetainv_list[i] = convertToMontMulR(thetainv_list[i],RP,MY_PRIME);
// 	// //	cout << x[i] << " ";
// 	// }

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		x[i] = MulModMont(x[i],
// 				theta_list_conv64_step2[i], MY_PRIME, INV_PRIME);
// 	}

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		y[i] = MulModMont(y[i],
// 				theta_list_conv64_step2[i], MY_PRIME, INV_PRIME);
// 	}
// 	//w = convertToMontMulR(w,RP,MY_PRIME);

// 	//	GFPF_DFT_64(&x[0], w_conv64_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_64(&x[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout << convertOut(x[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;
// 	//	GFPF_DFT_64(&y[0], w_conv64_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_64(&y[0], vec_pow_omega_conv_step2, MY_PRIME, INV_PRIME);
// 	// for (int i = 0; i < 8; ++i)
// 	// {
// 	// 	cout <<  convertOut(y[i],MY_PRIME,INV_PRIME) << " ";
// 	// }
// 	// cout << endl;
// 	// cout << endl;

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		x[i] = MulModMont(x[i], y[i], MY_PRIME,
// 				INV_PRIME);
// 	}

// 	//	winv = convertToMontMulR(winv,RP,MY_PRIME);
// 	//	GFPF_DFT_64(&x[0], winv_conv64_step2, MY_PRIME, INV_PRIME);
// 	GFPF_DFT_64(&x[0], vec_pow_omega_inv_conv_step2, MY_PRIME, INV_PRIME);

// 	//	ninv = convertToMontMulR(ninv,RP,MY_PRIME);

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int i = 0; i < 64; ++i)
// 	{
// 		x[i] = MulModMont(x[i],
// 				thetainv_list_conv64_step2[i], MY_PRIME, INV_PRIME);
// 		x[i] = MulModMont(x[i], ninv_conv64_step2,
// 				MY_PRIME, INV_PRIME);
// 	}
// 	//	return x;
// }

// /**************************************/
// /**************************************/

// void init_gfpf_mult_data(crt_u192_data* t_crt_data, usfixn64 r_inv_0,
// 		usfixn64 r_inv_1, usfixn64 r)
// {

// 	usfixn64 p2 = 4179340454199820289;
// 	usfixn64 m2 = 740506764262110493;

// 	usfixn64 p1 = 2485986994308513793;
// 	usfixn64 m1 = 1244909922527606046;

// 	u64_mod_u64(&m2, &p1);
// 	m2 = p1 - m2;

// 	//m1 == m1%p2
// 	u64_mod_u64(&m1, &p2);

// 	usfixn64 u64_mod_p1, u64_mod_p2;
// 	//u64_mod_p1 = (2^64)%p1
// 	//u64_mod_p2 = (2^64)%p2
// 	u64_mod_p1 = (U64_MASK);
// 	u64_mod_p2 = (U64_MASK);
// 	u64_mod_u64(&u64_mod_p1, &p1);
// 	u64_mod_u64(&u64_mod_p2, &p2);
// 	u64_mod_p1 += 1;
// 	u64_mod_p2 += 1;

// 	usfixn64 p1_inv_q, p1_inv_m, p2_inv_q, p2_inv_m;
// 	inv_p_mod_u128(&p1, &p1_inv_m, &p1_inv_q);
// 	inv_p_mod_u128(&p2, &p2_inv_m, &p2_inv_q);
// 	//	__int128 p1p2 = (__int128) p1 * (__int128) p2;
// 	t_crt_data->m1 = m1;
// 	t_crt_data->m2 = m2;
// 	t_crt_data->p1 = p1;
// 	t_crt_data->p2 = p2;
//   t_crt_data->m1_mont = convertIn(m1, R_11, CONVOLUTION_PRIME_1,
//           INV_CONVOLUTION_PRIME_1);
//   t_crt_data->m2_mont = convertIn(m2, R_22, CONVOLUTION_PRIME_2,
//           INV_CONVOLUTION_PRIME_2);

// 	t_crt_data->p1_inv_m = p1_inv_m;
// 	t_crt_data->p1_inv_q = p1_inv_q;
// 	t_crt_data->p2_inv_m = p2_inv_m;
// 	t_crt_data->p2_inv_q = p2_inv_q;

// 	//	t_crt_data->p1p2 = p1p2;
// 	//	t_crt_data->p1p2_q = p1p2 >> 64;
// 	//	t_crt_data->p1p2_m = p1p2 & (U64_MASK);

// 	mult_u64_u64(p1, p2, &t_crt_data->p1p2_m, &t_crt_data->p1p2_q);
// 	div_by_const_R(0, 1, r_inv_0, r_inv_1, &t_crt_data->qb, &t_crt_data->mb, r);

// 	t_crt_data->radix = r;
// 	t_crt_data->r_inv_0 = r_inv_0;
// 	t_crt_data->r_inv_1 = r_inv_1;
// }

// void __inline__
// crt_mult_sub_u192_with_reduction_test (usfixn64 *a1, usfixn64 *a2, int k)
// {

//   //// backup gmp implementation.

//   usfixn64 s0, s1;
//   usfixn64 t[2];
//   usfixn64 m0, m1;
//   usfixn64 data_m1, data_m2;
//   const usfixn64 one = 1;

// //  for (int j = 0; j < k; j++)
//   for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
//     for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
//       {
//   t[0] = a1[j];
//   t[1] = a2[j];

//   data_m1 = t_crt_data_global.m1_mont;
//   data_m2 = t_crt_data_global.m2_mont;


// // //a1m1%p1 and a2m2%p2
// //   mult_ab_mod_p_ptr (&t[0], &data_m1, global_P1);
// //   mult_ab_mod_p_ptr (&t[1], &data_m2, global_P2);
//   t[0] = MulModMont(t[0],data_m1,CONVOLUTION_PRIME_1,INV_CONVOLUTION_PRIME_1);
//   t[1] = MulModMont(t[1],data_m2,CONVOLUTION_PRIME_2,INV_CONVOLUTION_PRIME_2);
//   s0 = convertOut(t[0],  CONVOLUTION_PRIME_1,
//           INV_CONVOLUTION_PRIME_1);
//   s1 = convertOut(t[1],  CONVOLUTION_PRIME_2,
//           INV_CONVOLUTION_PRIME_2);


// //   convertOut_GLOBAL_ptr (&s0, global_P1);
// //   convertOut_GLOBAL_ptr (&s1, global_P2);

//   __asm__ __volatile__(
//       "movq  %2, %%rax;\n\t"          // rax = a
// //      "movq  %3, %%rdx;\n\t"// rax = a
//       "mulq  %3;\n\t"// rdx:rax = a * b
//       "movq  %%rax, %0;\n\t"// s0 = rax
//       "movq  %%rdx, %1;\n\t"// s1 = rdx

//       "movq  %4, %%rax;\n\t"// rax = a
// //      "movq  %5, %%rdx;\n\t"// rax = a
//       "mulq  %5;\n\t"// rdx:rax = a * b
// //      "addq  %0, %%rax;\n\t"// s0 = rax
// //      "adcq  %1, %%rdx;\n\t"// s1 = rdx
// //      "movq  %%rax, %0;\n\t"// s0 = rax
// //      "movq  %%rdx, %1;\n\t"// s1 = rdx
//       "addq  %%rax, %0;\n\t"// s0 = rax
//       "adcq  %%rdx, %1;\n\t"// s1 = rdx
//       : "+g" (t[0]),"+g"(t[1])
//       : "g"(s0), "g"(t_crt_data_global.p2),"g"(s1), "g"(t_crt_data_global.p1)
//       : "%rax", "%rdx");


//   if ((t[1] > t_crt_data_global.p1p2_q)
//       || ((t[1] == t_crt_data_global.p1p2_q)
//     && (t[0] > t_crt_data_global.p1p2_m)))
//     {
//       m0 = U64_MASK - t_crt_data_global.p1p2_m;
//       m1 = U64_MASK - t_crt_data_global.p1p2_q;
//       __asm__ __volatile__(
//     "addq %2, %0; \n\t"
//     "adcq %3, %1; \n\t"
//     "addq $0x1, %0; \n\t"
//     "adcq $0x0, %1; \n\t"
//     : "+g" (t[0]),"+g"(t[1])
//     : "g"(m0), "g"(m1)
//     : );

//     }

//   a1[j] = t[0];
//   a2[j] = t[1];

//       }

// }

// /**************************************/

// void  GFPFMultiplication(sfixn64 *x, sfixn64 *y, int k,
// 		const crt_u192_data* t_crt_data)
// {

// 	int is_zero = 1;
// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// //	for (int i = 0; i < k; i++)
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			if (x[i] != 0)
// 				is_zero = 0;
// 		}

// 	if (is_zero == 1)
// 	{
		
// 		return;
// 	}

// 	is_zero = 1;
// #pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int i = 0; i < k; i++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// 		//	for (int i = 0; i < k; i++)
// 		for (int i = b; i < b + CONVOLUTION_CACHE_SIZE; i++)
// 		{
// 			//		is_zero = is_zero | (y[i] == 0);
// 			if (y[i] != 0)
// 				is_zero = 0;
// 		}

// 	if (is_zero == 1)
// 	{
// 		memset(x, 0x00, k * sizeof(sfixn64));
// 		return;
// 	}

	
// 	int n_bytes = k * sizeof(usfixn64);
	
// 	sfixn64* x1, *x2, *y1, *y2, *l_vec, *h_vec, *c_vec;
// 	usfixn64* s0, *s1;
// 	x1 = global_x1;
// 	x2 = global_x2;
// 	y1 = global_y1;
// 	y2 = global_y2;
// 	s0 = global_s0;
// 	s1 = global_s1;
// 	l_vec = global_l_vec;
// 	h_vec = global_h_vec;
// 	c_vec = global_c_vec;

// 	usfixn64 l_temp, h_temp, c_temp;
// 	// r_inv = U128_mask/r;
// 	// r0,r1= r_inv%u64, r_inv/u64;
// 	usfixn64 r_inv_0;
// 	usfixn64 r_inv_1;
// 	sfixn64 r;

// 	r_inv_0 = t_crt_data->r_inv_0;
// 	r_inv_1 = t_crt_data->r_inv_1;
// 	r = t_crt_data->radix;
	
// 	memcpy(x1, x, n_bytes);
// 	memcpy(x2, x, n_bytes);

// 	memcpy(y1, y, n_bytes);
// 	memcpy(y2, y, n_bytes);

// #pragma unroll LOOP_UNROLLING_ENABLED
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 			x2[j] = convertIn(x2[j], R_22, CONVOLUTION_PRIME_2,
// 					INV_CONVOLUTION_PRIME_2);

// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 			x1[j] = convertIn(x1[j], R_11, CONVOLUTION_PRIME_1,
// 					INV_CONVOLUTION_PRIME_1);

// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 			y1[j] = convertIn(y1[j], R_11, CONVOLUTION_PRIME_1,
// 					INV_CONVOLUTION_PRIME_1);

// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 			y2[j] = convertIn(y2[j], R_22, CONVOLUTION_PRIME_2,
// 					INV_CONVOLUTION_PRIME_2);

	

// 	switch (k)
// 	{
// 	case 8:
// 		convolution8_step1(x1, y1, CONVOLUTION_PRIME_1,
// 				INV_CONVOLUTION_PRIME_1);
// 		convolution8_step2(x2, y2, CONVOLUTION_PRIME_2,
// 				INV_CONVOLUTION_PRIME_2);
// 		break;
// 	case 16:
// 		convolution16_step1(x1, y1, CONVOLUTION_PRIME_1,
// 				INV_CONVOLUTION_PRIME_1);
// 		convolution16_step2(x2, y2, CONVOLUTION_PRIME_2,
// 				INV_CONVOLUTION_PRIME_2);
// 		break;   //execution starts at this case label
// 	case 32:
// 		convolution32_step1(x1, y1, CONVOLUTION_PRIME_1,
// 				INV_CONVOLUTION_PRIME_1);
// 		convolution32_step2(x2, y2, CONVOLUTION_PRIME_2,
// 				INV_CONVOLUTION_PRIME_2);
// 		break;
// 	case 64:
// 		convolution64_step1(x1, y1, CONVOLUTION_PRIME_1,
// 				INV_CONVOLUTION_PRIME_1);
// 		convolution64_step2(x2, y2, CONVOLUTION_PRIME_2,
// 				INV_CONVOLUTION_PRIME_2);
// 		break;
// 	}

// 	//	timer_record_stop(t_convolution);
// 	//	timer_record_start(t_conv_out);

// 	//	memcpy(x, x1, k * sizeof(sfixn64));
// 	//	return;

// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 		{
// 			x1[j] = convertOut(x1[j], CONVOLUTION_PRIME_1,
// 					INV_CONVOLUTION_PRIME_1);
// 		}

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 		{
// 			x2[j] = convertOut(x2[j], CONVOLUTION_PRIME_2,
// 					INV_CONVOLUTION_PRIME_2);
// 		}

	
// 	char *post = global_post;
	

// 	usfixn64 qb, mb;
// 	qb = t_crt_data->qb;
// 	mb = t_crt_data->mb;

// 	usfixn64 a0 = 3332663724254167040;
// 	usfixn64 a1 = 281615714199011328;

// 	usfixn64 b0 = 6665327448508334080;
// 	usfixn64 b1 = 563231428398022656;

	

// 	usfixn64 x22 = 0, x11 = 0;

// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 		{
// 			//		post[j] = 1;

// 			x22 = (usfixn64) x2[j];
// 			x11 = (usfixn64) x1[j];
// 			//		crt_mult_sub_u192_with_reduction((usfixn64) x2[j], (usfixn64) x1[j],
// 			//				t_crt_data, s0[j], s1[j]);

// 			crt_mult_sub_u192_with_reduction(x22, x11, t_crt_data, &s0[j],
// 					&s1[j]);

// 			//		crt_mult_sub_u192_with_reduction(x22, x11, t_crt_data, t0 , t1 );
// 		}
// // crt_mult_sub_u192_with_reduction_test ((usfixn64*)x1, (usfixn64*)x2,(usfixn64)k);
// //   for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// // #pragma unroll LOOP_UNROLLING_ENABLED
// //     for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// //     {
// //       post[j] = 0;
// //       if ((x2[j] > a1) || (x2[j] == a1 && x1[j] > a0))
// //       {
      
// //         x2[j] = b1 - x2[j];


// //         if (x1[j] <= b0)
// //         {
// //           //        printf("--- if 1\n");
// //           x1[j] = b0 - x1[j] + 1;
// //         }
// //         else
// //         {
// //           //        printf("--- if 2\n");
// //           x2[j]--;
// //           x1[j] = U64_MASK - (x1[j] - b0 - 1) + 1;
// //         }
// //         post[j] = LHC_NEGATIVE_SIGN;
// //       }
// //       //    printf("===================\n");
// //     }

// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 		{
// 			post[j] = 0;
// 			if ((s1[j] > a1) || (s1[j] == a1 && s0[j] > a0))
// 			{
			
// 				s1[j] = b1 - s1[j];


// 				if (s0[j] <= b0)
// 				{
// 					//				printf("--- if 1\n");
// 					s0[j] = b0 - s0[j] + 1;
// 				}
// 				else
// 				{
// 					//				printf("--- if 2\n");
// 					s1[j]--;
// 					s0[j] = U64_MASK - (s0[j] - b0 - 1) + 1;
// 				}
// 				post[j] = LHC_NEGATIVE_SIGN;
// 			}
// 			//		printf("===================\n");
// 		}


// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// #pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 		{

// 			//		printf("lhc-step [%d]\n", j);
// 			l_temp = 0;
// 			h_temp = 0;
// 			c_temp = 0;
// 			lhc_by_R_u128(s0[j], s1[j], r_inv_0, r_inv_1, qb, mb, &l_temp,
// 					&h_temp, &c_temp, r);
			
// 			l_vec[j] = l_temp;
// 			h_vec[j] = h_temp;
// 			c_vec[j] = c_temp;
// 		}

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// //	#pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 			if (post[j])
// 				l_vec[j] = -l_vec[j];

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// 		//	#pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 			if (post[j])
// 				h_vec[j] = -h_vec[j];

// //#pragma unroll LOOP_UNROLLING_ENABLED
// //	for (int j = 0; j < k; j++)
// 	for (int b = 0; b < k; b += CONVOLUTION_CACHE_SIZE)
// 		//	#pragma unroll LOOP_UNROLLING_ENABLED
// 		for (int j = b; j < b + CONVOLUTION_CACHE_SIZE; j++)
// 			if (post[j])
// 				c_vec[j] = -c_vec[j];

	

// 	oneShiftRight(h_vec, k);
// 	twoShiftRight(c_vec, k);

// //	BEGIN_PROFILE
// 	Add(l_vec, h_vec, k, r);
// 	Add(l_vec, c_vec, k, r);


// 	memcpy(x, l_vec, k * sizeof(sfixn64));

// }
