
#ifndef _PACKED_MON_DIV_DEFS_
#define _PACKED_MON_DIV_DEFS_


typedef int (*DivTest_ptr)(degrees_t, degrees_t);

//monomial divide test for nvar 1
static inline int monomialDivideTest_1(degrees_t adegs, degrees_t bdegs) {
	return adegs >= bdegs;
}

//monomial divide test for nvar 2
static inline int monomialDivideTest_2(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V2) < (bdegs & EXP_1_V2) ||
			   (adegs & EXP_2_V2) < (bdegs & EXP_2_V2)  );
}

//monomial divide test for nvar 3
static inline int monomialDivideTest_3(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V3) < (bdegs & EXP_1_V3) ||
			   (adegs & EXP_2_V3) < (bdegs & EXP_2_V3) ||
			   (adegs & EXP_3_V3) < (bdegs & EXP_3_V3)  );
}

//monomial divide test for nvar 4
static inline int monomialDivideTest_4(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V4) < (bdegs & EXP_1_V4) ||
			   (adegs & EXP_2_V4) < (bdegs & EXP_2_V4) ||
			   (adegs & EXP_3_V4) < (bdegs & EXP_3_V4) ||
			   (adegs & EXP_4_V4) < (bdegs & EXP_4_V4)  );
}

//monomial divide test for nvar 5
static inline int monomialDivideTest_5(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V5) < (bdegs & EXP_1_V5) ||
			   (adegs & EXP_2_V5) < (bdegs & EXP_2_V5) ||
			   (adegs & EXP_3_V5) < (bdegs & EXP_3_V5) ||
			   (adegs & EXP_4_V5) < (bdegs & EXP_4_V5) ||
			   (adegs & EXP_5_V5) < (bdegs & EXP_5_V5)  );
}

//monomial divide test for nvar 6
static inline int monomialDivideTest_6(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V6) < (bdegs & EXP_1_V6) ||
			   (adegs & EXP_2_V6) < (bdegs & EXP_2_V6) ||
			   (adegs & EXP_3_V6) < (bdegs & EXP_3_V6) ||
			   (adegs & EXP_4_V6) < (bdegs & EXP_4_V6) ||
			   (adegs & EXP_5_V6) < (bdegs & EXP_5_V6) ||
			   (adegs & EXP_6_V6) < (bdegs & EXP_6_V6)  );
}

//monomial divide test for nvar 7
static inline int monomialDivideTest_7(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V7) < (bdegs & EXP_1_V7) ||
			   (adegs & EXP_2_V7) < (bdegs & EXP_2_V7) ||
			   (adegs & EXP_3_V7) < (bdegs & EXP_3_V7) ||
			   (adegs & EXP_4_V7) < (bdegs & EXP_4_V7) ||
			   (adegs & EXP_5_V7) < (bdegs & EXP_5_V7) ||
			   (adegs & EXP_6_V7) < (bdegs & EXP_6_V7) ||
			   (adegs & EXP_7_V7) < (bdegs & EXP_7_V7)  );
}

//monomial divide test for nvar 8
static inline int monomialDivideTest_8(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V8) < (bdegs & EXP_1_V8) ||
			   (adegs & EXP_2_V8) < (bdegs & EXP_2_V8) ||
			   (adegs & EXP_3_V8) < (bdegs & EXP_3_V8) ||
			   (adegs & EXP_4_V8) < (bdegs & EXP_4_V8) ||
			   (adegs & EXP_5_V8) < (bdegs & EXP_5_V8) ||
			   (adegs & EXP_6_V8) < (bdegs & EXP_6_V8) ||
			   (adegs & EXP_7_V8) < (bdegs & EXP_7_V8) ||
			   (adegs & EXP_8_V8) < (bdegs & EXP_8_V8)  );
}

//monomial divide test for nvar 9
static inline int monomialDivideTest_9(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V9) < (bdegs & EXP_1_V9) ||
			   (adegs & EXP_2_V9) < (bdegs & EXP_2_V9) ||
			   (adegs & EXP_3_V9) < (bdegs & EXP_3_V9) ||
			   (adegs & EXP_4_V9) < (bdegs & EXP_4_V9) ||
			   (adegs & EXP_5_V9) < (bdegs & EXP_5_V9) ||
			   (adegs & EXP_6_V9) < (bdegs & EXP_6_V9) ||
			   (adegs & EXP_7_V9) < (bdegs & EXP_7_V9) ||
			   (adegs & EXP_8_V9) < (bdegs & EXP_8_V9) ||
			   (adegs & EXP_9_V9) < (bdegs & EXP_9_V9)  );
}

//monomial divide test for nvar 10
static inline int monomialDivideTest_10(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V10) < (bdegs & EXP_1_V10) ||
			   (adegs & EXP_2_V10) < (bdegs & EXP_2_V10) ||
			   (adegs & EXP_3_V10) < (bdegs & EXP_3_V10) ||
			   (adegs & EXP_4_V10) < (bdegs & EXP_4_V10) ||
			   (adegs & EXP_5_V10) < (bdegs & EXP_5_V10) ||
			   (adegs & EXP_6_V10) < (bdegs & EXP_6_V10) ||
			   (adegs & EXP_7_V10) < (bdegs & EXP_7_V10) ||
			   (adegs & EXP_8_V10) < (bdegs & EXP_8_V10) ||
			   (adegs & EXP_9_V10) < (bdegs & EXP_9_V10) ||
			   (adegs & EXP_10_V10) < (bdegs & EXP_10_V10)  );
}

//monomial divide test for nvar 11
static inline int monomialDivideTest_11(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V11) < (bdegs & EXP_1_V11) ||
			   (adegs & EXP_2_V11) < (bdegs & EXP_2_V11) ||
			   (adegs & EXP_3_V11) < (bdegs & EXP_3_V11) ||
			   (adegs & EXP_4_V11) < (bdegs & EXP_4_V11) ||
			   (adegs & EXP_5_V11) < (bdegs & EXP_5_V11) ||
			   (adegs & EXP_6_V11) < (bdegs & EXP_6_V11) ||
			   (adegs & EXP_7_V11) < (bdegs & EXP_7_V11) ||
			   (adegs & EXP_8_V11) < (bdegs & EXP_8_V11) ||
			   (adegs & EXP_9_V11) < (bdegs & EXP_9_V11) ||
			   (adegs & EXP_10_V11) < (bdegs & EXP_10_V11) ||
			   (adegs & EXP_11_V11) < (bdegs & EXP_11_V11)  );
}

//monomial divide test for nvar 12
static inline int monomialDivideTest_12(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V12) < (bdegs & EXP_1_V12) ||
			   (adegs & EXP_2_V12) < (bdegs & EXP_2_V12) ||
			   (adegs & EXP_3_V12) < (bdegs & EXP_3_V12) ||
			   (adegs & EXP_4_V12) < (bdegs & EXP_4_V12) ||
			   (adegs & EXP_5_V12) < (bdegs & EXP_5_V12) ||
			   (adegs & EXP_6_V12) < (bdegs & EXP_6_V12) ||
			   (adegs & EXP_7_V12) < (bdegs & EXP_7_V12) ||
			   (adegs & EXP_8_V12) < (bdegs & EXP_8_V12) ||
			   (adegs & EXP_9_V12) < (bdegs & EXP_9_V12) ||
			   (adegs & EXP_10_V12) < (bdegs & EXP_10_V12) ||
			   (adegs & EXP_11_V12) < (bdegs & EXP_11_V12) ||
			   (adegs & EXP_12_V12) < (bdegs & EXP_12_V12)  );
}

//monomial divide test for nvar 13
static inline int monomialDivideTest_13(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V13) < (bdegs & EXP_1_V13) ||
			   (adegs & EXP_2_V13) < (bdegs & EXP_2_V13) ||
			   (adegs & EXP_3_V13) < (bdegs & EXP_3_V13) ||
			   (adegs & EXP_4_V13) < (bdegs & EXP_4_V13) ||
			   (adegs & EXP_5_V13) < (bdegs & EXP_5_V13) ||
			   (adegs & EXP_6_V13) < (bdegs & EXP_6_V13) ||
			   (adegs & EXP_7_V13) < (bdegs & EXP_7_V13) ||
			   (adegs & EXP_8_V13) < (bdegs & EXP_8_V13) ||
			   (adegs & EXP_9_V13) < (bdegs & EXP_9_V13) ||
			   (adegs & EXP_10_V13) < (bdegs & EXP_10_V13) ||
			   (adegs & EXP_11_V13) < (bdegs & EXP_11_V13) ||
			   (adegs & EXP_12_V13) < (bdegs & EXP_12_V13) ||
			   (adegs & EXP_13_V13) < (bdegs & EXP_13_V13)  );
}

//monomial divide test for nvar 14
static inline int monomialDivideTest_14(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V14) < (bdegs & EXP_1_V14) ||
			   (adegs & EXP_2_V14) < (bdegs & EXP_2_V14) ||
			   (adegs & EXP_3_V14) < (bdegs & EXP_3_V14) ||
			   (adegs & EXP_4_V14) < (bdegs & EXP_4_V14) ||
			   (adegs & EXP_5_V14) < (bdegs & EXP_5_V14) ||
			   (adegs & EXP_6_V14) < (bdegs & EXP_6_V14) ||
			   (adegs & EXP_7_V14) < (bdegs & EXP_7_V14) ||
			   (adegs & EXP_8_V14) < (bdegs & EXP_8_V14) ||
			   (adegs & EXP_9_V14) < (bdegs & EXP_9_V14) ||
			   (adegs & EXP_10_V14) < (bdegs & EXP_10_V14) ||
			   (adegs & EXP_11_V14) < (bdegs & EXP_11_V14) ||
			   (adegs & EXP_12_V14) < (bdegs & EXP_12_V14) ||
			   (adegs & EXP_13_V14) < (bdegs & EXP_13_V14) ||
			   (adegs & EXP_14_V14) < (bdegs & EXP_14_V14)  );
}

//monomial divide test for nvar 15
static inline int monomialDivideTest_15(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V15) < (bdegs & EXP_1_V15) ||
			   (adegs & EXP_2_V15) < (bdegs & EXP_2_V15) ||
			   (adegs & EXP_3_V15) < (bdegs & EXP_3_V15) ||
			   (adegs & EXP_4_V15) < (bdegs & EXP_4_V15) ||
			   (adegs & EXP_5_V15) < (bdegs & EXP_5_V15) ||
			   (adegs & EXP_6_V15) < (bdegs & EXP_6_V15) ||
			   (adegs & EXP_7_V15) < (bdegs & EXP_7_V15) ||
			   (adegs & EXP_8_V15) < (bdegs & EXP_8_V15) ||
			   (adegs & EXP_9_V15) < (bdegs & EXP_9_V15) ||
			   (adegs & EXP_10_V15) < (bdegs & EXP_10_V15) ||
			   (adegs & EXP_11_V15) < (bdegs & EXP_11_V15) ||
			   (adegs & EXP_12_V15) < (bdegs & EXP_12_V15) ||
			   (adegs & EXP_13_V15) < (bdegs & EXP_13_V15) ||
			   (adegs & EXP_14_V15) < (bdegs & EXP_14_V15) ||
			   (adegs & EXP_15_V15) < (bdegs & EXP_15_V15)  );
}

//monomial divide test for nvar 16
static inline int monomialDivideTest_16(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V16) < (bdegs & EXP_1_V16) ||
			   (adegs & EXP_2_V16) < (bdegs & EXP_2_V16) ||
			   (adegs & EXP_3_V16) < (bdegs & EXP_3_V16) ||
			   (adegs & EXP_4_V16) < (bdegs & EXP_4_V16) ||
			   (adegs & EXP_5_V16) < (bdegs & EXP_5_V16) ||
			   (adegs & EXP_6_V16) < (bdegs & EXP_6_V16) ||
			   (adegs & EXP_7_V16) < (bdegs & EXP_7_V16) ||
			   (adegs & EXP_8_V16) < (bdegs & EXP_8_V16) ||
			   (adegs & EXP_9_V16) < (bdegs & EXP_9_V16) ||
			   (adegs & EXP_10_V16) < (bdegs & EXP_10_V16) ||
			   (adegs & EXP_11_V16) < (bdegs & EXP_11_V16) ||
			   (adegs & EXP_12_V16) < (bdegs & EXP_12_V16) ||
			   (adegs & EXP_13_V16) < (bdegs & EXP_13_V16) ||
			   (adegs & EXP_14_V16) < (bdegs & EXP_14_V16) ||
			   (adegs & EXP_15_V16) < (bdegs & EXP_15_V16) ||
			   (adegs & EXP_16_V16) < (bdegs & EXP_16_V16)  );
}

//monomial divide test for nvar 17
static inline int monomialDivideTest_17(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V17) < (bdegs & EXP_1_V17) ||
			   (adegs & EXP_2_V17) < (bdegs & EXP_2_V17) ||
			   (adegs & EXP_3_V17) < (bdegs & EXP_3_V17) ||
			   (adegs & EXP_4_V17) < (bdegs & EXP_4_V17) ||
			   (adegs & EXP_5_V17) < (bdegs & EXP_5_V17) ||
			   (adegs & EXP_6_V17) < (bdegs & EXP_6_V17) ||
			   (adegs & EXP_7_V17) < (bdegs & EXP_7_V17) ||
			   (adegs & EXP_8_V17) < (bdegs & EXP_8_V17) ||
			   (adegs & EXP_9_V17) < (bdegs & EXP_9_V17) ||
			   (adegs & EXP_10_V17) < (bdegs & EXP_10_V17) ||
			   (adegs & EXP_11_V17) < (bdegs & EXP_11_V17) ||
			   (adegs & EXP_12_V17) < (bdegs & EXP_12_V17) ||
			   (adegs & EXP_13_V17) < (bdegs & EXP_13_V17) ||
			   (adegs & EXP_14_V17) < (bdegs & EXP_14_V17) ||
			   (adegs & EXP_15_V17) < (bdegs & EXP_15_V17) ||
			   (adegs & EXP_16_V17) < (bdegs & EXP_16_V17) ||
			   (adegs & EXP_17_V17) < (bdegs & EXP_17_V17)  );
}

//monomial divide test for nvar 18
static inline int monomialDivideTest_18(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V18) < (bdegs & EXP_1_V18) ||
			   (adegs & EXP_2_V18) < (bdegs & EXP_2_V18) ||
			   (adegs & EXP_3_V18) < (bdegs & EXP_3_V18) ||
			   (adegs & EXP_4_V18) < (bdegs & EXP_4_V18) ||
			   (adegs & EXP_5_V18) < (bdegs & EXP_5_V18) ||
			   (adegs & EXP_6_V18) < (bdegs & EXP_6_V18) ||
			   (adegs & EXP_7_V18) < (bdegs & EXP_7_V18) ||
			   (adegs & EXP_8_V18) < (bdegs & EXP_8_V18) ||
			   (adegs & EXP_9_V18) < (bdegs & EXP_9_V18) ||
			   (adegs & EXP_10_V18) < (bdegs & EXP_10_V18) ||
			   (adegs & EXP_11_V18) < (bdegs & EXP_11_V18) ||
			   (adegs & EXP_12_V18) < (bdegs & EXP_12_V18) ||
			   (adegs & EXP_13_V18) < (bdegs & EXP_13_V18) ||
			   (adegs & EXP_14_V18) < (bdegs & EXP_14_V18) ||
			   (adegs & EXP_15_V18) < (bdegs & EXP_15_V18) ||
			   (adegs & EXP_16_V18) < (bdegs & EXP_16_V18) ||
			   (adegs & EXP_17_V18) < (bdegs & EXP_17_V18) ||
			   (adegs & EXP_18_V18) < (bdegs & EXP_18_V18)  );
}

//monomial divide test for nvar 19
static inline int monomialDivideTest_19(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V19) < (bdegs & EXP_1_V19) ||
			   (adegs & EXP_2_V19) < (bdegs & EXP_2_V19) ||
			   (adegs & EXP_3_V19) < (bdegs & EXP_3_V19) ||
			   (adegs & EXP_4_V19) < (bdegs & EXP_4_V19) ||
			   (adegs & EXP_5_V19) < (bdegs & EXP_5_V19) ||
			   (adegs & EXP_6_V19) < (bdegs & EXP_6_V19) ||
			   (adegs & EXP_7_V19) < (bdegs & EXP_7_V19) ||
			   (adegs & EXP_8_V19) < (bdegs & EXP_8_V19) ||
			   (adegs & EXP_9_V19) < (bdegs & EXP_9_V19) ||
			   (adegs & EXP_10_V19) < (bdegs & EXP_10_V19) ||
			   (adegs & EXP_11_V19) < (bdegs & EXP_11_V19) ||
			   (adegs & EXP_12_V19) < (bdegs & EXP_12_V19) ||
			   (adegs & EXP_13_V19) < (bdegs & EXP_13_V19) ||
			   (adegs & EXP_14_V19) < (bdegs & EXP_14_V19) ||
			   (adegs & EXP_15_V19) < (bdegs & EXP_15_V19) ||
			   (adegs & EXP_16_V19) < (bdegs & EXP_16_V19) ||
			   (adegs & EXP_17_V19) < (bdegs & EXP_17_V19) ||
			   (adegs & EXP_18_V19) < (bdegs & EXP_18_V19) ||
			   (adegs & EXP_19_V19) < (bdegs & EXP_19_V19)  );
}

//monomial divide test for nvar 20
static inline int monomialDivideTest_20(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V20) < (bdegs & EXP_1_V20) ||
			   (adegs & EXP_2_V20) < (bdegs & EXP_2_V20) ||
			   (adegs & EXP_3_V20) < (bdegs & EXP_3_V20) ||
			   (adegs & EXP_4_V20) < (bdegs & EXP_4_V20) ||
			   (adegs & EXP_5_V20) < (bdegs & EXP_5_V20) ||
			   (adegs & EXP_6_V20) < (bdegs & EXP_6_V20) ||
			   (adegs & EXP_7_V20) < (bdegs & EXP_7_V20) ||
			   (adegs & EXP_8_V20) < (bdegs & EXP_8_V20) ||
			   (adegs & EXP_9_V20) < (bdegs & EXP_9_V20) ||
			   (adegs & EXP_10_V20) < (bdegs & EXP_10_V20) ||
			   (adegs & EXP_11_V20) < (bdegs & EXP_11_V20) ||
			   (adegs & EXP_12_V20) < (bdegs & EXP_12_V20) ||
			   (adegs & EXP_13_V20) < (bdegs & EXP_13_V20) ||
			   (adegs & EXP_14_V20) < (bdegs & EXP_14_V20) ||
			   (adegs & EXP_15_V20) < (bdegs & EXP_15_V20) ||
			   (adegs & EXP_16_V20) < (bdegs & EXP_16_V20) ||
			   (adegs & EXP_17_V20) < (bdegs & EXP_17_V20) ||
			   (adegs & EXP_18_V20) < (bdegs & EXP_18_V20) ||
			   (adegs & EXP_19_V20) < (bdegs & EXP_19_V20) ||
			   (adegs & EXP_20_V20) < (bdegs & EXP_20_V20)  );
}

//monomial divide test for nvar 21
static inline int monomialDivideTest_21(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V21) < (bdegs & EXP_1_V21) ||
			   (adegs & EXP_2_V21) < (bdegs & EXP_2_V21) ||
			   (adegs & EXP_3_V21) < (bdegs & EXP_3_V21) ||
			   (adegs & EXP_4_V21) < (bdegs & EXP_4_V21) ||
			   (adegs & EXP_5_V21) < (bdegs & EXP_5_V21) ||
			   (adegs & EXP_6_V21) < (bdegs & EXP_6_V21) ||
			   (adegs & EXP_7_V21) < (bdegs & EXP_7_V21) ||
			   (adegs & EXP_8_V21) < (bdegs & EXP_8_V21) ||
			   (adegs & EXP_9_V21) < (bdegs & EXP_9_V21) ||
			   (adegs & EXP_10_V21) < (bdegs & EXP_10_V21) ||
			   (adegs & EXP_11_V21) < (bdegs & EXP_11_V21) ||
			   (adegs & EXP_12_V21) < (bdegs & EXP_12_V21) ||
			   (adegs & EXP_13_V21) < (bdegs & EXP_13_V21) ||
			   (adegs & EXP_14_V21) < (bdegs & EXP_14_V21) ||
			   (adegs & EXP_15_V21) < (bdegs & EXP_15_V21) ||
			   (adegs & EXP_16_V21) < (bdegs & EXP_16_V21) ||
			   (adegs & EXP_17_V21) < (bdegs & EXP_17_V21) ||
			   (adegs & EXP_18_V21) < (bdegs & EXP_18_V21) ||
			   (adegs & EXP_19_V21) < (bdegs & EXP_19_V21) ||
			   (adegs & EXP_20_V21) < (bdegs & EXP_20_V21) ||
			   (adegs & EXP_21_V21) < (bdegs & EXP_21_V21)  );
}

//monomial divide test for nvar 22
static inline int monomialDivideTest_22(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V22) < (bdegs & EXP_1_V22) ||
			   (adegs & EXP_2_V22) < (bdegs & EXP_2_V22) ||
			   (adegs & EXP_3_V22) < (bdegs & EXP_3_V22) ||
			   (adegs & EXP_4_V22) < (bdegs & EXP_4_V22) ||
			   (adegs & EXP_5_V22) < (bdegs & EXP_5_V22) ||
			   (adegs & EXP_6_V22) < (bdegs & EXP_6_V22) ||
			   (adegs & EXP_7_V22) < (bdegs & EXP_7_V22) ||
			   (adegs & EXP_8_V22) < (bdegs & EXP_8_V22) ||
			   (adegs & EXP_9_V22) < (bdegs & EXP_9_V22) ||
			   (adegs & EXP_10_V22) < (bdegs & EXP_10_V22) ||
			   (adegs & EXP_11_V22) < (bdegs & EXP_11_V22) ||
			   (adegs & EXP_12_V22) < (bdegs & EXP_12_V22) ||
			   (adegs & EXP_13_V22) < (bdegs & EXP_13_V22) ||
			   (adegs & EXP_14_V22) < (bdegs & EXP_14_V22) ||
			   (adegs & EXP_15_V22) < (bdegs & EXP_15_V22) ||
			   (adegs & EXP_16_V22) < (bdegs & EXP_16_V22) ||
			   (adegs & EXP_17_V22) < (bdegs & EXP_17_V22) ||
			   (adegs & EXP_18_V22) < (bdegs & EXP_18_V22) ||
			   (adegs & EXP_19_V22) < (bdegs & EXP_19_V22) ||
			   (adegs & EXP_20_V22) < (bdegs & EXP_20_V22) ||
			   (adegs & EXP_21_V22) < (bdegs & EXP_21_V22) ||
			   (adegs & EXP_22_V22) < (bdegs & EXP_22_V22)  );
}

//monomial divide test for nvar 23
static inline int monomialDivideTest_23(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V23) < (bdegs & EXP_1_V23) ||
			   (adegs & EXP_2_V23) < (bdegs & EXP_2_V23) ||
			   (adegs & EXP_3_V23) < (bdegs & EXP_3_V23) ||
			   (adegs & EXP_4_V23) < (bdegs & EXP_4_V23) ||
			   (adegs & EXP_5_V23) < (bdegs & EXP_5_V23) ||
			   (adegs & EXP_6_V23) < (bdegs & EXP_6_V23) ||
			   (adegs & EXP_7_V23) < (bdegs & EXP_7_V23) ||
			   (adegs & EXP_8_V23) < (bdegs & EXP_8_V23) ||
			   (adegs & EXP_9_V23) < (bdegs & EXP_9_V23) ||
			   (adegs & EXP_10_V23) < (bdegs & EXP_10_V23) ||
			   (adegs & EXP_11_V23) < (bdegs & EXP_11_V23) ||
			   (adegs & EXP_12_V23) < (bdegs & EXP_12_V23) ||
			   (adegs & EXP_13_V23) < (bdegs & EXP_13_V23) ||
			   (adegs & EXP_14_V23) < (bdegs & EXP_14_V23) ||
			   (adegs & EXP_15_V23) < (bdegs & EXP_15_V23) ||
			   (adegs & EXP_16_V23) < (bdegs & EXP_16_V23) ||
			   (adegs & EXP_17_V23) < (bdegs & EXP_17_V23) ||
			   (adegs & EXP_18_V23) < (bdegs & EXP_18_V23) ||
			   (adegs & EXP_19_V23) < (bdegs & EXP_19_V23) ||
			   (adegs & EXP_20_V23) < (bdegs & EXP_20_V23) ||
			   (adegs & EXP_21_V23) < (bdegs & EXP_21_V23) ||
			   (adegs & EXP_22_V23) < (bdegs & EXP_22_V23) ||
			   (adegs & EXP_23_V23) < (bdegs & EXP_23_V23)  );
}

//monomial divide test for nvar 24
static inline int monomialDivideTest_24(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V24) < (bdegs & EXP_1_V24) ||
			   (adegs & EXP_2_V24) < (bdegs & EXP_2_V24) ||
			   (adegs & EXP_3_V24) < (bdegs & EXP_3_V24) ||
			   (adegs & EXP_4_V24) < (bdegs & EXP_4_V24) ||
			   (adegs & EXP_5_V24) < (bdegs & EXP_5_V24) ||
			   (adegs & EXP_6_V24) < (bdegs & EXP_6_V24) ||
			   (adegs & EXP_7_V24) < (bdegs & EXP_7_V24) ||
			   (adegs & EXP_8_V24) < (bdegs & EXP_8_V24) ||
			   (adegs & EXP_9_V24) < (bdegs & EXP_9_V24) ||
			   (adegs & EXP_10_V24) < (bdegs & EXP_10_V24) ||
			   (adegs & EXP_11_V24) < (bdegs & EXP_11_V24) ||
			   (adegs & EXP_12_V24) < (bdegs & EXP_12_V24) ||
			   (adegs & EXP_13_V24) < (bdegs & EXP_13_V24) ||
			   (adegs & EXP_14_V24) < (bdegs & EXP_14_V24) ||
			   (adegs & EXP_15_V24) < (bdegs & EXP_15_V24) ||
			   (adegs & EXP_16_V24) < (bdegs & EXP_16_V24) ||
			   (adegs & EXP_17_V24) < (bdegs & EXP_17_V24) ||
			   (adegs & EXP_18_V24) < (bdegs & EXP_18_V24) ||
			   (adegs & EXP_19_V24) < (bdegs & EXP_19_V24) ||
			   (adegs & EXP_20_V24) < (bdegs & EXP_20_V24) ||
			   (adegs & EXP_21_V24) < (bdegs & EXP_21_V24) ||
			   (adegs & EXP_22_V24) < (bdegs & EXP_22_V24) ||
			   (adegs & EXP_23_V24) < (bdegs & EXP_23_V24) ||
			   (adegs & EXP_24_V24) < (bdegs & EXP_24_V24)  );
}

//monomial divide test for nvar 25
static inline int monomialDivideTest_25(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V25) < (bdegs & EXP_1_V25) ||
			   (adegs & EXP_2_V25) < (bdegs & EXP_2_V25) ||
			   (adegs & EXP_3_V25) < (bdegs & EXP_3_V25) ||
			   (adegs & EXP_4_V25) < (bdegs & EXP_4_V25) ||
			   (adegs & EXP_5_V25) < (bdegs & EXP_5_V25) ||
			   (adegs & EXP_6_V25) < (bdegs & EXP_6_V25) ||
			   (adegs & EXP_7_V25) < (bdegs & EXP_7_V25) ||
			   (adegs & EXP_8_V25) < (bdegs & EXP_8_V25) ||
			   (adegs & EXP_9_V25) < (bdegs & EXP_9_V25) ||
			   (adegs & EXP_10_V25) < (bdegs & EXP_10_V25) ||
			   (adegs & EXP_11_V25) < (bdegs & EXP_11_V25) ||
			   (adegs & EXP_12_V25) < (bdegs & EXP_12_V25) ||
			   (adegs & EXP_13_V25) < (bdegs & EXP_13_V25) ||
			   (adegs & EXP_14_V25) < (bdegs & EXP_14_V25) ||
			   (adegs & EXP_15_V25) < (bdegs & EXP_15_V25) ||
			   (adegs & EXP_16_V25) < (bdegs & EXP_16_V25) ||
			   (adegs & EXP_17_V25) < (bdegs & EXP_17_V25) ||
			   (adegs & EXP_18_V25) < (bdegs & EXP_18_V25) ||
			   (adegs & EXP_19_V25) < (bdegs & EXP_19_V25) ||
			   (adegs & EXP_20_V25) < (bdegs & EXP_20_V25) ||
			   (adegs & EXP_21_V25) < (bdegs & EXP_21_V25) ||
			   (adegs & EXP_22_V25) < (bdegs & EXP_22_V25) ||
			   (adegs & EXP_23_V25) < (bdegs & EXP_23_V25) ||
			   (adegs & EXP_24_V25) < (bdegs & EXP_24_V25) ||
			   (adegs & EXP_25_V25) < (bdegs & EXP_25_V25)  );
}

//monomial divide test for nvar 26
static inline int monomialDivideTest_26(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V26) < (bdegs & EXP_1_V26) ||
			   (adegs & EXP_2_V26) < (bdegs & EXP_2_V26) ||
			   (adegs & EXP_3_V26) < (bdegs & EXP_3_V26) ||
			   (adegs & EXP_4_V26) < (bdegs & EXP_4_V26) ||
			   (adegs & EXP_5_V26) < (bdegs & EXP_5_V26) ||
			   (adegs & EXP_6_V26) < (bdegs & EXP_6_V26) ||
			   (adegs & EXP_7_V26) < (bdegs & EXP_7_V26) ||
			   (adegs & EXP_8_V26) < (bdegs & EXP_8_V26) ||
			   (adegs & EXP_9_V26) < (bdegs & EXP_9_V26) ||
			   (adegs & EXP_10_V26) < (bdegs & EXP_10_V26) ||
			   (adegs & EXP_11_V26) < (bdegs & EXP_11_V26) ||
			   (adegs & EXP_12_V26) < (bdegs & EXP_12_V26) ||
			   (adegs & EXP_13_V26) < (bdegs & EXP_13_V26) ||
			   (adegs & EXP_14_V26) < (bdegs & EXP_14_V26) ||
			   (adegs & EXP_15_V26) < (bdegs & EXP_15_V26) ||
			   (adegs & EXP_16_V26) < (bdegs & EXP_16_V26) ||
			   (adegs & EXP_17_V26) < (bdegs & EXP_17_V26) ||
			   (adegs & EXP_18_V26) < (bdegs & EXP_18_V26) ||
			   (adegs & EXP_19_V26) < (bdegs & EXP_19_V26) ||
			   (adegs & EXP_20_V26) < (bdegs & EXP_20_V26) ||
			   (adegs & EXP_21_V26) < (bdegs & EXP_21_V26) ||
			   (adegs & EXP_22_V26) < (bdegs & EXP_22_V26) ||
			   (adegs & EXP_23_V26) < (bdegs & EXP_23_V26) ||
			   (adegs & EXP_24_V26) < (bdegs & EXP_24_V26) ||
			   (adegs & EXP_25_V26) < (bdegs & EXP_25_V26) ||
			   (adegs & EXP_26_V26) < (bdegs & EXP_26_V26)  );
}

//monomial divide test for nvar 27
static inline int monomialDivideTest_27(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V27) < (bdegs & EXP_1_V27) ||
			   (adegs & EXP_2_V27) < (bdegs & EXP_2_V27) ||
			   (adegs & EXP_3_V27) < (bdegs & EXP_3_V27) ||
			   (adegs & EXP_4_V27) < (bdegs & EXP_4_V27) ||
			   (adegs & EXP_5_V27) < (bdegs & EXP_5_V27) ||
			   (adegs & EXP_6_V27) < (bdegs & EXP_6_V27) ||
			   (adegs & EXP_7_V27) < (bdegs & EXP_7_V27) ||
			   (adegs & EXP_8_V27) < (bdegs & EXP_8_V27) ||
			   (adegs & EXP_9_V27) < (bdegs & EXP_9_V27) ||
			   (adegs & EXP_10_V27) < (bdegs & EXP_10_V27) ||
			   (adegs & EXP_11_V27) < (bdegs & EXP_11_V27) ||
			   (adegs & EXP_12_V27) < (bdegs & EXP_12_V27) ||
			   (adegs & EXP_13_V27) < (bdegs & EXP_13_V27) ||
			   (adegs & EXP_14_V27) < (bdegs & EXP_14_V27) ||
			   (adegs & EXP_15_V27) < (bdegs & EXP_15_V27) ||
			   (adegs & EXP_16_V27) < (bdegs & EXP_16_V27) ||
			   (adegs & EXP_17_V27) < (bdegs & EXP_17_V27) ||
			   (adegs & EXP_18_V27) < (bdegs & EXP_18_V27) ||
			   (adegs & EXP_19_V27) < (bdegs & EXP_19_V27) ||
			   (adegs & EXP_20_V27) < (bdegs & EXP_20_V27) ||
			   (adegs & EXP_21_V27) < (bdegs & EXP_21_V27) ||
			   (adegs & EXP_22_V27) < (bdegs & EXP_22_V27) ||
			   (adegs & EXP_23_V27) < (bdegs & EXP_23_V27) ||
			   (adegs & EXP_24_V27) < (bdegs & EXP_24_V27) ||
			   (adegs & EXP_25_V27) < (bdegs & EXP_25_V27) ||
			   (adegs & EXP_26_V27) < (bdegs & EXP_26_V27) ||
			   (adegs & EXP_27_V27) < (bdegs & EXP_27_V27)  );
}

//monomial divide test for nvar 28
static inline int monomialDivideTest_28(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V28) < (bdegs & EXP_1_V28) ||
			   (adegs & EXP_2_V28) < (bdegs & EXP_2_V28) ||
			   (adegs & EXP_3_V28) < (bdegs & EXP_3_V28) ||
			   (adegs & EXP_4_V28) < (bdegs & EXP_4_V28) ||
			   (adegs & EXP_5_V28) < (bdegs & EXP_5_V28) ||
			   (adegs & EXP_6_V28) < (bdegs & EXP_6_V28) ||
			   (adegs & EXP_7_V28) < (bdegs & EXP_7_V28) ||
			   (adegs & EXP_8_V28) < (bdegs & EXP_8_V28) ||
			   (adegs & EXP_9_V28) < (bdegs & EXP_9_V28) ||
			   (adegs & EXP_10_V28) < (bdegs & EXP_10_V28) ||
			   (adegs & EXP_11_V28) < (bdegs & EXP_11_V28) ||
			   (adegs & EXP_12_V28) < (bdegs & EXP_12_V28) ||
			   (adegs & EXP_13_V28) < (bdegs & EXP_13_V28) ||
			   (adegs & EXP_14_V28) < (bdegs & EXP_14_V28) ||
			   (adegs & EXP_15_V28) < (bdegs & EXP_15_V28) ||
			   (adegs & EXP_16_V28) < (bdegs & EXP_16_V28) ||
			   (adegs & EXP_17_V28) < (bdegs & EXP_17_V28) ||
			   (adegs & EXP_18_V28) < (bdegs & EXP_18_V28) ||
			   (adegs & EXP_19_V28) < (bdegs & EXP_19_V28) ||
			   (adegs & EXP_20_V28) < (bdegs & EXP_20_V28) ||
			   (adegs & EXP_21_V28) < (bdegs & EXP_21_V28) ||
			   (adegs & EXP_22_V28) < (bdegs & EXP_22_V28) ||
			   (adegs & EXP_23_V28) < (bdegs & EXP_23_V28) ||
			   (adegs & EXP_24_V28) < (bdegs & EXP_24_V28) ||
			   (adegs & EXP_25_V28) < (bdegs & EXP_25_V28) ||
			   (adegs & EXP_26_V28) < (bdegs & EXP_26_V28) ||
			   (adegs & EXP_27_V28) < (bdegs & EXP_27_V28) ||
			   (adegs & EXP_28_V28) < (bdegs & EXP_28_V28)  );
}

//monomial divide test for nvar 29
static inline int monomialDivideTest_29(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V29) < (bdegs & EXP_1_V29) ||
			   (adegs & EXP_2_V29) < (bdegs & EXP_2_V29) ||
			   (adegs & EXP_3_V29) < (bdegs & EXP_3_V29) ||
			   (adegs & EXP_4_V29) < (bdegs & EXP_4_V29) ||
			   (adegs & EXP_5_V29) < (bdegs & EXP_5_V29) ||
			   (adegs & EXP_6_V29) < (bdegs & EXP_6_V29) ||
			   (adegs & EXP_7_V29) < (bdegs & EXP_7_V29) ||
			   (adegs & EXP_8_V29) < (bdegs & EXP_8_V29) ||
			   (adegs & EXP_9_V29) < (bdegs & EXP_9_V29) ||
			   (adegs & EXP_10_V29) < (bdegs & EXP_10_V29) ||
			   (adegs & EXP_11_V29) < (bdegs & EXP_11_V29) ||
			   (adegs & EXP_12_V29) < (bdegs & EXP_12_V29) ||
			   (adegs & EXP_13_V29) < (bdegs & EXP_13_V29) ||
			   (adegs & EXP_14_V29) < (bdegs & EXP_14_V29) ||
			   (adegs & EXP_15_V29) < (bdegs & EXP_15_V29) ||
			   (adegs & EXP_16_V29) < (bdegs & EXP_16_V29) ||
			   (adegs & EXP_17_V29) < (bdegs & EXP_17_V29) ||
			   (adegs & EXP_18_V29) < (bdegs & EXP_18_V29) ||
			   (adegs & EXP_19_V29) < (bdegs & EXP_19_V29) ||
			   (adegs & EXP_20_V29) < (bdegs & EXP_20_V29) ||
			   (adegs & EXP_21_V29) < (bdegs & EXP_21_V29) ||
			   (adegs & EXP_22_V29) < (bdegs & EXP_22_V29) ||
			   (adegs & EXP_23_V29) < (bdegs & EXP_23_V29) ||
			   (adegs & EXP_24_V29) < (bdegs & EXP_24_V29) ||
			   (adegs & EXP_25_V29) < (bdegs & EXP_25_V29) ||
			   (adegs & EXP_26_V29) < (bdegs & EXP_26_V29) ||
			   (adegs & EXP_27_V29) < (bdegs & EXP_27_V29) ||
			   (adegs & EXP_28_V29) < (bdegs & EXP_28_V29) ||
			   (adegs & EXP_29_V29) < (bdegs & EXP_29_V29)  );
}

//monomial divide test for nvar 30
static inline int monomialDivideTest_30(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V30) < (bdegs & EXP_1_V30) ||
			   (adegs & EXP_2_V30) < (bdegs & EXP_2_V30) ||
			   (adegs & EXP_3_V30) < (bdegs & EXP_3_V30) ||
			   (adegs & EXP_4_V30) < (bdegs & EXP_4_V30) ||
			   (adegs & EXP_5_V30) < (bdegs & EXP_5_V30) ||
			   (adegs & EXP_6_V30) < (bdegs & EXP_6_V30) ||
			   (adegs & EXP_7_V30) < (bdegs & EXP_7_V30) ||
			   (adegs & EXP_8_V30) < (bdegs & EXP_8_V30) ||
			   (adegs & EXP_9_V30) < (bdegs & EXP_9_V30) ||
			   (adegs & EXP_10_V30) < (bdegs & EXP_10_V30) ||
			   (adegs & EXP_11_V30) < (bdegs & EXP_11_V30) ||
			   (adegs & EXP_12_V30) < (bdegs & EXP_12_V30) ||
			   (adegs & EXP_13_V30) < (bdegs & EXP_13_V30) ||
			   (adegs & EXP_14_V30) < (bdegs & EXP_14_V30) ||
			   (adegs & EXP_15_V30) < (bdegs & EXP_15_V30) ||
			   (adegs & EXP_16_V30) < (bdegs & EXP_16_V30) ||
			   (adegs & EXP_17_V30) < (bdegs & EXP_17_V30) ||
			   (adegs & EXP_18_V30) < (bdegs & EXP_18_V30) ||
			   (adegs & EXP_19_V30) < (bdegs & EXP_19_V30) ||
			   (adegs & EXP_20_V30) < (bdegs & EXP_20_V30) ||
			   (adegs & EXP_21_V30) < (bdegs & EXP_21_V30) ||
			   (adegs & EXP_22_V30) < (bdegs & EXP_22_V30) ||
			   (adegs & EXP_23_V30) < (bdegs & EXP_23_V30) ||
			   (adegs & EXP_24_V30) < (bdegs & EXP_24_V30) ||
			   (adegs & EXP_25_V30) < (bdegs & EXP_25_V30) ||
			   (adegs & EXP_26_V30) < (bdegs & EXP_26_V30) ||
			   (adegs & EXP_27_V30) < (bdegs & EXP_27_V30) ||
			   (adegs & EXP_28_V30) < (bdegs & EXP_28_V30) ||
			   (adegs & EXP_29_V30) < (bdegs & EXP_29_V30) ||
			   (adegs & EXP_30_V30) < (bdegs & EXP_30_V30)  );
}

//monomial divide test for nvar 31
static inline int monomialDivideTest_31(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V31) < (bdegs & EXP_1_V31) ||
			   (adegs & EXP_2_V31) < (bdegs & EXP_2_V31) ||
			   (adegs & EXP_3_V31) < (bdegs & EXP_3_V31) ||
			   (adegs & EXP_4_V31) < (bdegs & EXP_4_V31) ||
			   (adegs & EXP_5_V31) < (bdegs & EXP_5_V31) ||
			   (adegs & EXP_6_V31) < (bdegs & EXP_6_V31) ||
			   (adegs & EXP_7_V31) < (bdegs & EXP_7_V31) ||
			   (adegs & EXP_8_V31) < (bdegs & EXP_8_V31) ||
			   (adegs & EXP_9_V31) < (bdegs & EXP_9_V31) ||
			   (adegs & EXP_10_V31) < (bdegs & EXP_10_V31) ||
			   (adegs & EXP_11_V31) < (bdegs & EXP_11_V31) ||
			   (adegs & EXP_12_V31) < (bdegs & EXP_12_V31) ||
			   (adegs & EXP_13_V31) < (bdegs & EXP_13_V31) ||
			   (adegs & EXP_14_V31) < (bdegs & EXP_14_V31) ||
			   (adegs & EXP_15_V31) < (bdegs & EXP_15_V31) ||
			   (adegs & EXP_16_V31) < (bdegs & EXP_16_V31) ||
			   (adegs & EXP_17_V31) < (bdegs & EXP_17_V31) ||
			   (adegs & EXP_18_V31) < (bdegs & EXP_18_V31) ||
			   (adegs & EXP_19_V31) < (bdegs & EXP_19_V31) ||
			   (adegs & EXP_20_V31) < (bdegs & EXP_20_V31) ||
			   (adegs & EXP_21_V31) < (bdegs & EXP_21_V31) ||
			   (adegs & EXP_22_V31) < (bdegs & EXP_22_V31) ||
			   (adegs & EXP_23_V31) < (bdegs & EXP_23_V31) ||
			   (adegs & EXP_24_V31) < (bdegs & EXP_24_V31) ||
			   (adegs & EXP_25_V31) < (bdegs & EXP_25_V31) ||
			   (adegs & EXP_26_V31) < (bdegs & EXP_26_V31) ||
			   (adegs & EXP_27_V31) < (bdegs & EXP_27_V31) ||
			   (adegs & EXP_28_V31) < (bdegs & EXP_28_V31) ||
			   (adegs & EXP_29_V31) < (bdegs & EXP_29_V31) ||
			   (adegs & EXP_30_V31) < (bdegs & EXP_30_V31) ||
			   (adegs & EXP_31_V31) < (bdegs & EXP_31_V31)  );
}

//monomial divide test for nvar 32
static inline int monomialDivideTest_32(degrees_t adegs, degrees_t bdegs) {
	return !(  (adegs & EXP_1_V32) < (bdegs & EXP_1_V32) ||
			   (adegs & EXP_2_V32) < (bdegs & EXP_2_V32) ||
			   (adegs & EXP_3_V32) < (bdegs & EXP_3_V32) ||
			   (adegs & EXP_4_V32) < (bdegs & EXP_4_V32) ||
			   (adegs & EXP_5_V32) < (bdegs & EXP_5_V32) ||
			   (adegs & EXP_6_V32) < (bdegs & EXP_6_V32) ||
			   (adegs & EXP_7_V32) < (bdegs & EXP_7_V32) ||
			   (adegs & EXP_8_V32) < (bdegs & EXP_8_V32) ||
			   (adegs & EXP_9_V32) < (bdegs & EXP_9_V32) ||
			   (adegs & EXP_10_V32) < (bdegs & EXP_10_V32) ||
			   (adegs & EXP_11_V32) < (bdegs & EXP_11_V32) ||
			   (adegs & EXP_12_V32) < (bdegs & EXP_12_V32) ||
			   (adegs & EXP_13_V32) < (bdegs & EXP_13_V32) ||
			   (adegs & EXP_14_V32) < (bdegs & EXP_14_V32) ||
			   (adegs & EXP_15_V32) < (bdegs & EXP_15_V32) ||
			   (adegs & EXP_16_V32) < (bdegs & EXP_16_V32) ||
			   (adegs & EXP_17_V32) < (bdegs & EXP_17_V32) ||
			   (adegs & EXP_18_V32) < (bdegs & EXP_18_V32) ||
			   (adegs & EXP_19_V32) < (bdegs & EXP_19_V32) ||
			   (adegs & EXP_20_V32) < (bdegs & EXP_20_V32) ||
			   (adegs & EXP_21_V32) < (bdegs & EXP_21_V32) ||
			   (adegs & EXP_22_V32) < (bdegs & EXP_22_V32) ||
			   (adegs & EXP_23_V32) < (bdegs & EXP_23_V32) ||
			   (adegs & EXP_24_V32) < (bdegs & EXP_24_V32) ||
			   (adegs & EXP_25_V32) < (bdegs & EXP_25_V32) ||
			   (adegs & EXP_26_V32) < (bdegs & EXP_26_V32) ||
			   (adegs & EXP_27_V32) < (bdegs & EXP_27_V32) ||
			   (adegs & EXP_28_V32) < (bdegs & EXP_28_V32) ||
			   (adegs & EXP_29_V32) < (bdegs & EXP_29_V32) ||
			   (adegs & EXP_30_V32) < (bdegs & EXP_30_V32) ||
			   (adegs & EXP_31_V32) < (bdegs & EXP_31_V32) ||
			   (adegs & EXP_32_V32) < (bdegs & EXP_32_V32)  );
}


static DivTest_ptr getMonomialDivideTestFuncPtr(int nvar) {
	switch (nvar)  {
		case 1: { return monomialDivideTest_1; }
		case 2: { return monomialDivideTest_2; }
		case 3: { return monomialDivideTest_3; }
		case 4: { return monomialDivideTest_4; }
		case 5: { return monomialDivideTest_5; }
		case 6: { return monomialDivideTest_6; }
		case 7: { return monomialDivideTest_7; }
		case 8: { return monomialDivideTest_8; }
		case 9: { return monomialDivideTest_9; }
		case 10: { return monomialDivideTest_10; }
		case 11: { return monomialDivideTest_11; }
		case 12: { return monomialDivideTest_12; }
		case 13: { return monomialDivideTest_13; }
		case 14: { return monomialDivideTest_14; }
		case 15: { return monomialDivideTest_15; }
		case 16: { return monomialDivideTest_16; }
		case 17: { return monomialDivideTest_17; }
		case 18: { return monomialDivideTest_18; }
		case 19: { return monomialDivideTest_19; }
		case 20: { return monomialDivideTest_20; }
		case 21: { return monomialDivideTest_21; }
		case 22: { return monomialDivideTest_22; }
		case 23: { return monomialDivideTest_23; }
		case 24: { return monomialDivideTest_24; }
		case 25: { return monomialDivideTest_25; }
		case 26: { return monomialDivideTest_26; }
		case 27: { return monomialDivideTest_27; }
		case 28: { return monomialDivideTest_28; }
		case 29: { return monomialDivideTest_29; }
		case 30: { return monomialDivideTest_30; }
		case 31: { return monomialDivideTest_31; }
		case 32: { return monomialDivideTest_32; }
	}
	return NULL;
}


/**
 * Determine if monomial b divides monomial a in lex ordering.
 * These exponent vectors are assumed to be unpacked.
 * returns non-zero iff monomial b divides monomial a.s
 */
static inline int monomialDivideTest_unpk(degrees_t adegs, degrees_t bdegs, int nvar) {
	degree_t* adegs_pt = (degree_t*) adegs;
	degree_t* bdegs_pt = (degree_t*) bdegs;
	for(int i = 0; i < nvar; ++i) {
		if (adegs_pt[i] < bdegs_pt[i]) {
			return 0;
		}
	}
	return 1;
}

/**
 * Helper to determine if monomial b divides monomial a.
 *
 * Note: since we assume working in Q, do not need to check coefficients 
 * nvar: number of variables of monomials a and b
 */
static inline int monomialDivideTest(degrees_t adegs, degrees_t bdegs, int nvar) {
	switch (nvar) {
		case 1: {
			return adegs >= bdegs;
		}
		case 2: {
			return monomialDivideTest_2(adegs,bdegs);
		}
		case 3: {
			return monomialDivideTest_3(adegs,bdegs);
		}
		case 4: {
			return monomialDivideTest_4(adegs,bdegs);
		}
		case 5: {
			return monomialDivideTest_5(adegs,bdegs);
		}
		case 6: {
			return monomialDivideTest_6(adegs,bdegs);
		}
		case 7: {
			return monomialDivideTest_7(adegs,bdegs);
		}
		case 8: {
			return monomialDivideTest_8(adegs,bdegs);
		}
		case 9: {
			return monomialDivideTest_9(adegs,bdegs);
		}
		case 10: {
			return monomialDivideTest_10(adegs,bdegs);
		}
		case 11: {
			return monomialDivideTest_11(adegs,bdegs);
		}
		case 12: {
			return monomialDivideTest_12(adegs,bdegs);
		}
		case 13: {
			return monomialDivideTest_13(adegs,bdegs);
		}
		case 14: {
			return monomialDivideTest_14(adegs,bdegs);
		}
		case 15: {
			return monomialDivideTest_15(adegs,bdegs);
		}
		case 16: {
			return monomialDivideTest_16(adegs,bdegs);
		}
		case 17: {
			return monomialDivideTest_17(adegs,bdegs);
		}
		case 18: {
			return monomialDivideTest_18(adegs,bdegs);
		}
		case 19: {
			return monomialDivideTest_19(adegs,bdegs);
		}
		case 20: {
			return monomialDivideTest_20(adegs,bdegs);
		}
		case 21: {
			return monomialDivideTest_21(adegs,bdegs);
		}
		case 22: {
			return monomialDivideTest_22(adegs,bdegs);
		}
		case 23: {
			return monomialDivideTest_23(adegs,bdegs);
		}
		case 24: {
			return monomialDivideTest_24(adegs,bdegs);
		}
		case 25: {
			return monomialDivideTest_25(adegs,bdegs);
		}
		case 26: {
			return monomialDivideTest_26(adegs,bdegs);
		}
		case 27: {
			return monomialDivideTest_27(adegs,bdegs);
		}
		case 28: {
			return monomialDivideTest_28(adegs,bdegs);
		}
		case 29: {
			return monomialDivideTest_29(adegs,bdegs);
		}
		case 30: {
			return monomialDivideTest_30(adegs,bdegs);
		}
		case 31: {
			return monomialDivideTest_31(adegs,bdegs);
		}
		case 32: {
			return monomialDivideTest_32(adegs,bdegs);
		}
		default : {
			fprintf(stderr, "MONOMIAL DIVIDE TEST NOT IMPLEMENTED FOR NVAR %d\n", nvar);
			exit(-1);
		}
	}

	return 0;
}


#endif