#ifndef BPAS_PARSER_H
#define BPAS_PARSER_H

// #define _GNU_SOURCE
// #define __STDC_WANT_LIB_EXT2__ 1
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "../../src/Parser/parser_type.h"
#include "../../src/Parser/parser_print.h"
#include "../RationalNumberPolynomial/SMQP_Support-AA.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    # define __BEGIN_DECLS extern "C" {
    # define __END_DECLS }
#else
    # define __BEGIN_DECLS /* empty */
    # define __END_DECLS /* empty */
#endif

__BEGIN_DECLS   

/**
* Creates Consecutive string array of size arrSize.
* Use the function to specify the variables of a polynomials.
* example: create_dyanmic_str_array(2, "x", "y");
* 	arrSize : size of the to be created
* 	... : arrSize number of string variables
* 	return: array of strings
**/
char **create_dynamic_str_array(int arrSize, ...);

/**
* Generates an alternating array of a polynomial from the string input poly_str.
* 	poly_str: string polynomial 
* 	return: AltArr_t representation of the string polynomial
**/
AltArr_t* generate_altarr(const char* poly_str);

/**
* Generates an alternating array of a polynomial from the string input poly_str.
* The variables and number of variables need to be specified to efficiently generate
* the alternating array format of the given string polynomial.
* 	poly_str: string polynomial 
*	variables: all the variables in the polynomial
	num_var: number of variables in the polynomial
* 	return: AltArr_t representation of the string polynomial
**/
AltArr_t* generate_altarr_var_defined(const char* poly_str,  const char** variables, int num_var);

/**
* Generates an alternating array pack containing the alternating array version of the polynomial, 
* all the variables, and the number of variables from the string input poly_str.
*       poly_str: string polynomial 
*       return: altarr_pack as stated above 
**/
altarr_pack* generate_altarr_pack(const char* poly_str);

__END_DECLS

#endif

