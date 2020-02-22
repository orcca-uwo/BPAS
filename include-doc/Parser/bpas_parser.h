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
AltArr_t* generate_altarr_var_defined(const char* poly_str,  char** variables, int num_var);

/**
* Generates an alternating array pack containing the alternating array version of the polynomial, 
* all the variables, and the number of variables from the string input poly_str.
*       poly_str: string polynomial 
*       return: altarr_pack as stated above 
**/
altarr_pack* generate_altarr_pack(const char* poly_str);

__END_DECLS

#endif

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


