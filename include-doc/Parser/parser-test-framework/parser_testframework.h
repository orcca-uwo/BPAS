#ifndef PARSER_TESTFRAMEWORK_H
#define PARSER_TESTFRAMEWORK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ASSERT(msg, expression) (parser_assert((msg), (#expression), (expression)?1:0)) 
#define ASSERT_STRING_EQUAL(expected, actual) ASSERT((#expected), strcmp((expected), (actual))==0)
#define ASSERT_EQUAL(expected, actual) ASSERT((#expected), (expected)==(actual))

#define TT_COLOR_CODE 0x1B
#define TT_COLOR_RED "[1;31m"
#define TT_COLOR_GREEN "[1;32m"
#define TT_COLOR_RESET "[0m"

int parser_pass = 0;
int parser_fail = 0;
int fail = 0;
int pass = 0;
int parser_current_test_fail = 0;
const char *parser_msg = NULL;
const char *parser_expression = NULL;

int parser_assert(const char *msg, const char *expression, int pass);
void parser_test(const char *msg, const char* expr);
void parser_global_var_reset();
void parser_report();

int parser_assert(const char *msg, const char *expression, int pass){
    parser_msg = msg;
    parser_expression = expression;
    parser_pass = pass;
    return pass;
}

void parser_test(const char *msg, const char* expr){
    printf("print test\n");
    if(parser_pass){
        pass++;
        printf("%s , %s, %s%c%s%s%c%s%s\n", msg, expr, "[", TT_COLOR_CODE, TT_COLOR_GREEN, "PASS", TT_COLOR_CODE, TT_COLOR_RESET, "]");
        parser_global_var_reset();
    }
    else{
        fail++;
        printf("%s , %s, %s%c%s%s%c%s%s\n", msg, expr, "[", TT_COLOR_CODE, TT_COLOR_RED, "FAIL", TT_COLOR_CODE, TT_COLOR_RESET, "]");
        parser_global_var_reset();
    }
}

void parser_global_var_reset(){
    parser_pass = 0;
    parser_fail = 0;
    parser_current_test_fail = 0;
    parser_msg = NULL;
    parser_expression = NULL;
}

void parser_report(){
    printf("Total %d: PASS and %d: FAIL\n", pass, fail);
}



#ifdef __cplusplus
}
#endif

#endif /* PARSER_TESTFRAMEWORK_H */

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


