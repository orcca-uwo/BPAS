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

