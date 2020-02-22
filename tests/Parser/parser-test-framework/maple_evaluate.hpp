#ifndef MAPLE_EVALUTE_H
#define MAPLE_EVALUTE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../MapleTestTool/MapleTestTool.hpp"

/**
 * @brief callback used for directing result output
 * 
 * @param data 
 * @param tag 
 * @param output 
 */
static void textCallBack( void *data, int tag, char *output );
/**
 * @brief Evaluate polynomial input to a Maple and return
 * result as a string.
 * 
 * @param in - maple input
 * @return char* - result output as string
 */
static std::string maple_evaluate(const char *in);
/**
 * @brief Evaluate and compare the result of two maple polynomial
 * expansion, and return 0 if they are equal. The comparison is 
 * done by subtracting the actual with expected input.
 * 
 * @param actual 
 * @param expected 
 * @return int 
 */
static int maple_evaluate_equal(const char *actual, const char *expected);
/**
 * @brief Handle fatal error, by exiting the program.
 * 
 * @param message 
 */
static void handle_parser_error(const char *message);





static void textCallBack( void *data, int tag, char *output )
{
    printf("%s\n",output);
}

static std::string maple_evaluate(const char *in){
	const char *pch = strchr(in, ']');
	char *input = (char*)malloc((strlen(in)+1)*sizeof(char));
	if(input == NULL){
		handle_parser_error("Unable to create heap memory.");
	}

	strcpy(input, in+(pch-in+1));
	const char *col = ":";
	strncat(input, col, strlen(col));

	MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
	MKernelVector kv = mapleTest->getMKernelVector();
    ALGEB r, l, m, n , p, q, s, t, e1, e2;  /* Maple data-structures */


	r = EvalMapleStatement(kv,"expand:");
	l = EvalMapleStatement(kv,input);
    e1 = EvalMapleProc(kv, r, 1, l);
	std::string msg = mapleTest->algebToString(kv, e1);
	return msg;
}

static int maple_evaluate_equal(const char *actual, const char *expected){
	char *input = (char*)malloc((strlen(actual)+1)*sizeof(char));
	if(input == NULL){
		handle_parser_error("Unable to create heap memory.");
	}
	char *input2 = (char*)malloc((strlen(expected)+1)*sizeof(char));
	if(input2 == NULL){
		handle_parser_error("Unable to create heap memory.");
	}

	if((strchr(actual, '[')) != NULL){
		const char *pch = strchr(actual, ']');
		strcpy(input, actual+(pch-actual+1));
	}else{
		strcpy(input, (actual));
	}
	strcpy(input2, expected);

const char *col = ":";
	strncat(input, col, strlen(col));
	strncat(input2, col, strlen(col));

	MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
	MKernelVector kv = mapleTest->getMKernelVector();
    ALGEB r, l, m, n , p, q, s, t, e1, e2;  /* Maple data-structures */

    /* initialize Maple */
	r = EvalMapleStatement(kv,"expand:");
	l = EvalMapleStatement(kv,input);
    e1 = EvalMapleProc(kv, r, 1, l);

	// MapleALGEB_Printf(kv, "e1:: %a\n\n\n", e1);
	n = EvalMapleStatement(kv, input2);
	e2 = EvalMapleProc(kv, r, 1, n);

	// MapleALGEB_Printf(kv, "e2:: %a\n\n\n", e2);
	
	free(input);
	free(input2);
	return !(mapleTest->testEquality(e1, e2)); //0 encodes same in this context
}

static void handle_parser_error(const char *message){
	fprintf(stderr, "%s :  [%s]: %s: @%d\n", message, __FILE__, __func__, __LINE__);
    exit(EXIT_FAILURE)	;
}



#endif  /*MAPLE_EVALUATE_H*/
