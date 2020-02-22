#ifndef MAPLE_EVALUTE_H
#define MAPLE_EVALUTE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <maplec.h>

/**
 * @brief callback used for directing result output
 * 
 * @param data 
 * @param tag 
 * @param output 
 */
static void M_DECL textCallBack( void *data, int tag, char *output );
/**
 * @brief Evaluate polynomial input to a Maple and return
 * result as a string.
 * 
 * @param in - maple input
 * @return char* - result output as string
 */
char *maple_evaluate(const char *in);
/**
 * @brief Evaluate and compare the result of two maple polynomial
 * expansion, and return 0 if they are equal. The comparison is 
 * done by subtracting the actual with expected input.
 * 
 * @param actual 
 * @param expected 
 * @return int 
 */
int maple_evaluate_equal(const char *actual, const char *expected);
/**
 * @brief Handle fatal error, by exiting the program.
 * 
 * @param message 
 */
void handle_parser_error(const char *message);





static void M_DECL textCallBack( void *data, int tag, char *output )
{
    printf("%s\n",output);
}

char *maple_evaluate(const char *in){
	char *pch = strchr(in, ']');
	char *input = (char*)malloc((strlen(in)+1)*sizeof(char));
	if(input == NULL){
		handle_parser_error("Unable to create heap memory.");
	}

	strcpy(input, in+(pch-in+1));
	char *col = ":";
	strncat(input, col, strlen(col));
	char err[2048];  /* command input and error string buffers */
    MKernelVector kv;  /* Maple kernel handle */
    MCallBackVectorDesc cb = {  (void*)textCallBack, 
				0,   /* errorCallBack not used */
				0,   /* statusCallBack not used */
				0,   /* readLineCallBack not used */
				0,   /* redirectCallBack not used */
				0,   /* streamCallBack not used */
			        0,   /* queryInterrupt not used */ 
				0    /* callBackCallBack not used */
			    };
    ALGEB r, l, m, n , p, q, s, t, e1, e2;  /* Maple data-structures */

    /* initialize Maple */
    if( (kv=StartMaple(0,NULL,&cb,NULL,NULL,err)) == NULL ) {
		printf("Fatal error, %s\n",err);
		exit(EXIT_FAILURE);
    }
    r = EvalMapleStatement(kv,"expand:");
	l = EvalMapleStatement(kv,input);
	m = ToMapleFunction(kv, r, 1, l);
	e1 = MapleEval(kv, m);
	ALGEB msg = MapleALGEB_SPrintf(kv, "%a", e1);
	return MapleToString(kv, msg);
}

int maple_evaluate_equal(const char *actual, const char *expected){
	char *input = (char*)malloc((strlen(actual)+1)*sizeof(char));
	if(input == NULL){
		handle_parser_error("Unable to create heap memory.");
	}
	char *input2 = (char*)malloc((strlen(expected)+1)*sizeof(char));
	if(input2 == NULL){
		handle_parser_error("Unable to create heap memory.");
	}

	if((strchr(actual, '[')) != NULL){
		char *pch = strchr(actual, ']');
		strcpy(input, actual+(pch-actual+1));
	}else{
		strcpy(input, (actual));
	}
	strcpy(input2, expected);

	char *col = ":";
	strncat(input, col, strlen(col));
	strncat(input2, col, strlen(col));

	// fprintf(stdout, "input: %s\n", input);
	// fprintf(stdout, "input2: %s\n", input2);

	char err[2048];  /* command input and error string buffers */
    MKernelVector kv;  /* Maple kernel handle */
    MCallBackVectorDesc cb = {  (void*)textCallBack, 
				0,   /* errorCallBack not used */
				0,   /* statusCallBack not used */
				0,   /* readLineCallBack not used */
				0,   /* redirectCallBack not used */
				0,   /* streamCallBack not used */
			    0,   /* queryInterrupt not used */ 
				0    /* callBackCallBack not used */
			    };
    ALGEB r, l, m, n , p, q, s, t, e1, e2;  /* Maple data-structures */

    /* initialize Maple */
    if( (kv=StartMaple(0,NULL,&cb,NULL,NULL,err)) == NULL ) {
		printf("Fatal error, %s\n",err);
		exit(EXIT_FAILURE);
    }
    r = EvalMapleStatement(kv,"expand:");
	l = EvalMapleStatement(kv,input);
	m = ToMapleFunction(kv, r, 1, l);
	e1 = MapleEval(kv, m);

	// MapleALGEB_Printf(kv, "e1:: %a", e1);

	n = EvalMapleStatement(kv, input2);
	q = ToMapleFunction(kv, r, 1, n);
	e2 = MapleEval(kv, n);

	// MapleALGEB_Printf(kv, "e2:: %a", e2);
	
	free(input);
	free(input2);
	return e1-e2;
}

void handle_parser_error(const char *message){
	fprintf(stderr, "%s :  [%s]: %s: @%d\n", message, __FILE__, __func__, __LINE__);
    exit(EXIT_FAILURE);
}



#ifdef __cplusplus
}
#endif


#endif  /*MAPLE_EVALUATE_H*/
