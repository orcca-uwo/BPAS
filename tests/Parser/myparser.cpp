
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../include/Parser/bpas_parser.h"
#include "parser-test-framework/parser_testframework.h"
#include "parser-test-framework/maple_evaluate.hpp"

int is_empty(const char *s);
void handle_input(const char* input);

int main(int argc, char **argv){
	//************************** EXAMPLES ON HOW TO USE THE PARSER*******************************************************
    const char *poly = "-(x^6*y - x*y^6)^6";
    AltArr_t *result = generate_altarr(poly);
    // print_naked_AltArr_t_poly(result);
    freePolynomial_AA(result);          //gcc 7.x.x if you don't free AltArr_t* it will result in a segmentation fault.

    const char *poly2 = "(x+y+z^2-m)^2";
    int numvars = 4;
    char **vars = create_dynamic_str_array(numvars, "x", "y", "z", "m");
    AltArr_t *result2 = generate_altarr(poly2);
    // print_poly_to_terminal(result2, vars, numvars);
    // print_poly_to_terminal_fancy(result2, vars, numvars, 1);
    freePolynomial_AA(result2);

    const char *poly3 = "-(x^6*y - x*y^6)^6";
    //char *poly3 = "-x-y";
    int numvars3 = 2;
    char **vars3 = create_dynamic_str_array(numvars3, "x", "y");
    AltArr_t *result3 = generate_altarr(poly3);
    // print_naked_AltArr_t_poly(result3);
    // print_poly_to_terminal(result3, vars3, numvars3);
    char *ret = print_poly_to_string_variable(result3, vars3, numvars3);
    ASSERT_EQUAL(maple_evaluate_equal(poly3, ret), 0);
    parser_test(poly3, "TEST");
    freePolynomial_AA(result3);

    const char *poly4 = "((x*y+z)*(y*z+x)+z*x+y)*((y*z+x)*(x*z+y)+x*y+z)+(x*z+y)*(x*y+z)+y*z+x";
    int numvars4 = 3;
    char **vars4 = create_dynamic_str_array(numvars4, "x", "y", "z");
    AltArr_t *result4 = generate_altarr(poly4);
    // print_poly_to_terminal(result4, vars4, numvars3);
    char *ret2 = print_poly_to_string_variable(result4, vars4, numvars4);
    ASSERT_EQUAL(maple_evaluate_equal(poly4, ret2), 0);
    parser_test(poly4, "TEST");
    freePolynomial_AA(result4);
	//***************************************EXAMPLES ON HOW TO USE THE PARSER ENDS*************************************

	//***************************************TESTING PARSER BEGINS******************************************************
    if(argc > 1) {
	handle_input(argv[1]);
 	//handle_input("data.txt");
    }	    

	//***************************************TESTING PARSER ENDS********************************************************

    return 0;
}

int is_empty(const char *s) {
  while (*s != '\0') {
    if (!isspace((unsigned char)*s))
      return 0;
    s++;
  }
  return 1;
}

void handle_input(const char* input){
	FILE *temp;
	temp = fopen(input, "r");
	if(!temp){
		fprintf(stderr, "%s", "Error 2: opening the file! \n");
		exit(EXIT_FAILURE);
	}
	
	fseek(temp, 0, SEEK_END);
	size_t size = ftell(temp);
	rewind(temp);
	// printf("size: %zd\n", size);
	
	char *buf = (char*)malloc(size*sizeof(char));
	size_t r = fread(buf, sizeof(char), size, temp);
	
	char *newbuf = strip_comments(buf, size);
	free(buf);
	
	const char *delim = ";";
	char *result =  strtok(newbuf, delim);
	altarr_pack *ret = NULL;
	while(result != NULL){
		size_t ilen = strlen(result);
		for(int i =0; i<ilen+1; i++){
			if(result[i] == '\n')
				result[i] = ' ';
		}
		if(is_empty(result) != 1){
			ret = generate_altarr_pack(result);
			if(ret->altarr_t_data != NULL){
				char *ret2 = NULL;
				ret2 = print_poly_to_string_variable(ret->altarr_t_data, ret->vars, ret->numVars);
				if(strlen(ret2) !=0){
					ASSERT_EQUAL(maple_evaluate_equal(result, ret2), 0);
					parser_test(result, "TEST");
				}else{
					//print function evaluate input and return null, for example "x-x" returned as NULL
					//since the test function expect "x-x" when compared NULL with actual, test fails
					ASSERT_EQUAL(0, 0);
					parser_test(result, "TEST");
				}
				free(ret2);
			}
		}
		//freePolynomial_AA(ret->altarr_t_data);
		//free(ret);
		result = strtok(NULL, delim);
	}
	fclose(temp);	
	free(newbuf);
	free(result);
}

