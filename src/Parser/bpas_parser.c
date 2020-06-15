#include "../include/Parser/bpas_parser.h"

#include "../../src/Parser/parser_helper.h"
#include "../../src/Parser/parser_print.h"
#include "../../src/Parser/parser_errno.h"
#include "../../src/Parser/parser_grammar.tab.h"


#include <pthread.h>

extern AltArr_t *altarr_data;
extern char** g_variables;
extern int g_num_variables;
extern int is_all_var_defined;
extern void yy_scan_string(char*);
extern int yylex_destroy  (void);

static pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

static int total_var_size(const char** vars, const int numVars){
    int size = 0;
    for(int i=0; i<numVars; i++){
        size += strlen(vars[i]);
    }
    return size;
}

static char *create_var_list(const char **vars, const int numVars){
    int size = total_var_size(vars, numVars);
    int t_size = size + (numVars-1) + 3;
    char *var_list = (char*)calloc(t_size, sizeof(char));
   var_list[0] = '\0';
    if(!var_list){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }

    strncat(var_list, "[", strlen("[")+1);
    for(int i=0; i<numVars; i++){
        strncat(var_list, vars[i], strlen(vars[i]));
        if(i != (numVars-1)){
            strncat(var_list, ",", strlen(",")+1);
        }
    }
    strncat(var_list, "]", strlen("]")+1);
    strncat(var_list, "\0", strlen("\0")+1);
    return var_list;
}

char **create_dynamic_str_array(int arrSize, ...){
    char **temp_str_arr = (char**)malloc(sizeof(char*)*arrSize);
    if(!temp_str_arr){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }
    va_list args;
    va_start(args, arrSize);
    for(int i=0; i<arrSize; i++){
        char *temp_var = strdup(va_arg(args, char*));
        if(!temp_var){
            parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
            exit(PARSER_NOALLOC);
        }
        temp_str_arr[i] = (char*)malloc(sizeof(char)*strlen(temp_var)+1);
        strncpy(temp_str_arr[i], temp_var, strlen(temp_var)+1);
    }
    va_end(args);

    return temp_str_arr;
}

AltArr_t* generate_altarr(const char* poly_str){
pthread_mutex_lock( &mutex1 );
    char *temp_polystr = strdup(poly_str);

    yy_scan_string(temp_polystr);
    yyparse();
    yylex_destroy();

    free(temp_polystr);

    AltArr_t* temp = altarr_data;
    altarr_data = NULL;
    for (int i = 0; i < g_num_variables; ++i) {
        free(g_variables[i]);
    }
    free(g_variables);
    g_variables = NULL;
    g_num_variables = 0;
    is_all_var_defined = 0;


pthread_mutex_unlock( &mutex1 );
    return temp;
}

altarr_pack* generate_altarr_pack(const char* poly_str){
pthread_mutex_lock( &mutex1 );
    char *temp_polystr = strdup(poly_str);
    altarr_pack *pack = (altarr_pack*)malloc(sizeof(altarr_pack));

    yy_scan_string(temp_polystr);
    yyparse();
    yylex_destroy();

    pack->numVars = g_num_variables;
    g_num_variables = 0;

    pack->altarr_t_data = altarr_data;
    altarr_data = NULL;

    pack->vars = g_variables;
    g_variables = NULL;

    is_all_var_defined = 0;

    free(temp_polystr);
pthread_mutex_unlock( &mutex1 );
    return pack;
}

AltArr_t* generate_altarr_var_defined(const char* poly_str, const char** variables, int num_var){
pthread_mutex_lock( &mutex1 );

    char *temp_var_list = create_var_list(variables, num_var);
    char *temp_poly = (char*)malloc(sizeof(char)*(strlen(temp_var_list) + strlen(poly_str) + 1));
    if(!temp_poly){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }
    //given var (x,y,z) are copied to temp_poly as [x,y,z] and concat the poly_str
    // oct 31,2018 bug fix
    strncpy(temp_poly, temp_var_list, strlen(temp_var_list)+1);
    strncat(temp_poly, poly_str, strlen(poly_str)+1); //suppose to be strncat not strcpy

    yy_scan_string(temp_poly);
    yyparse();
    yylex_destroy();

    free(temp_var_list);
    free(temp_poly);

    AltArr_t* temp = altarr_data;
    altarr_data = NULL;
    for (int i = 0; i < g_num_variables; ++i) {
        free(g_variables[i]);
    }
    free(g_variables);
    g_variables = NULL;
    g_num_variables = 0;
    is_all_var_defined = 0;

pthread_mutex_unlock( &mutex1 );

    return temp;
}









/* AltArr_t* generate_altarr_var_defined(const char* poly_str, const char** variables, int num_var){
    int var_size = 0;
    for(int i=0; i<num_var; i++){
        var_size += sizeof(variables[i]);
    }
    char *temp_var = (char*)calloc((num_var+var_size+2), sizeof(char));
    if(temp_var == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }
    strncat(temp_var, "[", sizeof("["));
    for(int i=0; i<num_var; i++){
        strncat(temp_var, variables[i], sizeof(variables[i]));
        if(i != (num_var-1)){
            strncat(temp_var, ",", sizeof(","));
        }
    }
    strncat(temp_var, "]", sizeof("]"));

    char *temp_str = (char*)malloc((sizeof(temp_var)+sizeof(poly_str)+1)*sizeof(char));
    if(temp_str == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }

    char* temp_polystr = strdup(poly_str);
    strncat(temp_str, temp_var, sizeof(temp_var));
    strncat(temp_str, temp_polystr, sizeof(poly_str));
    strncat(temp_str, "\0", sizeof("\0"));

    yy_scan_string(temp_str);
    yyparse();

    free(temp_polystr);
    free(temp_var);
    free(temp_str);

    return altarr_data;
} */
