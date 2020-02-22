#include "parser_helper.h"


powervar *create_power_var(char *var, degree_t exp){
    powervar *temp_powervar = (powervar*)malloc(sizeof(powervar));
    if(temp_powervar == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
	temp_powervar->var = (char*)malloc(strlen(var)*sizeof(char)+1);
    if(temp_powervar->var == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
	memcpy(temp_powervar->var, var, strlen(var)+1);
	temp_powervar->exp = exp;

	return temp_powervar;
}


term *create_term(int arraySize){
	term* temp_term = (term*)malloc(sizeof(term));
    if(temp_term == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
    temp_term->exp = calloc(arraySize, sizeof(degree_t));
    if(temp_term->exp == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
    mpq_init(temp_term->coef);
    mpq_set_str(temp_term->coef, "1", 10);
	return temp_term;
}

degree_t* generate_long_int_exponents(int arraySize){
    degree_t* temp_exponents = calloc(arraySize, sizeof(unsigned long int));
    if(!temp_exponents){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }

    return temp_exponents;
}

void free2DArray(char** arr, int size){
    if(arr != NULL){
        for(int i=0; i<size; i++)
            free(arr[i]);
        free(arr);
    }
}

void fill_term_exponent(char **orderedVarArray, degree_t *exponentToFill, char* currentVar, degree_t currentExp, int *arraySize){
    int pos = -1;
    //The for loop is a bottle neck, unfortunately it is required to reset the array because
    //of a wierd bug in the exponentToFill array. he problem is, instead of keeping exponentToFill
    //filled with zero, the program fills the certain index with random number if thier length is
    //greater than 3. This index (if they are greater than 3) will be filled with the value 1 or 33
    //instead of 0.
    // for(int i=0; i<numVariables; i++){
    //     exponentToFill[i] = 0;
    // }
    for(int i=0; i<(*arraySize); i++){
        if(strcmp(currentVar, orderedVarArray[i])==0){ //unsafe function
            pos = i;
            break;
        }
    }
    if(pos >= 0){
        if(exponentToFill[pos] > 0){//added may 4 , when input [x] x*x only returned x
            exponentToFill[pos] += currentExp;
        }else{
            exponentToFill[pos] = currentExp;
        }
    }else{
        // TODO: printf("new varaiable is discovered\n"); 
        //if the variable used unspecified, mostly the code resides here 
        //to enumirate through the array to update the exponents
        // printf("variable is not found ::)\n");
        // This condition is useful if the data structure linked-list
	}
}

int check_if_it_exists(char** vars, char* v, int numvar){
    if(vars){
        for(int i=0; i<numvar; ++i){
            if(strcmp(vars[i], v)==0){
                return 1;
            }
        }
        return 0;
    }
    return 0;
}


char** push_back_dynamic(char **varArray, int *numvar, char* var){
    if(!check_if_it_exists(varArray, var, *numvar)){
        (*numvar)++;
        if(!varArray){
            varArray = (char**)malloc((*numvar)*sizeof(char*));
        }else{
            varArray = realloc(varArray, (*numvar)*sizeof(char*));
        }
        varArray[(*numvar)-1] = (char*)malloc(strlen(var)+1*sizeof(char));
        memcpy(varArray[(*numvar)-1], var, strlen(var)+1);
    }

    return varArray;
}


void deep_degrees_cpy(degree_t *dst, degree_t *src, int size){
    for(int i=0; i<size; i++){
        dst[i] = src[i];
    }
}


void clear_term(term *tempTerm, int numVar){
	for(int i=0; i<numVar; i++){
        tempTerm->exp[i] = 0;
    }
}


void free_term(term* temp){
	free(temp->exp);
    mpq_clear(temp->coef);
	free(temp);
}


void add_unpacked_term_to_smqp_aa(AltArr_t* aa, const degree_t* deg, ratNum_t coef, int numVar) {
    addTerm_AA(aa, deg, numVar, coef);
}

// void add_packed_degree_term_to_smqp_aa(AltArr_t *aa, degree_t *deg, ratNum_t coef, int numvar){
    
//     if(numvar == 0){ //will dump error message in the console complaining EXP and NVARS are empty
//         numvar = 1;
//     }
//     int *size = getExpOffsetArray(numvar);
//     degree_t d = 0;

//     for(int i = 0; i < numvar; ++i){
//         d |= (deg[i] << size[i]);
//     }

//     addTermSafe_AA(aa, d, coef);
//     free(size);
// }


degree_t* unpacked_degs(degrees_t packed_degs, int numvars){
    if(numvars == 0){ //will dump error message in the console complaining EXP and NVARS are empty
        numvars = 1;
    }
    const degrees_t* __restrict__ oldmasks = getExpMaskArray(numvars);
    const int* __restrict__ oldsizes = getExpOffsetArray(numvars);
    
    degree_t* ret = (degree_t*)calloc(numvars, sizeof(degree_t));
    if(!ret){
        fprintf(stderr, "%s\n", "error callocing");
        return NULL;
    }
    for(int j=0; j<numvars; j++){
        ret[j] = GET_NTH_EXP(packed_degs, oldmasks[j], oldsizes[j]);
    } 
    return ret;
}

char* strip_comments(const char* buf, size_t size){
    char *newbuf = (char*)malloc(size*sizeof(char));
    int i =0;
    int j = 0;
    while(i<size){
        if(buf[i]=='/'&&buf[i+1]=='*'){
            i = i+2;
            while(i<size){
                if(buf[i]=='*' && buf[i+1]=='/'){
                    break;
                }
                i++;
            }                
            i = i+2;            
        }
        newbuf[j] = buf[i];
        j++;
        i++;
    }
    
    return newbuf;
}