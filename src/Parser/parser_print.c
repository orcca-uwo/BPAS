#include "parser_print.h"

void print_naked_AltArr_t_poly(AltArr_t* aa) {
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		gmp_fprintf(stdout, "%Qd*%llx + ", aa->elems[i].coef, aa->elems[i].degs);
	}
	fprintf(stdout, "\n");
}


void print_naked_AltArr_t_poly2(AltArr_t* aa, char **vars, int numvars, char *msg) {
	for (int i = 0; i < aa->size; ++i) {
        degree_t *deg = unpacked_degs(aa->elems[i].degs, aa->nvar);
		gmp_fprintf(stdout, "%s: %Qd ", msg, aa->elems[i].coef);
        for(int j=0; j<aa->nvar; j++){
            fprintf(stdout, "%s^%d", vars[j], deg[j]);
        }
        free(deg);
	}
	fprintf(stdout, "\n");
}


void print_poly_to_terminal(AltArr_t *l, char** var, int numVar ){
    AltArr_t *temp = l;
    fprintf(stdout, ">>> ");
    if(temp == NULL){
        fprintf(stdout, "%s", "0");
    }
    
    int first = 1;
    int needMult = 0;
    int isConst = 1;
    ratNum_t coef;
    mpq_init(coef);
    for (int i = 0; i < AA_SIZE(temp); i++){
        mpq_set(coef, temp->elems[i].coef);
        isConst = 1;
        if(mpq_cmp_si(temp->elems[i].coef, 0,1) < 0){
            // fprintf(stdout, "%s", "-");
        }else if(!first){
            fprintf(stdout, "%s", "+");
        }

        if(mpq_cmp_si(temp->elems[i].coef, 1,1) != 0 ){
            if(mpq_cmp_si(temp->elems[i].coef, -1,1) == 0){
                fprintf(stdout, "%s", "-");
                needMult = 0;
            }else{
                gmp_fprintf(stdout, "%Qd", temp->elems[i].coef);
                needMult = 1;
            }
        }
        degree_t *expon = unpacked_degs(temp->elems[i].degs, temp->nvar);
        for(int j=0; j<temp->nvar; j++){
            if(expon[j] == 0){
                continue;
            }
            isConst = 0;
            if(needMult){
                fprintf(stdout, "%s", "*");
            }
            fprintf(stdout, "%s", var[j]);
            if(expon[j] > 1){
                fprintf(stdout, "^%d", expon[j]);
            }
            needMult = 1;
        }
        // temp = temp->next;
        first = 0;
        needMult = 0;
    }
    if(isConst && (mpq_cmp_si(coef, 1,1) == 0)){
        gmp_fprintf(stdout, "%Qd", coef);
    }
    fprintf(stdout, "%s", "\n");
}


void print_poly_to_terminal_fancy(AltArr_t *l, char** var, int numVar, int term_vars_separated){
    AltArr_t *temp = l;
    fprintf(stdout, ">>> ");
    if(temp == NULL){
        fprintf(stdout, "%s", "0");
    }
    
    int first = 1;
    int needMult = 0;
    int isConst = 1;
    ratNum_t coef;
    mpq_init(coef);
    for (int i = 0; i < AA_SIZE(temp); i++){
        mpq_set(coef, temp->elems[i].coef);
        isConst = 1;
        if(mpq_cmp_si(temp->elems[i].coef, 0,1) < 0){
            // fprintf(stdout, "%s", "-");
        }else if(!first){
            if(term_vars_separated){
                fprintf(stdout, "%s", " + ");
            }else{
                fprintf(stdout, "%s", "+");
            }
        }

        if(mpq_cmp_si(temp->elems[i].coef, 1,1) != 0 ){
            if(mpq_cmp_si(temp->elems[i].coef, -1,1) == 0){
                if(term_vars_separated){
                    fprintf(stdout, "%s", " - ");
                }else{
                    fprintf(stdout, "%s", "-");
                }
                needMult = 0;
            }else{
                gmp_fprintf(stdout, "%Qd", temp->elems[i].coef);
                needMult = 1;
            }
        }
        degree_t *expon = unpacked_degs(temp->elems[i].degs, temp->nvar);
        for(int j=0; j<temp->nvar; j++){
            if(expon[j] == 0){
                continue;
            }
            isConst = 0;
            if(term_vars_separated){
                if(needMult){
                    fprintf(stdout, "%s", " ");
                }
            }
            fprintf(stdout, "%s", var[j]);
            if(expon[j] > 1){
                fprintf(stdout, "%s", generate_super(expon[j]));
            }
            needMult = 1;
        }
        // temp = temp->next;
        first = 0;
        needMult = 0;
    }
    if(isConst && (mpq_cmp_si(coef, 1,1) == 0)){
        gmp_fprintf(stdout, "%Qd", coef);
    }
    fprintf(stdout, "%s", "\n");
}


void print_poly_to_file(AltArr_t *l, char** var, int numVar, char *filename, int insertVarsAtTheBeginning){
    AltArr_t *temp = l;

    FILE *f = fopen(filename, "w");
    if(!f){
        fprintf(stderr, "%s  [%s]: %s: @%d\n", "File can not be opened for writing!: ", __FILE__, __func__, __LINE__);
        return;
    }
    if(insertVarsAtTheBeginning){
        fprintf(f, "[");
        for(int i=0; i<numVar; i++){
            fprintf(f, "%s",var[i]);
            if(i != numVar-1){
                fprintf(f, ",");
            }
        }
        fprintf(f, "]");
    }

    if(temp == NULL){
        fprintf(f, "%s", "0");
    }
    
    int first = 1;
    int needMult = 0;
    int isConst = 1;
    ratNum_t coef;
    mpq_init(coef);
    for (int i = 0; i < AA_SIZE(temp); i++){
        mpq_set(coef, temp->elems[i].coef);
        isConst = 1;
        if(mpq_cmp_si(temp->elems[i].coef, 0,1) < 0){
            // fprintf(f, "%s", "-");
        }else if(!first){
            fprintf(f, "%s", "+");
        }

        if(mpq_cmp_si(temp->elems[i].coef, 1,1) != 0 ){
            if(mpq_cmp_si(temp->elems[i].coef, -1,1) == 0){
                fprintf(f, "%s", "-");
                needMult = 0;
            }else{
                gmp_fprintf(f, "%Qd", temp->elems[i].coef);
                needMult = 1;
            }
        }
        degree_t *expon = unpacked_degs(temp->elems[i].degs, temp->nvar);
        for(int j=0; j<temp->nvar; j++){
            if(expon[j] == 0){
                continue;
            }
            isConst = 0;
            if(needMult){
                fprintf(f, "%s", "*");
            }
            fprintf(f, "%s", var[j]);
            if(expon[j] > 1){
                fprintf(f, "^%d", expon[j]);
            }
            needMult = 1;
        }
        // temp = temp->next;
        first = 0;
        needMult = 0;
    }
    if(isConst && (mpq_cmp_si(coef, 1,1) == 0)){
        gmp_fprintf(f, "%Qd", coef);
    }
    fprintf(f, "%c", '\0');
    fclose(f);
}


char* print_poly_to_string_variable(AltArr_t *altarr, char** var, int numVar ){
    char *filename = "tempFileToConverToStringVar.txt";
    print_poly_to_file(altarr, var, numVar, filename, 0);
    FILE *f = fopen(filename, "r");
    if(!f){
        fprintf(stderr, "%s  [%s]: %s: @%d\n", "File can not be opened for reading!: ", __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }

    fseek(f, 0, SEEK_END);
	size_t size = ftell(f);
	rewind(f);

    // char *buf = (char*)malloc(size*sizeof(char)+1);
    char *buf = (char*)calloc(size, sizeof(char));
	size_t r = fread(buf, sizeof(char), size, f);
    // buf[size+1] = '\0';
    // printf("printing buffer INSIDE : %s\n", buf);
    fclose(f);
    remove(filename);
    return buf;
}

void print_term(term* t, char **vars, int numvars, char* message){
    gmp_printf("%s: %Qd ", message, t->coef);
    for(int i=0; i<numvars; ++i){
        printf("%s^%d ", vars[i], t->exp[i]);
    }
    printf("\n");
}