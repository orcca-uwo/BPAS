/**
 * FME_Support_inequality.c
 * *************************
 * implementation of functions on inequality_t data-structre
 */
#include "PolyhedralSets/FME_Support_inequality.h"

int allocInequality(int varNum, inequality_t* newIneq)
{
	newIneq->dimension = varNum;

	newIneq->coeff = (mpz_t *) malloc(sizeof(mpz_t) * varNum);
	for(int i = 0 ; i < varNum ; i++)
		mpz_init(newIneq->coeff[i]);

	mpz_init(newIneq->constant);

	return(EXIT_SUCCESS);
}

void setInequality(inequality_t * newIneq , mpz_t * coeffData , mpz_t constantData)
{
	int varNum = newIneq->dimension;

	for(int i = 0 ; i < varNum ; i++)
		mpz_set(newIneq->coeff[i] , coeffData[i]);

	mpz_set(newIneq->constant , constantData);
}


void copyInequality(inequality_t * source , inequality_t * dest)
{
	for(int i = 0 ; i < source->dimension ; i++)
		mpz_set(dest->coeff[i] , source->coeff[i]);

	mpz_set(dest->constant , source->constant);
}

void printInequality(inequality_t * input , char var)
{
	for(int i = 0 ; i < input->dimension ; i++)
		gmp_printf("(%Zd) * %c%d +  ", input->coeff[i], var, i);
	gmp_printf(" <= %Zd\n", input->constant);
}


void getFromFile(mpz_t * data, char * fileName, int varNum, int ineqNum)
{
	FILE * file = fopen(fileName, "r");
	char * u = (char *) malloc(MAX_LEN * sizeof(char));

	rewind(file);
	char c;
	int x = 0;
    int counter = 0;

	while ((c = fgetc(file)) != EOF)
	{
		if (c != ' ')
			u[x++] = c;
		else
		{
            u[x++] = '\0';
			mpz_set_str(data[counter++], u , 10);
			x = 0;
		}
	}
	fclose(file);

	free(u);
}


void getInputArray(char * fileName , inequality_t * data , int varNum , int ineqNum)
{
	int len = (varNum + 1) * ineqNum;
	mpz_t * intData = (mpz_t *) malloc(sizeof(mpz_t) * len);

	for(int i = 0 ; i < len ; i++)
		mpz_init(intData[i]);

	for(int i = 0 ; i < ineqNum ; i++)
		allocInequality(varNum , &data[i]);

	getFromFile(intData , fileName , varNum , ineqNum);

    int counter = 0;
	for(int i = 0 ; i < ineqNum ; i++)
    {
		for(int j = 0 ; j < varNum ; j++)
			mpz_set(data[i].coeff[j] , intData[counter++]);
        mpz_set(data[i].constant , intData[counter++]);
    }

	for(int i = 0 ; i < len ; i++)
		mpz_clear(intData[i]);
	free(intData);
}

void freeInequality(inequality_t * a)
{
	for(int i = 0 ; i < a->dimension ; i++)
	{
		mpz_clear(a->coeff[i]);
	}
	free(a->coeff);
	mpz_clear(a->constant);
}
