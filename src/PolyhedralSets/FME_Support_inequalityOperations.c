/**
 * FME_Support_inequalityOperations.c
 * *************************
 * implementation of mathematical functions on inequality_t data-structre
 */
#include "PolyhedralSets/FME_Support_inequalityOperations.h"


int ineqSgn(inequality_t a, int varIndex)
{
	if (mpz_sgn(a.coeff[varIndex]) > 0)
		return (1);
	else if (mpz_sgn(a.coeff[varIndex]) < 0)
		return (-1);
	else
		return (0);
}

void combineTwoIneq(inequality_t ineq1, inequality_t ineq2,
		inequality_t * result, int varIndex)
{
	mpz_t mult1, mult2;
	mpz_init(mult1);
	mpz_init(mult2);
	
	findFMEMultiplier(ineq1 , ineq2 , varIndex , mult1 , mult2); //lcm

	inequality_t * result1 = (inequality_t *) malloc(sizeof(inequality_t));
	inequality_t * result2 = (inequality_t *) malloc(sizeof(inequality_t));
	inequality_t * result12 = (inequality_t *) malloc(sizeof(inequality_t));

	allocInequality(ineq1.dimension , result1);
	allocInequality(ineq2.dimension , result2);
	allocInequality(ineq2.dimension , result12);

	multScalarToIneq(&ineq1 , mult1, result1);
	multScalarToIneq(&ineq2 , mult2, result2);
	addTwoIneq(*result1 , *result2 , result12);

	simplifyIneq(result12 , result);

	/// cleaning

	mpz_clear(mult1);
	mpz_clear(mult2);
	freeInequality(result1);
	freeInequality(result2);
	freeInequality(result12);
	free(result1);
	free(result2);
	free(result12);
}

void findFMEMultiplier(inequality_t ineq1, inequality_t ineq2,
		int varIndex , mpz_t result1 , mpz_t result2)
{
	mpz_t temp1,temp2;
	mpz_init(temp1);
	mpz_init(temp2);
	
	mpz_abs(temp1 , ineq1.coeff[varIndex]);
	mpz_abs(temp2 , ineq2.coeff[varIndex]);
	
	findCommenMult(temp1, temp2 , result1 , result2);

	mpz_clear(temp1);
	mpz_clear(temp2);
}

void findCommenMult(mpz_t in1, mpz_t in2 , mpz_t mult1 , mpz_t mult2)
{
	mpz_t lcm;
	mpz_init(lcm);

	mpz_lcm(lcm, in1, in2);

	mpz_div(mult1 , lcm , in1);
	mpz_div(mult2 , lcm , in2);

	mpz_clear(lcm);
}

void multScalarToIneq(inequality_t * a, mpz_t mult, inequality_t * result)
{
	for (int i = 0; i < a->dimension; i++)
		mpz_mul(result->coeff[i], a->coeff[i], mult);
	mpz_mul(result->constant , a->constant , mult);
}


void addTwoIneq(inequality_t ineq1, inequality_t ineq2,
		inequality_t * result)
{
	mpz_t tmp;
	mpz_init(tmp);

	for (int i = 0; i < ineq1.dimension; i++)
	{
		mpz_add(tmp , ineq1.coeff[i] , ineq2.coeff[i]);
		mpz_set(result->coeff[i], tmp);
	}
	mpz_add(result->constant , ineq1.constant , ineq2.constant);

	mpz_clear(tmp);
}

void simplifyIneq(inequality_t * ineq1 , inequality_t * result)
{
	mpz_t gcd;
	mpz_init(gcd);

	mpz_t *tmp = (mpz_t*) malloc(sizeof(mpz_t) * (ineq1->dimension+1));
	for(int i = 0 ; i < ineq1->dimension+1 ; i++)
		mpz_init(tmp[i]);

	int c = 0;
	for(int i = 0 ; i < ineq1->dimension ; i++)
		if(mpz_sgn(ineq1->coeff[i]) != 0)
		{
			mpz_set(tmp[c] , ineq1->coeff[i]);
			mpz_abs(tmp[c] , tmp[c]);
			c++;
		}

	if(mpz_sgn(ineq1->constant) != 0)
	{
		mpz_set(tmp[c] , ineq1->constant);
		mpz_abs(tmp[c] , tmp[c]);
	}

	mpz_gcd(gcd , tmp[0] , tmp[1]);

	for(int i = 2 ; i < c+1 ; i++)
		mpz_gcd(gcd , gcd , tmp[i]);


	for(int i = 0 ; i < ineq1->dimension ; i++)
		mpz_div(result->coeff[i] , ineq1->coeff[i] , gcd);
	mpz_div(result->constant , ineq1->constant , gcd);

	mpz_clear(gcd);
	for(int i = 0 ; i < ineq1->dimension+1 ; i++)
		mpz_clear(tmp[i]);
	free(tmp);

}
