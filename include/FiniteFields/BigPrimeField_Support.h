#ifndef _BIGPRIMEFIELD_SUPPORT_H_
#define _BIGPRIMEFIELD_SUPPORT_H_

// #ifdef __cplusplus
// extern "C" { 
// #endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <unistd.h>

void BigPrimeField_Add_inplace(mpz_t &X, mpz_t &Y, mpz_t &P);

void BigPrimeField_Sub_inplace(mpz_t &X, mpz_t &Y, mpz_t &P);

void BigPrimeField_Multiplication_inplace(mpz_t &X, mpz_t &Y, mpz_t &P);

void BigPrimeField_Division_inplace(mpz_t &X, mpz_t &Y, mpz_t &P);

void BigPrimeField_Inverse_inplace(mpz_t &X, mpz_t &P);

void BigPrimeField_Add(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z);

void BigPrimeField_Sub(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z);

void BigPrimeField_Multiplication(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z);

void BigPrimeField_Division(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z);

void BigPrimeField_Inverse(mpz_t &X, mpz_t &P, mpz_t &Z);
#endif