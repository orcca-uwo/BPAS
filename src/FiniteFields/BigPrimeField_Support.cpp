#include "FiniteFields/BigPrimeField_Support.h"

void BigPrimeField_Add_inplace(mpz_t &X, mpz_t &Y, mpz_t &P)
  {
    mpz_add(X, X, Y);
    mpz_mod(X, X, P);
  }
void BigPrimeField_Sub_inplace(mpz_t &X, mpz_t &Y, mpz_t &P)
  {
    mpz_sub(X, X, Y);
    mpz_mod(X, X, P);
  }
void BigPrimeField_Multiplication_inplace(mpz_t &X, mpz_t &Y, mpz_t &P)
  {
    mpz_mul(X, X, Y);
    mpz_mod(X, X, P);
  }

void BigPrimeField_Division_inplace(mpz_t &X, mpz_t &Y, mpz_t &P)
  {
    mpz_invert (Y, Y, P);
    mpz_mul(X, X, Y);
    mpz_mod(X, X, P);
  }

void BigPrimeField_Inverse_inplace(mpz_t &X, mpz_t &P)
  {
    mpz_invert(X,X,P);
  }

void BigPrimeField_Add(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z)
  {
    mpz_add(Z, X, Y);
    mpz_mod(Z, Z, P);
  }
void BigPrimeField_Sub(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z)
  {
    mpz_sub(Z, X, Y);
    mpz_mod(Z, Z, P);
  }
void BigPrimeField_Multiplication(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z)
  {
    mpz_mul(Z, X, Y);
    mpz_mod(Z, Z, P);
  }

void BigPrimeField_Division(mpz_t &X, mpz_t &Y, mpz_t &P, mpz_t &Z)
  {
    mpz_invert (Z, Y, P);
    mpz_mul(Z, X, Z);
    mpz_mod(Z, Z, P);
  }

void BigPrimeField_Inverse(mpz_t &X, mpz_t &P, mpz_t &Z)
  {
    mpz_invert(Z,X,P);
  }
