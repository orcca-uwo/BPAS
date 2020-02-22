#ifndef __DFT32_H
#define __DFT32_H

template <class FiniteField>
extern FiniteField* DFT_32(FiniteField* a ,FiniteField omega);

template <class FiniteField>
extern  FiniteField* DFT_32(FiniteField* a,FiniteField* omegas);

//long int* DFT_32(long int* a,long int* omegas, long int prime,long long int R,long long int pP);
//GeneralizedFermatPrimeField* DFT_32(GeneralizedFermatPrimeField* A, GeneralizedFermatPrimeField* omegas);

#endif
