#ifndef __DFT16_H
#define __DFT16_H

template <class FiniteField>
extern FiniteField* DFT_16(FiniteField* a ,FiniteField omega);

template <class FiniteField>
extern  FiniteField* DFT_16(FiniteField* a,FiniteField* omegas);

//long int* DFT_16(long int* a,long int* omegas, long int prime,long long int R,long long int pP);
//GeneralizedFermatPrimeField* DFT_16(GeneralizedFermatPrimeField* a, GeneralizedFermatPrimeField* omegas);

#endif
