#ifndef __DFT8_H
#define __DFT8_H

//#include "../../../include/FiniteFields/GeneralizedFermatPrimeField.hpp"

template <class FiniteField>
extern  FiniteField* DFT_8(FiniteField* a,FiniteField omega);

template <class FiniteField>
extern  FiniteField* DFT_8(FiniteField* a,FiniteField* omegas);

//GeneralizedFermatPrimeField* DFT_8(GeneralizedFermatPrimeField* a,GeneralizedFermatPrimeField* omegas);
//long int* DFT_8(long int* a,long int* omegas, long int prime,long long int R,long long int pP);

#endif
