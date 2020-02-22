#ifndef __DFT_general_H
#define __DFT_general_H

template <class FiniteField>
extern void DFT_general(FiniteField* vector, int K, int e,FiniteField omega_w);
template <class FiniteField>
extern void inverse_DFT(FiniteField* vector, int K, int e,FiniteField omega_w);

template<class FiniteField>
extern void precomputed_DFT_general(FiniteField* vector, int K, int e, FiniteField omega_w, FiniteField** powersofomega_outer, FiniteField* powersofomega_inner);

extern long int FFTPrime(int exponent);


#endif
