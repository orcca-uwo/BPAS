#ifndef __DFT_general_fermat_H
#define __DFT_general_fermat_H


extern void dft_general_fermat(GeneralizedFermatPrimeField* vector, int K, int e, GeneralizedFermatPrimeField omega_w);
extern void inverse_fermat_DFT(GeneralizedFermatPrimeField* vector, int K, int e, GeneralizedFermatPrimeField omega_w);
#endif
