

#include "malloc.h"
#include <string.h>

#include "FFT/src/fft_iter_main.h"


//forward declare from auto-generated file
void DFT_eff_H1024(int n, int r,long long int *A,long long int *W,long long int *B, const Prime_ptr* Pptr);
void InvDFT_eff_H1024(int n, int r,long long int *A,long long int *W,long long int *B, long long int invn, const Prime_ptr* Pptr);


void getNthRoots_inForm_spf(long long int n, long long int* omega, const Prime_ptr* Pptr) {
	long long int root = smallprimefield_PrimitiveRootofUnity(n, Pptr);
	root = smallprimefield_convert_in(root, Pptr);
	omega[0] = smallprimefield_convert_in(1ll, Pptr);
	omega[1] = root;
	for (int j = 2; j < n; ++j) {
		omega[j] = smallprimefield_mul(omega[j-1], root, Pptr);
	}
}


/*
*n = 2^r
*e.g. n = 16, contents of rootsTable:
* 1, w, w^2, w^3, w^4, w^5, w^6, w^7, w^8, w^9, w^10, w^11, w^12, w^13, w^14, w^15
* 1, w^2, w^4, w^6, w^8, w^10, w^12, w^14
* 1, w^4, w^8, w^12
* 1, w^8
*/
void RootsTable_spf(int n, int r,
	long long int *T, const Prime_ptr* Pptr)
{//---------------------------------------------------------------

	getNthRoots_inForm_spf(n, T, Pptr);

	long long int *T_ptr = T;
	int ni = n;
	for (int i=1; i<r; i++){
	//std::cout<<"i: "<<i<<"\n";
		T_ptr = T_ptr + ni;
		ni = n>>i; //n/2^i
		for (int j=0; j<ni; j++){
			T_ptr[j] = T[(j<<i)]; //j*2^i
		}
	}
}


// Ti is the reversal of T except the first column
/*
*n = 2^r
*e.g. n = 16, contents of inverserootsTable:
*  m, m^2, m^3, m^4, m^5, m^6, m^7, m^8, m^9, m^10, m^11, m^12, m^13, m^14, m^15
*  m^2, m^4, m^6, m^8, m^10, m^12, m^14
*  m^4, m^8, m^12
*  m^8
*/
void InverseRootsTable(int n, long long int *T, long long int *Ti)
{//-------------------------------------------------
//std::cout<<"--InverseRootsTable, n: "<<n<<std::endl;
	Ti[0] = T[0];
	int m=n;
//int end=n;
	int i=0;
	long long int *Ti_p = Ti;
	long long int *T_p = T;
	while(m>=2){
		Ti_p[0] = T_p[0];
		for (i=1; i<m; i++){
			Ti_p[m-i] = T_p[i];
		}
		Ti_p+=m;
		T_p +=m;
		m=(m>>1);
//end+=m;
	}
//std::cout<<end<<"\n";
}


void DFT_eff(int n, int r,
	       long long int *A,
	       long long int *W,
	       const Prime_ptr* Pptr,
	       int H,
	       long long int *B)
{
	switch (H) {
		case 1024 :
			DFT_eff_H1024(n, r, A, W, B, Pptr);
			break;
		default :
			fprintf(stderr, "DFT_eff not valid for H=%d\n", H);
			exit(1);
			break;
	}
}



void InvDFT_eff(int n, int r,
	  long long int *A,
	  long long int *W,
   	  const Prime_ptr* Pptr,
	  int H,
	  long long int *B,
	  long long int invn)
{
	switch (H) {
		case 1024 :
			InvDFT_eff_H1024(n, r, A, W, B, invn, Pptr);
			break;
		default :
			fprintf(stderr, "InvDFT_eff not valid for H=%d\n", H);
			exit(1);
			break;
	}
}


