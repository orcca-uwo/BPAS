/**
 * FME_Support_blockElimination.c
 * *******************************
 * implementation of block elimination to get the
 * projection of the redundancy test cone, W0
 */
#include "PolyhedralSets/FME_Support_blockElimination.h"
/**
 * find p' <= p
 * mat: coefficient matrix of input polyhedron
 * p: number of variables to eliminate
 * m: number of inequalities in the input system
 */
int pPrime(mpq_t * mat, int p, int m)
{
	int prime = p;
	int flag = 1;
	while (1)
	{
		flag = 1;
		for (int i = 0; i < m; i++)
			if (mpq_sgn(mat[(prime - 1) + i * m]) != 0)
				flag = 0;
		if (flag == 0)
			prime--;
		else
			break;
	}
	return prime;
}
/**
 * decompose the input matrix to get coefficient vector of the projection cone
 * mat: input matrix
 * out1: coefficient matrix of projection cone
 * out2: coefficient matrix of other variables
 * pprime: number less that p, computed in pPrime function
 * q: number of variables that remain
 * m: number of inequalities
 */
void projectionCone(mpq_t * mat, mpq_t * out1, mpq_t * out2, int pprime, int q,
		int m)
{

	int c1 = 0;
	int c2 = 0;
	int qprime = q + pprime;

	//1- get the first chunck of matrix

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m - qprime; j++) ///<<<....
			mpq_set(out1[c1++], mat[i * m + pprime + j]);
		for (int j = m - q; j < m; j++) ////<<<.....
			mpq_set(out2[c2++], mat[(i * m) + j]);
	}

	//2- get the second chunck of matrix

	for (int j = 0; j < q; j++)
	{
		mpq_set(out2[q * (qprime - 1) + j], out2[((m - 1) * q + j)]);
		mpq_set_d(out2[((m - 1) * q + j)], 0);
	}

}

/**
 * use the structure of the cone and prepare date for finding extreme rays
 * mat: coefficient of the input cone
 * extrIn: input data prepared for finding extreme rays
 * m: number of inequalities
 * qprime: integer number related to pprime and q
 */

void makeExtrInput(mpq_t * mat/*out1*/, mpq_t * extrIn, int m, int qprime)
{
	mpq_t * extrIn1 = (mpq_t *) malloc(m * (m - qprime) * sizeof(mpq_t));
	for (int i = 0; i < m * (m - qprime); i++)
		mpq_init(extrIn1[i]);

	mpq_t * extrIn1t = (mpq_t *) malloc(m * (m - qprime) * sizeof(mpq_t));
	for (int i = 0; i < m * (m - qprime); i++)
		mpq_init(extrIn1t[i]);

	int c = 0;
	for (int i = 0; i < (qprime - 1) * (m - qprime); i++)
	{
		mpq_set(extrIn1[c], mat[c]);
		c++;
	}
	for (int i = (m - 1) * (m - qprime); i < m * (m - qprime); i++)
	{
		mpq_set(extrIn1[c], mat[i]);
		c++;
	}

	matrixTranspose(extrIn1, extrIn, qprime, m - qprime);

	c = 0;
	for (int i = 0; i < m * (m - qprime); i++)
		mpq_neg(extrIn[i], extrIn[i]);
}

/**
 * put everythings together to find and return projection of W0 cone, using block elimination
 * w0: the cone W0
 * p: number of variables to eliminate
 * q: number of variables that remain
 * m: number of inequalities
 * size: store number of inequalities in the projection
 */

mpq_t * blockElimination(inequalityULL_t w0, int p, int q, int m, int * size)
{
	mpq_t * coefMat = (mpq_t *) malloc((m * m) * sizeof(mpq_t));


	//1- decompose inpute coefficient matrix and find coefficient
	//matrix of the projection cone

	for (int i = 0; i < (m * m); i++)
		mpq_init(coefMat[i]);

	coefMatrix(w0, coefMat, m, m);

	int pprime = pPrime(coefMat, p, m);

	q = q + 1;
	int qprime = q + pprime;

	mpq_t * out1 = (mpq_t *) malloc(((m - qprime) * m) * sizeof(mpq_t));
	mpq_t * out2 = (mpq_t *) malloc(((q) * m) * sizeof(mpq_t));

	for (int i = 0; i < ((m - qprime) * m); i++)
		mpq_init(out1[i]);

	for (int i = 0; i < ((q) * m); i++)
		mpq_init(out2[i]);

	projectionCone(coefMat, out1, out2, pprime, q, m);
	//

	//2- make input ready for block elimination function

	mpq_t * extrIn = (mpq_t *) malloc((m - qprime) * qprime * sizeof(mpq_t));
	for (int i = 0; i < (m - qprime) * qprime; i++)
		mpq_init(extrIn[i]);
	for (int i = 0; i < (m - qprime) * qprime; i++)
		mpq_set_d(extrIn[i], 0);

	makeExtrInput(out1, extrIn, m, qprime);
    //
	mpq_t * extrs = NULL;

	int grow;
	int gcol;


	//3- use the function which uses CDD library to find extreme rays

	extr(extrIn, &extrs, m /*row*/, qprime /*col*/, &grow, &gcol); //<<--- CDD function

	*size = grow;

	mpq_t * mat = (mpq_t *) malloc(grow * q * sizeof(mpq_t));
	for (int i = 0; i < grow * q; i++)
		mpq_init(mat[i]);

	mpq_t * extrVec = (mpq_t *) malloc(m * sizeof(mpq_t));
	for (int i = 0; i < m; i++)
		mpq_init(extrVec[i]);

	mpq_t * rowResult = (mpq_t *) malloc(q * sizeof(mpq_t));
	for (int i = 0; i < q; i++)
		mpq_init(rowResult[i]);

	int v = 0;
	int c = 0;
	int cmat = 0;

	//3- expand extreme rays, using structure of the cone

	for (int i = 0; i < grow; i++)
	{
		for (int j = 0; j < m; j++)
			mpq_set(extrVec[j], extrs[c++]);
		matrixMatrixMult(extrVec, out2, rowResult, 1, m, q); //A: N * M , B: M * L, C: N * L

		for (int j = 0; j < q; j++)
			mpq_set(mat[cmat++], rowResult[j]);
	}

	return (mat);
}

