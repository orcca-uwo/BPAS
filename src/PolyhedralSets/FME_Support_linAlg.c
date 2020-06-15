/**
 * FME_Support_linAlg.c
 * **********************
 * interface to Flint
 */

#include "PolyhedralSets/FME_Support_linAlg.h"

void allocMatrix(int r , int c , matrix_t * m)
{
	m->rowNum = r;
	m->colNum = c;

	m->entries = (mpq_t *) malloc(r * c * sizeof(mpq_t));
	for(int i = 0 ; i < r * c ; i++)
		mpq_init(m->entries[i]);
}

void setMatrix(matrix_t * m ,  mpq_t * data)
{
	int num = m->rowNum * m->colNum;

	for(int i = 0 ; i < num ; i++)
		mpq_set(m->entries[i] , data[i]);

}

void printMatrix(matrix_t * m)
{
	printf("\n Matrix is:\n");
	for(int i = 0 ; i < m->rowNum ; i++)
	{
		for(int j = 0 ; j < m->colNum ; j++)
			gmp_printf("%Qd ", m->entries[i * m->colNum + j]);
		printf("\n");
	}
}

void matrixToFlint( matrix_t source , fmpq_mat_t * dest)
{
	int num = 0;
	int rowNumber = source.rowNum;
	int colNumber = source.colNum;

	for(int i = 0 ; i < rowNumber ; i++ )
		for(int j = 0 ; j < colNumber ; j++)
			fmpq_set_mpq(fmpq_mat_entry(*dest, i, j), source.entries[num++]);
}

void flintToMatrix( fmpq_mat_t source , matrix_t * dest)
{
	int s = 0;
	int rowNumber = fmpq_mat_nrows(source);
	int colNumber = fmpq_mat_ncols(source);

	allocMatrix(rowNumber , colNumber , dest);

	for(int i = 0 ; i < rowNumber ; i++)
		for(int j = 0 ; j < colNumber ; j++)
			fmpq_get_mpq(dest->entries[s++], fmpq_mat_entry(source, i, j));


}

void matrixMatrixMult(matrix_t * mult1 , matrix_t * mult2 , matrix_t * out)
{
	if(mult1->colNum != mult2->rowNum)
	{
		printf("Unable to multiply matrices!\n");
		return;
	}

	fmpq_mat_t flintMult1;
	fmpq_mat_t flintMult2;
	fmpq_mat_t flintResult;

	fmpq_mat_init(flintMult1 , mult1->rowNum , mult1->colNum);
	fmpq_mat_init(flintMult2 , mult2->rowNum , mult2->colNum);
	fmpq_mat_init(flintResult , mult1->rowNum , mult2->colNum);

	matrixToFlint(*mult1 , &flintMult1);
	matrixToFlint(*mult2 , &flintMult2);

	fmpq_mat_mul(flintResult , flintMult1 , flintMult2);

	flintToMatrix(flintResult , out);

	fmpq_mat_clear(flintMult1);
	fmpq_mat_clear(flintMult2);
	fmpq_mat_clear(flintResult);

}

int matrixRank(matrix_t * input)
{
	fmpq_mat_t flintInput;
	fmpq_mat_t flintTemp;

	fmpq_mat_init(flintInput, input->rowNum, input->colNum);
	fmpq_mat_init(flintTemp , input->rowNum , input->colNum);

	matrixToFlint(*input , &flintInput);

	int rank = fmpq_mat_rref(flintTemp , flintInput);

	fmpq_mat_clear(flintInput);
	fmpq_mat_clear(flintTemp);

	return(rank);
}

void matrixTranspose(matrix_t * source , matrix_t * dest)
{	
	fmpq_mat_t flintSource;
	fmpq_mat_t flintDest;

	fmpq_mat_init(flintSource , source->rowNum , source->colNum);
	fmpq_mat_init(flintDest , source->colNum , source->rowNum);

	matrixToFlint(*source , &flintSource);

	fmpq_mat_transpose(flintDest , flintSource);

	flintToMatrix(flintDest , dest);


	fmpq_mat_clear(flintSource);
	fmpq_mat_clear(flintDest);
}

void matrixInverse(matrix_t * mat , matrix_t * inverse)
{
	fmpq_mat_t flintMat;
	fmpq_mat_t flintInverse;

	fmpq_mat_init(flintMat , mat->rowNum , mat->colNum);
	fmpq_mat_init(flintInverse , mat->rowNum , mat->colNum);

	matrixToFlint(*mat , &flintMat);
	fmpq_mat_inv(flintInverse , flintMat);
	flintToMatrix(flintInverse , inverse); //<<<<<----

	fmpq_mat_clear(flintMat);
	fmpq_mat_clear(flintInverse);
}

void subMat( matrix_t * mat , int * index , int size , matrix_t * subMat)
{
	int colNum = mat->colNum;
	allocMatrix(size , colNum , subMat);

	int count = 0;
	for(int i = 0 ; i < size ; i++)
	{
		int row = index[i];
		int p = row * colNum;
		for(int j = 0 ; j < colNum ; j++)
			mpq_set(subMat->entries[count++] , mat->entries[p+j]);
	}

}

void freeMatrix(matrix_t * m)
{
	for(int i = 0 ; i < m->rowNum * m->colNum ; i++)
		mpq_clear(m->entries[i]);
	free(m->entries);
}

