/**
 * FME_Support_extreme.c
 * ****************************
 * implementation of interface to CDD library
 */

#include "PolyhedralSets/FME_Support_extreme.h"


void findExtremeRays(const matrix_t in , matrix_t * out)
{
	printf("Enter finding extrem rays\n");
	dd_MatrixPtr A , G;
	dd_PolyhedraPtr poly;
	dd_ErrorType err;

	int col = in.colNum;
	int row = in.rowNum + col;

	dd_set_global_constants();

	A = dd_CreateMatrix(row , col);

	for(int i = 0 ; i < col ; i++)
		dd_set_d(A->matrix[i][i] , 1);

	int counter = 0;
	for(int i = col ; i < row ; i++)
		for(int j = 0 ; j < col ; j++)
				dd_set(A->matrix[i][j], in.entries[counter++]);

	A->representation = dd_Inequality;
//	dd_WriteMatrix(stdout,A);

	poly = dd_DDMatrix2Poly(A, &err); //<<<---- Main function

	G = dd_CopyGenerators(poly);


//	dd_WriteMatrix(stdout,G);

	int grow = G->rowsize;
	int gcol = G->colsize;

//	printf("grow is %d and gcol is %d\n",grow,gcol);
	allocMatrix(grow, gcol , out);


	counter = 0;
	for(int i = 0 ; i <  grow ; i++)
		for(int j = 0 ; j < gcol ; j++)
			mpq_set(out->entries[counter++] , G->matrix[i][j]);

	dd_FreeMatrix(A);
	dd_FreeMatrix(G);
	dd_FreePolyhedra(poly);

}


