/*
 * FME_Support_redundancyTest.c
 * ****************************
 * implementation of funcitons to find the
 * redundancy test cone (W^0) of an input polyhedron.
 */

#include "PolyhedralSets/FME_Support_redundancyTest.h"


void getCoeffMatrix(inequality_t * l , matrix_t * mat , int ineqNum)
{
    int varNum = l[0].dimension;
    allocMatrix(ineqNum , varNum , mat);

    int counter = 0;
    for(int i = 0 ; i < ineqNum ; i++)
        for(int j = 0 ; j < varNum ; j++)
            mpq_set_z(mat->entries[counter++] , l[i].coeff[j]);
}

void getConstVector(inequality_t * l , matrix_t * vec , int ineqNum)
{
	allocMatrix(ineqNum , 1 , vec);
	for(int i = 0 ; i < ineqNum ; i++)
		mpq_set_z(vec->entries[i] , l[i].constant);
}
//assumption: matrix is full column rank
void extendMatrix(matrix_t source , matrix_t * dest)
{
    int r = source.rowNum;
    int c = source.colNum;
    int extNum = r - c;

    allocMatrix(r , r , dest);

    int counter = 0;
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < r; j++)
		{
			if (j < c)
				mpq_set(dest->entries[i*r+j], source.entries[i*c+j]);
			else
			{
				if (i < c)
					mpq_set_d(dest->entries[i*r + j], 0);
				else
				{
					if (j-c == i-c)
						mpq_set_d(dest->entries[i*r+j], 1);
					else
						mpq_set_d(dest->entries[i*r +j], 0);
				}
			}
		}
	}
}
//A_0^{-1} , c are the inputs
void remVarMat(matrix_t inputMat , matrix_t inputVec , matrix_t * outMat , int varNum)
{
	int m = inputMat.rowNum;
	int n = varNum;

//	allocMatrix(n+1 , n+1 , outMat);

	matrix_t * temp = (matrix_t *) malloc(sizeof(matrix_t));
//	allocMatrix(m, 1, temp);
	matrixMatrixMult(&inputMat, &inputVec, temp);

	matrix_t * matTrans = (matrix_t *) malloc(sizeof(matrix_t));
	allocMatrix(n+1 , n+1 , matTrans);

	int counter = 0;
	int index = 0;
	int tindex = 0;

	for(int i = 0 ; i < n ; i++)
	{
		for(int j = 0 ; j < n ; j++)
			mpq_set(matTrans->entries[counter++] , inputMat.entries[index++]);

		mpq_neg(temp->entries[tindex],temp->entries[tindex]);
		mpq_set(matTrans->entries[counter++] , temp->entries[tindex++]);
		index = index + m - n;
	}

	for(int i = 0 ; i < n ; i++)
		mpq_set_d(matTrans->entries[counter++] , 0);
	mpq_set_d(matTrans->entries[counter++] , 1);

	matrixTranspose(matTrans , outMat);

	freeMatrix(temp);
	freeMatrix(matTrans);
	free(temp);
	free(matTrans);
}
//A_0^{-1} , c are the inputs
void elimVarMat(matrix_t inputMat , matrix_t inputVec , matrix_t * outMat , int varNum)
{
	int m = inputMat.rowNum;
	int n = varNum;

	allocMatrix(m-n , n+1 , outMat);

	matrix_t * temp = (matrix_t *) malloc(sizeof(matrix_t));
	matrixMatrixMult(& inputMat , & inputVec , temp);
//
	int counter = 0;
	int index = (m*n);
	int tIndex = n;
	mpq_t tmpRn;
	mpq_init(tmpRn);

	for(int i = 0 ; i < m-n ; i++)
	{
		for(int j = 0 ; j < n ; j++)
		{
			mpq_neg(tmpRn , inputMat.entries[index++]);
			mpq_set(outMat->entries[counter++] , tmpRn);
		}
		mpq_set(outMat->entries[counter++] , temp->entries[tIndex++]);
		index = index + m - n;
	}
	printf("enter free temp\n");
	freeMatrix(temp);
	free(temp);
	mpq_clear(tmpRn);
}

void blocElimination(matrix_t inputMat , matrix_t * extrMat , matrix_t * outMat)
{
	int extrNum = extrMat->rowNum;
	int varNum = extrMat->colNum;

	matrixMatrixMult(extrMat , &inputMat , outMat); //this is the initial test cone
}

void initialTestCone(inequality_t * inputSys , int ineqNum , matrix_t * initTestCone)
{
	int varNum = inputSys[0].dimension;

	matrix_t * inputCoeffMat = (matrix_t *) malloc(sizeof(matrix_t));
	getCoeffMatrix(inputSys, inputCoeffMat, ineqNum);

	matrix_t * inputConsVec = (matrix_t *) malloc(sizeof(matrix_t));
	getConstVector(inputSys, inputConsVec, ineqNum);

	matrix_t * extendedMat = (matrix_t *) malloc(sizeof(matrix_t));
	extendMatrix(*inputCoeffMat, extendedMat);

	matrix_t * inv = (matrix_t *) malloc(sizeof(matrix_t));
	matrixInverse(extendedMat, inv);

	matrix_t * elim = (matrix_t *) malloc(sizeof(matrix_t));
	elimVarMat(*inv, *inputConsVec, elim, varNum);


	matrix_t * rem = (matrix_t *) malloc(sizeof(matrix_t)); //<<-- problem
	remVarMat(*inv, *inputConsVec, rem, varNum);

	matrix_t * extremeRays = (matrix_t *) malloc(sizeof(matrix_t));
	findExtremeRays(*elim, extremeRays);

	blocElimination(*rem, extremeRays, initTestCone);

	//cleaning
	freeMatrix(inputCoeffMat);
	freeMatrix(inputConsVec);
	freeMatrix(extendedMat);
	freeMatrix(inv);
	freeMatrix(elim);
	freeMatrix(rem);
	freeMatrix(extremeRays);
	free(inputCoeffMat);
	free(inputConsVec);
	free(extendedMat);
	free(inv);
	free(elim);
	free(rem);
	free(extremeRays);
}

void testCone(matrix_t * initialCone , int step , matrix_t * output)
{
	int rowNum = initialCone->rowNum;
	int colNum = initialCone->colNum;

	int count = 0;
	for(int i = 0 ; i < rowNum ; i++)
		for(int j = step ; j < colNum  ; j++)
			mpq_set(output->entries[count++] , initialCone->entries[i*colNum + j]);

}


int extremeRayTest(matrix_t * coeffMat , matrix_t * vec)
{
	matrix_t * mulOutput = (matrix_t *) malloc(sizeof(matrix_t));

//	allocMatrix(coeffMat->rowNum , 1 , mulOutput);
	matrixMatrixMult(coeffMat , vec , mulOutput);

	int n = 0;
	for(int i = 0 ; i < mulOutput->rowNum ; i++)
		if(mpq_sgn(mulOutput->entries[i]) == 0)
			n++;
	int * idx = (int *) malloc(sizeof(int) * n);
	int c = 0;
	for(int i = 0 ; i < mulOutput->rowNum ; i++)
		if(mpq_sgn(mulOutput->entries[i]) == 0)
			idx[c++] = i;

	matrix_t * subMatrix = (matrix_t *) malloc(sizeof(matrix_t));
	subMat(coeffMat , idx , n , subMatrix);

	int r = matrixRank(subMatrix);

	//cleaning
	freeMatrix(mulOutput);
	freeMatrix(subMatrix);
	free(mulOutput);
	free(subMatrix);
	free(idx);

	if(r == coeffMat->colNum-1)
		return(1);
	else
		return(0);
}

int redundancyTest(matrix_t * testCone , inequality_t * toBeTestedIneq , int idx)
{

	matrix_t * toBeTestedVec = (matrix_t *) malloc(sizeof(matrix_t));
	allocMatrix(testCone->colNum , 1 , toBeTestedVec);

	for(int i = 0 ; i < testCone->colNum-1 ; i++)
		mpq_set_z(toBeTestedVec->entries[i], toBeTestedIneq->coeff[idx+i+1]);

	mpq_set_z(toBeTestedVec->entries[testCone->colNum-1] , toBeTestedIneq->constant);
	int r = extremeRayTest(testCone , toBeTestedVec);

	//cleaning
	freeMatrix(toBeTestedVec);
	free(toBeTestedVec);

	if(r == 1)
		return(1);
	else
		return(0);
}


