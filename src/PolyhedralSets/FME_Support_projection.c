/**
 * FME_Support_projection.c
 * *************************
 * implementation of minimal representation projection
 */
#include "PolyhedralSets/FME_Support_projection.h"

void project(char * fileName, int varNum, int ineqNum)
{

	inequality_t * inputData = (inequality_t *) malloc(
			ineqNum * sizeof(inequality_t));
	getInputArray(fileName, inputData, varNum, ineqNum);

	matrix_t * initTestCone = (matrix_t *) malloc(sizeof(matrix_t));
	initialTestCone(inputData, ineqNum, initTestCone);

	FMEDS_t * inputFMEDS = (FMEDS_t *) malloc(sizeof(FMEDS_t));
	FMEDS_t * outputFMEDS = (FMEDS_t *) malloc(sizeof(FMEDS_t));

	allocFMEDS(varNum, inputFMEDS);

	allocFMEDS(varNum, outputFMEDS);


	for (int i = 0; i < ineqNum; i++)
	{
		addToFMEDS(inputFMEDS, &inputData[i], 0);
	}


	for (int i = 0; i < ineqNum; i++)
		freeInequality(&(inputData[i]));
	free(inputData);

	matrix_t * theTestCone = (matrix_t *) malloc(sizeof(matrix_t));

	printf("********************FME data struct allocation Done\n");

	for (int i = 0; i < varNum - 1; i++)
	{
		printf("Step %d Statrs\n", i);

		allocMatrix(initTestCone->rowNum, initTestCone->colNum - i - 1,
				theTestCone);

		testCone(initTestCone, i + 1, theTestCone);

		oneStepFME(inputFMEDS, i, outputFMEDS, theTestCone);
//
		printf("Step %d finishes\n", i);

		int ppSize = getFMEDSSize(outputFMEDS);
		printf("***** %d\n", ppSize);
//
		freeMatrix(theTestCone);

		freeFMEDS(inputFMEDS);

		inputFMEDS = outputFMEDS;

		outputFMEDS = (FMEDS_t *) malloc(sizeof(FMEDS_t));
		allocFMEDS(varNum , outputFMEDS);

	}

	freeMatrix(initTestCone);
	free(initTestCone);
	free(theTestCone);
	freeFMEDS(inputFMEDS);
	free(inputFMEDS);
	freeFMEDS(outputFMEDS);
	free(outputFMEDS);
}
