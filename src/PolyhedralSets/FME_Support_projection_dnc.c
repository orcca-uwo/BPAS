/**
 * FME_Support_projection.c
 * *************************
 * implementation of minimal representation projection
 */
#include "PolyhedralSets/FME_Support_projection_dnc.h"

void project_dnc(char * fileName, int varNum, int ineqNum, int thr)
{

	inequality_t * inputData = (inequality_t *) malloc(
			ineqNum * sizeof(inequality_t));
	getInputArray(fileName, inputData, varNum, ineqNum);

	matrix_t * initTestCone = (matrix_t *) malloc(sizeof(matrix_t));
	initialTestCone(inputData, ineqNum, initTestCone);

	FMEDS_dnc_t * inputFMEDS = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));
	FMEDS_dnc_t * outputFMEDS = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));

	allocFMEDS_dnc(varNum, inputFMEDS);

	allocFMEDS_dnc(varNum, outputFMEDS);


	for (int i = 0; i < ineqNum; i++)
		addToFMEDS_dnc(inputFMEDS, &inputData[i], 0);
	


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

	        oneStepFME_dnc(inputFMEDS, i, outputFMEDS, theTestCone, thr);
//
		printf("Step %d finishes\n", i);

		int ppSize = getFMEDSSize_dnc(outputFMEDS);
		printf("***** %d\n", ppSize);
//
		freeMatrix(theTestCone);

		freeFMEDS_dnc(inputFMEDS);

		inputFMEDS = outputFMEDS;

		outputFMEDS = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));
		allocFMEDS_dnc(varNum , outputFMEDS);

	}

	freeMatrix(initTestCone);
	free(initTestCone);
	free(theTestCone);
	freeFMEDS_dnc(inputFMEDS);
	free(inputFMEDS);
	freeFMEDS_dnc(outputFMEDS);
	free(outputFMEDS);
}
