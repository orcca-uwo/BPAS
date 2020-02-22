/**
 * FME_Support_fme.c
 * *******************
 * implementation of Fourier-Motzkin elimination
 * using Balas method to eliminate 
 * redundant inequalities
 */


#include "PolyhedralSets/FME_Support_fme.h"

void oneStepFME(FMEDS_t * input, int varIndex, FMEDS_t * output,
		matrix_t * testCone)
{
	printf("Enter fme \n");
	int dim = input->dimension;
	inequality_t ** newPositiveIneq = (inequality_t **) malloc(
			sizeof(inequality_t *));
	inequality_t ** newNegativeIneq = (inequality_t **) malloc(
			sizeof(inequality_t *));
	inequality_t ** newZeroIneq = (inequality_t **) malloc(
			sizeof(inequality_t *));


	FMEDSRewind(input);

	
////////////////////////////////////////////////////////////////////////////////
//
	inequality_t * newIneq = (inequality_t *) malloc(sizeof(inequality_t));
	allocInequality(dim, newIneq);
//
	while (input->currentPos != NULL && input->posSubSys->number != 0)
	{
		int flagp = getNext((input->currentPos), input->posIndex,
				newPositiveIneq);
		if (flagp == 0)
			break;

		while (input->currentNeg != NULL && input->negSubSys->number != 0)
		{
//
			int flagn = getNext((input->currentNeg), input->negIndex,
					newNegativeIneq);
			if (flagn == 0)
				break;
//
			combineTwoIneq((**newPositiveIneq), (**newNegativeIneq), newIneq,
					varIndex);
//
////			printf("********** new ineqalite is : \n");
////			printInequality(newIneq, 'x');
//
			if (redundancyTest(testCone, newIneq, varIndex) == 1)
			{
				printInequality(newIneq, 'x');
				addToFMEDS(output, newIneq, varIndex + 1);
			}
//
		}
		FMEDSListRewind(input, 'n');
	}
//
	while (input->currentZer != NULL && input->zerSubSys->number != 0)
	{
		int flagz = getNext((input->currentZer), input->zerIndex, newZeroIneq);
		if (flagz == 0)
			break;
		if (redundancyTest(testCone, (*newZeroIneq) , varIndex) == 1)
		{
//			printInequality((*newZeroIneq), 'x');
			addToFMEDS(output, (*newZeroIneq), varIndex + 1);
		}
	}

	freeInequality(newIneq);
	free(newIneq);
	free(newPositiveIneq);
	free(newNegativeIneq);
	free(newZeroIneq);
}





