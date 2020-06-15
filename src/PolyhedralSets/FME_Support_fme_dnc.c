#include "PolyhedralSets/FME_Support_fme_dnc.h"

void oneStepFME_dncSupp(FMEDS_dnc_t * input, int varIndex, FMEDS_dnc_t * output,
		matrix_t * testCone)
{
	int dim = input->dimension;
	inequality_t ** newPositiveIneq = (inequality_t **) malloc(
			sizeof(inequality_t *));
	inequality_t ** newNegativeIneq = (inequality_t **) malloc(
			sizeof(inequality_t *));

	FMEDSRewind_dnc(input);

	inequality_t * newIneq = (inequality_t *) malloc(sizeof(inequality_t));
	allocInequality(dim, newIneq);
//
	int countP = 0;
	int countN = 0;

	while (/*input->currentPos != NULL*/countP < input->posSubSys->number
			&& input->posSubSys->number != 0)
	{
		int flagp = getNext_dnc((input->currentPos), input->posIndex,
				newPositiveIneq);
		if (flagp == 0)
			break;

		countP++;
		countN = 0;
		while (/*input->currentNeg != NULL*/countN < input->negSubSys->number
				&& input->negSubSys->number != 0)
		{
//
			int flagn = getNext_dnc((input->currentNeg), input->negIndex,
					newNegativeIneq);
			if (flagn == 0)
				break;
//
			countN++;
			combineTwoIneq((**newPositiveIneq), (**newNegativeIneq), newIneq,
					varIndex);
//
			if (redundancyTest(testCone, newIneq, varIndex) == 1)
			{
		//		printf("Added\n");
			//	printInequality(newIneq, 'x');
				addToFMEDS_dnc(output, newIneq, varIndex + 1);
			}
//
		}
		FMEDSListRewind_dnc(input, 'n');
	}

	freeInequality(newIneq);
	free(newIneq);
	free(newPositiveIneq);
	free(newNegativeIneq);
}

void oneStepFME_dnc(FMEDS_dnc_t * input, int varIndex, FMEDS_dnc_t * output,
		matrix_t * testCone, int thr)
{

	int dim = input->dimension;
	
	if (input->posSubSys->number <= thr || input->negSubSys->number <= thr)
	{
		oneStepFME_dncSupp(input, varIndex, output, testCone);
		return;
	}
	else
	{

		FMEDS_dnc_t * input1 = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));
		FMEDS_dnc_t * input2 = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));

		FMEDS_dnc_t * output1 = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));
		FMEDS_dnc_t * output2 = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));

		FMEDS_dnc_t * output3 = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));
		FMEDS_dnc_t * output4 = (FMEDS_dnc_t *) malloc(sizeof(FMEDS_dnc_t));

		allocFMEDS_dnc(dim, input1);
		allocFMEDS_dnc(dim, input2);

		allocFMEDS_dnc(dim, output1);
		allocFMEDS_dnc(dim, output2);
		allocFMEDS_dnc(dim, output3);
		allocFMEDS_dnc(dim, output4);

		divideFMEDS_dnc(input, input1, input2);


	 	//__cilkrts_set_param("nworkers","4");

//		cilk_spawn
		oneStepFME_dnc(input1, varIndex, output1, testCone, thr);

//		cilk_spawn
		oneStepFME_dnc(input2, varIndex, output2, testCone, thr);
//		cilk_sync;

		inequalityNode_dnc_t * keep = input1->negSubSys->head->next;
		int keepNum = input1->negSubSys->number;
		input1->negSubSys->head->next = input2->negSubSys->head->next;
		input1->negSubSys->number = input2->negSubSys->number;
		input2->negSubSys->head->next = keep;
		input2->negSubSys->number = keepNum;

//		cilk_spawn
		oneStepFME_dnc(input1, varIndex, output3, testCone, thr);
		


//         	cilk_spawn
		oneStepFME_dnc(input2, varIndex, output4, testCone, thr);
//		cilk_sync;

		mergeFMEDS_dnc(output1, output2, output3, output4, output);
	

		addZeros_dnc(input, varIndex, output, testCone);	
//		return;

	}

}

