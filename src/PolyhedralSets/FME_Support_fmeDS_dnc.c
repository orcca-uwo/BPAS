#include "PolyhedralSets/FME_Support_fmeDS_dnc.h"

void allocFMEDS_dnc(int varNum, FMEDS_dnc_t * newds)
{
	newds->dimension = varNum;

	newds->posSubSys = (unrolledLl_dnc_t *) malloc(sizeof(unrolledLl_dnc_t));
	allocList_dnc(varNum, newds->posSubSys);

	newds->negSubSys = (unrolledLl_dnc_t *) malloc(sizeof(unrolledLl_dnc_t));
	allocList_dnc(varNum, newds->negSubSys);

	newds->zerSubSys = (unrolledLl_dnc_t *) malloc(sizeof(unrolledLl_dnc_t));
	allocList_dnc(varNum, newds->zerSubSys);

	newds->currentPos = (inequalityNode_dnc_t **) malloc(
			sizeof(inequalityNode_dnc_t *));
	newds->currentNeg = (inequalityNode_dnc_t **) malloc(
			sizeof(inequalityNode_dnc_t *));
	newds->currentZer = (inequalityNode_dnc_t **) malloc(
			sizeof(inequalityNode_dnc_t *));

	*(newds->currentPos) = newds->posSubSys->head;
	*(newds->currentNeg) = newds->negSubSys->head;
	*(newds->currentZer) = newds->zerSubSys->head;


	newds->posIndex = (int *) malloc(sizeof(int));
	newds->negIndex = (int *) malloc(sizeof(int));
	newds->zerIndex = (int *) malloc(sizeof(int));

	*(newds->posIndex) = 1;
	*(newds->negIndex) = 1;
	*(newds->zerIndex) = 1;

}

void addZeros_dnc(FMEDS_dnc_t * input, int varIndex, FMEDS_dnc_t * output, 
		matrix_t * testCone)
{

//	printf("Enter add zeros\n");
	int countZ = 0;
	input->currentZer = &(input->zerSubSys->head);
	inequality_t ** newZeroIneq = (inequality_t **) malloc(
			sizeof(inequality_t *));

	while (/*input->currentZer != NULL*/countZ <= input->zerSubSys->number
			&& input->zerSubSys->number != 0)
	{
		int flagz = getNext_dnc((input->currentZer), input->zerIndex, newZeroIneq);
		if (flagz == 0)
			break;
		countZ++;
		if (redundancyTest(testCone, (*newZeroIneq) , varIndex) == 1)

			addToFMEDS_dnc(output, (*newZeroIneq), varIndex + 1);
	}

	free(newZeroIneq);
}

void divideFMEDS_dnc(FMEDS_dnc_t * input, FMEDS_dnc_t * out1, FMEDS_dnc_t * out2)
{

//	printf("******* Enter divid\n");
//	printList(input->posSubSys);
	inequalityNode_dnc_t * halfHeadPos;
	halfHeadPos = input->posSubSys->head;
	int halfNumPos1;
	int halfNumPos2;
	if (input->posSubSys->number % 2 == 0)
	{
		halfNumPos1 = input->posSubSys->number / 2;
		halfNumPos2 = input->posSubSys->number - halfNumPos1;
	}
	else
	{
		halfNumPos1 = (input->posSubSys->number + 1) / 2;
		halfNumPos2 = input->posSubSys->number - halfNumPos1;
	}

	for (int i = 0; i < halfNumPos1; i++)
		halfHeadPos = halfHeadPos->next;

	out1->posSubSys->head->next = input->posSubSys->head->next;
	out1->posSubSys->number = halfNumPos1;

	out2->posSubSys->number = halfNumPos2;
	out2->posSubSys->head->next = halfHeadPos->next; /////<<<<=====


//////////////////////////////////////

	inequalityNode_dnc_t * halfHeadNeg;
	halfHeadNeg = input->negSubSys->head;
	int halfNumNeg1;
	int halfNumNeg2;
	if (input->negSubSys->number % 2 == 0)
	{
		halfNumNeg1 = input->negSubSys->number / 2;
		halfNumNeg2 = input->negSubSys->number - halfNumNeg1;
	}
	else
	{
		halfNumNeg1 = (input->negSubSys->number + 1) / 2;
		halfNumNeg2 = input->negSubSys->number - halfNumNeg1;
	}

	for (int i = 0; i < halfNumNeg1; i++)
		halfHeadNeg = halfHeadNeg->next;

	out1->negSubSys->head->next = input->negSubSys->head->next;
	out1->negSubSys->number = halfNumNeg1;

	out2->negSubSys->head->next = halfHeadNeg->next;
	out2->negSubSys->number = halfNumNeg2;


}

void mergeFMEDS_dnc(FMEDS_dnc_t * out1, FMEDS_dnc_t * out2, FMEDS_dnc_t * out3, FMEDS_dnc_t * out4,
		FMEDS_dnc_t * output)
{
	//printf("Enter merge \n");
	output->posSubSys->number = out1->posSubSys->number
			+ out2->posSubSys->number + out3->posSubSys->number
			+ out4->posSubSys->number;

	output->negSubSys->number = out1->negSubSys->number
			+ out2->negSubSys->number + out3->negSubSys->number
			+ out4->negSubSys->number;

	inequalityNode_dnc_t * endPos1 = out1->posSubSys->head;
	if (out1->posSubSys->number != 0)
		for (int i = 0; i < out1->posSubSys->number; i++)
			endPos1 = endPos1->next;

	inequalityNode_dnc_t * endPos2 = out2->posSubSys->head;
	if (out2->posSubSys->number != 0)
		for (int i = 0; i < out2->posSubSys->number; i++)
			endPos2 = endPos2->next;
//
	inequalityNode_dnc_t * endPos3 = out3->posSubSys->head;
	if (out3->posSubSys->number != 0)
		for (int i = 0; i < out3->posSubSys->number; i++)
			endPos3 = endPos3->next;

	inequalityNode_dnc_t * endPos4 = out4->posSubSys->head;
	if (out4->posSubSys->number != 0)
		for (int i = 0; i < out4->posSubSys->number; i++)
			endPos4 = endPos4->next;

	inequalityNode_dnc_t * posAttachPoint = output->posSubSys->head;
//
	if (out1->posSubSys->number != 0)
	{
		posAttachPoint->next = out1->posSubSys->head->next;
		posAttachPoint = endPos1;
	}
	if (out2->posSubSys->number != 0)
	{
//		printf("Enter 2 if\n");
		posAttachPoint->next = out2->posSubSys->head->next;
		posAttachPoint = endPos2;
	}

	if (out3->posSubSys->number != 0)
	{
//		printf("Enter 3 if\n");
		posAttachPoint->next = out3->posSubSys->head->next;
		posAttachPoint = endPos3;
	}

	if (out4->posSubSys->number != 0)
	{
		posAttachPoint->next = out4->posSubSys->head->next;
		posAttachPoint = endPos4;
	}

	posAttachPoint->next = NULL;

/////////////////////////////////////////////////////negatives

	output->negSubSys->number = out1->negSubSys->number
			+ out2->negSubSys->number + out3->negSubSys->number
			+ out4->negSubSys->number;

	output->negSubSys->number = out1->negSubSys->number
			+ out2->negSubSys->number + out3->negSubSys->number
			+ out4->negSubSys->number;

	inequalityNode_dnc_t * endneg1 = out1->negSubSys->head;
	if (out1->negSubSys->number != 0)
		for (int i = 0; i < out1->negSubSys->number; i++)
			endneg1 = endneg1->next;

	inequalityNode_dnc_t * endneg2 = out2->negSubSys->head;
	if (out2->negSubSys->number != 0)
		for (int i = 0; i < out2->negSubSys->number; i++)
			endneg2 = endneg2->next;
	//
	inequalityNode_dnc_t * endneg3 = out3->negSubSys->head;
	if (out3->negSubSys->number != 0)
		for (int i = 0; i < out3->negSubSys->number; i++)
			endneg3 = endneg3->next;

	inequalityNode_dnc_t * endneg4 = out4->negSubSys->head;
	if (out4->negSubSys->number != 0)
		for (int i = 0; i < out4->negSubSys->number; i++)
			endneg4 = endneg4->next;

	inequalityNode_dnc_t * negAttachPoint = output->negSubSys->head;
	//
	if (out1->negSubSys->number != 0)
	{
		negAttachPoint->next = out1->negSubSys->head->next;
		negAttachPoint = endneg1;
	}
	if (out2->negSubSys->number != 0)
	{
		//printf("Enter 2 if\n");
		negAttachPoint->next = out2->negSubSys->head->next;
		negAttachPoint = endneg2;
	}

	if (out3->negSubSys->number != 0)
	{
	//	printf("Enter 3 if\n");
		negAttachPoint->next = out3->negSubSys->head->next;
		negAttachPoint = endneg3;
	}

	if (out4->negSubSys->number != 0)
	{
		negAttachPoint->next = out4->negSubSys->head->next;
		negAttachPoint = endneg4;
	}

	negAttachPoint->next = NULL;


	//////////////////////////////////////////////***zeros
	
	output->zerSubSys->number = out1->zerSubSys->number
			+ out2->zerSubSys->number + out3->zerSubSys->number
			+ out4->zerSubSys->number;

	inequalityNode_dnc_t * endzer1 = out1->zerSubSys->head;
	if (out1->zerSubSys->number != 0)
		for (int i = 0; i < out1->zerSubSys->number; i++)
			endzer1 = endzer1->next;

	inequalityNode_dnc_t * endzer2 = out2->zerSubSys->head;
	if (out2->zerSubSys->number != 0)
		for (int i = 0; i < out2->zerSubSys->number; i++)
			endzer2 = endzer2->next;
	//
	inequalityNode_dnc_t * endzer3 = out3->zerSubSys->head;
	if (out3->zerSubSys->number != 0)
		for (int i = 0; i < out3->zerSubSys->number; i++)
			endzer3 = endzer3->next;

	inequalityNode_dnc_t * endzer4 = out4->zerSubSys->head;
	if (out4->zerSubSys->number != 0)
		for (int i = 0; i < out4->zerSubSys->number; i++)
			endzer4 = endzer4->next;

	inequalityNode_dnc_t * zerAttachPoint = output->zerSubSys->head;
	//
	if (out1->zerSubSys->number != 0)
	{
		zerAttachPoint->next = out1->zerSubSys->head->next;
		zerAttachPoint = endzer1;
	}
	if (out2->zerSubSys->number != 0)
	{
	//	printf("Enter 2 if\n");
		zerAttachPoint->next = out2->zerSubSys->head->next;
		zerAttachPoint = endzer2;
	}

	if (out3->zerSubSys->number != 0)
	{
		//printf("Enter 3 if\n");
		zerAttachPoint->next = out3->zerSubSys->head->next;
		zerAttachPoint = endzer3;
	}

	if (out4->zerSubSys->number != 0)
	{
		zerAttachPoint->next = out4->zerSubSys->head->next;
		zerAttachPoint = endzer4;
	}

	zerAttachPoint->next = NULL;

	///////////////////////////////////////////mange pointeres
	
	*(output->currentPos) = (output->posSubSys->head);
        *(output->currentNeg) = (output->negSubSys->head);
	*(output->currentZer) = (output->zerSubSys->head);

	

	for(int i = 0 ; i < output->posSubSys->number; i++)
		*(output->currentPos) = (*(output->currentPos))->next;

	for(int i = 0 ; i < output->negSubSys->number; i++)
		*(output->currentNeg) = (*(output->currentNeg))->next;

	
	for(int i = 0 ; i < output->zerSubSys->number; i++)
		*(output->currentZer) = (*(output->currentZer))->next;

}


void FMEDSRewind_dnc(FMEDS_dnc_t * a)
{
	*(a->currentPos) = (a->posSubSys->head);
	*(a->currentNeg) = a->negSubSys->head;
	*(a->currentZer) = a->zerSubSys->head;

	*(a->posIndex) = 1;
	*(a->negIndex) = 1;
	*(a->zerIndex) = 1;
}

void FMEDSListRewind_dnc(FMEDS_dnc_t * a, char n)
{
	if (n == 'n')
	{
		*(a->currentNeg) = a->negSubSys->head;
		*(a->negIndex) = UNROLLED_LL_SIZE;
	}
	if (n == 'p')
	{
		*(a->currentPos) = a->posSubSys->head;
		*(a->posIndex) = UNROLLED_LL_SIZE;
	}
	if (n == 'z')
	{
		*(a->currentZer) = a->zerSubSys->head;
		*(a->zerIndex) = UNROLLED_LL_SIZE;
	}

}

void addToFMEDS_dnc(FMEDS_dnc_t * a, inequality_t * newIneq, int varIndex)
{

	if (ineqSgn(*newIneq, varIndex) == 1)
		specificAdd_dnc(a->posSubSys, a->currentPos, a->posIndex, newIneq);
	else if (ineqSgn(*newIneq, varIndex) == -1)
		specificAdd_dnc(a->negSubSys, a->currentNeg, a->negIndex, newIneq);
	else
		specificAdd_dnc(a->zerSubSys, a->currentZer, a->zerIndex, newIneq);

}

int getFMEDSSize_dnc(FMEDS_dnc_t * a)
{
	int pSize = getListSize_dnc((a->posSubSys));
	int nSize = getListSize_dnc((a->negSubSys));
	int zSize = getListSize_dnc((a->zerSubSys));

	return (pSize + nSize + zSize);
	return (nSize);
}



void freeFMEDS_dnc(FMEDS_dnc_t * inputFMEDS)
{
//	printf("from free head address is %p\n", inputFMEDS->negSubSys->head);
	freeList_dnc(inputFMEDS->negSubSys);
	freeList_dnc(inputFMEDS->posSubSys);
	freeList_dnc(inputFMEDS->zerSubSys);

//	if(inputFMEDS->currentNeg != NULL)
	free((inputFMEDS->currentNeg));
//	if(inputFMEDS->currentPos != NULL)
	free(inputFMEDS->currentPos);
//	if(inputFMEDS->currentZer != NULL)
//	free(inputFMEDS->currentZer); <<<---

	free(inputFMEDS->negIndex);
	free(inputFMEDS->posIndex);
	free(inputFMEDS->zerIndex);
}
