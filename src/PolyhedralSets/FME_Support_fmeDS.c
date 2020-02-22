/**
 * FME_Support_fmeDS.c
 * *******************
 * implementation of Fourier-Motzkin elimination data structure
 * and functions
 */


#include "PolyhedralSets/FME_Support_fmeDS.h"
void allocFMEDS(int varNum, FMEDS_t * newds)
{
	newds->dimension = varNum;

	newds->posSubSys = (unrolledLl_t *) malloc(sizeof(unrolledLl_t));
	allocList(varNum, newds->posSubSys);

	newds->negSubSys = (unrolledLl_t *) malloc(sizeof(unrolledLl_t));
	allocList(varNum, newds->negSubSys);

	newds->zerSubSys = (unrolledLl_t *) malloc(sizeof(unrolledLl_t));
	allocList(varNum, newds->zerSubSys);

	newds->currentPos = (inequalityNode_t **) malloc(
			sizeof(inequalityNode_t *));
	newds->currentNeg = (inequalityNode_t **) malloc(
			sizeof(inequalityNode_t *));
	newds->currentZer = (inequalityNode_t **) malloc(
			sizeof(inequalityNode_t *));

	*(newds->currentPos) = newds->posSubSys->head;
	*(newds->currentNeg) = newds->negSubSys->head;
	*(newds->currentZer) = newds->zerSubSys->head;

//	printf("from alloc head is %p and current is %p\n", newds->zerSubSys->head,
//			*(newds->currentZer));

	newds->posIndex = (int *) malloc(sizeof(int));
	newds->negIndex = (int *) malloc(sizeof(int));
	newds->zerIndex = (int *) malloc(sizeof(int));

	*(newds->posIndex) = UNROLLED_LL_SIZE;
	*(newds->negIndex) = UNROLLED_LL_SIZE;
	*(newds->zerIndex) = UNROLLED_LL_SIZE;

}

void FMEDSRewind(FMEDS_t * a)
{
	*(a->currentPos) = (a->posSubSys->head);
	*(a->currentNeg) = a->negSubSys->head;
	*(a->currentZer) = a->zerSubSys->head;

	*(a->posIndex) = UNROLLED_LL_SIZE;
	*(a->negIndex) = UNROLLED_LL_SIZE;
	*(a->zerIndex) = UNROLLED_LL_SIZE;
}

void FMEDSListRewind(FMEDS_t * a, char n)
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

void addToFMEDS(FMEDS_t * a, inequality_t * newIneq, int varIndex)
{

	if (ineqSgn(*newIneq, varIndex) == 1)
		specificAdd(a->posSubSys, a->currentPos, a->posIndex, newIneq);
	else if (ineqSgn(*newIneq, varIndex) == -1)
		specificAdd(a->negSubSys, a->currentNeg, a->negIndex, newIneq);
	else
	{
		specificAdd(a->zerSubSys, a->currentZer, a->zerIndex, newIneq);
	}

}

int getFMEDSSize(FMEDS_t * a)
{
	printf("Enter get size\n");
	int pSize = getListSize((a->posSubSys));
	int nSize = getListSize((a->negSubSys));
	int zSize = getListSize((a->zerSubSys));

	return (pSize + nSize + zSize);
	return (nSize);
}

void freeFMEDS(FMEDS_t * inputFMEDS)
{
//	printf("from free head address is %p\n", inputFMEDS->negSubSys->head);
	freeList(inputFMEDS->negSubSys);
	freeList(inputFMEDS->posSubSys);
	freeList(inputFMEDS->zerSubSys);

//	if(inputFMEDS->currentNeg != NULL)
	free((inputFMEDS->currentNeg));
//	if(inputFMEDS->currentPos != NULL)
	free(inputFMEDS->currentPos);
//	if(inputFMEDS->currentZer != NULL)
	free(inputFMEDS->currentZer);

	free(inputFMEDS->negIndex);
	free(inputFMEDS->posIndex);
	free(inputFMEDS->zerIndex);
}

