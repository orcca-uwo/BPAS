#include "PolyhedralSets/FME_Support_unrolledll_dnc.h"

void allocNode_dnc(int varNum, inequalityNode_dnc_t * new)
{
	new->data = (inequality_t *) malloc(
			sizeof(inequality_t));

	for (int i = 0; i < 1; i++)
		allocInequality(varNum, &(new->data[i]));

	new->fill = 0;
	new->next = NULL;
}

void copyNode_dnc(inequalityNode_dnc_t source, inequalityNode_dnc_t * destination)
{
	destination->fill = source.fill;
	destination->next = source.next;

	for (int i = 0; i < destination->fill; i++)
		copyInequality(&source.data[i], &destination->data[i]);
}

/*
 * note: output of this allocation function is the object, not a pointer
 * to the object.
 */

void allocList_dnc(int varNum, unrolledLl_dnc_t * new)
{
	inequalityNode_dnc_t * headNode = (inequalityNode_dnc_t *) malloc(
			sizeof(inequalityNode_dnc_t));
	allocNode_dnc(varNum, headNode);

	new->head = headNode;
	new->head->fill = 1; //we do not add any data to the head node of a list
	new->number = 0; //index number of the head node is zero
}

void specificAdd_dnc(unrolledLl_dnc_t * l, inequalityNode_dnc_t ** current, int * index,
		inequality_t * newIneq)
{


	{
//		printf("Enter else\n");
		(*current)->next = (inequalityNode_dnc_t *) malloc(sizeof(inequalityNode_dnc_t));

		allocNode_dnc(newIneq->dimension, (*current)->next);
//		printf("alloc node\n");

		*current = (*current)->next;

		(l->number)++;
		(*index) = 0;
		copyInequality(newIneq, &((*current)->data[*index]));
		((*current)->fill)++;
		(*index) = 1;
	}
}



int getListSize_dnc(unrolledLl_dnc_t * d)
{

	inequalityNode_dnc_t * current = d->head->next;

	int size = 0;

	for (int i = 0; i < d->number; i++)
	{
		size = size + 1;
		current = current->next;
	}

	return (size);
}

void freeNode_dnc(inequalityNode_dnc_t * a)
{
	for (int i = 0; i < UNROLLED_LL_SIZE ; i++)
		freeInequality(&(a->data[i]));

	free(a->data);
}

void freeList_dnc(unrolledLl_dnc_t * l)
{
	inequalityNode_dnc_t * tmp1 = l->head;
	inequalityNode_dnc_t * tmp2 = tmp1;

	while (tmp1 != NULL)
	{
		tmp2 = tmp1;
		tmp1 = tmp1->next;
		freeNode_dnc(tmp2);
		free(tmp2);
	}
	free(l);
}

int getNext_dnc(inequalityNode_dnc_t ** currentNode, int * index,
		inequality_t ** nextIneq)
{
	if ((*index) < ((*currentNode)->fill) - 1)
		(*index)++;

	else
	{
		if ((*currentNode)->next == NULL)
			return (0);
		*currentNode = (*currentNode)->next;
		(*index) = 0;
	}

	*nextIneq = &((*currentNode)->data[*index]);


	return (1);
}

void printList_dnc(unrolledLl_dnc_t * d)
{

	inequalityNode_dnc_t * current = d->head->next;
	for(int i = 0 ; i < d->number ; i++)
	{
		printInequality(&(current->data[0]),'x');
		current = current->next;
	}
}
