/**
 * FEM_Support_unrolledll.c
 * ***************************
 * implementation of inequality_t unrolled linked list
 */

#include "PolyhedralSets/FME_Support_unrolledll.h"

void allocNode(int varNum, inequalityNode_t * new)
{
	new->data = (inequality_t *) malloc(
			sizeof(inequality_t) * UNROLLED_LL_SIZE);

	for (int i = 0; i < UNROLLED_LL_SIZE; i++)
	{
		allocInequality(varNum, &(new->data[i]));
	}

	new->fill = 0;
	new->next = NULL;
}

void copyNode(const inequalityNode_t source, inequalityNode_t * destination)
{
	destination->fill = source.fill;
	destination->next = source.next;

	for (int i = 0; i < destination->fill; i++)
		copyInequality(&source.data[i], &destination->data[i]);
}


void allocList(int varNum, unrolledLl_t * new)
{
	inequalityNode_t * headNode = (inequalityNode_t *) malloc(
			sizeof(inequalityNode_t));
	allocNode(varNum, headNode);

	new->head = headNode;
	new->head->fill = UNROLLED_LL_SIZE; //we do not add any data to the head node of a list
	new->number = 0; //index number of the head node is zero
}

void specificAdd(unrolledLl_t * l, inequalityNode_t ** current, int * index,
		const inequality_t * newIneq)
{
	if ((*index) < UNROLLED_LL_SIZE)
	{
//		printf("Enter if\n");
		copyInequality(newIneq, &((*current)->data[*index]));
		(*index)++;
		((*current)->fill)++;
	}
	else
	{
//		printf("Enter else\n");

		(*current)->next = (inequalityNode_t *) malloc(sizeof(inequalityNode_t));

		allocNode(newIneq->dimension, (*current)->next);

		*current = (*current)->next;

		(l->number)++;
		(*index) = 0;
		copyInequality(newIneq, &((*current)->data[*index]));
		((*current)->fill)++;
		(*index) = 1;
	}
}



int getListSize(unrolledLl_t * d)
{

	inequalityNode_t * current = d->head->next;
	int size = 0;

	for (int i = 0; i < d->number; i++)
	{
		size = size + current->fill;
		current = current->next;
	}

	return (size);
}

void freeNode(inequalityNode_t * a)
{
	for (int i = 0; i < UNROLLED_LL_SIZE ; i++)
		freeInequality(&(a->data[i]));

	free(a->data);
}

void freeList(unrolledLl_t * l)
{
	inequalityNode_t * tmp1 = l->head;
	inequalityNode_t * tmp2 = tmp1;

	while (tmp1 != NULL)
	{
		tmp2 = tmp1;
		tmp1 = tmp1->next;
		freeNode(tmp2);
		free(tmp2);
	}
	free(l);
}

int getNext(inequalityNode_t ** currentNode, int * index,
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

//    printf("from getnext function, value of the pointer is %p\n", *nextIneq);

	return (1);
}
