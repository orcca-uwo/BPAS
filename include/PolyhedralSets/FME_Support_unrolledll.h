/*
 * Implementation of unrolled linked list strucutre and its related functions
 */

#ifndef UNROLLEDLL_H
#define UNROLLEDLL_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "FME_Support_inequality.h"

#define UNROLLED_LL_SIZE 10

/*
 * This data strucuter is a linkedlist node for saving inequalities.
 */

typedef struct inequalityNode_t
{
	inequality_t * data;
	int fill;
	struct inequalityNode_t * next;
}inequalityNode_t;


typedef struct unrolledLl_t
{
	inequalityNode_t * head;
	int number;
}unrolledLl_t;


void allocNode(int varNum , inequalityNode_t * newIneq);

void copyNode(const inequalityNode_t source  , inequalityNode_t * destination);

void freeNode(inequalityNode_t * a);

void allocList(int varNum , unrolledLl_t * newList);


void specificAdd(unrolledLl_t * l, inequalityNode_t ** current, int * index,
		const inequality_t * newIneq);

int getListSize(unrolledLl_t * d);

void freeList(unrolledLl_t * l);


int getNext(inequalityNode_t ** currentNode, int * index,
		inequality_t ** nextIneq);

void printList(unrolledLl_t * d);

#endif
