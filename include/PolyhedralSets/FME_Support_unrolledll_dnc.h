/*
 * Implementation of unrolled linked list strucutre and its related functions
 */

#ifndef UNROLLEDLL_dnc_H
#define UNROLLEDLL_dnc_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "FME_Support_inequality.h"

#define UNROLLED_LL_SIZE 1

/*
 * This data strucuter is a linkedlist node for saving inequalities.
 */

typedef struct inequalityNode_dnc_t
{
	inequality_t * data;
	int fill;
	struct inequalityNode_dnc_t * next;
}inequalityNode_dnc_t;


typedef struct unrolledLl_dnc_t
{
	inequalityNode_dnc_t * head;
	int number;
}unrolledLl_dnc_t;


void allocNode_dnc(int varNum , inequalityNode_dnc_t * newIneq);

void copyNode_dnc(inequalityNode_dnc_t source  , inequalityNode_dnc_t * destination);

void freeNode_dnc(inequalityNode_dnc_t * a);

void allocList_dnc(int varNum , unrolledLl_dnc_t * newList);


void specificAdd_dnc(unrolledLl_dnc_t * l, inequalityNode_dnc_t ** current, int * index,
		inequality_t * newIneq);

int getListSize_dnc(unrolledLl_dnc_t * d);

void freeList_dnc(unrolledLl_dnc_t * l);


int getNext_dnc(inequalityNode_dnc_t ** currentNode, int * index,
		inequality_t ** nextIneq);

void printList_dnc(unrolledLl_dnc_t * d);

#endif
