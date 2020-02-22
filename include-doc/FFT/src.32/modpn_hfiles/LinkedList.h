#ifndef __LinkedList_h
#define __LinkedList_h

#include "Types.h"
#include "MPMMTS.h"

RegularPair *EX_RegularPair_Init(preFFTRep * poly, TriSet *ts);

void EX_RegularPair_Print(void *element);

void EX_RegularPair_Free(void *element);

RegularListPair *EX_RegularListPair_Init(sfixn no, preFFTRep **polyList, TriSet *ts);

void EX_RegularListPair_Print(void *element);

void EX_RegularListPair_Free(void *element);

void EX_Poly_Print(void *element);

void EX_Poly_Free(void *element);

TaskPair *EX_TaskPair_Init(int index, TriSet *ts);

void EX_TaskPair_Print(void *element);

void EX_TaskPair_Free(void *element);

LinearNode *EX_LinearNode_Init(void *element);

void EX_LinearNode_Print(LinearNode *node, void (*printELement)(void *));

void EX_LinearNode_Free(LinearNode *node, void (*freeElement)(void *));

LinkedQueue *EX_LinkedQueue_Init();

int EX_LinkedQueue_IsEmpty(LinkedQueue *queue);

void EX_LinkedQueue_Print(LinkedQueue *queue, void (*printELement)(void *));

void EX_LinkedQueue_Free(LinkedQueue *queue, void (*freeElement)(void *));

void EX_LinkedQueue_Enqeue(LinkedQueue *queue, void *element);

void *EX_LinkedQueue_Deqeue(LinkedQueue *queue);

LinkedQueue *EX_LinkedQueue_Copy(LinkedQueue *queue, void *(*copyElement)(void *) );

int EX_LinkedQueue_Size(LinkedQueue *queue);

LinkedQueue * EX_LinkedQueue_Concat_1(LinkedQueue *queue1, LinkedQueue *queue2);

void **LinkedQueue2Array(LinkedQueue *queue, void *(*copyElement)(void *));

void *EX_CopyRegularPair(void *pair);

void *EX_CopyRegularListPair(void *pair);

void *EX_CopyPoly(void *element);

void EX_RegularListPair_List_Free(void *element);

void EX_RegularPair_List_Free(void *element);

#endif
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


