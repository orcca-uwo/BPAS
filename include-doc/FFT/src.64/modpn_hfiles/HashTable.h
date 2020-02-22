/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __HashTable_h
#define __HashTable_h 

#include "Types.h"
#include "UniHensel.h"

#define  HashModuli 101

typedef struct node_struct NODE; 
typedef struct hashtable_struct HASHTABLE;

struct node_struct {
  operand opnd;
  struct node_struct *next;	
};


// array index is the key value.
struct hashtable_struct
{
        NODE ** table;
	uint32 length;
};


NODE *list_create(operand opnd);

NODE *list_insert_after(NODE *node, operand opnd);

NODE *list_insert_beginning(NODE *list, operand opnd);

int32 list_remove(NODE *list, NODE *node);

int32 compareOperand(operand opnd1, operand opnd2);

NODE *list_find(NODE *node, operand opnd);

void list_free(NODE *list);

int32 hashFunction(operand opnd);


// yes, return the existing node.
// no, return a NULL;
operand searchNodeinHashTable(operand opnd, HASHTABLE *hash_table);



HASHTABLE *newHashTable(uint32 length);


void freeHashTable(HASHTABLE *hash_table);


operand  searchNodeinHashTable(operand opnd, HASHTABLE *hash_table);

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


