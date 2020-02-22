/**
 * FME_Support_kohler.c
 * ***********************
 * implementation of checks for redundancies
 */
#include "PolyhedralSets/FME_Support_kohler.h"


/**
 * checking if an inequality_t pass kohler's check w.r.t. a system
 * returns 0 if does not pass, 1 otherwise
 * a: inequality_t to check
 * d: inequality_t system
 * varNum: number of variables
 * ineqNum: number of inequalities
 */

int kohlerCheck(inequality_t a, inequalityULL_t * d, int v , int varNum , int ineqNum)
{
	int k = 0;
	int row = 0;
	mpq_t temp;
	mpq_init(temp);

	int col = v + 1  ;
	mpq_t * mat = (mpq_t *) malloc(ineqNum * varNum * sizeof(mpq_t));

	for (int i = 0; i < ineqNum * varNum; i++)
		mpq_init(mat[i]);

	node * current = d->head;

	for (int i = 0; i < ineqNum; i++) //n: number of initial inequalities
		if (a.history[i] == 1)
		{
			int no = (i / UNROLLED_LL_SIZE) + 1;
			int da = i % UNROLLED_LL_SIZE; //same

			current = d->head;

			for(int i1 = 0 ; i1 < no ; i1++)
				current = current->next;

			for (int j = 0; j <= v; j++)
				mpq_set(mat[k++], current->data[da].coef[j]);

			row++;
		}

	int rank = rankOfMatrix(mat, row, col);

	for(int i = 0 ; i < ineqNum * varNum ; i++)
		mpq_clear(mat[i]);
	free(mat);

	//if rank is row-1 then, inequality_t is can be redundant or non-redundant
	if (rank == row - 1)
		return 1;
	else
		//inequality_t is redundant
		return 0;
}

/**
 * balas check function returns 1 if inequality_t is not redundant, 0 otherwise
 * a: inequality_t to check
 * pw0: projection of redundancy test cone
 * varNum: number of variables
 * ineqNum: number of inequalities
 * sizepw0: number of inequalities in pw0
 * elimNum: number of eliminated variables
 */

int balasCheck(inequality_t a, mpq_t * pw0, int varNum, int ineqNum , int sizepw0 , int elimNum)
{

	int vcheck = varNum - elimNum;

	mpq_t * v = (mpq_t *) malloc((varNum - elimNum) * sizeof(mpq_t));
	for(int i = 0 ; i < varNum - elimNum; i++)
		mpq_init(v[i]);


	int c = 0;
	for(int i = elimNum + 1 ; i < varNum  ; i++)
		mpq_set(v[c++] , a.coef[i]);

	mpq_set(v[varNum - elimNum - 1] , a.constant);

	//check if the coefficient vector of the inequality_t is an extreme ray
	//of the projection of the redundancy test cone

	int r = checkBoundary(pw0, v, varNum - elimNum , sizepw0 , vcheck);
	return(r);
}

