/*
 * FME_Support_balas.c
 * ****************************
 * implementation of funcitons to find the
 * redundancy test cone (W^0) of an input polyhedron.
 */

#include "PolyhedralSets/FME_Support_balas.h"

/*
 * find coefficient matrix of eliminating variables.
 * input: the input inequality_t system
 * result: one-dimensional array containing coefficients
 * elimVarsNum: number of varialbes to be eliminated
 */

void whatIsA(inequalityULL_t input, mpq_t * result, int elimVarsNum)
{
	int c = 0;

	node * current = input.head->next;

	for (int i = 0; i < input.number; i++)
	{
		for (int j = 0; j < current->fill; j++)
			for (int k = 0; k < elimVarsNum; k++)
				mpq_set(result[c++], current->data[j].coef[k]);

		current = current->next;
	}

}

/*
 * find remained variables coefficient matrix
 * input: the input inequality_t system
 * result: one-dimensional array containing coefficients
 * elimVarsNum: number of varialbes to be eliminated
 * ineqNum: number of inequalities in the input system
 * varNum: dimension of the polyhedron (number of varialbes in the input system)
 */

void whatIsB(inequalityULL_t input, mpq_t * result, int elimVarsNum, int ineqNum,
		int varNum)
{
	int c = 0;

	node * current = input.head->next;

	for (int i = 0; i < input.number; i++)
	{
		for (int j = 0; j < current->fill; j++)
			for (int k = elimVarsNum; k < varNum; k++)
				mpq_set(result[c++], current->data[j].coef[k]);

		current = current->next;
	}

}

/*
 * find right hand side vector
 * input: the input inequality_t system
 * result: one-dimensional array containing coefficients
 * ineqNum: number of inequalities in the input system.
 */

void whatIsd(inequalityULL_t input, mpq_t * result, int ineqNum)
{
	int c = 0;

	node * current = input.head->next;

	for (int i = 0; i < input.number; i++)
	{
		for (int j = 0; j < current->fill; j++)
			mpq_set(result[c++], current->data[j].constant);

		current = current->next;
	}
}

/*
 * find B0 matrix from B.
 * The B0 matrix is used to find minimal representation of the projected polyhedron.
 * B is m * q and it is full columns rank
 * B: m*q, full-column rank matrix of remaining variables
 * B0: computed B0 matrix
 * p: number of eliminated variables
 */

void whatIsB0(mpq_t * B, mpq_t * B0, int m, int q, int p) //
{
	int c = 0;
	int r = rankOfMatrix(B, m, q); //this would be q, but we write it for generality.
	int n = m + q - r;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j < q)
				mpq_set(B0[i * n + j], B[i * q + j]);
			else
			{
				if (i < r)
					mpq_set_d(B0[i * n + j], 0);
				else
				{
					if (j - q == i - r)
						mpq_set_d(B0[i * n + j], 1);
					else
						mpq_set_d(B0[i * n + j], 0);
				}
			}
		}
	}
}

/*
 * find a variable indexed by v from an equation and substitute it in an
 * inequality_t system.
 * s1: the equation used for evaluating
 * s2: system of inequalities
 * v: index of a variable to be substitute
 * ineqNum: number of inequalities in the input system
 * varNum: number of variables in the input system
 */

void eval(inequality_t s1, inequalityULL_t s2, int v, int ineqNum, int varNum)
{
	int pprime;

	mpq_t mul, temp1, temp2, temp3;
	mpq_t * coef = (mpq_t *) malloc(varNum * sizeof(mpq_t));
	for (int i = 0; i < varNum; i++)
	{
		mpq_init(coef[i]);
		mpq_set(coef[i], s1.coef[i]);
	}

	mpq_init(mul);
	mpq_init(temp1);
	mpq_init(temp2);
	mpq_init(temp3);

	mpq_set(mul, coef[v]);

	node * current = s2.head->next;

	for (int i = 0; i < s2.number; i++)
	{
		for (int j = 0; j < current->fill; j++)
		{
			for (int k = 0; k < varNum; k++)

				if (k != v)
				{
					mpq_neg(temp3, coef[k]);
					mpq_div(temp1, temp3, mul);
					mpq_mul(temp2, current->data[j].coef[v], temp1);
					mpq_add(current->data[j].coef[k], current->data[j].coef[k],
							temp2);

					mpq_div(temp1, s1.constant, mul);
					mpq_mul(temp2, current->data[j].coef[v], temp1);
					mpq_sub(current->data[j].constant,
							current->data[j].constant, temp2);
				}
			mpq_set_d(current->data[j].coef[v], 0);
		}
		current = current->next;
	}

	for (int i = 0; i < varNum; i++)
		mpq_clear(coef[i]);
	free(coef);
	mpq_clear(mul);
	mpq_clear(temp1);
	mpq_clear(temp2);
	mpq_clear(temp3);

}
/*
 * check if coefficient of a variable is zero in all inequalities
 * returns 1 if this is true, zero otherwise.
 * l: system of linear inequalities
 * c: index of the variable to check
 */
int zeroCol(inequalityULL_t l, int c)
{
	mpq_t zero;
	mpq_init(zero);
	mpq_set_d(zero, 0);

	node * current = l.head->next;

	for (int i = 0; i < l.number; i++)
	{
		for (int j = 0; j < current->fill; j++)
			if (mpq_equal(current->data[j].coef[c], zero) == 0)
				return (0);
		current = current->next;
	}
	return (1);
}

int isUnit(node * current, int varNum, int j)
{
	mpq_t zero;
	mpq_init(zero);
	mpq_set_d(zero, 0);

	mpq_t one;
	mpq_init(one);
	mpq_set_d(one, 1);

	int s = 0;

	for (int i = 0; i < varNum; i++)
	{
		if (mpq_equal(current->data[j].coef[i], zero) == 0)
			if (mpq_equal(current->data[j].coef[i], one) == 0)
				return (0);
		if (mpq_equal(current->data[j].coef[i], one) != 0)
			s++;
		if (s > 1)
			return (0);

	}
	return (1);
}

/*
 * swap two variables in a system
 * l: inequality_t system
 * c1: index of the first variable
 * c2: index of the second variable
 */

void swapVar(inequalityULL_t * l, int c1, int c2)
{
	node * currentNode = l->head->next;

	mpq_t temp;
	mpq_init(temp);

	for (int i = 0; i < l->number; i++)
	{
		for (int j = 0; j < currentNode->fill; j++)
		{
			mpq_set(temp, currentNode->data[j].coef[c1]);
			mpq_set(currentNode->data[j].coef[c1],
					currentNode->data[j].coef[c2]);
			mpq_set(currentNode->data[j].coef[c2], temp);
		}
		currentNode = currentNode->next;
	}
}

//
/*find and return the redundancy test cone, W0.
 * W0 = {omega (B0 ^ -1)A = 0, omega (B0 ^ -1) >= 0, omega (B0 ^ -1) d >=0 }
 * A0: coefficient vector of variables to eliminate (find by "whatIsA" function )
 * A0: coefficient vector of variables to remain (find by "whatIsB" function )
 * d0: right-hand-side vector (find by "whatIsd" function )
 * m: number of inequalities in the input system
 * q: number of remaining variables
 * p: number of eliminating variables
 * r: is equal to q (for future improvements)
 */
inequalityULL_t whatIsW0(mpq_t * A0, mpq_t * B0, mpq_t * d0, int m, int q, int p,
		int r)
{
	int n = m + q - r; //r = q => n = m, this part is for possible generalization and it is not important
	int c;

	inequalityULL_t setEq = makeList(n + 1, n + 1);
	inequalityULL_t setInq = makeList(n + 1, n + 1);


	//1- find invers of the B0 matrix

	mpq_t * invB0 = (mpq_t *) malloc(n * n * sizeof(mpq_t));
	for (int i = 0; i < n * n; i++)
		mpq_init(invB0[i]);

	mpq_t * invB01 = (mpq_t *) malloc(n * n * sizeof(mpq_t));
	for (int i = 0; i < n * n; i++)
		mpq_init(invB01[i]);

	inverse(B0, invB01, n);

	c = 0;

	//1.1- change order of columns to have the desired structure
	for (int i = q * n; i < n * n; i++)
		mpq_set(invB0[c++], invB01[i]);
	for (int i = 0; i < q * n; i++)
		mpq_set(invB0[c++], invB01[i]);

	//1.2- make inequalities using B0 ^ -1
	c = 0;
	inequality_t * newp = (inequality_t *) malloc(sizeof(inequality_t));

	makeNewInequality(newp, n + 1, n + 1);

	for (int i = 0; i < n; i++) //iterate over inequalities
	{
		for (int j = 0; j < n; j++) //iterate over variables
			mpq_set(newp->coef[j], invB0[j * n + i]);

		mpq_set_d(newp->coef[n], 0);

		for (int k = 0; k < n + 1; k++)
			newp->history[k] = 0;

		newp->history[i] = 1;

		mpq_set_d(newp->constant, 0);

		addList(&setInq, *newp, n + 1, n + 1);
	}

	//2.1- find A0 * B0 ^ -1 matrix
	c = 0;
	mpq_t * C = (mpq_t *) malloc(n * p * sizeof(mpq_t));
	for (int i = 0; i < n * p; i++)
		mpq_init(C[c++]);

	c = 0;

	mpq_t * C1 = (mpq_t *) malloc(n * p * sizeof(mpq_t));
	for (int i = 0; i < n * p; i++)
		mpq_init(C1[c++]);

	matrixMatrixMult(invB01, A0, C1, n, n, p);

	//2.2- putting things in order
	c = 0;
	for (int i = p * q; i < p * n; i++)
		mpq_set(C[c++], C1[i]);
	for (int i = 0; i < p * q; i++)
		mpq_set(C[c++], C1[i]);

	//2.3- Make p equalities
	for (int i = 0; i < p; i++) //iterate over equalities
	{
		for (int j = 0; j < n; j++) //iterate over variables
			mpq_set((*newp).coef[j], C[j * p + i]);

		for (int k = 0; k < n + 1; k++)
			newp->history[k] = 0;

		(*newp).history[i] = 1;

		mpq_set_d((*newp).constant, 0);

		addList(&setEq, *newp, n + 1, n + 1);
	}

    //3- find the last inequality_t

	mpq_t * E = (mpq_t *) malloc(n * sizeof(mpq_t));
	mpq_t * E1 = (mpq_t *) malloc(n * sizeof(mpq_t));
	for (int i = 0; i < n; i++)
		mpq_init(E[i]);

	for (int i = 0; i < n; i++)
		mpq_init(E1[i]);

	matrixMatrixMult(invB01, d0, E1, n, n, 1);

	c = 0;
	for (int i = q; i < n; i++)
		mpq_set(E[c++], E1[i]);
	for (int i = 0; i < q; i++)
		mpq_set(E[c++], E1[i]);

	for (int j = 0; j < n; j++) //iterate over variables
	{
		mpq_neg(E[j], E[j]);
		mpq_set((*newp).coef[j], E[j]);
	}
	for (int k = 0; k < n + 1; k++)
		newp->history[k] = 0;

	(*newp).history[n] = 1;

	mpq_set_d((*newp).coef[n], 1);

	addList(&setInq, *newp, n + 1, n + 1);

	//4- make the cone to be full-dimensional
	// ** begin evaluating ** \\

	node * current = setEq.head->next;

	int pprime = 0;
	int evalflag = 1;
	int vflag = 0;

	for (int i = 0; i < setEq.number; i++)
	{
		for (int j = 0; j < current->fill; j++)
		{
			int v = 0;
			while (1)
			{
				if (v > n)
				{
					evalflag = 0;
					break;
				}
				if (mpq_sgn(current->data[j].coef[v]) == 0)
					v++;
				else
					break;
			}

			if (v != pprime)
				vflag = 1;
			if (evalflag == 1)
			{
				pprime += 1;
				eval(current->data[j], setInq, v, n + 1, n + 1);
				eval(current->data[j], setEq, v, p, n + 1);
			}
		}

		current = current->next;
	}

   //5- change order of some rows and some columns to make the desired structure

	if (vflag == 1)
	{
		int qprime = q + pprime;
		for (int i1 = 0; i1 < pprime; i1++)
			if (zeroCol(setInq, i1) == 0) //it is not zero
			{
				for (int j1 = i1 + 1; j1 < m - qprime; j1++)
					if (zeroCol(setInq, j1) == 1)
					{
						swapVar(&setInq, i1, j1);
						break;
					}
			}

		node * currentSub = setInq.head->next;
		node * currentNext = setInq.head->next;

		inequalityULL_t units = makeList(n + 1, n + 1);
		int unitsInx = 0;

		inequalityULL_t nounits = makeList(n + 1, n + 1);
		int nounitsInx = 0;

		int p1 = -1;
		int p2 = 0;
		mpq_t subtemp;
		mpq_init(subtemp);

		int histemp;

		for (int i1 = 0; i1 < setInq.number; i1++)
		{
			for (int j1 = 0; j1 < currentSub->fill; j1++)
			{
				p1++;
				if (isUnit(currentSub, n, j1) == 1 && p1 < qprime)
					addList(&units, currentSub->data[j1], n + 1, n + 1);

				if (isUnit(currentSub, n, j1) == 0 && p1 >= qprime
						&& p1 < m + 1)
					addList(&nounits, currentSub->data[j1], n + 1, n + 1);
			}
			currentSub = currentSub->next;
		}

		currentSub = setInq.head->next;
		p1 = 0;

		node * up = units.head->next;
		node * nup = nounits.head->next;

		for (int i1 = 0; i1 < setInq.number; i1++)
		{
			for (int j1 = 0; j1 < currentSub->fill; j1++)
			{
				p1++;
				if (isUnit(currentSub, n, j1) == 1 && p1 <= qprime)
				{

					for (int k1 = 0; k1 < n + 1; k1++)
						mpq_set(currentSub->data[j1].coef[k1],
								nup->data[nounitsInx].coef[k1]);
					mpq_set(currentSub->data[j1].constant,
							nup->data[nounitsInx].constant);
					if (nounitsInx == nup->fill)
						if (nup->next != NULL)
							nup = nup->next;
					if (nounitsInx != nup->fill)
						nounitsInx++;
				}

			}
		}

	}
////////////////////////////////////////////////////
	//6- cleaning section

	for (int i = 0; i < n * n; i++)
		mpq_clear(invB0[i]);
	free(invB0);

	for (int i = 0; i < n * n; i++)
		mpq_clear(invB01[i]);
	free(invB01);

	for (int i = 0; i < n * p; i++)
		mpq_clear(C[i]);
	free(C);

	for (int i = 0; i < n * p; i++)
		mpq_clear(C1[i]);
	free(C1);

	for (int i = 0; i < n; i++)
		mpq_clear(E[i]);
	free(E);

	for (int i = 0; i < n; i++)
		mpq_clear(E1[i]);
	free(E1);

	freeIneq(newp, n + 1, n + 1);
	return (setInq);
//
}
/*
 * find coefficient matrix of an inequality_t system
 */
void coefMatrix(inequalityULL_t input, mpq_t * mat, int varNum, int ineqNum)
{
	int c = 0;

	node * current = input.head->next;

	for (int i = 0; i < input.number; i++)
	{
		for (int j = 0; j < current->fill; j++)
			for (int k = 0; k < varNum; k++)
			{
				//if(isUnit(current , varNum , j) == 0)
				mpq_set(mat[c++], current->data[j].coef[k]);
			}

		current = current->next;
	}
}

/*
 * put everything together to find redundancy test cone from input system
 * data: input polyhedron
 * elimNumber: number of variables to eliminate
 * varNum: number of variables in the input system
 * ineqNum: number of inequalities in the input system
 * W0ll: the resulting cone
 */

void balasW0(inequalityULL_t data, int elimNumber, int varNum, int ineqNum,
		inequalityULL_t * W0ll)
{
	//1- find matrix A of coefficients of variables to eliminate
	mpq_t * A1 = (mpq_t *) malloc(elimNumber * ineqNum * sizeof(mpq_t));
	for (int k = 0; k < elimNumber * ineqNum; k++)
		mpq_init(A1[k]);

	mpq_t * A = (mpq_t *) malloc(elimNumber * ineqNum * sizeof(mpq_t));
	for (int k = 0; k < elimNumber * ineqNum; k++)
		mpq_init(A[k]);

	whatIsA(data, A, elimNumber);

	//

	//2- find matrix B of coefficients of variables that remain

	mpq_t * B1 = (mpq_t *) malloc(
			(varNum - elimNumber) * ineqNum * sizeof(mpq_t));
	for (int k = 0; k < (varNum - elimNumber) * ineqNum; k++)
		mpq_init(B1[k]);

	mpq_t * B = (mpq_t *) malloc(
			(varNum - elimNumber) * ineqNum * sizeof(mpq_t));
	for (int k = 0; k < (varNum - elimNumber) * ineqNum; k++)
		mpq_init(B[k]);

	whatIsB(data, B, elimNumber, ineqNum, varNum);
	//

   //3- find the right-hand-side vector
	mpq_t * d1 = (mpq_t *) malloc(ineqNum * sizeof(mpq_t));
	for (int k = 0; k < ineqNum; k++)
		mpq_init(d1[k]);

	mpq_t * d = (mpq_t *) malloc(ineqNum * sizeof(mpq_t));
	for (int k = 0; k < ineqNum; k++)
		mpq_init(d[k]);

	whatIsd(data, d, ineqNum);

	//

	//4- swap some rows if needed
	orderMatrix(A, B, d, A1, B1, d1, ineqNum, (varNum - elimNumber),
			elimNumber);

	//

	//5- find B0 matrix
	mpq_t * B0 = (mpq_t *) malloc(ineqNum * ineqNum * sizeof(mpq_t));
	for (int k = 0; k < ineqNum * ineqNum; k++)
		mpq_init(B0[k]);

	whatIsB0(B1, B0, ineqNum, varNum - elimNumber, elimNumber);

	//

	//6- using A,B0,d and their dimensions to find W0
	*W0ll = whatIsW0(A1, B0, d1, ineqNum, varNum - elimNumber, elimNumber,
			varNum - elimNumber);

	//7- cleaning
	for (int i = 0; i < elimNumber * ineqNum; i++)
	{
		mpq_clear(A1[i]);
		mpq_clear(A[i]);
	}

	free(A1);
	free(A);

	for (int i = 0; i < ineqNum; i++)
	{
		mpq_clear(d1[i]);
		mpq_clear(d[i]);
	}

	free(d1);
	free(d);

	for (int i = 0; i < ineqNum * ineqNum; i++)
		mpq_clear(B0[i]);

	free(B0);
	return;
}
/*
 * check if a vector is an extreme ray of a cone, using algebraic test
 * returns 1 if it is, 0 if it is not.
 * A: coefficient matrix of the cone
 * v: the vector to check
 * varNum: number of variables
 * ineqNum: number of inequalities
 * vcheck: an integer computed in another function
 */
int checkBoundary(mpq_t * A, mpq_t * v, int varNum, int ineqNum, int vcheck)
{
	mpq_t * d = (mpq_t *) malloc(ineqNum * sizeof(mpq_t));
	for (int i = 0; i < ineqNum; i++)
	{
		mpq_init(d[i]);
		mpq_set_d(d[i], 21);
	}
	matrixMatrixMult(A, v, d, ineqNum, varNum, 1);

	int e = 0;
	for (int i = 0; i < ineqNum; i++)
		if (mpq_sgn(d[i]) == 0)
			e++;

	mpq_t * w = (mpq_t *) malloc(e * varNum * sizeof(mpq_t));
	for (int i = 0; i < e * varNum; i++)
		mpq_init(w[i]);

	int c = 0;
	for (int i = 0; i < ineqNum; i++)
		if (mpq_sgn(d[i]) == 0)
			for (int j = 0; j < varNum; j++)
				mpq_set(w[c++], A[i * varNum + j]);

	c = 0;
	int r = rankOfMatrix(w, e, varNum);

	if (r == vcheck - 1)
		return (1);
	else
		return (0);
}

