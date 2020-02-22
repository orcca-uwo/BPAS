

#include "RationalNumberPolynomial/SMQP_Support_Recursive.h"
#include "RationalNumberPolynomial/SMQP_Support_Test.h"
#include <time.h>

/** 
 * Build a random recursive polynomial. 
 * Nterms is the number of terms wrt the main variable. 
 * coefNTerms is the number of terms in the coefficients (when viewed recursively, and nvar > 1);
 * coefBound is the absolute upper limit on the rational numbers in the polynomial.
 * sparsity is the difference between degrees (wrt to main variable) of successive terms - 1. Therefore 2 is dense. 
 * sparsity is also passed when constructing coefficient polynomials.
 */
RecNode_t* buildRandomRecPoly(int nvar, int nterms, int coefNterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
	
	static int initRand = 0;
	if (!initRand) {
		time_t seed = time(NULL); 
		srand(seed);
		initRand = 1;
	}

	RecNode_t* poly = NULL, *tail = NULL;
	degree_t nextDeg = 0;
	// if (sparsity == 2) {
	// 	//if dense always include the constant
	// 	nextDeg = 0;
	// } else {
	// 	nextDeg = rand() % sparsity; //50/50 chance of first term being a constant;
	// }
	for (int i = 0; i < nterms; ++i) {
		RecNode_t* node = (RecNode_t*) malloc(sizeof(RecNode_t));
		node->next = NULL;
		node->exp = nextDeg;

		Node* coefNode = NULL;
		if (nvar == 1) {
			coefNode = (Node*) malloc(sizeof(Node));
			coefNode->next = NULL;
			long int coef_l = rand() % coefBound;
			if (coef_l == 0) {
				coef_l = 1l;
			}
			if (includeNeg && rand() % 2) {
				//50/50 chance of being negative
				coef_l *= -1;
			}
			coefNode->degs = (degrees_t) calloc(1, sizeof(degree_t));
			mpq_init(coefNode->coef);
			mpq_set_si(coefNode->coef, coef_l, 1ul);	
		} else {
			coefNode = buildRandomPoly(nvar-1, coefNterms, coefBound, sparsity, includeNeg);
			degrees_t degs, tempDegs;
			Node* curNode = coefNode;
			while (curNode != NULL) {
				//expand degs by 1. 
				degs = curNode->degs;
				tempDegs = (degrees_t) calloc(nvar, sizeof(degree_t));
				memcpy(&tempDegs[1], degs, (nvar-1)*sizeof(degree_t));
				curNode->degs = tempDegs;
				free(degs);
				curNode = curNode->next;
			}
		}
		node->coef = coefNode;

		if (poly == NULL) {
			poly = node;
			tail = node;
		} else {
			tail->next = node;
			tail = node;
		}

		//get next degree
		degree_t step = sparsity - 1;// rand() % sparsity;
		// while (step == 0) {
			// step = rand() % sparsity;
		// }
		nextDeg += step;
	}

	//Now reverse the poly so it is in decreasing order;
	RecNode_t* cur = poly;
	RecNode_t* prev = poly;
	RecNode_t* next = poly->next;
	cur->next = NULL;
	while(next != NULL) {
		cur = next;
		next = cur->next;
		cur->next = prev;
		prev = cur;
	}
	poly = cur;

	return poly;
} 


RecNode_t* convertToRecursiveNodeAtIdx(Node* poly, int idx) {

	if (poly == NULL) {
		return NULL;
	}

	Node* node = poly;
	RecNode_t* head = NULL;
	RecNode_t* tail = NULL;
	Node* coefHead = poly;
	Node* coefTail = poly;
	degree_t curExp = poly->degs[idx];
	degree_t lastExp = curExp;
	while (node != NULL) {
		curExp = node->degs[idx];
		if (curExp != lastExp) {
			//then we need to commit our current coef chain.
			tail = addRecTerm(tail, coefHead, lastExp);
			if (head == NULL) {
				head = tail;
			}

			coefTail->next = NULL;
			coefHead = node;
			lastExp = curExp;
		}

		node->degs[idx] = 0;
		coefTail = node;

		node = node->next;
	}
	//one last "commit" needed.
	tail = addRecTerm(tail, coefHead, curExp);
	if (head == NULL) {
		head = tail;
	}

	return head;
}

RecNode_t* convertToRecursiveNode(Node* poly) {
	if (poly == NULL) {
		return NULL;
	}

	Node* node = poly;
	RecNode_t* head = NULL;
	RecNode_t* tail = NULL;
	Node* coefHead = poly;
	Node* coefTail = poly;
	degree_t curExp = poly->degs[0];
	degree_t lastExp = curExp;
	while (node != NULL) {
		curExp = node->degs[0];
		if (curExp != lastExp) {
			//then we need to commit our current coef chain.
			tail = addRecTerm(tail, coefHead, lastExp);
			if (head == NULL) {
				head = tail;
			}
			
			coefTail->next = NULL;
			coefHead = node;
			lastExp = curExp;
		}

		node->degs[0] = 0;
		coefTail = node;

		node = node->next;
	}
	//one last "commit" needed.
	tail = addRecTerm(tail,coefHead,curExp);
	if (head == NULL) {
		head = tail;
	}

	return head;

}

RecNode_t* convertToRecursiveNodeCopy(Node* poly, int nvar) {
	if (poly == NULL) {
		return NULL;
	}

	Node* node = poly;
	RecNode_t* head = NULL;
	RecNode_t* tail = NULL;
	Node* coefHead = NULL;
	Node* coefTail = NULL;
	degree_t curExp = poly->degs[0];
	degree_t lastExp = curExp;
	while (node != NULL) {
		curExp = node->degs[0];
		if (curExp != lastExp) {
			//then we need to commit our current coef chain.
			tail = addRecTerm(tail, coefHead, lastExp);
			if (head == NULL) {
				head = tail;
			}
			coefHead = NULL;
			coefTail = NULL;
			lastExp = curExp;
		}

		degrees_t newDegs = (degrees_t) malloc(sizeof(degree_t) * nvar-1);
		for (int i = 0; i < nvar-1; ++i) {
			newDegs[i] = node->degs[i+1];
		}
		coefTail = addTerm(coefTail, newDegs, node->coef);
		if (coefHead == NULL) {
			coefHead = coefTail;
		}

		node = node->next;
	}
	//one last "commit" needed.
	tail = addRecTerm(tail,coefHead,curExp);
	if (head == NULL) {
		head = tail;
	}

	return head;

}

Node* convertFromRecursiveNodeAtIdx(RecNode_t* poly, int idx) {
	if (poly == NULL) {
		return NULL;
	}

	Node* coefHead = poly->coef;
	Node* coefNode = NULL;
	Node* coefTail = NULL;
	RecNode_t* recNode = poly, *recNext = NULL;
	degree_t exp;
	while (recNode != NULL) {
		exp = recNode->exp;
		coefNode = recNode->coef;

		//Link last coef to new coef.
		if (coefTail != NULL) {
			coefTail->next = coefNode;
		}

		while (coefNode != NULL) {
			coefNode->degs[idx] = exp;
			coefTail = coefNode;
			coefNode = coefNode->next;
		}

		recNext = recNode->next;
		free(recNode);
		recNode = recNext;
	}

	return coefHead;
}


Node* convertFromRecursiveNode(RecNode_t* poly) {
	if (poly == NULL) {
		return NULL;
	}

	Node* coefHead = poly->coef;
	Node* coefNode = NULL;
	Node* coefTail = NULL;
	RecNode_t* recNode = poly, *recNext = NULL;
	degree_t exp;
	while (recNode != NULL) {
		exp = recNode->exp;
		coefNode = recNode->coef;

		//Link last coef to new coef.
		if (coefTail != NULL) {
			coefTail->next = coefNode;
		}

		while (coefNode != NULL) {
			coefNode->degs[0] = exp;
			coefTail = coefNode;
			coefNode = coefNode->next;
		}

		recNext = recNode->next;
		free(recNode);
		recNode = recNext;
	}

	return coefHead;
}

Node* convertFromRecursiveNodeCopy(RecNode_t* poly, int nvar) {
	if (poly == NULL) {
		return NULL;
	}

	Node* head = NULL;
	Node* tail = NULL;
	Node* coefNode = NULL;
	RecNode_t* recNode = poly;
	degree_t exp;
	degrees_t degs;
	while (recNode != NULL) {
		exp = recNode->exp;
		coefNode = recNode->coef;

		while (coefNode != NULL) {
			degs = coefNode->degs;
			degrees_t newDegs = (degrees_t) malloc(sizeof(degree_t)*(nvar+1));
			newDegs[0] = exp;
			for (int i = 1; i < nvar+1; ++i) {
				newDegs[i] = degs[i-1]; 
			}
			tail = addTerm(tail, newDegs, coefNode->coef);
			if (head == NULL) {
				head = tail;
			}

			coefNode = coefNode->next;
		}

		recNode = recNode->next;
	}

	return head;
}

RecNode_t* deepCopyRecPolynomial(RecNode_t* poly, int nvar) {
	RecNode_t* recNode = poly;
	RecNode_t* head = NULL;
	RecNode_t* tail = NULL;
	Node* coefCopy;
	while (recNode != NULL) {
		coefCopy = deepCopyPolynomial(recNode->coef, nvar);
		tail = addRecTerm(tail, coefCopy, recNode->exp);
		if (head == NULL) {
			head = tail;
		}
		recNode = recNode->next;
	}
	return head;
}


char* vars[] = {"x", "y", "z", "t", "u", "v", "w"};

void printPolyNoNewline(Node* node, int nvar) {
	while (node != NULL) {
		gmp_fprintf(stderr, "%Qd*", node->coef);
		for(int i = 0; i < nvar; ++i) {
			fprintf(stderr, "%s^%lu", vars[i], node->degs[i]);
			if (i < nvar -1) {
				fprintf(stderr, "*");
			}
		}
		fprintf(stderr, " + ");
		node = node->next;
	}
}


void printPoly(Node* node, int nvar) {
	printPolyNoNewline(node, nvar);
	fprintf(stderr, "\n");
}

void printRecPoly(RecNode_t* recNode, int nvar) {
	while(recNode != NULL) {
		fprintf(stderr, "(");
		printPoly(recNode->coef, nvar);
		fprintf(stderr, ")");
		fprintf(stderr, "*x^%lu + ",recNode->exp);
		recNode = recNode->next;
	}
	fprintf(stderr, "\n");
}



/***********
 * Pesudo Division
 ***********/


/**
 * Print the product degrees_t currently in the heap.
 */
void recProdHeapPrint(RecProdHeap_t* h) {
	for (int i = 0; i < h->heapSize; ++i){
		fprintf(stderr, "( %lx, (", (unsigned long int) h->elements[i]);
		Node* prodCoef = multiplyPolynomials(h->elements[i]->a_i->coef, h->elements[i]->b->coef, h->nvar);
		printPolyNoNewline(prodCoef, h->nvar);
		free(prodCoef);
		// fprintf(stderr, "%lx, ", (unsigned long int) h->elements[i]->prodCoef);
		fprintf(stderr, "), %lu", (h->elements[i]->prodExp));
		RecProdHeapElem_t* next = h->elements[i]->next;
		while (next != NULL) {
			fprintf(stderr, "; ");
			fprintf(stderr, "%lx, (", (unsigned long int) next);
			prodCoef = multiplyPolynomials(next->a_i->coef, next->b->coef, h->nvar);
			printPolyNoNewline(prodCoef, h->nvar);
			free(prodCoef);
			// fprintf(stderr, "%lx, ", (unsigned long int) next->prodCoef);
			fprintf(stderr, "), %lu", (next->prodExp));
			next = next->next;
		}
		fprintf(stderr, ")");
		fprintf(stderr," ");
	}
	fprintf(stderr, "\n");
}

RecProdHeapElem_t* recProdHeapMakeElement(RecNode_t* a_i, RecNode_t* b, int nvar) {
	RecProdHeapElem_t* elem = (RecProdHeapElem_t*) malloc(sizeof(RecProdHeapElem_t));
	elem->a_i = a_i;
	elem->b = b;
	elem->prodExp = a_i->exp + b->exp;
	// elem->prodCoef = multiplyPolynomials(a_i->coef, b->coef, nvar);
	elem->next = NULL;
	return elem;
}

RecProdHeap_t* recProdHeapCreate(int nvar) {
	RecProdHeap_t* h = (RecProdHeap_t*) malloc(sizeof(ProductHeap));
	h->elements = NULL;
	h->heapSize = 0;
	h->maxHeapSize = 0;
	h->nvar = nvar;
	return h;
}

void recProdHeapInsert(RecProdHeap_t* h, RecProdHeapElem_t* elem) {
	polysize_t s = h->heapSize;
	polysize_t N = h->maxHeapSize;
	if (s == 0) {
		if (s >= N) {
			N = (s == 0) ? 1 : s * 2;
			h->elements = (RecProdHeapElem_t**) realloc(h->elements, N*sizeof(RecProdHeapElem_t*));
			h->maxHeapSize = N;
		}
		h->elements[0] = elem;
		h->heapSize = 1;
		return;
	}

	RecProdHeapElem_t** elems = h->elements;
	
	degree_t insertExp = elem->prodExp;	

	//first check if we can chain off the root
	if (insertExp == elems[0]->prodExp) {
		elem->next = elems[0];
		elems[0] = elem;
		return;
	}

	//otherwise, we must search the heap to find the new product's insertion point
	//note that since we are looking for chains we cannot use the simple swim method
	//we sort of fake the swimming, looking for a chain to be made or the eventual
	//place where the swim would stop. At this point, we insert the new elem
	//in that spot, and "push" the entire path we took down a level. Assuming 
	//that we insert e and it ends up at the root, we push down the 'x' path
	//
	//     x     --->    e
	//    / \           / \ 
	//   x   o         x   o
	//                /
 	//               x

	polysize_t i = (s-1)/2; //i is parent 
	polysize_t j = s;       //j is current insertion point
	long long unsigned int path = 1;
	degree_t cmpExp;
	while (j > 0) {
		cmpExp = elems[i]->prodExp;
		if (cmpExp == insertExp) {
			elem->next = elems[i];
			elems[i] = elem;
			return;
		} else if (cmpExp < insertExp) {
			path <<= 1;
			if (!(j & 1)) {
				//set trailing bit to 1 to distinguish left/right path
				path += 1;
			}
			j = i;
			i = (i - 1) >> 1;
		} else {
			break; //found the insert point
		}
	}

	//if we made it this far, then we are adding a new element to the heap. 
	//resize if necessary
	if (s >= N) {
		N = (s == 0) ? 1 : s * 2;
		h->elements = (RecProdHeapElem_t**) realloc(h->elements, N*sizeof(RecProdHeapElem_t*));
		elems = h->elements;
		h->maxHeapSize = N;
	}

	//then j is now the place we need to insert elem;
	//do so, and then push all others down the path, inserting the last
	//as the new element in elems[s];
	RecProdHeapElem_t* temp;
	while (j <= s) {
		temp = elems[j];
		elems[j] = elem;
		elem = temp;
		j = 2*j + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
}

degree_t recProdHeapPeek(RecProdHeap_t* h) {
	if (h->heapSize > 0) {
		return h->elements[0]->prodExp;
	}

	return -1;
}

RecProdHeapElem_t* recProdHeapRemoveMax(RecProdHeap_t* h) {
	//TODO
	if (h->heapSize == 0) {
		return NULL;
	}

	RecProdHeapElem_t** elems = h->elements;
	RecProdHeapElem_t* maxElem = elems[0];
	int i = 0, j = 1;
	int s = --(h->heapSize);
	
	//promote largest children
	while (j < s) {
		if (j+1 < s && elems[j]->prodExp < elems[j+1]->prodExp) {
			++j;
		}
		elems[i] = elems[j];
		i = j;
		j = (2*j) + 1;
	}
	//now place last element into i and swim up to make tree complete 
	j = (i-1)/2;
	while(i > 0) {
		if (elems[s]->prodExp < elems[j]->prodExp) {
			break;
		}
		elems[i] = elems[j];
		i = j;
		j = (j-1)/2;
	}
	elems[i] = elems[s]; 
	return maxElem;
}

Node* recProdHeapGetNextCoef(RecProdHeap_t* h) {
	if (h->heapSize == 0) {
		return NULL;
	}


	RecProdHeapElem_t** elems = h->elements;
	polysize_t size = h->heapSize;
	int nvar = h->nvar;

	degree_t maxExp = elems[0]->prodExp;
	Node* ret = NULL;
	RecProdHeapElem_t* insertChain = NULL;
	RecProdHeapElem_t* maxElem, * nextMaxElem;

	//It's possible that not all elemtns of same degree are chained. So much loop like this.
	while(size > 0 && elems[0]->prodExp == maxExp) {
		maxElem = recProdHeapRemoveMax(h);
		--size;

		//go through the chain now;
		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;

			//Notice here we delay the calculation of the actual product coef.
			//This is because a_i will change over the course of the algorithm.
			if (ret == NULL) {
				ret = multiplyPolynomials(maxElem->a_i->coef, maxElem->b->coef, nvar);
				// ret = maxElem->prodCoef;
			} else {
				Node* prod = multiplyPolynomials(maxElem->a_i->coef, maxElem->b->coef, nvar);
				ret = addPolynomials_inp(ret, prod, nvar);
				freeNode(prod);
			}
			maxElem->b = maxElem->b->next;
			if (maxElem->b != NULL) {
				maxElem->prodExp = maxElem->a_i->exp + maxElem->b->exp;
				// maxElem->prodCoef = multiplyPolynomials(maxElem->a_i->coef, maxElem->b->coef, nvar);
				maxElem->next = insertChain;
				insertChain = maxElem;
			} else {
				free(maxElem);
			}
			maxElem = nextMaxElem;
		}
	}

	while(insertChain != NULL) {
		nextMaxElem = insertChain->next;
		insertChain->next = NULL;
		recProdHeapInsert(h, insertChain);
		insertChain = nextMaxElem;
	}

	return ret;
}

// Deprecated. Left here for future reference.
// void recProdHeapUpdateQuoByH (RecProdHeap_t* h, Node* hNode) {
// 	RecProdHeapElem_t** elems = h->elements;
// 	polysize_t s = h->heapSize;
// 	int nvar = h->nvar;
// 	RecProdHeapElem_t* elem;
// 	for (int i = 0; i < s; ++i) {
// 		elem = elems[i];
// 		while (elem != NULL) {
// 			// elem->prodCoef = multiplyPolynomials_inp(elem->prodCoef, hNode, nvar);
// 			elem = elem->next;
// 		}
// 	}
// }


/**
 * Multiply a recursively-viewed polynomial by a "coefficient" polynomial, in place.
 * That is, the coefficient supplied is a polynomial whose variables do not include
 * the mvar of the supplied RecNode_t a. 
 *
 * For optimization purposes, it may be that there is a slot in the coefficient poly's
 * exponent vector for the mvar of a, but it should all be 0. 
 *
 * Returns a pointer to the new head of a after multiplication. 
 */
RecNode_t* recNodeMultiplyByCoef_inp(RecNode_t* a, Node* coef, int nvar) {
	RecNode_t* curA = a; 
	while (curA != NULL) {
		curA->coef = multiplyPolynomials_inp(curA->coef, coef, nvar);
		curA = curA->next;
	}
	return a;
}





void pesudoDivideOneTerm(RecNode_t* c, RecNode_t* b, RecNode_t** res_a, RecNode_t** res_r, int* e, Node** hPow, int nvar, int lazy) {
	if (b == NULL) {
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (b->exp > c->exp) {
		*res_r = deepCopyRecPolynomial(c, nvar);

		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		if (hPow != NULL) {
			mpq_t mpqCoef;
			mpq_init(mpqCoef);
			mpq_set_si(mpqCoef, 1l, 1l);
			degrees_t degs = calloc(nvar, sizeof(degree_t));
			*hPow = addTerm(NULL,degs,mpqCoef);
			mpq_clear(mpqCoef);
		}
		return;
	}

	int i = 1;
	Node* h = b->coef;
	degree_t bDeg = b->exp, eps = c->exp - bDeg;
	Node* multerm = deepCopyPolynomial(c->coef, nvar);
	RecNode_t* heada = addRecTerm(NULL, multerm, eps); 
	RecNode_t* a = heada;
	RecNode_t* r = NULL, *rhead = NULL;
	RecNode_t* k = c->next;

	Node* hI = deepCopyPolynomial(h, nvar);

	while (k != NULL) {
		eps = k->exp;
		multerm = multiplyPolynomials(k->coef, hI, nvar);
		k = k->next;

		if (eps >= bDeg) {
			heada = recNodeMultiplyByCoef_inp(heada, h, nvar);
			a = addRecTerm(a, multerm, eps - bDeg);
			++i;
			hI = multiplyPolynomials_inp(hI, h, nvar); //in-place wrt first arg.
		} else {
			//accumulate remainder.
			r = addRecTerm(r, multerm, eps);
			if (rhead == NULL) {
				rhead = r;
			}
		}
	}

	if (!lazy) {
		int d = c->exp - b->exp + 1 - i;
		i += d;
		for (int j = 0; j < d; ++j) {
			heada = recNodeMultiplyByCoef_inp(heada, h, nvar);
			rhead = recNodeMultiplyByCoef_inp(rhead, h, nvar);
			hI = multiplyPolynomials_inp(hI, h, nvar);
		}
	}

	*res_a = heada;
	*res_r = rhead;

	//Return number of division steps;
	if (e != NULL) {
		*e = i;
	}

	//Return the initial to power of e in hPow;
	if (hPow != NULL) {
		*hPow = hI;
	}

}


/** 
 * Do the pesudo division of c by b. Quotient is returned in res_a, and remainder
 * in res_r. 
 *
 * If e is not NULL then it returns the exact number of division steps which occurred.
 *
 * If hPow is not NULL then *hPow is set to the initial of b to the power of e.
 *
 * If lazy is 0 then e will be equal to deg(c) - deg(b) + 1.
 */
void pesudoDivide(RecNode_t* c, RecNode_t* b, RecNode_t** res_a, RecNode_t** res_r, int* e, Node** hPow, int nvar, int lazy) {

	if (b == NULL) {
		fprintf(stderr, "Division by zero! Exiting...\n");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (b->exp > c->exp) {
		*res_r = deepCopyRecPolynomial(c, nvar);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		if (hPow != NULL) {
			mpq_t mpqCoef;
			mpq_init(mpqCoef);
			mpq_set_si(mpqCoef, 1l, 1l);
			degrees_t degs = calloc(nvar, sizeof(degree_t));
			*hPow = addTerm(NULL,degs,mpqCoef);
			mpq_clear(mpqCoef);
		}
		return;
	}

	if (b->next == NULL) {
		pesudoDivideOneTerm(c, b, res_a, res_r, e, hPow, nvar, lazy);
		return;
	}

	//manually do first div as we know it comes from first term of c;
	int i = 1;
	Node* h = b->coef;
	degree_t bDeg = b->exp, eps = c->exp - bDeg;
	long long int delta = 1;
	Node* multerm = deepCopyPolynomial(c->coef, nvar);
	Node* multerm2 = NULL;
	RecNode_t* heada = addRecTerm(NULL, multerm, eps); 
	RecNode_t* a = heada;
	RecNode_t* r = NULL, *rhead = NULL;
	RecNode_t* k = c->next;
	long long int kDeg = k->exp;

	//TODO find a better way to do this accumulation.. 
	Node* hI = deepCopyPolynomial(h, nvar);

	RecProdHeap_t* prodHeap = recProdHeapCreate(nvar);
	recProdHeapInsert(prodHeap, recProdHeapMakeElement(a, b->next, nvar));


	while (delta > -1 || kDeg > -1) {
		//get the leading term of dividend quotient-product difference
		delta = recProdHeapPeek(prodHeap);

		// fprintf(stderr, "af peek delta: %ld , kdeg: %ld\n", delta, kDeg);

		if (delta == -1 && kDeg == -1) {
			break;
		}

		if (delta > kDeg) {
			eps = delta;
			//then ctilde is from quotient product
			multerm = recProdHeapGetNextCoef(prodHeap);
			
			if (multerm == NULL) {
				//the elements in the product with degree delta
				// ended up canceling out coefs
				continue; 
			}
			negatePolynomial(multerm);
		} else if (delta < kDeg) {
			eps = kDeg;
			//then ctilde is from dividend
			multerm = multiplyPolynomials(k->coef, hI, nvar);
			k = k->next;
			kDeg = k == NULL ? -1 : k->exp;
			// fprintf(stderr, "delta < kdeg. Eps = %u\n", eps);
		} else {
			eps = delta;
			//combine both
			multerm2 = recProdHeapGetNextCoef(prodHeap);
			if (multerm2 == NULL) {
				continue;
			}
			multerm = multiplyPolynomials(k->coef, hI, nvar);
			// fprintf(stderr,"Multerm in choosing: \n") ;
			// printPoly(multerm, nvar);
			// fprintf(stderr, "\nMulterm2: \n" );
			// printPoly(multerm2, nvar);
			//actually getting the next term might result in 0
			// fprintf(stderr, "doing subtraction\n");
			// This is in-place wrt to first variable.
			multerm = subPolynomials_inp(multerm, multerm2, nvar);
			freePolynomial(multerm2);
			// fprintf(stderr,"Multerm in choosing: \n") ;
			// printPoly(multerm, nvar);
			
			k = k->next;
			kDeg = k == NULL ? -1 : k->exp;
			// fprintf(stderr, "delta == kdeg. Eps = %u\n", eps);
			
			if (multerm == NULL) {
				//if sub resulted in a zero then we we must get a new multerm
				continue;
			}
		}

		// fprintf(stderr, "\n\nAfter extract: \n\n");
		// recProdHeapPrint(prodHeap);

		//multerm is now the leading coef (eps the degree) of the current difference

		// fprintf(stderr, "\n\nBefore insert: \n\n");
		// recProdHeapPrint(prodHeap);		

		if (eps >= bDeg) {
			// fprintf(stderr, "Adding to quotient: %lu\n", eps - bDeg);
			// fprintf(stderr, "multerm: ");
			// gmp_fprintf(stderr, "%Qd\n", multerm->coef);
			// printPoly(multerm, nvar);
			//do division of ctilde and leading term of b
			//since we multiply all by h, h cancels in division of 
			//ctilde by b, leaving only exponent subtraction necessary.

			//update quotient, first multiply cur quotient by h, then append new term
			//this is done in-place for heada. 
			// gmp_fprintf(stderr, "h: %Qd\n", h->coef);
			heada = recNodeMultiplyByCoef_inp(heada, h, nvar);
			a = addRecTerm(a, multerm, eps - bDeg);

			// IF the above rec node mult was not done in place, then
			// we would need to do the below. Since we delay the actual
			// calculation of the coef product from heap elements until just 
			// before extraction.
			// recProdHeapUpdateQuoByH(prodHeap, h);

			//insert b_2! Since we constructed a to exactly cancel b_1, useless.
			recProdHeapInsert(prodHeap, recProdHeapMakeElement(a, b->next, nvar));

			// fprintf(stderr, "\nCurrent quotient: \n");
			// printRecPoly(heada, nvar);

			//update h counter
			//this is only done when we actually add a new element to the quotient.
			++i;
			hI = multiplyPolynomials_inp(hI, h, nvar); //in-place wrt first arg.
		} else {
			// fprintf(stderr, "Adding to remainder: %lu\n\n", eps);
			// fprintf(stderr, "multerm: ");
			// printPoly(multerm, nvar);
			//accumulate remainder.
			r = addRecTerm(r, multerm, eps);
			if (rhead == NULL) {
				rhead = r;
			}
		}

		// fprintf(stderr, "\n\nAfter insert: \n\n");
		// recProdHeapPrint(prodHeap);		
		
		// fprintf(stderr, "end loop\n\n");

	}


	if (!lazy) {
		int d = c->exp - b->exp + 1 - i;
		i += d;
		for (int j = 0; j < d; ++j) {
			heada = recNodeMultiplyByCoef_inp(heada, h, nvar);
			rhead = recNodeMultiplyByCoef_inp(rhead, h, nvar);
			hI = multiplyPolynomials_inp(hI, h, nvar);
		}

	}

	*res_a = heada;
	*res_r = rhead;


	//Return number of division steps;
	if (e != NULL) {
		*e = i;
	}

	//Return the initial to power of e in hPow;
	if (hPow != NULL) {
		*hPow = hI;
	}


}


/*
 * Addition
 */

/**
 * Add two recursive polynomials given their head nodes, a and b.
 * nvar: number of variables in the polynomial.
 * returns a pointer to the head Node of the sum.
 */
RecNode_t* addRecPolynomials(RecNode_t *a, RecNode_t *b, int nvar) {
	RecNode_t* c = NULL;
	RecNode_t* trailC = c;
	RecNode_t* curA = a;
	RecNode_t* curB = b;

	Node* ccoef = (Node*) malloc(sizeof(Node*));
	degree_t exp;
	while (curA != NULL && curB != NULL) {

		if (curA->exp < curB->exp) {
			// a < b:
			exp = curB->exp;
			trailC = addRecTerm(trailC, deepCopyPolynomial(curB->coef, nvar), exp);
			curB = curB->next;
		} else if (curA->exp == curB->exp) {
			// a == b:
			ccoef = addPolynomials(curA->coef, curB->coef, nvar);
			if (ccoef != NULL) {
				exp = curB->exp;
				trailC = addRecTerm(trailC, ccoef, exp);
			}
			curA = curA->next;
			curB = curB->next;
		} else {
			exp = curA->exp;
			trailC = addRecTerm(trailC, deepCopyPolynomial(curA->coef, nvar),
					exp);
			curA = curA->next;
		}
		if (c == NULL) {
			c = trailC;
		}
	}
	if (curA != NULL) {
		RecNode_t* tail = deepCopyRecPolynomial(curA, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
	}
	if (curB != NULL) {
		RecNode_t* tail = deepCopyRecPolynomial(curB, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
	}
	return c;
}

void multiDivisorPseudoDivision(RecNode_t* f, RecNode_t** T, RecNode_t** quoSet,
		RecNode_t** rem, Node** H, int nvar, int nSet, int lazy) {
	if (nSet < 1) {
		fprintf(stderr, "nSet should be a positive integer!");
		return;
	}
	int e = 1;
	// base case is when nSet = 1:
	if (nSet == 1) {
		pesudoDivide(f, T[0], &quoSet[0], rem, &e, H, nvar, lazy);
		return;
	}

	// nSet > 1:

	// // // // // // // // // // // // // // // // // //

	Node* totalH = NULL; // At the end, *H = totalH;
	Node* loopH = NULL; // loopH is the init generated in first loop
	RecNode_t* rIn = deepCopyRecPolynomial(f, nvar);
	degree_t loopSize = (rIn->exp) + 1;
	RecNode_t** r_i = (RecNode_t**) calloc(loopSize, sizeof(RecNode_t*));
	Node** h_i = (Node**) calloc(loopSize, sizeof(Node*));
	RecNode_t** loopQ[loopSize];
	for (unsigned long i = 0; i < loopSize; i++)
		loopQ[i] = (RecNode_t**) calloc(nSet - 1, sizeof(RecNode_t*));

	// TODO: delete this tmpT (it is not necessary!):
	RecNode_t* tmpT[nSet - 1]; //copy the first nSet-1 member of the triangular set in tmpT for using in loops
	for (int i = 0; i < nSet - 1; ++i)
		tmpT[i] = T[i];

	RecNode_t* rCur = rIn;
	for (degree_t i = 0; i < loopSize; ++i) {
		if (rCur->exp == loopSize - i - 1) {
			multiDivisorPseudoDivision(
					convertToRecursiveNodeAtIdx(rCur->coef, nvar - nSet + 1),
					tmpT, loopQ[i], &r_i[i], &h_i[i], nvar, nSet - 1, lazy);
			rCur = rCur->next;
		} else {
			degrees_t hDegs = (degrees_t) calloc(nvar, sizeof(degree_t));
			ratNum_t hCoef;
			mpq_init(hCoef);
			mpq_set_ui(hCoef, 1ul, 1ul);
			h_i[i] = addPolynomials_inp(h_i[i], addTerm(NULL, hDegs, hCoef),
					nvar);
		}
		if (i == 0)
			loopH = addPolynomials_inp(loopH, h_i[i], nvar);
		else
			loopH = multiplyPolynomials(loopH, h_i[i], nvar);
		//rCur = rCur->next;
	}

	RecNode_t* rHead = NULL; // head of new remainder s.t. rHead = \sum\limits_{i} loopH/h_i[i] * r_i[i] * x^i;

	for (degree_t i = 0; i < loopSize; ++i) {
		// compute loopH/h_i[i] :
		Node* quo_H_div_hi = NULL;
		Node* rem_H_div_hi = NULL;
		dividePolynomials(loopH, h_i[i], &quo_H_div_hi, &rem_H_div_hi, nvar);
		freePolynomial(rem_H_div_hi); // don't need this part for the rest of algorithm!

		//make new remainder:
		rHead = addRecPolynomials(rHead,
				addRecTerm(NULL,
						multiplyPolynomials(quo_H_div_hi,
								convertFromRecursiveNodeAtIdx(r_i[i],
										nvar - nSet + 1), nvar),
						loopSize - i - 1), nvar);

		// update quoSet:
		for (int j = 0; j < nSet - 1; ++j) {
			quoSet[j] = addRecPolynomials(quoSet[j],
					addRecTerm(NULL,
							multiplyPolynomials(quo_H_div_hi,
									convertFromRecursiveNodeAtIdx(loopQ[i][j],
											nvar - nSet + 1), nvar),
							loopSize - i - 1), nvar); // TODO: Does it need addPolynomials in Future?!
		}
	}

	// // // // // // // // // // // // // // // // // //

	RecNode_t* newR = NULL; // it is the new remainder generated for the last part of algorithm (second loop)
	Node* hPow = NULL;
	RecNode_t* tmp_q = NULL;

	pesudoDivide(rHead, T[nSet - 1], &tmp_q, &newR, &e, &hPow, nvar, lazy);
	totalH = multiplyPolynomials(loopH, hPow, nvar); //update totalH
	for (int j = 0; j < nSet; ++j) {
		RecNode_t* curQj = quoSet[j];
		while (curQj != NULL) {
			curQj->coef = multiplyPolynomials(curQj->coef, hPow, nvar);
			curQj = curQj->next;
		}
	}

	quoSet[nSet - 1] = addRecPolynomials(quoSet[nSet - 1], tmp_q, nvar);

	// // // // // // // // // // // // // // // // // //

	loopSize = (newR->exp) + 1; //update the size w.r.t newR (remainder)
	RecNode_t** rr_i = (RecNode_t**) calloc(loopSize, sizeof(RecNode_t*));
	Node* loopHH = NULL;
	Node** hh_i = (Node**) calloc(loopSize, sizeof(Node*));
	RecNode_t** loopQQ[loopSize]; //loopQ[loopSize][nSet-1]
	for (unsigned long i = 0; i < loopSize; i++)
		loopQQ[i] = (RecNode_t**) calloc(nSet - 1, sizeof(RecNode_t*));

	RecNode_t* rrCur = newR;

	for (degree_t i = 0; i < loopSize; ++i) {
		if (rrCur->exp == loopSize - i - 1) {
			multiDivisorPseudoDivision(
					convertToRecursiveNodeAtIdx(rrCur->coef, nvar - nSet + 1),
					tmpT, loopQQ[i], &rr_i[i], &hh_i[i], nvar, nSet - 1, lazy);
			rrCur = rrCur->next;
		} else {
			degrees_t hhDegs = (degrees_t) calloc(nvar, sizeof(degree_t));
			ratNum_t hhCoef;
			mpq_init(hhCoef);
			mpq_set_ui(hhCoef, 1ul, 1ul);
			hh_i[i] = addPolynomials_inp(hh_i[i], addTerm(NULL, hhDegs, hhCoef),
					nvar);
		}

		if (i == 0)
			loopHH = addPolynomials(loopHH, hh_i[i], nvar);
		else
			loopHH = multiplyPolynomials(loopHH, hh_i[i], nvar);
	}
	// update quoSet[i] s.t.  quoSet[i] = loopHH * quoSet[i]
	for (int i = 0; i < nSet - 1; ++i) {
		RecNode_t* curQi = quoSet[i];
		while (curQi != NULL) {
			curQi->coef = multiplyPolynomials(curQi->coef, loopHH, nvar);
			curQi = curQi->next;
		}
	}

	RecNode_t* rrHead = NULL;
	for (degree_t i = 0; i < loopSize; ++i) {
		Node* quo_H_div_hhi = NULL;
		Node* rem_H_div_hhi = NULL;
		dividePolynomials(loopHH, hh_i[i], &quo_H_div_hhi, &rem_H_div_hhi,
				nvar);
		freePolynomial(rem_H_div_hhi);
		//make new remainder
		rrHead = addRecPolynomials(rrHead,
				addRecTerm(NULL,
						multiplyPolynomials(quo_H_div_hhi,
								convertFromRecursiveNodeAtIdx(rr_i[i],
										nvar - nSet + 1), nvar),
						loopSize - i - 1), nvar);
		for (int j = 0; j < nSet - 1; ++j) {
			quoSet[j] = addRecPolynomials(quoSet[j],
					addRecTerm(NULL,
							multiplyPolynomials(quo_H_div_hhi,
									convertFromRecursiveNodeAtIdx(loopQQ[i][j],
											nvar - nSet + 1), nvar),
							loopSize - i - 1), nvar);
		}
	}
	totalH = multiplyPolynomials_inp(totalH, loopHH, nvar);
	*H = deepCopyPolynomial(totalH, nvar);
	*rem = deepCopyRecPolynomial(rrHead, nvar);
	// // // // // // // // // // // // // // // // // //
	return;
}
