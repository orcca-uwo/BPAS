
#include "RationalNumberPolynomial/SMQP_Support.h"


/*****************
 * Exponent Vector functions.
 *****************/

/**
 * Compare two exponent vectors for lexicographical term order.
 * a: int array representing exponent vector
 * b: int array representing exponent vector
 * v: size of the arrays
 * returns -1 if a < b, 0 if equal, 1 if a > b
 */
int compareExponentVectors(degrees_t a, degrees_t b, int v) {

#if SMQP_COUNT_COMPARISONS
	++smqpCompareCount;
#endif
	degree_t diff = a[0] - b[0];
	int i = 1;
	while (diff == 0 && i < v) {
		diff = a[i] - b[i];
		++i;
	}
	return diff;
	// for (int i = 0; i < v; ++i) {
	// 	if (a[i] < b[i])
	// 	    return -1;
	// 	else if (a[i] > b[i])
	// 	    return 1;
	// }
	// return 0;
}

/**
 * Check if an exponent vector has all 0 exponents.
 * a: the exponent vector to check
 * v: size of a
 */ 
int isZeroExponentVector(degrees_t a, int v) {
	if (a == NULL) {
		return 1;
	}
	for (int i = 0; i < v; ++i) {
		if (a[i] > 0) {
			return 0;
		}
	}
	return 1;
}

#if SMQP_SUPPORT_DEBUG
/**
 * Print the contents of an exponent vector in a comma-separated list
 * d: the exponent vector
 * nvar: the number of elements in the vector
 */
void printDegs(degrees_t d, int nvar) {
	for(int i = 0; i < nvar-1; ++i) {
		fprintf(stderr, "%d, ", d[i]);
	}
	fprintf(stderr, "%d", d[nvar-1]);
}
#endif




/*****************
 * Node defintion and helper functions.
 *****************/

/**
 * Free a node, and all it's sucessors. 
 */
void freePolynomial(Node* node) {
	Node* nextNode;
	while (node != NULL) {
		nextNode = node->next;
		freeNode(node);
		node = nextNode;
	}
}

/**
 * Deep copy the polynomail whose head in the input Node*.
 */
Node* deepCopyPolynomial(Node* node, int nvar) {
	Node* head = NULL;
	Node* tail = NULL;
	
	degrees_t degs;
	degrees_t nodeDegs;
	while (node != NULL) {
		nodeDegs = node->degs;
		degs = (degrees_t) malloc(sizeof(degree_t) * nvar);
		memcpy(degs, nodeDegs, sizeof(degree_t) * nvar);
		tail = addTerm(tail, degs, node->coef);
		if (head == NULL) {
			head = tail;
		}
		node = node->next;
	}

	return head;
}

/**
 * Determine the number of terms in a polynomial given the head node 
 * of the linked-list representation of the polynomial.
 * returns the number of terms including and following the supplied node 
 */
 polysize_t numberOfTermsNode(Node* node) {	
	polysize_t val = 0;
	Node *h = node;
	while(h != NULL){
		h = h->next;
		val++;
	}				
	return val;
}


/**
 * Given a linked-list of nodes, whose head is poly, 
 * sort the list based on the ordering of compareExponentVectors.
 * Updates the polynomial in place by modifying *poly to be the new head node.
 * Also returns a pointer which is the new head node. 
 */
Node* sortPolynomial(Node** poly, int nvar) {
	//TODO not insertion sort.
	Node* sorted = NULL;

	Node* node = *poly;
	Node* next;
	while(node != NULL) {
		next = node->next;

		if (sorted == NULL || compareExponentVectors(node->degs, sorted->degs, nvar) > 0) {
			//if inserted node is greated than current head
			node->next = sorted;
			sorted = node;
		} else {
			//else traverse and find insertion point
			Node* cur = sorted;
			while (cur->next != NULL && compareExponentVectors(cur->next->degs, node->degs, nvar) > 0) {
				cur = cur->next;
			}
			node->next = cur->next;
			cur->next = node;
		}

		node = next;
	}

	condensePolyomial(sorted, nvar);

	*poly = sorted;
	return sorted;
}

/**
 * Given a polynomial in sorted order but with possible duplicates, 
 * condense the list such that like-terms are combined into a single node.
 */
void condensePolyomial(Node* poly, int nvar) {
	Node* prev = poly;
	if (prev == NULL ) {
		return;
	}

	Node* node = prev->next;
	while (node != NULL) {
		if (compareExponentVectors(prev->degs, node->degs, nvar) == 0) {
			mpq_add(prev->coef, prev->coef, node->coef);
			prev->next = node->next;
			node->next = NULL;
			freeNode(node);
			node = prev->next;
		} else {
			prev = node;
			node = node->next;
		}
	}
}

/**
 * Add a term to the end of the polynomial linked list given the tail node,
 * trailingTerm, the exponent vector for the new term, d, and the coefficient, 
 * ratNum_t. A copy of the input coef is made and stored in the node. But the
 * degrees_t is stored directly. 
 * 
 * Note that the input coef can be null. This signifies to create a Node with a 
 * 0 coef.
 *
 * returns a pointer to the new tail of the list (the node created)
 */
Node* addTerm(Node* trailingTerm, degrees_t d, const ratNum_t coef){
	Node* temp = (Node*) malloc(sizeof(Node));

	temp->degs = d;
	mpq_init(temp->coef);
	mpq_set(temp->coef, coef);
	// temp->coef = coef;
	temp->next = NULL;
	
	if (trailingTerm != NULL) {
		trailingTerm->next = temp;
	}
	return temp;
}

/**
 * Add a term to the end of linked list with tail trailingTerm. The new term
 * will have degs d and a coef of 0.
 * 
 * returns a pointer to the new tail of the list (the node created)
 */
Node* addZeroTerm(Node* trailingTerm, degrees_t d) {
	mpq_t mpqZero;
	mpq_init(mpqZero);
	Node* ret = addTerm(trailingTerm, d, mpqZero);
	mpq_clear(mpqZero);
	return ret;
}

/** 
 * Given two terms, as nodes a and b, multiply these two and return the
 * single result as a node.
 */
Node* multiplyTerms(Node* a, Node* b, int nvar) {
	if (a == NULL || b == NULL) {
		return NULL;
	}

	ratNum_t coef;
	mpq_init(coef);
	mpq_mul(coef, a->coef, b->coef);
	// ratNum_t coef = a->coef * b->coef;
	degrees_t degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
	addExponentVectors(a->degs, b->degs, degs, nvar);
	Node* c = addTerm(NULL, degs, coef);
	mpq_clear(coef);
	return c;
}

/**
 * Negate a polynomial. 
 * Given the head of a polynomial, negate all coefficients in place.
 */
void negatePolynomial(Node* a) {
	Node* node = a;
	while (node != NULL) {
		mpq_neg(node->coef, node->coef);
		node = node->next;
	}
}

void evalPolyToVal(Node* a, mpq_t* vals, int nvar, mpq_t res) {
	if (a == NULL) {
		mpq_set_ui(res, 0ul, 1ul);
		return;
	}

	//TODO multivariate horner's method should replace the below
	//TODO there is a 1/10000 (ish) bug in the below.. 

	int i;
	mpq_set_ui(res, 0ul, 1ul);
	Node* current = a;
	mpq_t acc[nvar];
	for (i = 0; i < nvar; ++i) {
		mpq_init(acc[i]);
	}
	degrees_t degs;
	while (current != NULL) {
		degs = current->degs;
		//exponentiate each vals[i] to degs[i]
		for (i = 0; i < nvar; ++i) {
			if (degs[i] == 0) {
				continue;
			}
			mpz_pow_ui(mpq_numref(acc[i]), mpq_numref(vals[i]), degs[i]);
			mpz_pow_ui(mpq_denref(acc[i]), mpq_denref(vals[i]), degs[i]);
			mpq_canonicalize(acc[i]);
		}

		//accumulate and reset all variable values into acc[0]
		if (mpq_sgn(acc[0]) == 0) {
			mpq_set_ui(acc[0], 1ul, 1ul);
		}
		for (i = 1; i < nvar; ++i) {
			if (degs[i] == 0) {
				continue;
			}
			mpq_mul(acc[0], acc[0], acc[i]);
			mpq_set_ui(acc[i], 0ul, 1ul);
		}
		mpq_mul(acc[0], current->coef, acc[0]);
		mpq_add(res, res, acc[0]);
		mpq_set_ui(acc[0], 0ul, 1ul);

		current = current->next;
	}

	for (i = 0; i < nvar; ++i) {
		mpq_clear(acc[i]);
	}
}

Node* evaluatePoly(Node* a, int* active, mpq_t* vals, int nvar) {
	if (a == NULL) {
		return NULL;
	}

	int i;
	int newNvar = nvar;
	for (i = 0; i < nvar; ++i) {
		if (active[i]) {
			--newNvar;
		}
	}

	if (newNvar == 0) {
		mpq_t res;
		mpq_init(res);
		evalPolyToVal(a, vals, nvar, res);
		degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
		Node* ret = addTerm(NULL, degs, res);
		mpq_clear(res);
		return ret;
	} 

	Node* newPoly = NULL;
	Node* newTail = NULL;
	Node* current = a;
	degrees_t curDegs;
	degrees_t newDegs;
	mpq_t acc[nvar];
	for (i = 0; i < nvar; ++i) {
		mpq_init(acc[i]);
	}
	mpq_t newCoef;
	mpq_init(newCoef);
	int j;

	while(current != NULL) {
		j = 0;
		curDegs = current->degs;
		newDegs = (degrees_t) calloc(newNvar, sizeof(degrees_t));
		for (i = 0; i < nvar; ++i) {
			if (curDegs[i] == 0) {
				continue;
			}
			if (!active[i]) {
				newDegs[j] = curDegs[i];
				++j;
				continue;
			}

			mpz_pow_ui(mpq_numref(acc[i]), mpq_numref(vals[i]), curDegs[i]);
			mpz_pow_ui(mpq_denref(acc[i]), mpq_denref(vals[i]), curDegs[i]);
			mpq_canonicalize(acc[i]);
		}

		if (mpq_sgn(acc[0]) == 0) {
			mpq_set_ui(acc[0], 1ul, 1ul);
		}
		for (i = 1; i < nvar; ++i) {
			if (active[i] && curDegs[i] != 0) {
				mpq_mul(acc[0], acc[0], acc[i]);
				mpq_set_ui(acc[i], 0ul, 1ul);
			}
		}
		mpq_mul(newCoef, acc[0], current->coef); 
		mpq_set_ui(acc[0], 0ul, 1ul);

		newTail = addTerm(newTail, newDegs, newCoef);
		if (newPoly == NULL) {
			newPoly = newTail;
		}

		current = current->next;
	}


	for (i = 0; i < nvar; ++i) {
		mpq_clear(acc[i]);
	}
	mpq_clear(newCoef);

	condensePolyomial(newPoly, nvar);

	return newPoly;
}

Node* leadingTerm (Node* a){
  Node* n = addTerm(NULL, a->degs, a->coef);
  return n;
}


/*****************
 * SMQP Addition & Subraction
 *****************/

/**
 * Add two polynomials given their head nodes, a and b.
 * nvar: number of variables in the polynomials.
 * returns a pointer to the head Node of the sum.
 * Note that Nodes and exponent vectors are reused as much as possible.
 */
Node* addPolynomials(Node *a, Node *b, int nvar) {
	Node* c = NULL;
	Node* trailC = c;
	Node* curA = a;
	Node* curB = b;

	ratNum_t ccoef;
	mpq_init(ccoef);
	degrees_t degs;
	while (curA != NULL && curB != NULL) {
		int cmp = compareExponentVectors(curA->degs, curB->degs, nvar);
		if (cmp < 0) {
			//a < b
			degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			memcpy(degs, curB->degs, sizeof(degree_t)*nvar);
			trailC = addTerm (trailC, degs, curB->coef);
			curB = curB->next;
		} else if (cmp == 0) {
			// a==b
			mpq_add(ccoef, curA->coef, curB->coef);
			// ccoef = curA->coef + curB->coef;
			if (mpq_sgn(ccoef) != 0) {
				degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
				memcpy(degs, curA->degs, sizeof(degree_t)*nvar);
				trailC = addTerm (trailC, degs, ccoef);
			}
			curA = curA->next;
			curB = curB->next;
		} else {
			//a > b
			degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			memcpy(degs, curA->degs, sizeof(degree_t)*nvar);
			trailC = addTerm (trailC, degs, curA->coef);
			curA = curA->next;
		}
		if (c == NULL) {
			c = trailC;
		}
	}

	if (curA != NULL) {
		Node* tail = deepCopyPolynomial(curA, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
		// trailC = addTerm(trailC, curA->degs, curA->coef);
		// if (c == NULL) {
		// 	c = trailC;
		// }
		// curA = curA->next;
	}
	if (curB != NULL) {
		Node* tail = deepCopyPolynomial(curB, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
		// trailC = addTerm(trailC, curB->degs, curB->coef);
		// if (c == NULL) {
		// 	c = trailC;
		// }
		// curB = curB->next;
	}

	mpq_clear(ccoef);
	return c;
}

Node* addPolynomials_inp(Node* a, Node* b, int nvar) {
	Node* c = NULL;
	Node* trailC = c;
	Node* curA = a;
	Node* curB = b;

	while (curA != NULL && curB != NULL) {
		int cmp = compareExponentVectors(curA->degs, curB->degs, nvar);
		if (cmp < 0) {
			//a < b
			if (trailC != NULL) {
				trailC->next = deepCopyNode(curB, nvar);
				trailC = trailC->next;
			} else {
				trailC = deepCopyNode(curB, nvar);
			}
			mpq_set(trailC->coef, trailC->coef); 
			curB = curB->next;
		} else if (cmp == 0) {
			// a==b
			mpq_add(curA->coef, curA->coef, curB->coef);
			if (mpq_sgn(curA->coef) == 0) {
				Node* del = curA;
				curA = curA->next;
				freeNode(del);
			} else {
				if(trailC != NULL) {
					trailC->next = curA;
					trailC = trailC->next;
				} else {
					trailC = curA;
				}
				curA = curA->next;
			}
			curB = curB->next;
		} else {
			if (trailC == NULL) {
				trailC = curA;
			} else {
				trailC->next = curA;
				trailC = curA;
			}
			curA = curA->next;
		}

		if (c == NULL) {
			c = trailC;
		}
	}

	if (curA != NULL) {
		if (trailC != NULL) {
			trailC->next = curA;
		} else {
			c = curA;
		}
	}
	if (curB != NULL) {
		Node* tail = deepCopyPolynomial(curB, nvar);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
	}

	return c;
}

Node* subPolynomials_inp(Node* a, Node* b, int nvar) {
	Node* c = NULL;
	Node* trailC = c;
	Node* curA = a;
	Node* curB = b;

	while (curA != NULL && curB != NULL) {
		int cmp = compareExponentVectors(curA->degs, curB->degs, nvar);
		if (cmp < 0) {
			//a < b
			if (trailC != NULL) {
				trailC->next = deepCopyNode(curB, nvar);
				trailC = trailC->next;
			} else {
				trailC = deepCopyNode(curB, nvar);
			}
			//negate since from b
			mpq_neg(trailC->coef, trailC->coef); 
			
			curB = curB->next;
		} else if (cmp == 0) {
			// a==b
			mpq_sub(curA->coef, curA->coef, curB->coef);
			if (mpq_sgn(curA->coef) == 0) {
				Node* del = curA;
				curA = curA->next;
				freeNode(del);
			} else {
				if(trailC != NULL) {
					trailC->next = curA;
					trailC = trailC->next;
				} else {
					trailC = curA;
				}
				curA = curA->next;
			}
			curB = curB->next;
		} else {
			if (trailC == NULL) {
				trailC = curA;
			} else {
				trailC->next = curA;
				trailC = curA;
			}
			curA = curA->next;
		}

		if (c == NULL) {
			c = trailC;
		}
	}

	if (curA != NULL) {
		if (trailC != NULL) {
			trailC->next = curA;
		} else {
			c = curA;
		}
	}
	if (curB != NULL) {
		Node* tail = deepCopyPolynomial(curB, nvar);
		negatePolynomial(tail);
		if (trailC != NULL) {
			trailC->next = tail;
		} else {
			c = tail;
		}
	}

	return c;
}



/*****************
 * SMQP Multiplication & Helpers
 *****************/

#if SMQP_SUPPORT_DEBUG
/**
 * Print the product degrees_t currently in the heap.
 */
void prodheapPrint(ProductHeap* h, int nvar) {
	for (int i = 0; i < h->heapSize; ++i){
		fprintf(stderr, "( %lx,", (unsigned long int) h->elements[i]);
		printDegs(h->elements[i]->product->degs, nvar);
#if SMQP_SUPPORT_CHAINED_HEAP
		productHeapElem* next = h->elements[i]->next;
		while (next != NULL) {
			fprintf(stderr, "; ");
			fprintf(stderr, "%lx,", (unsigned long int) next);
			printDegs(next->product->degs, nvar);
			next = next->next;
		}
#endif
		fprintf(stderr, ")");
		fprintf(stderr," ");
	}
	fprintf(stderr, "\n");
}
#endif

/**
 * Make an element for the product heap, combining nodes a and b as the 
 * element's product.
 */
productHeapElem* prodheapMakeElement(ProductHeap* h, Node* a, Node* b) {
	productHeapElem* elem = (productHeapElem*) malloc(sizeof(productHeapElem));
	elem->a_i = a;
	elem->b = b;
	elem->product = multiplyTerms(a, b, h->nvar);
#if SMQP_SUPPORT_CHAINED_HEAP
	elem->next = NULL;
#endif
	return elem;
}

/**
 * Cleanup memory initialized for the heap
 */
void prodheapFree(ProductHeap* h) {
	productHeapElem** elems = h->elements;
	for (int i = 0; i < h->heapSize; ++i) {
		prodheapFreeElement(elems[i]);
	}
	free(h->elements);
	free(h); 
}

/**
 * Create an empty product heap. 
 */
ProductHeap* prodheapCreate(int nvar) {
	ProductHeap* h = (ProductHeap*) malloc(sizeof(ProductHeap));
	h->elements = NULL;
	h->heapSize = 0;
	h->maxHeapSize = 0;
	h->nvar = nvar;
	return h;
}

/**
 * Initialize the product heap with the two polynomials to multiply, a and b.
 *
 * We know the maximum heap size is numTerms(a) as at most one of a_i*b
 * is in the heap at once. 
 */
ProductHeap* prodheapInit(Node* a, Node* b, int nvar) {
	ProductHeap* h = prodheapCreate(nvar);

#if SMQP_SUPPORT_CHAINED_HEAP 
	polysize_t num_a = numberOfTermsNode(a);
	productHeapElem** elems = (productHeapElem**) malloc(sizeof(productHeapElem*) * num_a);
	elems[0] = prodheapMakeElement(h, a, b);

	int heapCurSize = 1;
	int heapMaxSize = num_a;
#else 
	Node* curA = a;
	degrees_t d;
	degrees_t bdegs = b->degs;
	ratNum_t prodCoef;
	mpq_init(prodCoef);
	//Since terms of a are ordered in decreasing order, elems[] is already
	//heapified based on the product term.
	int heapCurSize = 0;
	polysize_t num_a = numberOfTermsNode(a);
	int heapMaxSize = num_a;
	productHeapElem** elems = (productHeapElem**) malloc(sizeof(productHeapElem*) * (num_a));
	for (heapCurSize = 0; heapCurSize < num_a; ++heapCurSize) {
		elems[heapCurSize] = (productHeapElem*) malloc(sizeof(productHeapElem));
		elems[heapCurSize]->a_i = curA;
		elems[heapCurSize]->b = b;
		d = (degrees_t) malloc(sizeof(degree_t)*nvar);
		addExponentVectors(curA->degs, bdegs, d, nvar);
		mpq_mul(prodCoef, curA->coef, b->coef);
		elems[heapCurSize]->product = addTerm(NULL, d, prodCoef);
		curA = curA->next;
	}
	mpq_clear(prodCoef);
#endif

	h->elements = elems;
	h->heapSize = heapCurSize;
	h->maxHeapSize = heapMaxSize;
	return h;
}

/**
 * Given an index, i, into the heap, swim that node up so that the heap 
 * is properly ordered. 
 */
void prodheapSwim(ProductHeap* h, int index) {
	productHeapElem** elems = h->elements;
	productHeapElem* temp;
	int nvar = h->nvar;
	int cmp;
	int i = index;
	int j;
	while (i > 0) {
		j = (i-1) >> 1;
		cmp = compareExponentVectors(elems[j]->product->degs, elems[i]->product->degs, nvar);
		if (cmp >= 0) {
			break;
		}
		temp = elems[j];
		elems[j] = elems[i];
		elems[i] = temp;
		i = j;
	}
}

/**
 * Given an index, i, into the heap, sink that node down so that the heap 
 * is properly ordered. 
 * Note: This is not used in the chained heap case
 */
void prodheapSink(ProductHeap* h, int i) {
	productHeapElem** elems = h->elements;
	productHeapElem* temp;
	int nvar = h->nvar;
	int N = h->heapSize;

	int j = (i << 1) + 1;
	int cmp;

	while(j < N) {
		if (j < N - 1 && compareExponentVectors(elems[j]->product->degs, elems[j+1]->product->degs, nvar) < 0) {
			++j;
		}

		if (compareExponentVectors(elems[i]->product->degs, elems[j]->product->degs, nvar) >= 0) {
			break;
		}

		temp = elems[j];
		elems[j] = elems[i];
		elems[i] = temp;
		i = j;
		j = (i << 1) + 1;
	}
}

/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
void prodheapInsert(ProductHeap* h, productHeapElem* elem) {

	int nvar = h->nvar;
	int s = h->heapSize;
	int N = h->maxHeapSize;
	productHeapElem** elems = h->elements;
#if SMQP_COUNT_CHAINS
	SMQP_INSERTS++;
#endif

	if (s == 0) {
		if (s >= N) {
			N = (s == 0) ? 1 : s * 2;
			h->elements = (productHeapElem**) realloc(h->elements, N*sizeof(productHeapElem*));
			elems = h->elements;
			h->maxHeapSize = N;
		}
		elems[0] = elem;
		h->heapSize = 1;
#if SMQP_COUNT_CHAINS
		SMQP_CHAINS++;
#endif
		return;		
	}

#if SMQP_SUPPORT_CHAINED_HEAP
	
	degrees_t prodDegs = elem->product->degs;

	//first check if we can chain off the root
	if (compareExponentVectors(elems[0]->product->degs, prodDegs, nvar) == 0) {
		elem->next = elems[0];
		elems[0] = elem;
#if SMQP_COUNT_CHAINS
		SMQP_CHAINS++;
#endif
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

	int i = (s-1)/2; //i is parent 
	int j = s;       //j is current insertion point
	long long unsigned int path = 1;
	int cmp;
	while (j > 0) {
		cmp = compareExponentVectors(elems[i]->product->degs, prodDegs, nvar);
		if (cmp == 0) {
#if SMQP_COUNT_CHAINS
			SMQP_CHAINS++;
#endif
			elem->next = elems[i];
			elems[i] = elem;
			return;
		} else if (cmp < 0) {
			path <<= 1;
			if (!(j & 1)) {
				//set the trailing bit to 1 to distinguish left/right of path
				path += 1; 
			}
			j = i;
			i = (i-1) / 2;
		} else { //cmp > 0
			break;
		}
	}

	//if we made it this far, then we are adding a new element to the heap. 
	//resize if necessary
	if (s >= N) {
		N = (s == 0) ? 1 : s * 2;
		h->elements = (productHeapElem**) realloc(h->elements, N*sizeof(productHeapElem*));
		elems = h->elements;
		h->maxHeapSize = N;
	}

	//then j is now the place we need to insert elem;
	//do so, and then push all others down the path, inserting the last
	//as the new element in elems[s];
	productHeapElem* temp;
	while (j <= s) {
		temp = elems[j];
		elems[j] = elem;
		elem = temp;
		j = 2*j + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
#else 

	if (s >= N) {
		N = (s == 0) ? 10 : s * 2;
		h->elements = (productHeapElem**) realloc(h->elements, N*sizeof(productHeapElem*));
		elems = h->elements;
		h->maxHeapSize = N;
	}
	elems[s] = elem;
	++(h->heapSize);
	prodheapSwim(h, s);
	
#endif
}

/**
 * Peek into the heap, h, to get the exponent vector of the product
 * of the max element.
 */
degrees_t prodheapPeek(ProductHeap* h) {
	if (h->heapSize > 0) {
		return h->elements[0]->product->degs;
	}
	return NULL;
}

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * returns NULL if no such element in the heap.
 */
productHeapElem* prodheapRemoveMax(ProductHeap* h) {
	if (h->heapSize == 0) {
		return NULL;
	}

	productHeapElem** elems = h->elements;
	productHeapElem* maxElem = elems[0];
	int nvar = h->nvar;
	int i = 0, j = 1;
	int s = --(h->heapSize);
	
	//promote largest children
	while (j < s) {
		if (j+1 < s && compareExponentVectors(elems[j]->product->degs, elems[j+1]->product->degs, nvar) < 0) {
			++j;
		}
		elems[i] = elems[j];
		i = j;
		j = (2*j) + 1;
	}
	//now place last element into i and swim up to make tree complete 
	j = (i-1)/2;
	while(i > 0) {
		if (compareExponentVectors(elems[s]->product->degs, elems[j]->product->degs, nvar) < 0) {
			break;
		}
		elems[i] = elems[j];
		i = j;
		j = (j-1)/2;
	}
	elems[i] = elems[s]; 
	return maxElem;
}

/**
 * Extract the maximum product from the heap, returned as the Node*. 
 * This automatically adds the next term in the stream of the chosen maximum
 * if such a term exists.
 */
Node* prodheapExtract(ProductHeap* h) {
	if (h->heapSize == 0) {
		return NULL;
	}

	productHeapElem** elems = h->elements;
	productHeapElem* maxElem = elems[0];
	Node* maxProd = maxElem->product;
	int nvar = h->nvar;
	int i,j;
	int N = h->heapSize;


	maxElem->b = (maxElem->b)->next;
	if(maxElem->b != NULL) {
		degrees_t d = (degrees_t) malloc(sizeof(degree_t)*nvar);
		addExponentVectors(maxElem->a_i->degs, maxElem->b->degs, d, nvar);
		ratNum_t coef;
		mpq_init(coef);
		mpq_mul(coef, maxElem->a_i->coef, maxElem->b->coef);
		maxElem->product = addTerm(NULL,d,coef);
		mpq_clear(coef);
		prodheapSink(h,0);
	} else {
		h->heapSize = --N;
		//there is a hole in elems[0] and we will push up the larger of the children.
		for(i = 0, j = 1; j < N;) {
			if (compareExponentVectors(elems[j]->product->degs, elems[j+1]->product->degs, nvar) < 0) {
				++j;
			}
			elems[i] = elems[j];
			i = j;
			j = 2*j + 1;
		}
		//hole is now at i, so fill it with the least element and swim it up;
		//while elems[N] > parent(i) switch the hole with its parent and compare
		//elems[N] to the hole's new parent
		j = (i-1)/2;
		while (i > 0 && compareExponentVectors(elems[N]->product->degs, elems[j]->product->degs, nvar) > 0) {
			elems[i] = elems[j];
			i = j;
			j = (j-1)/2;
		}
		elems[i] = elems[N];
	}
	return maxProd;
}

/**
 * Find the next product term of highest degree.
 * This searches through each i in I..n streams representing the expanded multiplcation
 * of node a_i with b. Picking the maximal a_i*b_fi.
 * 
 * returns the terms to multiply as a_s and b_fs and the value s.
 *
 * This search is done linearly, for use in the multiplication algorithm which does
 * not use heaps.
 */
void findNextProductTerm(Node* a_I, Node** b_fi, polysize_t numTermsA, polysize_t I, int* s, Node** a_s, Node** b_fs, int nvar) {
	degrees_t maxDegs = NULL; //(degrees_t) calloc(nvar, sizeof(degree_t));
	Node* curA = a_I;
	Node* maxA = curA;
	Node* maxB = b_fi[0];
	degrees_t curDegs = (degrees_t) malloc(sizeof(degree_t)*nvar);

	for (int i = I; i < numTermsA; ++i) {
		addExponentVectors(curA->degs, (b_fi[i])->degs, curDegs, nvar);
		if (maxDegs ==  NULL || 
				compareExponentVectors(maxDegs, curDegs, nvar) < 0) {
			free(maxDegs);
			maxDegs = curDegs;
			curDegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			*s = i;
			maxA = curA;
			maxB = b_fi[i];
		}
		curA = curA->next;
	}

	free(curDegs);
	free(maxDegs);
	*a_s = maxA;
	*b_fs = maxB;
}

/**
 * Multiply two polynomials given their head Nodes, a and b.
 * This algorithm makes use of heaps as an efficient search data structure. 
 * It is assumed that both a and b have compatible exponent vectors.
 * 
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the head Node of the product polynomial.
 */
Node* multiplyPolynomials(Node* a, Node* b, int nvar) {
	if (a == NULL || b == NULL) {
		return NULL;
	}

	ProductHeap* h = prodheapInit(a,b,nvar);
	degrees_t cdegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
	ratNum_t ccoef;
	mpq_init(ccoef);

#if SMQP_SUPPORT_CHAINED_HEAP
	Node* c = NULL;
	Node* prevC = NULL;
	Node* trailC = NULL;
#else
	Node* c = (Node*) malloc(sizeof(Node));
	Node* prevC = c;
	Node* trailC = c;
	mpq_init(c->coef);// coef = 0;
	addExponentVectors(a->degs, b->degs, cdegs, nvar);
	c->degs = cdegs;
	// c->coef = ccoef; //it is accumulated in the while loop
#endif

//Note that in this chained heap version of the algorithm the heap
//resuses its product nodes as elements get inserted back into the heap
//and we create additional nodes as needed for the result.
//This differs from the regular heap implementation (below) where the product nodes
//extracted from the heap are used in the result.  
#if SMQP_SUPPORT_CHAINED_HEAP
	productHeapElem* maxElem = NULL;
	productHeapElem* nextMaxElem = NULL;
	degrees_t nextDegs;
	while ( (nextDegs = prodheapPeek(h)) != NULL) {
		memcpy(cdegs, nextDegs, sizeof(degree_t)*nvar);
		
		while (nextDegs != NULL && compareExponentVectors(cdegs, nextDegs, nvar) == 0) {
			//we will extract and accumulate the coefficents 
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			productHeapElem* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax(h);
			while (maxElem != NULL) {
				mpq_add(ccoef, ccoef, maxElem->product->coef);

				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == b && maxElem->a_i->next != NULL) {
					productHeapElem* nextAHeap = prodheapMakeElement(h, maxElem->a_i->next, b);
					//push it to the head of the oldMaxElem
					nextAHeap->next = oldMaxElem;
					oldMaxElem = nextAHeap;
				}

				//cache next before freeing or overwriting 
				nextMaxElem = maxElem->next;

				//If the extracted term has another product in the stream, 
				//update the product and push onto the oldMaxElem chain
				maxElem->b = (maxElem->b)->next;
				if(maxElem->b != NULL) {
					// ratNum_t tempCoef;
					// mpq_init(tempCoef);
					// degrees_t d = (degrees_t) malloc(sizeof(degree_t)*nvar);
					addExponentVectors(maxElem->a_i->degs, maxElem->b->degs, maxElem->product->degs, nvar);
					mpq_mul(maxElem->product->coef, maxElem->a_i->coef, maxElem->b->coef);
					// mpq_clear(tempCoef);
					maxElem->next = oldMaxElem;
					oldMaxElem = maxElem;
				} else {
					//we are done with the maxElem productHeapElem
					prodheapFreeElement(maxElem);
				}

				maxElem = nextMaxElem;
			}

			//reset head of maxElem list
			maxElem = oldMaxElem;	

			nextDegs = prodheapPeek(h);		
		}

		//Commit new term to the product.
		if (mpq_sgn(ccoef) != 0) {
			trailC = addTerm(trailC, cdegs, ccoef);
			if (c == NULL) {
				c = trailC;
			}

			//reset accumulator variables
			cdegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			mpq_clear(ccoef);
			mpq_init(ccoef);
		}
		
		//Insert all successors of previously extracted products
		while(maxElem != NULL) {
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;
			prodheapInsert(h,maxElem);
			maxElem = nextMaxElem;

		}
	}
	free(cdegs);

#else

	Node* prod;
	degrees_t nextDegs;
	while ( (prod = prodheapExtract(h)) != NULL ) {
		nextDegs = prod->degs;
		if( compareExponentVectors(cdegs, nextDegs, nvar) != 0) {
			if (mpq_sgn(ccoef) == 0) {
				freeNode(trailC);
				prevC->next = prod;
			} else {
				//set the prod node's coef to the value of our accumulator ccoef
				mpq_set(trailC->coef,ccoef);
				trailC->next = prod;
				prevC = trailC;
			}
			trailC = prod;
			mpq_set(ccoef,prod->coef);
			cdegs = nextDegs;
		} else {
			mpq_add(ccoef, ccoef, prod->coef);
			// ccoef += prod->coef;
			freeNode(prod);
		}
	}

#endif

	mpq_clear(ccoef);
	prodheapFree(h);

	return c;
}

Node* multiplyPolynomials_inp(Node* a, Node* b, int nvar) {
	//TODO actually do this in-place.

	Node* prod = multiplyPolynomials(a, b, nvar);
	freePolynomial(a);
	return prod;
}

/** 
 * Initialize the multiplication producer's state given the 
 * two polynomials to multiply and the number of variables of
 * the polynomials.
 */
MultState* initProducerState(Node* a, Node* b, int nvar) {

	degrees_t cdegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
	addExponentVectors(a->degs, b->degs, cdegs, nvar);
	
	polysize_t n = numberOfTermsNode(a);
	Node** f = (Node**) malloc(sizeof(Node*)*n);
	for (int i = 0; i < n; ++i) {
		f[i] = b;
	}

	MultState* state = (MultState*) malloc(sizeof(MultState));
	state->a_I = a;
	state->f = f;
	state->I = 0;
	state->nTerms = n;
	state->nvar = nvar;
	mpq_init(state->lastCoef);
	// state->lastCoef = 0;
	state->lastDeg = cdegs;

	return state;
}

/** 
 * Cleanup the producer state memory allocation.
 */
void releaseProducerState(MultState* state) {
	free(state->f);
	if(state->lastDeg != NULL) {
		free(state->lastDeg);
	}
	mpq_clear(state->lastCoef);
	free(state);
}

/**
 * A 'producer' for the multiplication of two polynomials. 
 * This producer, using the input state, will produce a single
 * term of the product and then return. 
 *
 * Drivers must repeatedly call this to obtain terms, one after the other.
 * Terms are produced in descending order with like-terms already combined.
 * returns a NULL result if no more terms can be produced.
 *
 * Note the underlying multiplication algorithm is not efficient but this can 
 * be handy when only a few leading terms of the product are needed.
 */
MultResult multPolyProducer(MultState* state) {

	Node* a_I = state->a_I;
	Node** f = state->f;
	polysize_t numTerms = state->nTerms;
	int nvar = state->nvar;
	degrees_t lastDeg = state->lastDeg;
	polysize_t I = state->I;
	int s = I;

	Node* nextC = NULL;
	ratNum_t ccoef;
	mpq_init(ccoef);
	mpq_set(ccoef,state->lastCoef);

	Node* a_s;
	Node* b_fs; 

	while (a_I != NULL) {
		findNextProductTerm(a_I, f, numTerms, I, &s, &a_s, &b_fs, nvar);
		degrees_t nextDegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
		addExponentVectors(a_s->degs, b_fs->degs, nextDegs, nvar);
		int cmp = compareExponentVectors(lastDeg, nextDegs, nvar);
		if (cmp != 0) {
			//Then we are done accumulating into terms of exponent lastDeg
			//And create a new node of the accumulated term, if ccoef != 0
			if (mpq_sgn(ccoef) != 0) {
				nextC = addTerm(NULL, lastDeg, ccoef);
				mpq_clear(ccoef);
				mpq_init(ccoef);
				// ccoef = 0;
			} else {
				free(lastDeg);
				lastDeg = NULL;
			}
			lastDeg = nextDegs;
		} else {
			//Othweise, just reset our created nextDegs because it's not really next
			free(nextDegs);
			nextDegs = NULL;
		}

		ratNum_t prodCoef;
		mpq_init(prodCoef);
		mpq_mul(prodCoef, a_s->coef, b_fs->coef);
		mpq_add(ccoef, ccoef, prodCoef);
		mpq_clear(prodCoef);
		// ccoef += a_s->coef*b_fs->coef; //Accumulate the active term's coef
		f[s] = (f[s])->next; //Iterate the s'th stream's cursor
		if (f[s] == NULL) {
			a_I = a_I->next;
			++I;
		}				

		//If we have finished with a term, and it is non-zero, return it
		if (nextC != NULL) {
			mpq_set(state->lastCoef,ccoef);
			mpq_clear(ccoef);
			state->lastDeg = lastDeg;
			state->a_I = a_I;
			state->I = I;
			MultResult mr;
			mr.state = state;
			mr.nextRes = nextC;
			return mr;
		}
	}

	//The final term doesnt not get committed to a Node inside the while loop
	if (lastDeg != NULL) {
		nextC = addTerm(NULL, lastDeg, ccoef);
		state->lastDeg = NULL;
		mpq_clear(ccoef);
		mpq_set_ui(state->lastCoef, 0ul, 1ul);
		// state->lastCoef = 0;
		MultResult mr;
		mr.state = state;
		mr.nextRes = nextC;
		return mr;
	}

	mpq_clear(ccoef);
	MultResult mr;
	mr.state = state;
	mr.nextRes = NULL;
	return mr;
} 




/*****************
 * Polynomial exponentiation 
 *****************/

/**
 * Given a polynomial, a, compute a^n.
 * 
 * n: a positive integer
 * nvar: the number of variables of a
 * 
 */
Node* exponentiatePoly(Node* a, unsigned int n, int nvar) {

	if (n == 0) {
		degrees_t degs = (degrees_t) calloc(nvar, sizeof(degree_t));
		ratNum_t coef;
		mpq_init(coef);
		mpq_set_ui(coef, 1ul, 1ul);
		Node* ret = addTerm(NULL, degs, coef);
		mpq_clear(coef);
		return ret;
	} else if (n == 1) {
		return deepCopyPolynomial(a, nvar);
	}

	Node* r = NULL;
	Node* b = a;
	while (n > 1) {
		if (n & 1) {
			r = (r == NULL) ? b : multiplyPolynomials(r,b,nvar);
		}
		b = multiplyPolynomials(b,b,nvar);
		n >>= 1;
	}
	r = (r == NULL) ? b : multiplyPolynomials(r,b,nvar);

	return r;
}




/*****************
 * Polynomial division
 *****************/

/**
 * Helper to determine if monomial b divides monomial a.
 *
 * Note: since we assume working in Q, do not need to check coefficients 
 * nvar: number of variables of monomials a and b
 */
int monomialDivideTest(degrees_t adegs, degrees_t bdegs, int nvar) {
	for(int i = 0; i < nvar; ++i) {
		if (adegs[i] < bdegs[i]) {
			return 0;
		}
	}
	return 1;
}

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */
Node* divisionGetNextTerm(ProductHeap* h) {
	productHeapElem** elems = h->elements;
	int size = h->heapSize;
	int nvar = h->nvar;

#if SMQP_SUPPORT_CHAINED_HEAP
	if (size == 0) {
		return NULL;
	}

	degrees_t maxDegs = (degrees_t) malloc(sizeof(degree_t) * nvar);
	memcpy(maxDegs, elems[0]->product->degs, sizeof(degree_t)*nvar);
	Node* ret = addZeroTerm(NULL, maxDegs);
	productHeapElem* insertChain = NULL;
	productHeapElem* nextMaxElem;

	while(size > 0 && compareExponentVectors(elems[0]->product->degs, maxDegs, nvar) == 0) {
		productHeapElem* maxElem = prodheapRemoveMax(h);
		--size;

		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;
			mpq_add(ret->coef, ret->coef, maxElem->product->coef);
			maxElem->b = (maxElem->b)->next;
			if (maxElem->b != NULL) {
				addExponentVectors(maxElem->a_i->degs, maxElem->b->degs, maxElem->product->degs, nvar);
				mpq_mul(maxElem->product->coef, maxElem->a_i->coef, maxElem->b->coef);
				maxElem->next = insertChain;
				insertChain = maxElem;
			} else {
				prodheapFreeElement(maxElem);
			}
			maxElem = nextMaxElem;
		}
	}

	while(insertChain != NULL) {
		nextMaxElem = insertChain->next;
		insertChain->next = NULL;
		prodheapInsert(h,insertChain);
		insertChain = nextMaxElem;
	}

	if (mpq_sgn(ret->coef) == 0) {
		freeNode(ret);
		ret = NULL;
	}

#else

	Node* ret = NULL;
	Node* next;
	degrees_t	maxDegs = elems[0]->product->degs;

	while(size > 0 && compareExponentVectors(elems[0]->product->degs, maxDegs, nvar) == 0) {
		if (ret == NULL) {
			ret = elems[0]->product;
		} else {
			mpq_add(ret->coef, ret->coef, elems[0]->product->coef);
			freeNode(elems[0]->product);
			//either elems[0] is freed below or the product overwritten
		}

		next = elems[0]->b->next;
		if (next == NULL) {
			free(elems[0]);
			if (size > 1) {
				elems[0] = elems[--size];
				h->heapSize = size;
				prodheapSink(h,0);
			} else {
				elems[0] = NULL;
				--size;
				h->heapSize = size;
			}

		} else {
			elems[0]->b = next;
			elems[0]->product = multiplyTerms(elems[0]->a_i, next, nvar);	
			prodheapSink(h,0);
		}
	}

	if (ret != NULL && mpq_sgn(ret->coef) == 0) {
		freeNode(ret);
		ret = NULL;
	}


	/**
	 * The below continues to extract terms until there is a non-zerio
	 * coef on the product term to be extracted. This causes inconsistencies 
	 * when comparing delta (obtained from quoheapPeekMax) in dividePolynomials();
	 * One could use the below if quoheapPeekMax was modified to always return
	 * the max degree of the product term with a non-zero coef.
	 */
	// while (ret == NULL && size > 0) {
	// 	maxDegs = elems[0]->product->degs;

	// 	while(size > 0 && compareExponentVectors(elems[0]->product->degs, maxDegs, nvar) == 0) {
	// 		if (ret == NULL) {
	// 			ret = elems[0]->product;
	// 		} else {
	// 			mpq_add(ret->coef, ret->coef, elems[0]->product->coef);
	// 			freeNode(elems[0]->product);
	// 			//either elems[0] is freed below or the product overwritten
	// 		}

	// 		next = elems[0]->b->next;
	// 		if (next == NULL) {
	// 			free(elems[0]);
	// 			if (size > 1) {
	// 				elems[0] = elems[--size];
	// 				h->heapSize = size;
	// 				quoheapSink(h,0);
	// 			} else {
	// 				elems[0] = NULL;
	// 				--size;
	// 				h->heapSize = size;
	// 			}

	// 		} else {
	// 			elems[0]->b = next;
	// 			elems[0]->product = multiplyTerms(elems[0]->a_i, next, nvar);	
	// 			quoheapSink(h,0);
	// 		}
	// 	}

	// 	if (mpq_sgn(ret->coef) == 0) {
	// 		freeNode(ret);
	// 		ret = NULL;
	// 	}
	// }
#endif

	return ret;
}

/** 
 * Given a polynomial, c, and a term, b, determine polynomials a and r
 * such that c = b*a + r.
 * a and r are returned in res_a and res_r, respectively. 
 */
void divideBySingleTerm(Node* c, Node* b, Node** res_a, Node** res_r, int nvar) {
	if (b == NULL) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	Node* curC = c, *a = NULL, *r = NULL;
	degrees_t degs = NULL;
	ratNum_t coef;
	mpq_init(coef);

	*res_a = NULL;
	*res_r = NULL;

	while (curC != NULL) {
		if (monomialDivideTest(curC->degs,b->degs,nvar)) {
			mpq_div(coef, curC->coef, b->coef);
			// coef = curC->coef / b->coef;
			if (mpq_sgn(coef) != 0) {
				degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
				subtractExponentVectors(curC->degs, b->degs, degs, nvar);
				a = addTerm(a,degs,coef);
				if (*res_a == NULL) {
					*res_a = a;
				}
			}
		} else {
			degs = (degrees_t) malloc(sizeof(degree_t)*nvar);
			memcpy(degs, curC->degs, sizeof(degree_t)*nvar);
			// degs = curC->degs;
			mpq_set(coef, curC->coef);
			// coef = curC->coef;
			r = addTerm(r,degs,coef);
			if (*res_r == NULL) {
				*res_r = r;
			}
		}
		curC = curC->next;
	}
	mpq_clear(coef);
}

/** 
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r 
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials(Node* c, Node* b, Node** res_a, Node** res_r, int nvar) {

	if (b == NULL) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	// b does not divide c
	if(!monomialDivideTest(c->degs,b->degs,nvar)) {
		*res_a = NULL;
		*res_r = deepCopyPolynomial(c, nvar);
		return;
	}

	// b is a monomial so we can do a simple divide
	if (b->next == NULL) {
		divideBySingleTerm(c, b, res_a, res_r, nvar);
		return;
	}

	//init a with lt(c)/lt(b);
	ratNum_t bcoef;
	mpq_init(bcoef);
	mpq_set(bcoef, b->coef);
	// ratNum_t bcoef = b->coef;
	ratNum_t coef;
	mpq_init(coef);
	mpq_div(coef, c->coef, bcoef);
	degrees_t beta = b->degs;
	degrees_t eps = (degrees_t) malloc(sizeof(degree_t)*nvar);
	subtractExponentVectors(c->degs, beta, eps, nvar);
	Node* a = addTerm(NULL, eps, coef);
	Node* r = NULL;
	*res_a = a; //set head of resulting a to this leading term of a;
	*res_r = r;
	Node* k = c->next;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap* h = prodheapCreate(nvar);
	prodheapInsert(h, prodheapMakeElement(h, a, b->next));

	//loop variables
	degrees_t delta = NULL;
	degrees_t adegs = NULL;
	Node* multerm = NULL;
	int cmp;
	while(k != NULL || h->heapSize > 0) {
		ratNum_t e;
		mpq_init(e);

		delta = prodheapPeek(h);
		if (k == NULL) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors(delta, k->degs, nvar);
		}

		if (cmp > 0) {
			multerm = divisionGetNextTerm(h);
			if (multerm == NULL) {
				//in this case, the term with degree delta ended up 
				//having its coffeicient cancelled out (i.e. 0)
				continue;
			} else {
				eps = multerm->degs;
				mpq_neg(e, multerm->coef);
				// e = mpq_mul-1*multerm->coef;
			}
		} else if (cmp == 0) {
			multerm = divisionGetNextTerm(h);
			if (multerm == NULL) {
				continue;
			} else {
				eps = multerm->degs;
				mpq_sub(e, k->coef, multerm->coef);
				// e = k->coef - multerm->coef;
				k = k->next;
			}
		} else {
			multerm = NULL;
			eps = k->degs;
			mpq_set(e,k->coef);
			// e = k->coef;
			k = k->next;
		}

		if (mpq_sgn(e) != 0) {
			if (monomialDivideTest(eps, b->degs, nvar)) {
				adegs = (degrees_t) malloc(sizeof(degree_t)*nvar);
				subtractExponentVectors(eps, beta, adegs, nvar);
				mpq_div(e, e, bcoef);
				// e /= bcoef;
				a = addTerm(a, adegs, e);
				prodheapInsert(h, prodheapMakeElement(h, a, b->next));
				freeNode(multerm); //node and its degs
			} else {
				if (multerm != NULL) {
					r = addTerm(r,eps,e);
					free(multerm); //just the node, r holds multerm->degs
				} else {
					//else, eps is equal to k->degs, so make a copy.
					degrees_t temp = (degrees_t) malloc(sizeof(degree_t)*nvar);
					memcpy(temp, eps, sizeof(degree_t)*nvar);
					r = addTerm(r,temp,e);
				}
				if (*res_r == NULL) {
					*res_r = r;
				}
			}
		} else {
			freeNode(multerm);
		}
		mpq_clear(e);
	}
} 

/****************
* Primitive Factorization
*****************/

/**
 * @param[in] nvar The number of variable
 * @param[out] seq The recursive loop
 * \brief Return content of a where a = content(a) * primitivePart(a)
 */
int* recursiveLoop (int nvar)
{
  int size = (1 << nvar);
  int* seq = (int*) malloc( sizeof (int)* (size-1));
  for (int i = 0; i < (size-1); i++){ seq[i] = 0;}

  if (nvar == 1)
      return seq;
  else
    {
      int* seq_part = recursiveLoop (nvar-1);
      for (int i = 0; i < (size/2)-1; i++)
      {
    	  seq[i] = seq_part[i];
    	  seq[(size/2) +i] = seq_part[i];
      }
      free(seq_part);
      seq[(size/2)-1] = nvar-1;
      return seq;
    }
}

/** 
 * @param[in] a The polynomial
 * @param[out] cont The content of a
 * \brief Return content of a where a = content(a) * primitivePart(a)
 */
content_t* content (Node* a)
{
  Node* cur = a;
  mpzCoefs_t* n = (mpzCoefs_t*) malloc (sizeof (mpzCoefs_t));
  mpzCoefs_t* head;
  mpzCoefs_t* tail;
  mpzCoefs_t* cur_coef;
  content_t* gcd_res = (content_t*) malloc (sizeof (content_t));

  //make a mpzCoefs_t linked list which is coefs converted from mpq_t to mpz_t:
  mpz_init (n->mpzCoef);
  mpz_set (n->mpzCoef, mpq_numref (cur->coef));

  head = n;
  tail = n;

  cur = cur->next;
  while (cur != NULL)
    {
      n = (mpzCoefs_t*) malloc ( sizeof (mpzCoefs_t));
      mpz_init (n->mpzCoef);
      mpz_set (n->mpzCoef, mpq_numref(cur->coef));

      tail->next = n;
      tail = tail->next;

      cur = cur->next;
    }
  tail->next = NULL;

  // compute the content of mpz_t form of the coefs.
  cur_coef = head;
  mpz_init (gcd_res->cont);

  if (cur_coef->next != NULL)
    {
      mpz_gcd (gcd_res->cont, cur_coef->mpzCoef, cur_coef->next->mpzCoef);
      cur_coef = cur_coef->next;

      while (cur_coef->next != NULL)
	{
	  mpz_gcd (gcd_res->cont, gcd_res->cont, cur_coef->next->mpzCoef);
	  cur_coef = cur_coef->next;
	}
      return gcd_res;
    }
  else
    {
      mpz_set_ui (gcd_res->cont, 1ul);
      return gcd_res;
    }
}

/** 
 * @param[in] a The polynomial
 * @param[in] nvarThe number of variables
 * @param[out] returnContent The content of a
 * @param[out] head The primitive part of a
 * \brief Compute the primitive part of a and relevant content in returnContent.
 * For example, let  a = 4*x^2y^3 + 8*x*y^2 + 16 the primitivePart(a, returnContent, 2) 
 * returns x^2y^3 + 2*x*y^2 + 4 and the returnContent becomes 4.  
 */
Node* primitivePart (Node* a, content_t** returnContent,  int nvar)
{
  Node* cur = a;
  content_t* cont_tmp = content (a);
  Node* n = (Node*) malloc (sizeof (Node));
  Node* head;
  Node* tail;

  n->degs = cur->degs;

  mpq_init(n->coef);
  ratNum_t set_num;
  mpq_init (set_num);
  mpq_set_num (set_num, cont_tmp->cont);
  mpq_div(n->coef, cur->coef, set_num);

  head = n;
  tail = n;
  cur = cur->next;

  while (cur != NULL)
    {
      n = (Node*) malloc (sizeof (Node));
      n->degs = cur->degs;

      mpq_init(n->coef);
      ratNum_t set_num;
      mpq_init (set_num);
      mpq_set_num (set_num, cont_tmp->cont);
      mpq_div(n->coef, cur->coef, set_num);

      tail->next = n;
      tail = tail->next;
      cur = cur->next;
    }

  tail->next = NULL;

  *returnContent = cont_tmp;

  return head;
}

/****************
* Multi-Divisor Division
*****************/

/** Multi-Divisor Division (MDD) using Heap (Johnson's) Division
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void heapMDD(Node* f, Node** G, Node** Q, Node** r, int s, int nvar) {
	if (s < 0)
	{
		fprintf(stderr, "BPAS Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0)
	{
		*r = deepCopyPolynomial(f, nvar);
		return;
	}

	if (s == 1)  // simple multivariate polynomial division
	{
		dividePolynomials(f, *G, Q, r, nvar);
		return;
	}

	int i; // i= 0 ... (s-1)
	Node* h = deepCopyPolynomial(f, nvar);  // this algorithm runs until h != NULL
	Node* rem = NULL;
	Node* lt;
	Node* tmp_q = NULL;
	Node* tmp_r = NULL;

	while (h != NULL)
	{
		i = 0; // in each step i should start from the first divisor
		while (i < s)
		{
			if ( h != NULL && G[i] != NULL && monomialDivideTest(h->degs, G[i]->degs, nvar))
			{
				dividePolynomials(h, G[i], &tmp_q, &tmp_r, nvar);
				Q[i] = addPolynomials(Q[i], tmp_q, nvar);
				h = tmp_r;
				i = 0;
			} else
				i = i + 1;
		}
		if (h != NULL){
			lt = leadingTerm(h);
		}
		else
			lt = NULL;

		if (rem == NULL)
			rem = lt;
		else
		{
			if (lt != NULL) {
				Node* cur = rem;
				while (cur->next != NULL)
					cur = cur->next;
				cur->next = lt;
			}
		}
		Node* h_tmp = h;
		if (h_tmp != NULL){
			h = h_tmp->next;
			mpq_clear(h_tmp->coef);
			free(h_tmp);
		}
	}

	*r = rem;
	return;
}

/** Multi-Divisor Division (MDD) for Triangular Sets
 * @param[in] f The dividend
 * @param[in] G The triangular-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void triangularSetMDD(Node* f, Node** G, Node** Q, Node** r, int s, int nvar)
{
	if (s < 0)
	{
		fprintf(stderr, "BPAS Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0)
	{
		*r = deepCopyPolynomial(f, nvar);
		return;
	}

	if (s == 1)  // simple multivariate polynomial division
	{
		dividePolynomials(f, *G, Q, r, nvar);
		return;
	}

	int index = 0;
	int i = 0;
	int* orderList = recursiveLoop(s);
	int orderSize = (1 << s) - 1;
	Node* h =deepCopyPolynomial(f, nvar);
	Node* rem = NULL;
	Node* lt;
	Node* tmp_q = NULL;
	Node* tmp_r = NULL;

	while (h != NULL)
	{
		for (i = 0; i < orderSize; i++)
		{
			index = orderList[i];
			if (h != NULL && G[index] != NULL && monomialDivideTest(h->degs, G[index]->degs, nvar))
			{
				dividePolynomials(h, G[index], &tmp_q, &tmp_r, nvar);
				Q[index] = addPolynomials(Q[index], tmp_q, nvar);
				h = tmp_r;
			}
		}
		if (h != NULL)
			lt = leadingTerm(h);
		else
			lt = NULL;

		if (rem == NULL)
			rem = lt;
		else {
			if (lt != NULL){
				Node* cur = rem;
				while (cur->next != NULL)
					cur = cur->next;
				cur->next = lt;
			}
		}
		Node* h_tmp = h;
		if (h_tmp != NULL){
			h = h_tmp->next;
			mpq_clear(h_tmp->coef);
			free(h_tmp);
		}
	}

	*r = rem;
	free(orderList);
	return;
}

/** Multi-Divisor Division (MDD) for Triangular Sets using Primitive Factorization
 * @param[in] f The dividend
 * @param[in] G The triangular-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void primitiveFactorTriangularSetMDD(Node* f, Node** G, Node** Q, Node** r, int s, int nvar) {
	if (s < 0)
	{
		fprintf(stderr, "BPAS Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0)
	{
		*r = deepCopyPolynomial(f, nvar);
		return;
	}

	if (s == 1)  // simple multivariate polynomial division
			{
		dividePolynomials(f, *G, Q, r, nvar);
		return;
	}

	int index = 0;
	int i = 0;
	int* orderList = recursiveLoop(s);
	int orderSize = (1 << s) - 1;
	Node* rem = NULL;
	Node* lt;
	Node* tmp_q = NULL;
	Node* tmp_r = NULL;
	Node* pp = NULL;
	content_t* cont;

	pp = primitivePart(deepCopyPolynomial(f, nvar), &cont, nvar);
	Node* h = pp;

	while (h != NULL)
	{
		for (i = 0; i < orderSize; i++)
		{
			index = orderList[i];
			if (h != NULL && G[index] != NULL && monomialDivideTest(h->degs, G[index]->degs, nvar))
			{
				dividePolynomials(h, G[index], &tmp_q, &tmp_r, nvar);
				Q[index] = addPolynomials(Q[index], tmp_q, nvar);
				h = tmp_r;
			}
		}
		if (h != NULL)
			lt = leadingTerm(h);
		else
			lt = NULL;

		if (rem == NULL)
			rem = lt;
		else
		{
			if (lt != NULL){
				Node* cur = rem;
				while (cur->next != NULL)
					cur = cur->next;
				cur->next = lt;
			}
		}
		Node* h_tmp = h;
		if (h_tmp != NULL){
			h = h_tmp->next;
			mpq_clear(h_tmp->coef);
			free(h_tmp);
		}
	}

	mpz_t one;
	mpz_init(one);
	mpz_set_ui(one, 1ul);
	if (!mpz_cmp(cont->cont, one))
		*r = rem;
	else
	{
		ratNum_t set_num;
		mpq_init(set_num);
		mpq_set_num(set_num, cont->cont);
		Node* cur_r = rem;
		while (cur_r != NULL)
		{
		mpq_mul(cur_r->coef, cur_r->coef, set_num);
		cur_r = cur_r->next;
		}
		*r = rem;
	}

	free(orderList);
	return;
}
 
/** Multi-Divisor Division (MDD)
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[in] type The type of MDD (type = 0 ? HeapMDD : (type = 1 ? triangularSetMDD : primitiveFactorTriangularSetMDD))
 * @param[out] r The remainder
 * @param[out] Q The quotient-set
 * \brief This algorithm runs by f and G[s], and NULL polynomials r and Q[s],
 * Computes r, Q[0], ... and Q[s] such that f = Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r.
 */
void multiDivisorDivision (Node* f, Node** G, Node** Q, Node** r, int s, int nvar, int type)
{
	if (type == 0){
//		printf("Running HeapMDD ...\n");
		heapMDD (f, G, Q, r, s, nvar);
		return;
	} else if (type == 1){
//		printf("Running triangularSetMDD ...\n");
		triangularSetMDD (f, G, Q, r, s, nvar);
		return;
	} else{
//		printf("Running primitiveFactorTriangularSetMDD ...\n");
		primitiveFactorTriangularSetMDD (f, G, Q, r, s, nvar);
	}
	return;
}

/***************
* MDD Test and Verification
***************/

/**
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[out] 0-1 Return 1, unless there is a NULL polynomial
 */
int multiDivisorDivisionIsEmpty (Node* f, Node** G, int s)
{
  	if (f == NULL)
		return 0;
  	for (int i=0; i < s; ++i)
    {
    	if (G[i] == NULL)
	  		return 0;
    }
  	return 1;
}

/** Multi-Divisor Division Verification of Definition:
 * @param[in] f The dividend
 * @param[in] G The divisor-set
 * @param[in] s The size of G
 * @param[in] nvar The number of variables
 * @param[in] rr The remainder
 * @param[in] Q The quotient-set
 * @param[out] 0-1 Return 1, unless Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r != f
 * \brief This algorithm tests the Q[0]*G[0] + ... + Q[s-1]*G[s-1] + r -f = 0.
 */
int multiDivisorDivisionVerif (Node* ff, Node** G, Node** Q, Node* rr, int s, int nvar)
{
	Node* f = deepCopyPolynomial(ff, nvar);
	Node* r = deepCopyPolynomial(rr, nvar);
  	// generated f: gen_f
  	Node* gen_f = NULL;
  	for (int i = 0; i < s; i++)
    	gen_f = addPolynomials (gen_f, multiplyPolynomials (G[i], Q[i], nvar), nvar);

	//gen_f =  = G[0]*Q[0] + ... + G[s-1]*Q[s-1] + r
  	gen_f = addPolynomials (gen_f, r, nvar);
  	int fsize = numberOfTermsNode (f);
  	int gen_fsize = numberOfTermsNode (gen_f);
  	// Compare number of terms:
  	if (fsize == gen_fsize)
    {
      	Node* cur_f = f;
      	Node* cur_gen_f = gen_f;
      	for (int i = 0; i < fsize; i++)
		{
	  		// Compare Exponents:
	  		if (compareExponentVectors (cur_f->degs, cur_gen_f->degs, nvar) != 0)
	    	{
	    		fprintf(stderr, "Verification Error: Exponents are different!");
	      		return 0;
	    	}
	  		else
	    	{
	      		// Compare Coefficients:
	      		if (mpq_cmp (cur_f->coef, cur_gen_f->coef) != 0)
				{
					fprintf(stderr, "Verification Error: Coefficients are different!");
		  			return 0;
				}
	      		else
				{
		  			cur_f = cur_f->next;
		  			cur_gen_f = cur_gen_f->next;
				}
	    	}
		}
    	return 1;
    }
  	else
    {
    	fprintf(stderr, "Verification Error: Sizes are different!");
      	return 0;
    }
}



