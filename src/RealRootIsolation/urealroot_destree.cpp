#include "RationalNumberPolynomial/urpolynomial.h"

DesTree::DesTree() {
	root = new DesNode;

	root->k = 0;
	root->c = 0;

	root->parent = NULL;
        root->prev = NULL;
        root->left = NULL;
        root->right = NULL;
}

DesTree::DesTree(int _k, lfixz _c) {
	root = new DesNode;

        root->k = _k;
        root->c = _c;

        root->parent = NULL;
        root->prev = NULL;
        root->left = NULL;
        root->right = NULL;
}

DesTree::~DesTree() {
	delete root;
}

DesNode* DesTree::getNode() {
	return root;
}

void DesTree::addChildren(int k, lfixz c) {
	DesNode* leftChild = new DesNode;
	leftChild->k = k+1;
	leftChild->c = 2*c;
	leftChild->parent = root;
	leftChild->left = NULL;
	leftChild->right = NULL;
	leftChild->prev = NULL;

	root->left = leftChild;

	DesNode* rightChild = new DesNode;
	rightChild->k = k+1;
	rightChild->c = 2*c+1;
	rightChild->parent = root;
	rightChild->left = NULL;
	rightChild->right = NULL;
	rightChild->prev = NULL;

	root->right = rightChild;
}

DesNode* DesTree::nextNode(){
	if (root == NULL)
		return NULL;

	DesNode* prev = root;

	if (root->left != NULL) {
		root = root->left;
		root->prev = prev;
	}
	else {
		while (root->parent != NULL && root == root->parent->right)
			root = root->parent;
		root = root->parent;
		if (root != NULL) {
			root = root->right;
			if (root != NULL)
				root->prev = prev;
		}
	}

	return root;
}


long DenseUnivariateRationalPolynomial::taylorConstant(int k, lfixz c, int k1, lfixz c1) {
        int m = k - k1;

        lfixz elem;
        if (m > 0)
                elem = c1 << m;
        else if (m < 0)
                elem = c1 >> (-m);
        else
                elem = c1;

	elem -= c;
        return elem.get_si();
}


void DenseUnivariateRationalPolynomial::genDescartes(Intervals* pIs, DenseUnivariateRationalPolynomial* p, int ts) {
	DenseUnivariateRationalPolynomial A(*p);

        DesTree dT;
        DesNode* dN = dT.getNode();

        do {
                DenseUnivariateRationalPolynomial B(A);
                B.reciprocal();
                B.taylorShift(ts);

		int k = dN->k;
		lfixz c = dN->c;

                int s = B.numberOfSignChanges();
                if (s == 1) {
			lfixq left(c);
			left >>= k;
			lfixq right(c + 1);
			right >>= k;
                        pIs->pushInterval(left, right);
                }
                else if (s > 1) {
                        dT.addChildren(k, c);
                }

                dN = dT.nextNode();

                if (dN != NULL) {
			k = dN->k;
			c = dN->c;
                        int k1 = dN->prev->k;
                        lfixz c1 = dN->prev->c;

                        int m = taylorConstant(k1, c1, k, c);
                        if (m == 1)
                                A.taylorShift(ts);
                        else if (m) {
                                std::cout << "ERROR: Taylor constant " << m << "!" << std::endl;
                                exit(1);
                        }

			if (k1 < k)
                        	A.homothetic(k-k1);
			else if (k1 > k)
				A.scaleTransform(k1-k);

			// If current node is a right child, we need to check 
			// if there is an exact solution
                        if (dN->parent != NULL && dN == dN->parent->right && A.isConstantTermZero()) {
                                lfixq elem(c);
				elem >>= k;
				pIs->pushInterval(elem, elem);
                        }
                }

        } while (dN != NULL);
}
