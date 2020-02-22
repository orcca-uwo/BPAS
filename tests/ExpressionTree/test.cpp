

#include <bpas.h>

void testRNArray() {
	std::vector<RationalNumber> vec;
	vec.push_back(1);
	vec.push_back(2);
	vec.push_back(3);
	vec.push_back(4);
	vec.push_back(5);

	ExpressionTree tree;
	tree.fromVector<RationalNumber>(vec);
	
	std::string expected = "[1, 2, 3, 4, 5]";
	if (expected != tree.toString()) {
		std::cerr << "ExpressionTree: RNArray test FAILED" << std::endl;
		std::cerr << "Expected: " << expected << std::endl;
		std::cerr << "Got: " << tree.toString() << std::endl;
		exit(1);
	} else {
		std::cerr << "ExpressionTree: RNArray test PASSED" << std::endl;
	}
}

void testSMQPArray() {
	int nvar = 3;
	int numTerms = 10;
	int coefBound = 20;
	degree_t sparsity = 5;
	bool includeNeg = 1;

	SparseMultivariateRationalPolynomial p;
	p.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);
	SparseMultivariateRationalPolynomial q;
	q.randomPolynomial(nvar, numTerms, coefBound, sparsity, includeNeg);

	std::vector<SparseMultivariateRationalPolynomial> vec;
	vec.push_back(p);
	vec.push_back(q);

	ExpressionTree tree;
	tree.fromVector<SparseMultivariateRationalPolynomial>(vec);

	std::string expected = "[";
	expected += p.convertToExpressionTree().toString();
	expected += ", ";
	expected += q.convertToExpressionTree().toString();
	expected += "]";
	if (expected != tree.toString()) {
		std::cerr << "ExpressionTree: SMQP test FAILED" << std::endl;
		std::cerr << "Expected: " << expected << std::endl;
		std::cerr << "Got: " << tree.toString() << std::endl;
		exit(1);
	} else {
		std::cerr << "ExpressionTree: SMQP array test PASSED" << std::endl;
	}
}


int main() {
	testRNArray();
	testSMQPArray();
}