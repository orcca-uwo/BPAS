

#include <gmpxx.h>
#include <sstream>

#include "LinearAlgebra/SystemSolving-IML.h"
#include "../MapleTestTool/MapleTestTool.hpp"

void testRandomSystem(int n, int bitSize) {

	time_t t = time(NULL);
	// fprintf(stderr, "time: %ld\n", t);

	mpq_t* A = randomMPQMat(n, n, bitSize, t);
	mpq_t* b = randomMPQMat(n, 1, bitSize, t);
	// fprintf(stderr, "A: \n");
	// printMPQMat(A, n, n);
	// fprintf(stderr, "\n");
	// fprintf(stderr, "b: \n");
	// printMPQMat(b, n, 1);

	mpq_t X[n];
	for (int i = 0; i < n; ++i) {
		mpq_init(X[i]);
	}

  	long* rowList;
  	long rank = solveMPQSystem(A, b, n, X, &rowList);

  	if (rank != n) {
  		fprintf(stderr, "We have a rank defficient matrix!\n");
		fprintf(stderr, "These rows are good:\n");
  		for (int i = 0; i < rank; ++i) {
  			fprintf(stderr, "%ld\n", rowList[i]);
  		}
		fprintf(stderr, "These rows are bad:\n");
  		for (int i = rank; i < n; ++i) {
  			fprintf(stderr, "%ld\n", rowList[i]);
  		}
  		return;
  	}

	std::vector<std::string> inputs;

	std::stringstream ss;
	ss << "Matrix([";
	for (int i = 0; i < n; ++i) {
		ss << "[";
		for (int j = 0; j < n-1; ++j) {
			ss << mpq_class(A[i*n+j]) << ","; 
		}
		ss << mpq_class(A[i*n+(n-1)]);
		ss << "]";
		if (i < n-1) {
			ss << ",";
		}
	}
	ss << "])";
	inputs.push_back(ss.str());

	ss.str("");
	ss << "Matrix([";
	for (int i = 0; i < n; ++i) {
		ss << "[";
		ss << mpq_class(b[i]);
		ss << "]";
		if (i < n-1) {
			ss << ",";
		}
	}
	ss << "])";
	inputs.push_back(ss.str());

	ss.str("");
	ss << "Matrix([";
	for (int i = 0; i < n; ++i) {
		ss << "[";
		ss << mpq_class(X[i]);
		ss << "]";
		if (i < n-1) {
			ss << ",";
		}
	}
	ss << "])";

	std::string expected = ss.str();

    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
	mapleTest->restartMapleKernel();

	std::string retErr;
	if(mapleTest->testProcReturn("LinearAlgebra:-LinearSolve", inputs, expected, "Matrix", &retErr) == 0) {
		std::cerr << "Linear Algebra test FAILED!" << std::endl;
		std::cerr << "A:= " << inputs[0] << ";" << std::endl;
		std::cerr << "b:= " << inputs[1] << ";" << std::endl;
		std::cerr << "x:= " << expected << ";" << std::endl;
		std::cerr << "got:= " << retErr << ";" << std::endl;
		exit(1);
	}

	for (int i = 0; i < n; ++i) {
		mpq_clear(X[i]);
	}

	freeMPQMat(A, n, n);
	freeMPQMat(b, n, 1);

	std::cerr << "Linear Algebra test PASSED!\n" << std::endl;
}

int main(int argc, char *argv[]) {

	int n = 10;
	int bitSize = 8;

	testRandomSystem(n, bitSize);

	return 0;
}

