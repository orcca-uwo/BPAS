
#include <sstream>
#include <bpas.h>

void testDefaultConstructor() {
	Factors<RationalNumber> f;

	std::stringstream ss;
	ss << f;
	if (ss.str() != "[1, []]") {
		std::cerr << "Factors default constructor test: FAILED" << std::endl;
		std::cerr << "Got: " << f << std::endl;
		exit(1);
	} 

	std::cerr << "Factors default constructor test: PASSED" << std::endl;
}

void testRingConstructor() {
	Factors<RationalNumber> f(14);

	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(14, 1);

	std::vector<Factor<RationalNumber>> retV = f.factors();

	if (v.size() != retV.size()) {
		std::cerr << "Factors ring constructor test: FAILED" << std::endl;
		std::cerr << "Should only be 1 factor be got " << retV.size() << std::endl;
		exit(1);
	}

	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors ring constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	} 

	std::cerr << "Factors ring constructor test: PASSED" << std::endl;
}

void testVectorConstructor() {
	std::vector<RationalNumber> v;
	v.push_back(3);
	v.push_back(5);
	v.push_back(7);

	Factors<RationalNumber> f(v);

	std::vector<Factor<RationalNumber>> retV = f.factors();

	if (f.ringElement() != 1) {
		std::cerr << "Factors vector constructor test: FAILED" << std::endl;
		std::cerr << "Expected 1 but got " << f.ringElement() << std::endl;
		exit(1);
	}
	for (int i = 0; i < v.size(); ++i) {
		if (v[i] != retV[i].first || 1 != retV[i].second) {
			std::cerr << "Factors vector constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	} 

	std::vector<int> e;
	e.push_back(2);
	e.push_back(3);
	e.push_back(1);

	Factors<RationalNumber> f2(v, e);
	
	if (f.ringElement() != 1) {
		std::cerr << "Factors vector constructor test: FAILED" << std::endl;
		std::cerr << "Expected 1 but got " << f.ringElement() << std::endl;
		exit(1);
	}
	retV = f2.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i] != retV[i].first || e[i] != retV[i].second) {
			std::cerr << "Factors vector constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}

	RationalNumber u(17);
	Factors<RationalNumber> f3(v, e, u);

	if (u != f3.ringElement()) {	
		std::cerr << "Factors vector constructor test: FAILED" << std::endl;
		std::cerr << "Expected " << u << " but got " << f3.ringElement() << std::endl;
		exit(1);
	}
	retV = f3.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i] != retV[i].first || e[i] != retV[i].second) {
			std::cerr << "Factors vector constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}

	std::cerr << "Factors vector constructor test: PASSED" << std::endl;
}

void testFactorConstructor() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	Factors<RationalNumber> f(v);

	if (f.ringElement() != 1) {
		std::cerr << "Factors factor constructor test: FAILED" << std::endl;
		std::cerr << "Expected 1 but got " << f.ringElement() << std::endl;
		exit(1);
	}

	std::vector<Factor<RationalNumber>> retV = f.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors factor constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	} 

	RationalNumber u(17);

	Factors<RationalNumber> f2(v, u);

	if (f2.ringElement() != u) {
		std::cerr << "Factors factor constructor test: FAILED" << std::endl;
		std::cerr << "Expected " << u << " but got " << f.ringElement() << std::endl;
		exit(1);
	}

	retV = f2.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors factor constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	} 

	std::cerr << "Factors factor constructor test: PASSED" << std::endl;
}

void testCopyConstructor() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	Factors<RationalNumber> f(v);

	Factors<RationalNumber> f2(f);

	if (f.ringElement() != f2.ringElement()) {
		std::cerr << "Factors copy consturctor test: FAILED" << std::endl;
		std::cerr << "Expected " << f.ringElement() << " but got " << f2.ringElement() << std::endl;
		exit(1);
	}

	std::vector<Factor<RationalNumber>> v1 = f.factors();
	std::vector<Factor<RationalNumber>> v2 = f2.factors();
	if (v1.size() != v2.size()) {
		std::cerr << "Factors copy consturctor test: FAILED" << std::endl;
		std::cerr << "Vector of factors of different sizes" << std::endl;
		exit(1);
	};

	for (int i = 0; i < v1.size(); ++i) {
		if (v1[i].first != v2[i].first || v1[i].second != v2[i].second) {
			std::cerr << "Factors copy constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v1[i] << " but got " << v2[i] << std::endl;
			exit(1);
		}	
	}

	std::cerr << "Factors copy constructor test: PASSED" << std::endl;
}

void testMoveConstructor() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);

	Factors<RationalNumber> f2 = std::move(f);

	if (f.ringElement() != 1) {
		std::cerr << "Factors move consturctor test: FAILED" << std::endl;
		std::cerr << "F unit should be 1 but got " << f.ringElement() << std::endl;
		exit(1);	
	}
	if (f.factors().size() != 0) {
		std::cerr << "Factors move consturctor test: FAILED" << std::endl;
		std::cerr << "F should have no factors but has " << f.factors().size() << std::endl;
		exit(1);	
	}

	if (f2.ringElement() != u) {
		std::cerr << "Factors move constructor test: FAILED" << std::endl;
		std::cerr << "Expected " << u << "but got " << f2.ringElement() << std::endl;
		exit(1);
	}

	std::vector<Factor<RationalNumber>> retV = f2.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors move constructor test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}	

	std::cerr << "Factors move constructor test: PASSED" << std::endl;
}

void testCopyAssignment() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	Factors<RationalNumber> f(v);

	Factors<RationalNumber> f2;
	RationalNumber temp = f2.ringElement();

	f2 = f;

	if (f.ringElement() != f2.ringElement()) {
		std::cerr << "Factors copy assignment test: FAILED" << std::endl;
		std::cerr << "Expected " << f.ringElement() << " but got " << f2.ringElement() << std::endl;
		exit(1);
	}

	std::vector<Factor<RationalNumber>> v1 = f.factors();
	std::vector<Factor<RationalNumber>> v2 = f2.factors();
	if (v1.size() != v2.size()) {
		std::cerr << "Factors copy assignment test: FAILED" << std::endl;
		std::cerr << "Vector of factors of different sizes" << std::endl;
		exit(1);
	};

	for (int i = 0; i < v1.size(); ++i) {
		if (v1[i].first != v2[i].first || v1[i].second != v2[i].second) {
			std::cerr << "Factors copy assignment test: FAILED" << std::endl;
			std::cerr << "Expected " << v1[i] << " but got " << v2[i] << std::endl;
			exit(1);
		}	
	}

	std::cerr << "Factors copy assignment test: PASSED" << std::endl;
}

void testMoveAssignment() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);

	Factors<RationalNumber> f2;
	RationalNumber temp = f.ringElement();

	f2 = std::move(f);

	if (f.ringElement() != 1) {
		std::cerr << "Factors move assignment test: FAILED" << std::endl;
		std::cerr << "F unit should be 1 but got " << f.ringElement() << std::endl;
		exit(1);	
	}
	if (f.factors().size() != 0) {
		std::cerr << "Factors move assignment test: FAILED" << std::endl;
		std::cerr << "F should have no factors but has " << f.factors().size() << std::endl;
		exit(1);	
	}

	if (f2.ringElement() != u) {
		std::cerr << "Factors move assignment test: FAILED" << std::endl;
		std::cerr << "Expected " << u << "but got " << f2.ringElement() << std::endl;
		exit(1);
	}

	std::vector<Factor<RationalNumber>> retV = f2.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors move assignment test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}	

	std::cerr << "Factors move assignment test: PASSED" << std::endl;
}

void testRingElement() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);

	RationalNumber test(18);

	if (f.ringElement() != u) {
		std::cerr << "Factors ringElement() test: FAILED" << std::endl;
		std::cerr << "Expected " << u << " but got " << f.ringElement() << std::endl;
		exit(1);
	}
	f.setRingElement(test);
	if (f.ringElement() != test) {
		std::cerr << "Factors setRingElement() test: FAILED" << std::endl;
		std::cerr << "Expected " << test << " but got " << f.ringElement() << std::endl;
		exit(1);
	}

	std::cerr << "Factors ringElement() test: PASSED" << std::endl;
	std::cerr << "Factors setRingElement() test: PASSED" << std::endl;
}

void testMultiplyRingElement() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);

	RationalNumber u2(12);

	f.multiplyRingElement(u2);

	if (f.ringElement() != (u*u2)) {
		std::cerr << "Factors multiplyRingElement() test: FAILED" << std::endl;
		std::cerr << "Expected " << u*u2 << " but got " << f.ringElement() << std::endl;
		exit(1);
	}

	std::cerr << "Factors multiplyRingElement() test: PASSED" << std::endl;

}

void testFactors() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);

	std::vector<Factor<RationalNumber>> retV = f.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors factors() test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}	

	Factor<RationalNumber> test;
	test.first = 11;
	test.second = 2;
	v.insert(v.begin()+1, test);
	f.setFactors(v);
	retV = f.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors setFactors() test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}	

	std::cerr << "Factors factors() test: PASSED" << std::endl;
	std::cerr << "Factors setFactors() test: PASSED" << std::endl;

}

void testAddFactor() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	// v.emplace_back(3, 2);
	// v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);
	Factor<RationalNumber> fact(3, 2);
	v.push_back(fact);
	f.addFactor(fact);

	std::vector<Factor<RationalNumber>> retV = f.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors addFactor() test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}	

	v.emplace_back(5,3);
	f.addFactor(5,3);
	retV = f.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors addFactor() test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}	
	
	std::cerr << "Factors addFactor() test: PASSED" << std::endl;
}

void testAddFactors() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);
	Factors<RationalNumber> f2(v, u);
	f.addFactors(f2);

	for (int i = 0; i < v.size(); ++i) {
		v[i].second *= 2;
	}
	// v.insert(v.begin(), v.begin(), v.end());

	std::vector<Factor<RationalNumber>> retV = f.factors();
	for (int i = 0; i < v.size(); ++i) {
		if (v[i].first != retV[i].first || v[i].second != retV[i].second) {
			std::cerr << "Factors addFactors() test: FAILED" << std::endl;
			std::cerr << "Expected " << v[i] << " but got " << retV[i] << std::endl;
			exit(1);
		}
	}	

	std::cerr << "Factors addFactors() test: PASSED" << std::endl;	
}

void testSize() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);

	if (f.size() != v.size()) {
		std::cerr << "Factors size() test: FAILED" << std::endl;
		std::cerr << "Expected a size of " << v.size() << " but got " << f.size() << std::endl;
		exit(1);
	}

	Factors<RationalNumber> f2;
	if (f2.size() != 0) {
		std::cerr << "Factors size() test: FAILED" << std::endl;
		std::cerr << "Expected a size of 0 but got " << f.size() << std::endl;
		exit(1);
	}	

	std::cerr << "Factors size() test: PASSED" << std::endl;
}

void testArrayIndexing() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);

	for (int i = 0; i < v.size(); ++i) {
		if (v[i] != f[i]) {
			std::cerr << "Factors operator[] test: FAILED" << std::endl;
			std::cerr << "At index " << i << " expected " << v[i] << " but got " << f[i] << std::endl;
			exit(1);
		}
	}

	std::cerr << "Factors operator[] test: PASSED" << std::endl;
}

void testEquality() {
	std::vector<Factor<RationalNumber>> v;
	v.emplace_back(7, 5);
	v.emplace_back(3, 2);
	v.emplace_back(5, 3);

	RationalNumber u(17);

	Factors<RationalNumber> f(v, u);
	Factors<RationalNumber> f2(v, u);
	
	Factors<RationalNumber> f3(v, RationalNumber(18));


	Factor<RationalNumber> test;
	test.first = 11;
	test.second = 2;
	v.insert(v.begin()+1, test);
	Factors<RationalNumber> f4(v, u);

	if (!(f == f2) || f == f3 || f == f4 || f3 == f4) {
		std::cerr << "Factors operator== test: FAILED" << std::endl;
		exit(1);
	}
	
	if (f != f2 || !(f != f3) || !(f != f4) || !(f3 != f4)) {
		std::cerr << "Factors operator!= test: FAILED" << std::endl;
		exit(1);
	}

	std::cerr << "Factors operator== test: PASSED" << std::endl;
	std::cerr << "Factors operator!= test: PASSED" << std::endl;
}

int main(void) {
	testDefaultConstructor();
	testRingConstructor();
	testVectorConstructor();
	testFactorConstructor();
	testCopyConstructor();
	testMoveConstructor();
	testCopyAssignment();
	testMoveAssignment();
	testRingElement();
	testMultiplyRingElement();
	testFactors();
	testAddFactor();
	testAddFactors();
	testSize();
	testArrayIndexing();
	testEquality();

	return 0;
}

