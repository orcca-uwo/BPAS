#include <bpas.h>

using namespace std;

extern bool testDUZP();
extern bool testDUQP();
extern bool testSUP();

// a := 108+324*x+387*x^2+238*x^3+80*x^4+14*x^5+x^6
// having squrefree: (x+1)*(x+2)*(x+3)
// b :=108+108*x+108*x^2+216*x^3+135*x^4+135*x^5+144*x^6+63*x^7+63*x^8+40*x^9+13*x^10+13*x^11+4*x^12+x^13+x^14
// having squarefree: (x+1)*(x^2+2)*(x^3+3)
int main(int argc, char *argv[]) {
	bool isDUZP = testDUZP();
	bool isDUQP = testDUQP();
	bool isSUP = testSUP();
	if (isDUZP && isDUQP && isSUP)
		cout << "pass." << endl;
	else {
		if (!isDUZP)
			cout << "DUZP error." << endl;
		if (!isDUQP)
			cout << "DUQP error." << endl;
		if (!isSUP)
			cout << "SUP error." << endl;
	}
        return 0;
}

bool testDUZP() {
	bool isDUZP = 1;

	DUZP p (4), q (7);
	p.setVariableName(Symbol("x"));
	p.setCoefficient(0, 6);
	p.setCoefficient(1, 11);
	p.setCoefficient(2, 6);
	p.setCoefficient(3, 1);
	q.setVariableName(Symbol("x"));
	q.setCoefficient(0, 6);
	q.setCoefficient(1, 6);
	q.setCoefficient(2, 3);
	q.setCoefficient(3, 5);
	q.setCoefficient(4, 2);
	q.setCoefficient(5, 1);
	q.setCoefficient(6, 1);

	DUZP a (7), b (15);
	a.setVariableName(Symbol("x"));
	a.setCoefficient(0, 108);
	a.setCoefficient(1, 324);
	a.setCoefficient(2, 387);
	a.setCoefficient(3, 238);
	a.setCoefficient(4, 80);
	a.setCoefficient(5, 14);
	a.setCoefficient(6, 1);
	b.setVariableName(Symbol("x"));
	b.setCoefficient(0, 108);
	b.setCoefficient(1, 108);
	b.setCoefficient(2, 108);
	b.setCoefficient(3, 216);
	b.setCoefficient(4, 135);
	b.setCoefficient(5, 135);
	b.setCoefficient(6, 144);
	b.setCoefficient(7, 63);
	b.setCoefficient(8, 63);
	b.setCoefficient(9, 40);
	b.setCoefficient(10, 13);
	b.setCoefficient(11, 13);
	b.setCoefficient(12, 4);
	b.setCoefficient(13, 1);
	b.setCoefficient(14, 1);

	vector<Factor<DUZP>> c = a.squareFree().factors();
	vector<Factor<DUZP>> d = b.squareFree().factors();

	DUZP x, y;
	x.one(), y.one();
	for (int i = 0; i < c.size(); ++i)
		if (!c[i].first.isZero())
			x *= c[i].first;
	for (int i = 0; i < d.size(); ++i)
		if (!d[i].first.isZero())
			y *= d[i].first;

	if (p != x || q != y) {
		x = -x, y = -y;
		if (p != x || q != y)
			isDUZP = 0;
	}

	return isDUZP;
}

bool testDUQP() {
	bool isDUQP = 1;

	DUQP p (4), q (7);
        p.setVariableName(Symbol("x"));
        p.setCoefficient(0, 6);
        p.setCoefficient(1, 11);
        p.setCoefficient(2, 6);
        p.setCoefficient(3, 1);
        q.setVariableName(Symbol("x"));
        q.setCoefficient(0, 6);
        q.setCoefficient(1, 6);
        q.setCoefficient(2, 3);
        q.setCoefficient(3, 5);
        q.setCoefficient(4, 2);
        q.setCoefficient(5, 1);
        q.setCoefficient(6, 1);

	DUQP a (7), b (15);
        a.setVariableName(Symbol("x"));
        a.setCoefficient(0, 108);
        a.setCoefficient(1, 324);
        a.setCoefficient(2, 387);
        a.setCoefficient(3, 238);
        a.setCoefficient(4, 80);
        a.setCoefficient(5, 14);
        a.setCoefficient(6, 1);
        b.setVariableName(Symbol("x"));
        b.setCoefficient(0, 108);
        b.setCoefficient(1, 108);
        b.setCoefficient(2, 108);
        b.setCoefficient(3, 216);
        b.setCoefficient(4, 135);
        b.setCoefficient(5, 135);
        b.setCoefficient(6, 144);
        b.setCoefficient(7, 63);
        b.setCoefficient(8, 63);
        b.setCoefficient(9, 40);
        b.setCoefficient(10, 13);
        b.setCoefficient(11, 13);
        b.setCoefficient(12, 4);
        b.setCoefficient(13, 1);
        b.setCoefficient(14, 1);

	std::vector<Factor<DUQP>> c = a.squareFree().factors();	
	std::vector<Factor<DUQP>> d = b.squareFree().factors();

	DUQP x, y;
        x.one(), y.one();
        for (int i = 0; i < c.size(); ++i)
                if (!c[i].first.isZero())
                        x *= c[i].first;
        for (int i = 0; i < d.size(); ++i)
                if (!d[i].first.isZero())
                        y *= d[i].first;

        if (p != x || q != y)
		isDUQP = 0;

	return isDUQP;
}

bool testSUP() {
	bool isSUP = 1;

	SparseUnivariatePolynomial<Integer> p, q;
	p.setVariableName(Symbol("x"));
        p.setCoefficient(0, 6);
        p.setCoefficient(1, 11);
        p.setCoefficient(2, 6);
        p.setCoefficient(3, 1);
        q.setVariableName(Symbol("x"));
        q.setCoefficient(0, 6);
        q.setCoefficient(1, 6);
        q.setCoefficient(2, 3);
        q.setCoefficient(3, 5);
        q.setCoefficient(4, 2);
        q.setCoefficient(5, 1);
        q.setCoefficient(6, 1);

	SparseUnivariatePolynomial<Integer> a, b;
	a.setVariableName(Symbol("x"));
        a.setCoefficient(0, 108);
        a.setCoefficient(1, 324);
        a.setCoefficient(2, 387);
        a.setCoefficient(3, 238);
        a.setCoefficient(4, 80);
        a.setCoefficient(5, 14);
        a.setCoefficient(6, 1);
        b.setVariableName(Symbol("x"));
        b.setCoefficient(0, 108);
        b.setCoefficient(1, 108);
        b.setCoefficient(2, 108);
        b.setCoefficient(3, 216);
        b.setCoefficient(4, 135);
        b.setCoefficient(5, 135);
        b.setCoefficient(6, 144);
        b.setCoefficient(7, 63);
        b.setCoefficient(8, 63);
        b.setCoefficient(9, 40);
        b.setCoefficient(10, 13);
        b.setCoefficient(11, 13);
        b.setCoefficient(12, 4);
        b.setCoefficient(13, 1);
        b.setCoefficient(14, 1);

        Factors<SparseUnivariatePolynomial<Integer>> cf = a.squareFree();
        Factors<SparseUnivariatePolynomial<Integer>> df = b.squareFree();
	vector< Factor<SparseUnivariatePolynomial<Integer>>> c = cf.factors();
    vector< Factor<SparseUnivariatePolynomial<Integer>> > d = df.factors();

	SparseUnivariatePolynomial<Integer> x, y;
	if (c.size()) { x = cf.ringElement(); }
	if (d.size()) { y = df.ringElement(); }
        for (int i = 0; i < c.size(); ++i)
                if (!c[i].first.isZero())
                        x *= c[i].first;
        for (int i = 0; i < d.size(); ++i)
                if (!d[i].first.isZero())
                        y *= d[i].first;

        if (p != x || q != y) {
		x = -x, y = -y;
		if (p != x || q != y)
			isSUP = 0;
	}

	return isSUP;
}
