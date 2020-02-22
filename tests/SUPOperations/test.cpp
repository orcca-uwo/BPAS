#include <bpas.h>

using namespace std;

extern void testRationalNumber();
extern void testComplexRationalNumber();
extern void testUQPoperations(int);
extern void testSubresultant();
extern void testGCD();

void foo(Symbol s) {
	std::cerr << "foo" << std::endl;
}

int main(int argc, char *argv[]) {

	int s = atoi(argv[1]);
	//Ring<Integer>::SUP x;
	/*
	Integer x;
	SparseUnivariatePolynomial<Integer> f;
	x = 1;
	f = ((f - x << 2) + x + x << 1) - x;
	cout << "f: " << f << endl;
	
	DUQP g(4);
	//g.setCoefficient(0, 1);
	g.setCoefficient(1, 2);
	g.setCoefficient(2, -1);
	g.setCoefficient(3, 1);
	
	RationalNumber c;
	c = 0.5;
	g += c;
	g = g << 2 >> 1;
	g += mpq_class(1);
	cout << "g: " << g << endl;
	*/
	testUQPoperations(s);
	testRationalNumber();
	testComplexRationalNumber();
	return 0;
}

void testRationalNumber() {
	cout << "-- RationalNumber --" << endl;
	RationalNumber z1(1,2);
	RN z2(2,128);

	for (int i = 1; i <= 3; i++) {
		z1 = z1^i;
		//cout << z1 << endl;
	}
	if (z1 == z2)
		cout << "exponentiation pass." << endl;
		
	cout << "-- SUP<RationalNumber> --" << endl;
	z1.set(1,2);
	SparseUnivariatePolynomial<RN> P;
	SparseUnivariatePolynomial<RN> Q;
	P.setVariableName(Symbol("x"));
	Q.setVariableName(Symbol("x"));
	P.setCoefficient(0,z1);
	P.setCoefficient(2,z1);
	P *= P;
	Q.setCoefficient(2,z1);
	z1 /= 2;
	Q.setCoefficient(0,z1);
	Q.setCoefficient(4,z1);
	if (P == Q)
		cout << "multiplication equality pass: " << P << " = " << Q  << endl;
}


void testComplexRationalNumber() {
	cout << "-- ComplexRationalNumber --" << endl;
	ComplexRationalNumber z1(1,1);
	ComplexRationalNumber z2;

	for (int i = 1; i <= 8; i++) {
		z2 = z1^i;
		//cout << z2 << endl;
	}
	//z1 = 16;
	if (z2 == 16)
		cout << "exponentiation pass." << endl;

	z2.set(1,1);
	if (z1 == z2)
		cout << "set equality pass: " << z1 << " = " << z2 << endl;
		
	cout << "-- SUP<ComplexRationalNumber> --" << endl;
	CRN z(0,2);
	SparseUnivariatePolynomial<CRN> P;
	SparseUnivariatePolynomial<CRN> Q;
	P.setVariableName(Symbol("x"));
	Q.setVariableName(Symbol("x"));
	P.setCoefficient(0,z1);
	P.setCoefficient(2,z1);
	P *= P;
	Q.setCoefficient(0,z);
	Q.setCoefficient(4,z);
	z *= 2;
	Q.setCoefficient(2,z);
	if (P == Q)
		cout << "multiplication equality pass: " << P << " = " << Q  << endl;
}

void testGCD() {
	//cout << "---- GCD ----" << endl;
	Integer a, b;
	a = 40, b = 24;
	a.gcd(b);

	DUZP f(3), g(4);
	f.setVariableName(Symbol("t"));
	g.setVariableName(Symbol("t"));
	f.setCoefficient(1, -1);
	f.setCoefficient(2, 2);
	//f.content();
	g.setCoefficient(1, -3);
	g.setCoefficient(2, 5);
	g.setCoefficient(3, 2);
	//cout << f << " & " << g << endl;
	DUZP h = f.gcd(g);
	//cout << "DUZP gcd: " << h << endl;
	//h = f.gcd(g, 1);
	//cout << "modular gcd: " << h << endl;

	SparseUnivariatePolynomial<Integer> p, q;
	p.setVariableName(Symbol("t"));
	q.setVariableName(Symbol("t"));
	a = -1;
	p.setCoefficient(1, a);
	a = 2;
	p.setCoefficient(2, a);
	q.setCoefficient(3, a);
	a = -3;
	q.setCoefficient(1, a);
	a = 5;
	q.setCoefficient(2, a);
	//cout << p << " & " << q << endl;
	SparseUnivariatePolynomial<Integer> r = p.gcd(q);
	//cout << "SUP<Integer> gcd: " << r << endl;

	SparseUnivariatePolynomial< SparseUnivariatePolynomial<Integer> > x, y;
	x.setCoefficient(1, p);
	y.setCoefficient(2, q);

	SparseUnivariatePolynomial< SparseUnivariatePolynomial<Integer> > z = x.gcd(y);
	//cout << "SUP< SUP<Integer> > gcd: " << z << endl;

	SparseUnivariatePolynomial<DUZP> m, n;
	m.setCoefficient(1, f);
	n.setCoefficient(2, g);

	//cout << m << " & " << n << endl;
	SparseUnivariatePolynomial<DUZP> l = m.gcd(n);
	//cout << "SUP<DUZP> gcd: " << l << endl;

	cout << "gcd pass." << endl;
}


void testSubresultant() {
	DUZP a(1);
	vector< SparseUnivariatePolynomial<DUZP> > R;
	SparseUnivariatePolynomial<DUZP> p, q;
	//a.setCoefficient(0, 1);
	a = mpz_class(1) + a;
	p.setCoefficient(6, a);
	q.setCoefficient(0, a);
	q.setCoefficient(5, a);
	a.setCoefficient(0, 2);
	p.setCoefficient(2, a);
	p.setCoefficient(5, a);
	p.setCoefficient(7, a);
	p.setCoefficient(8, a);
	q.setCoefficient(2, a);
	q.setCoefficient(6, a);
	a.setCoefficient(0, 3);
	p.setCoefficient(1, a);
	q.setCoefficient(1, a);
	q.setCoefficient(4, a);
	a.setCoefficient(0, 4);
	p.setCoefficient(4, a);
	q.setCoefficient(7, a);

	//cout << "p = " << p << endl;
	//cout << "q = " << q << endl;

	R = p.subresultantChain(q);

	//for (int i = R.size()-1; i > -1; --i)
	//	cout << "R[" << i << "] = " << R[i] << endl;

	cout << "subresultant pass." << endl;

	Integer e;
	SparseUnivariatePolynomial<Integer> f, g;
	vector< SparseUnivariatePolynomial<Integer> > S;
	e = 1;
	f.setCoefficient(1, e);
	f.setCoefficient(5, e);
	g.setCoefficient(0, e);
	g.setCoefficient(1, e);
	g.setCoefficient(4, e);
	e = 2;
	f.setCoefficient(0, e);
	e = 3;
	g.setCoefficient(2, e);
	g.setCoefficient(3, e);
	e = 4;
	f.setCoefficient(2, e);
	f.setCoefficient(3, e);
	f.setCoefficient(4, e);

	//cout << "f = " << f << endl;
	//cout << "g = " << g << endl;

	S = f.subresultantChain(g);

	//for (int i = S.size()-1; i > -1; --i)
	//	cout << "S[" << i << "] = " << S[i] << endl;
}

void testUQPoperations(int s) {
	Integer temp, temp2;
	temp = mpz_class(0), temp2 = 40;
	temp += 2; 
	temp *= 7; 
	temp2 = temp + temp2;
	temp2 -= 12;
	temp2 /= temp;
	Integer temp3(temp);
	temp3 = temp^2;
	temp3 += temp * temp2;

	SparseUnivariatePolynomial<Integer> a, b, quo, rem;
	gmp_randclass rr (gmp_randinit_default);
	rr.seed(time(NULL));
	Integer elem;
	for (int i = 0; i < s; ++i) {
                elem = Integer((-1*rr.get_z_bits(s-1))-1);
                a.setCoefficient(i, elem);
                elem = Integer(rr.get_z_bits(s+1)+1);
		b.setCoefficient(i, elem);
	}
	elem = 1;
	a.setCoefficient(s, elem);
	b.setCoefficient(s, elem);

	DUZP qa = a.convertToDUZP();
	DUZP qb = b.convertToDUZP();

	DUZP qc = qa + qb;
	DUZP qd = qa - qb;
	DUZP qe = qa * qb;

        SparseUnivariatePolynomial<Integer> c = a + b;
        SparseUnivariatePolynomial<Integer> d = a - b;
        SparseUnivariatePolynomial<Integer> e = a * b;

	cout << "-- SUP<Integer> --" << endl;
	if (!(c == qc)) { cout << "addition error!" << endl; }
	else { cout << "addition pass." << endl; }
	if (!(d == qd)) { cout << "subtraction error!" << endl; }
	else { cout << "subtraction pass." << endl; }
	if (!(e == qe)) { cout << "multiplication error!" << endl; }
	else { cout << "multiplication pass." << endl; }

	DUZP qf = qe / qa;
	e /= a;
	if (!(e == qf)) { cout << "exact division error!" << endl; }
	else { cout << "exact division pass." << endl; }

	DUZP qr, qm = qc.monicDivide(qb, &qr);
	quo = c.monicDivide(b, &rem);
	if (!(quo == qm) || !(rem == qr)) { cout << "monic division error!" << endl; }
	else { cout << "monic division pass." << endl; }

	Integer qx, qy;
	DUZP qp = qc.lazyPseudoDivide(qd, &qr, &qx, &qy);
	Integer x, y;
	quo = c.lazyPseudoDivide(d, &rem, &x, &y);
	if (!(quo == qp) || !(rem == qr)) { cout << "lazy pseudo division error!" << endl; }
	else { cout << "lazy pseudo division pass." << endl; }

	//qc.pseudoDivide(qd, &qr, &qx);
	qc.pseudoDivide(qd);

	rem = rem ^ 2;
	qr = qr ^ 2;
	if (!(rem == qr)) { cout << "exponentiation error!" << endl; }
	else { cout << "exponentiation pass" << endl; }

	SparseUnivariatePolynomial<Integer> a2 = a;
	DUZP qa2 = qa;
	a2.differentiate(2);
	qa2.differentiate(2);
	if (!(a2 == qa2)) { cout << "differentiate error!" << endl; }
	else { cout << "differentiate pass." << endl; }

	temp = a.evaluate(temp2);
	temp3 = qa.evaluate(temp2);
	if (temp != temp3) { cout << "evaluate error!" << endl; }
	else { cout << "evaluate pass." << endl; }

	testGCD();

	cout << "-- SUP<DUZP> --" << endl;

	//SparseUnivariatePolynomial< SparseUnivariatePolynomial<Integer> > p, q, z, l;

	DUZP za = a.convertToDUZP();
	DUZP zb = b.convertToDUZP();
	za.setVariableName(Symbol("x"));
	zb.setVariableName(Symbol("x"));

	SparseUnivariatePolynomial< DUZP > p, q, z;
	p.setVariableName(Symbol("y"));
	q.setVariableName(Symbol("y"));
	p.setCoefficient(0, za);
	q.setCoefficient(0, za);
	p.setCoefficient(1, zb);

	//cout << "p := " << p << endl;
	//cout << "q := " << q << endl;

	z = p + q;
	//cout << "z := " << z << endl;
	cout << "addition pass." << endl;

	z *= q;
	//cout << "z := " << z << endl;
	cout << "multiplication pass." << endl;

	z /= q;
	//cout << "z := " << z << endl;
	cout << "exact division pass." << endl;

	testSubresultant();

	cout << "-- SUP< SUP< DUZP > > --" << endl;
	//SparseUnivariatePolynomial< SparseUnivariatePolynomial< SparseUnivariatePolynomial<Integer> > > f, g, h;
	SparseUnivariatePolynomial< SparseUnivariatePolynomial< DUZP > > f, g, h;
	f.setCoefficient(0, q);
	g.setCoefficient(0, q);
	f.setCoefficient(1, p);
	f.setVariableName(Symbol("z"));
	g.setVariableName(Symbol("z"));

	//cout << "f := " << f << endl;
	//cout << "g := " << g << endl;

	h = f + g;
	//cout << "h := " << h << endl;
	cout << "addition pass." << endl;

	h = f * g;
	//cout << "h := " << h << endl;
	cout << "multiplication pass." << endl;

	h /= f;
	//cout << "h := " << h << endl;
	cout << "exact division pass." << endl;
}
