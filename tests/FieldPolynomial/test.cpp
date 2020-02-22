#define MontSmallPrimeField 0
#include <bpas.h>
#include <time.h>
#include <gmp.h>
#include <typeinfo>


using namespace std;
#if (MontSmallPrimeField)
	long long int SmallPrimeField::R = 4294967296;
	long long int SmallPrimeField::Pp = 883949567;
#endif

template <class Field> class DenseUnivariatePolynomial;

extern void testUQPoperations();
extern void printPoly(DenseUnivariatePolynomial<SmallPrimeField> &p);
template <typename Field> void generatRandPoly(DenseUnivariatePolynomial<Field> &f, int n, int bitcnt, bool monic);
template <typename Field> void generatRandPolySmall(DenseUnivariatePolynomial<Field> &f, int n, bool monic);
template <typename Field> void divisionTest(DenseUnivariatePolynomial<SmallPrimeField>  &a, DenseUnivariatePolynomial<SmallPrimeField>  &b);

int main(int argc, char *argv[]) {
	sfixn p = 17;//469762049; // 4179340454199820289;
	int var = 2;
	if (argc > 1) { var = atoi(argv[1]); }

	srand (time(NULL));

	testUQPoperations();

	return 0;
	}

	/**
	* The following function generate a polynomail type DenseUnivariatePolynomial with long long int coefficient
	* f: DenseUnivariatePolynomial polynomail
	* n: number of terms,
	* bitcnt: bit count, size of the coefficient randomly genrated by gmp random number generator
	* monic: true if leading coefficient is 1
	**/
template <typename Field>
void generatRandPolySmall(DenseUnivariatePolynomial<Field> &f, int n, int bitcnt, bool monic){
	mpz_t rand_Num;
	unsigned long int i, seed;
	gmp_randstate_t r_state;

	seed = time(NULL);

	gmp_randinit_default (r_state);
	gmp_randseed_ui(r_state, seed);

	mpz_init(rand_Num);

	for(i = 0; i < n-2; ++i) {
		mpz_urandomb(rand_Num,r_state,bitcnt);
		long long int res;
		mpz_export(&res, 0, -1, sizeof res, 0, 0, rand_Num);
		f.setCoefficient(i, res);
		//gmp_printf("%Zd\n", rand_Num);
	}

	if(monic)
		f.setCoefficient(n-1, 1);
	else{
		mpz_urandomb(rand_Num,r_state,bitcnt);
		long long int res;
		mpz_export(&res, 0, -1, sizeof res, 0, 0, rand_Num);
		f.setCoefficient(n-1, res);
	}

	gmp_randclear(r_state);
	mpz_clear(rand_Num);
}

/**
* The following function generate a polynomail type DenseUnivariatePolynomial with mpz_class coefficient
* f: DenseUnivariatePolynomial polynomail
* n: number of terms,
* bitcnt: bit count, size of the coefficient randomly genrated by gmp random number generator
* monic: true if leading coefficient is 1
**/
template <typename Field>
void generatRandPoly(DenseUnivariatePolynomial<Field> &f, int n, int bitcnt, bool monic){
	mpz_t rand_Num;
	unsigned long int i;
	gmp_randstate_t r_state;

	gmp_randinit_default (r_state);
	gmp_randseed_ui(r_state, time(NULL));

	mpz_init(rand_Num);

	//std::cout << "not filled " << f << std::endl;

	for(i = 0; i < n-2; ++i) {
		//std::cout << i << "  "<< mpz_class(rand_Num) << std::endl;
		mpz_urandomb(rand_Num,r_state,bitcnt);
		f.setCoefficient(i, mpz_class(rand_Num));
		//gmp_printf("%Zd\n", rand_Num);
	}

	if(monic){
		f.setCoefficient(n-1, mpz_class(1));
	}
	else{
		mpz_urandomb(rand_Num,r_state,bitcnt);
		f.setCoefficient(n-1, mpz_class(rand_Num));
	}

	gmp_randclear(r_state);
	mpz_clear(rand_Num);
}

/**
*	The followning function perform a division test on two given polynomails
*
**/
template<class Field>
void divisionTest(DenseUnivariatePolynomial<Field>  &a, DenseUnivariatePolynomial<Field>  &b){
	//cout << "a ---> " << a << endl;
	//cout << "b ---> " << b << endl;

	//DenseUnivariatePolynomial<Field>  bb(b.degree()+1);
	string variable("x");
	//bb = b.Reverse();

	//cout << "Testing Reverse() call on b : " << bb << endl;

	//DenseUnivariatePolynomial<Field>  bbb(bb.degree()+1);
	//bbb = bb;
	//cout << "bbb(iter is destructive) copy of bb copy of b.Reverse : " << bbb << endl;
	//bbb.NewtonIterationInversion(/*a.degree()bb.degree()*/a.Reverse().degree()-bb.degree()+1);

	//cout << "Testing Newton Iteration call on b: "  << bbb << endl;

	cout << "division to calculate q" << endl;
	DenseUnivariatePolynomial<Field> q(a.degree());
	q.setVariableName(variable);
	DenseUnivariatePolynomial<Field> r(b.degree());
	r.setVariableName(variable);
	//timing

	int start_s=clock();
	a.NewtonDivisionQuotient(b, q,  r, a.Reverse().degree()-b.Reverse().degree()+1);
	int stop_s=clock();
	//cout << "SPF: q = " << q << endl;
	//cout << "SPF: r = " << r << endl;
	DenseUnivariatePolynomial<Field> ret(a.degree());
	ret = (b*q) + r;
	//cout << "r degree " << rtest.degree() << endl;
	//cout << "test2b degree " << test2b.degree() << endl;
	if((a==ret)  && (r.degree() < b.degree()))
		cout << "Division PASS " << endl;
	else
		cout << "Division Fail " << endl;

	cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec."<< endl;
	

}


void testUQPoperations() {
	//SmallPrimeField::setPrime(65537);
	//SmallPrimeField::setPrime(257);
	//mpz_class c;
	//c = "52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217";
	//BigPrimeField::setPrime(c);
	//BigPrimeField::setPrime(mpz_class("257", 10));
	//BigPrimeField::setPrime(mpz_class("883949569"));
	//SmallPrimeField::setPrime(2147483647);
	//SmallPrimeField::setPrime(7);
	string variable("x");

	//SmallPrimeField a;
	/*
	cout << "+++++++++++Testing Fast Division on Small Prime Field +++++++++++++++" << endl;
	//5x^5 + 4x^4 + 3x3 + 2x2 + x
	DenseUnivariatePolynomial<SmallPrimeField>  fa(20);
	fa.setVariableName(variable);
	a = 1;
	fa.setCoefficient(1, a);
	a = 2;
	fa.setCoefficient(2, a);
	a = 3;
	fa.setCoefficient(3, a);
	a = 4;
	fa.setCoefficient(4, a);
	a = 5;
	fa.setCoefficient(5, a);
	//cout << "fa * fa  =  " << fa*fa << endl;
	cout << "fa  =  " << fa << endl;


	//x2 + 2x + 3
	DenseUnivariatePolynomial<SmallPrimeField>  fb(4);
	fb.setVariableName(variable);
	a = -3;
	fb.setCoefficient(0, a);
	a = 2;
	fb.setCoefficient(1, a);
	a = 1;
	fb.setCoefficient(2, a);
	std::cout << fb << std::endl;


	divisionTest(fa, fb);

	//1+201*x+208*x^2+50*x^3+184*x^4
	DenseUnivariatePolynomial<GeneralizedFermatPrimeField>  fb(5);
	fb.setVariableName(variable);

	fb.setCoefficient(4, -184);
	fb.setCoefficient(3, -50);
	fb.setCoefficient(2, -208);
	fb.setCoefficient(1, -300);
	fb.setCoefficient(0, 1);

	std::cout << fb << std::endl;

	//a = 151+73*x+139*x^2+101*x^3+232*x^4+66*x^5+139*x^6+255*x^7+99*x^9
	//b = 151+73*x+139*x^2+x^4
	DenseUnivariatePolynomial<BigPrimeField>  biga(10);
	biga.setVariableName(variable);

	biga.setCoefficient(9, 99);
	biga.setCoefficient(8, 0);
	biga.setCoefficient(7, 255);
	biga.setCoefficient(6, 139);
	biga.setCoefficient(5, 66);
	biga.setCoefficient(4, 232);
	biga.setCoefficient(3, 101);
	biga.setCoefficient(2, 139);
	biga.setCoefficient(1, 73);
	biga.setCoefficient(0, 151);

	DenseUnivariatePolynomial<BigPrimeField>  bigb(5);
	bigb.setVariableName(variable);


	bigb.setCoefficient(4, 1);
	bigb.setCoefficient(3, 0);
	bigb.setCoefficient(2, 139);
	bigb.setCoefficient(1, 73);
	bigb.setCoefficient(0, 151);

	divisionTest(biga, bigb);

	DenseUnivariatePolynomial<SmallPrimeField>  smalla(10);
	smalla.setVariableName(variable);

	smalla.setCoefficient(9, 99);
	smalla.setCoefficient(8, 0);
	smalla.setCoefficient(7, 255);
	smalla.setCoefficient(6, 139);
	smalla.setCoefficient(5, 66);
	smalla.setCoefficient(4, 232);
	smalla.setCoefficient(3, 101);
	smalla.setCoefficient(2, 139);
	smalla.setCoefficient(1, 73);
	smalla.setCoefficient(0, 151);

	DenseUnivariatePolynomial<SmallPrimeField>  smallb(5);
	smallb.setVariableName(variable);


	smallb.setCoefficient(4, 1);
	smallb.setCoefficient(3, 0);
	smallb.setCoefficient(2, 139);
	smallb.setCoefficient(1, 73);
	smallb.setCoefficient(0, 151);

	//divisionTest(smalla, smallb);

	DenseUnivariatePolynomial<SmallPrimeField> test2a(100);
	test2a.setVariableName(variable);
	generatRandPolySmall(test2a, 100);
	//cout << "a : ------------->" << test2a << endl;
	//SmallPrimeField::setPrime(7);
	DenseUnivariatePolynomial<SmallPrimeField> test2b(50);
	test2b.setVariableName(variable);
	generatRandPolySmall(test2b, 50);
	//cout << "b : ------------->" << test2b << endl;


	//divisionTest(test2a, test2b);


	DenseUnivariatePolynomial<SmallPrimeField> test3a(10000);
	test3a.setVariableName(variable);
	generatRandPolySmall(test3a, 10000, 10, false);
	//cout << "a : ------------->" << test3a << endl;

	//SmallPrimeField::setPrime(7);
	DenseUnivariatePolynomial<SmallPrimeField> test3b(5000);
	test3b.setVariableName(variable);
	generatRandPolySmall(test3b, 5000, 10, true);
	//cout << "b : ------------->" << test2b << endl;
	divisionTest(test3a, test3b);
	*/

	DenseUnivariatePolynomial<BigPrimeField> test3c(100000);
	test3c.setVariableName(variable);
	generatRandPoly(test3c, 100000, 100, false);
	//cout << "a : ------------->" << test3c << endl;
	DenseUnivariatePolynomial<BigPrimeField> test3d(5000);
	test3d.setVariableName(variable);
	generatRandPoly(test3d, 5000, 100, true);
	//cout << "b : ------------->" << test2b << endl;


	divisionTest(test3c, test3d);


}
