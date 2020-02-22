#include <bpas.h>
#include "../../include/ring.h"
#include "../../include/RingPolynomial/upolynomial.h"
#include "../../include/RationalFunction/urationalfunction.h"
#include "../../include/Ring/Fraction.hpp"
#include "../../include/Ring/Integer.hpp"
#include "../../include/Ring/SmartFraction.hpp"

#include "../../include/Ring/RationalNumber.hpp"
#include "../../include/Ring/FactorRefinement.hpp"


#include <complex>
#include <mps/mps.h>
#include <mps/mpc.h>
#include <ctime>
#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"

using namespace std;



DUQP randomDUQP(int n, double sparsity, size_t bits, gmp_randclass &rr);
vector<DUQP> randomRF(int degA, int degD, double sparsity, mpz_class vMax, string variable);

extern void testFractionInteger();
extern void testFractionDUQP();
extern void testFactorRefinementInteger();
extern void testFactorRefinementDUQP();
extern void testSmartFractionDUQP();
//extern void experimentSmartFraction();
//extern void experimentSmartFraction2();
//extern void experimentSmartFraction3();
//extern void experimentSmartFraction4();


int main(int argc, char *argv[]) {
	testFractionInteger();

	cout  << "\n---------------------------------" << endl;
	testFactorRefinementInteger();

	testFractionDUQP();
	testFactorRefinementDUQP();
	testSmartFractionDUQP();
	//experimentSmartFraction();	
	//experimentSmartFraction2();
	//experimentSmartFraction3();
	//experimentSmartFraction4();
	return 0;
}

void testFractionInteger() {
	cout << "--test Fraction<Integer> --" << endl;

	int num1 = rand();
	int den1 = rand();
	Integer num_1(num1);
	Integer den_1(den1);
	Fraction<Integer> f1(num_1,den_1);

	int num2 = rand();
	int den2 = rand();
	Integer num_2(num2);
	Integer den_2(den2);
	Fraction<Integer> f2(num_2,den_2);

	RationalNumber r1(num1,den1);
	RationalNumber r2(num2,den2);


	if(f1.numerator()==r1.get_num()&&f1.denominator() == r1.get_den()){
		cout << "Fraction<Integer> constructor test:\tPASS"<<endl;
	}

	Fraction<Integer> f_sum = f1+f2;
	RationalNumber r_sum = r1+r2;

	if(f_sum.numerator()==r_sum.get_num() &&f_sum.denominator() ==r_sum.get_den() ){
		cout << "Fraction<Integer> +operator test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> +operator test:\tFAIL"<<endl;
	}

	Fraction<Integer> f_mul = f1*f2;
	RationalNumber r_mul = r1*r2;

	if(f_mul.numerator()==r_mul.get_num() &&f_mul.denominator() ==r_mul.get_den() ){
		cout << "Fraction<Integer> *operator test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> *operator test:\tFAIL"<<endl;
	}

	Fraction<Integer> f_sub = f2-f1;
	RationalNumber r_sub = r2-r1;
	if(f_sub.numerator()==r_sub.get_num() &&f_sub.denominator() ==r_sub.get_den() ){
		cout << "Fraction<Integer> -operator test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> -operator test:\tFAIL"<<endl;
	}

	Fraction<Integer> f_div = f1/f2;
	RationalNumber r_div = r1/r2;
	if(f_div.numerator()==r_div.get_num() &&f_div.denominator() ==r_div.get_den()){
		cout << "Fraction<Integer> /operator test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> /operator test:\tFAIL"<<endl;
	}

	Fraction<Integer> f_pow = f1^2;
	RationalNumber r_pow = r1^2;

	if(f_pow.numerator()==r_pow.get_num() &&f_pow.denominator() ==r_pow.get_den()){
		cout << "Fraction<Integer> ^operator test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> ^operator test:\tFAIL"<<endl;
	}

	Fraction<Integer> f_inv = f1.inverse();
	RationalNumber r_inv = r1.inverse();
	if(f_inv.numerator()==r_inv.get_num() &&f_inv.denominator() ==r_inv.get_den() ){
		cout << "Fraction<Integer> inverse test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> inverse test:\tFAIL"<<endl;
	}

	int num3 = 0;
	int den3 = 79;
	Integer num_3(num3);
	Integer den_3(den3);
	Fraction<Integer> f3(num_3,den_3);
	if(f3.isZero()==true){
		cout << "Fraction<Integer> isZero test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> isZero test:\tFAIL"<<endl;
	}
	f3.one();
	if(f3.isOne()==true){
		cout << "Fraction<Integer> isOne test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> isOne test:\tFAIL"<<endl;
	}

	f3.zero();
	if (f3.numerator().get_ui() == 0)
	{
		cout << "Fraction<Integer> zero test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> zero test:\tFAIL"<<endl;
	}
	f3.one();
	if (f3.numerator().get_ui() == 1&&f3.denominator().get_ui())
	{
		cout << "Fraction<Integer> one test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> one test:\tFAIL"<<endl;
	}

	int num4 = 9;
	int den4 = 3;
	Integer num_4(num4);
	Integer den_4(den4);
	Fraction<Integer> f4(num_4,den_4);

	Fraction<Integer> f_unitCanon = f4.unitCanonical();

	if(f_unitCanon.numerator().get_ui()==3&&f_unitCanon.denominator().get_ui()==1)
	{
		cout << "Fraction<Integer> unitCanonical test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> unitCanonical test:\tFAIL"<<endl;
	}

	Fraction<Integer> f_gcd = f1.gcd(f2);
	if(f_gcd.isOne()==true){
		cout << "Fraction<Integer> gcd test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> gcd test:\tFAIL"<<endl;
	}

	int num5 = 9;
	int den5 = 3;
	Integer num_5(num5);
	Integer den_5(den5);
	Fraction<Integer> f5(num5,den5);

	Fraction<Integer> f_euclideandiv = f1.euclideanDivision(f2);

	if (f_euclideandiv.isZero())
	{
		cout << "Fraction<Integer> euclideanDivision test:\tPASS"<<endl;
	}
	else{
		cout << "Fraction<Integer> euclideanDivision test:\tFAIL"<<endl;
	}
	return;
}






void testFactorRefinementInteger(){

	//SmartFraction<Integer> sf1;
	//Fraction<Integer> f1 = convertToStandardFraction(sf1);

	cout << "---test Factor Refinement: Algorithm1 PolyRefine---" <<endl;
	int inputa = 22;
	int inputb = 26;
	Integer a(inputa);
	int e = 1;
	Integer b(inputb);
	int f = 1;
	Integer ret_l1;
	int ret_e1;
	std::vector<Factor<Integer>> ret_G1;
	Integer ret_r1;
	int ret_f1;
	FactorRefinement::PolyRefine<Integer>(a,e,b,f,&ret_l1,&ret_e1,&ret_G1,&ret_r1,&ret_f1);
	cout << "input: " << inputa << " " << inputb << endl;
	cout << "output: "; 
	cout <<ret_l1.get_ui() << "^" << ret_e1 <<"  ";
	for (int i = 0; i < ret_G1.size(); ++i)
	{
		cout << ret_G1[i].first.get_ui() << "^" << ret_G1[i].second <<"  ";
	}
	cout <<ret_r1.get_ui() << "^" << ret_f1  <<endl;


	cout << "---test Factor Refinement: Algorithm2 MergeRefinePolySeq---" <<endl;

	int inputa2 = 5;
	Integer a2(inputa2);
	int inputb2 = 7;
	Integer b2(inputb2);	
	int inputc2 = 3;
	Integer c2(inputc2);
	int one_2 = 1;
	std::vector<Factor<Integer>> B2;
	
	B2.push_back(std::make_pair(b2,1));
	B2.push_back(std::make_pair(c2,1));

	cout << "input: " << inputa2 <<"^1 " <<inputb2 << "^1 " << inputc2 << "^1 " <<endl;
 

	Integer ret_l2;
	int ret_m2;
	std::vector<Factor<Integer>> ret_Q2;
	std::vector<Factor<Integer>> ret_S2;


	FactorRefinement::MergeRefinePolySeq<Integer>(a2,1,B2,&ret_l2,&ret_m2,&ret_Q2,&ret_S2);

	cout <<"output: ";

	cout << ret_l2.get_ui() << "^" << ret_m2 << "  ";
	for (int i = 0; i < ret_Q2.size(); ++i)
	{
		cout << ret_Q2[i].first.get_ui() << "^" << ret_Q2[i].second <<"  ";
	}
	for (int i = 0; i < ret_S2.size(); ++i)
	{
		cout << ret_S2[i].first.get_ui() << "^" << ret_S2[i].second <<"  ";
	}


	cout << endl;
	cout << "---test Factor Refinement: Algorithm3 MergeRefineTwoSeq---" <<endl;


	int inputa3 = 3;
	Integer a3(inputa3);
	int inputb3 = 2;
	Integer b3(inputb3);

	int inputc3 = 7;
	Integer c3(inputc3);
	int inputd3 = 10;
	Integer d3(inputd3);

	std::vector<Factor<Integer>> A3;

	std::vector<Factor<Integer>> B3;

	A3.push_back(std::make_pair(a3,1));
	A3.push_back(std::make_pair(b3,2));

	B3.push_back(std::make_pair(c3,1));
	B3.push_back(std::make_pair(d3,1));

	std::vector<Factor<Integer>> L;
	std::vector<Factor<Integer>> Q;
	std::vector<Factor<Integer>> S;

	cout << "input: " <<endl;

	cout << inputa3 << "^1 " << inputb3 << "^2 " << inputc3 << "^1 " << inputd3 << "^1" <<endl; 


	FactorRefinement::MergeRefineTwoSeq<Integer>(A3,B3,&L,&Q,&S);

	cout << "Output: " << endl;

	for (int i = 0; i < L.size(); ++i)
	{
		cout << L[i].first.get_ui() << "^" << L[i].second <<"  ";
	}


	for (int i = 0; i < Q.size(); ++i)
	{
		cout << Q[i].first.get_ui() << "^" << Q[i].second <<"  ";
	}

	for (int i = 0; i < S.size(); ++i)
	{
		cout << S[i].first.get_ui() << "^" << S[i].second <<"  ";
	}
	cout << endl;


	DUQP num;



	return;
}



DUQP randomDUQP(int n, double sparsity, size_t bits, gmp_randclass &rr){
	int k;
	mpz_class elem;
	DUQP P(n+1);
	for(int i = 0; i < n; i++) {
		elem = rr.get_z_bits(bits);
		k = (int) rand() % 2;
		// Positive or negative coefficients
		if (k)
			elem = -elem;
		// Set random coefficients with sparsity
		if (rand() > sparsity)
			P.setCoefficient(i, elem);
		else
			P.setCoefficient(i, 0);
	}
	// ensure P.degree() == n
	elem = rr.get_z_bits(bits);
	k = (int) rand() % 2;
	// Positive or negative coefficients
	if (k)
		elem = -elem;
	// Set random coefficients with sparsity
	if (abs(elem) > 1)
		P.setCoefficient(n, elem);
	else
		P.setCoefficient(n, 1);
	
	return(P);
}

vector<DUQP> randomRF(int degA, int degD, double sparsity, mpz_class vMax, Symbol variable){
	
	assert(degA < degD);
	assert((0 <= sparsity) && (sparsity < 1));
	assert(0 < vMax);
	
	DUQP one(2);
	one.one();
	size_t bits = mpz_sizeinbase(vMax.get_mpz_t(), 2);
	gmp_randclass rr (gmp_randinit_default);
	rr.seed(time(NULL));
	srand (time(NULL));
	
	DUQP A;
	DUQP D;
	vector<DUQP> RF;
	
	// Generate random numerator
	A = randomDUQP(degA,sparsity,bits,rr);
	// Attempt to generate coprime random denominator
	for (int i=0; i<10; i++){
		D = randomDUQP(degD,sparsity,bits,rr);
		if (A.gcd(D) == one)
			break;
	}
	if (!(A.gcd(D) == one)){
		mpz_class elem;
		elem = rr.get_z_bits(bits);
		if (elem == 0)
			D.setCoefficient(0,1);
		else
			D.setCoefficient(0,elem);
		if (!(A.gcd(D) == one))
			cout << "Failed to create relatively prime numerator and denominator" << endl;
	}
	
	A.setVariableName(variable);
	D.setVariableName(variable);
	
	RF.push_back(A);
	RF.push_back(D);
	
	return(RF);
}




void testFractionDUQP(){
	//vector<DUQP> randomRF(int degA, int degD, double sparsity, mpz_class vMax, string variable){
	vector<DUQP> rf = randomRF(5,10,0.5,7,'x');
	vector<DUQP> rf2 = randomRF(5,8,0.5,7,'x');
	Fraction<DUQP> f1(rf[0],rf[1]);
	Fraction<DUQP> f2(rf2[0],rf2[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> r1(rf[0],rf[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> r2(rf2[0],rf2[1]);

	if(f1.numerator() == r1.numerator()&&f1.denominator() == r1.denominator()){
		cout <<"Fraction<DUQP> constructor test:\tPASS" <<endl;
	}else{
		cout <<"Fraction<DUQP> constructor test:\tFAIL" <<endl;
	}
	Fraction<DUQP> f_sum = f1+f2;
	UnivariateRationalFunction<DUQP,RationalNumber> r_sum = r1+r2;
	if(f_sum.numerator() == r_sum.numerator()&&f_sum.denominator() == r_sum.denominator()){
		cout <<"Fraction<DUQP> +operator test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> +operator test:\tFAIL" << endl;
		exit(1);
	}

	Fraction<DUQP> f_sub = f1-f2;
	UnivariateRationalFunction<DUQP,RationalNumber> r_sub = r1-r2;
	if(f_sub.numerator() == r_sub.numerator()&&f_sub.denominator() == r_sub.denominator()){
		cout <<"Fraction<DUQP> -operator test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> -operator test:\tFAIL" << endl;
		exit(1);
	}


	Fraction<DUQP> f_mul = f1*f2;
	UnivariateRationalFunction<DUQP,RationalNumber> r_mul = r1*r2;
	//f_mul.print(std::cout);
	//r_mul.print(std::cout);
	if(f_mul.numerator() == r_mul.numerator()&&f_mul.denominator() == r_mul.denominator()){
		cout <<"Fraction<DUQP> *operator test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> *operator test:\tFAIL" << endl;
		exit(1);
	}


	Fraction<DUQP> f_div = f1/f2;
	UnivariateRationalFunction<DUQP,RationalNumber> r_div = r1/r2;
	if(f_div.numerator() == r_div.numerator()&&f_div.denominator() == r_div.denominator()){
		cout <<"Fraction<DUQP> /operator test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> /operator test:\tFAIL" << endl;
		exit(1);
	}


	Fraction<DUQP> f_pow = f1^2;
	cout << "Fraction exponential fine." << endl;
	UnivariateRationalFunction<DUQP,RationalNumber> r_pow = r1^2;
	if(f_pow.numerator() == r_pow.numerator()&&f_pow.denominator() == r_pow.denominator()){
		cout <<"Fraction<DUQP> /operator test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> /operator test:\tFAIL" << endl;
		exit(1);
	}

	Fraction<DUQP> f_inv = f1.inverse();
	UnivariateRationalFunction<DUQP,RationalNumber> r_inv = r1.inverse();
	if(f_inv.numerator() == r_inv.numerator()&&f_inv.denominator() == r_inv.denominator()){
		cout <<"Fraction<DUQP> inverse test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> inverse test:\tFAIL" << endl;
		exit(1);
	}

	f1.one();
	if(f1.isOne()){
		cout <<"Fraction<DUQP> one test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> one test:\tFAIL" << endl;
		exit(1);
	}

	f1.zero();
	if(f1.isZero()){
		cout <<"Fraction<DUQP> zero test:\tPASS" <<endl;
	}
	else{
		cout << "Fraction<DUQP> zero test:\tFAIL" << endl;
		exit(1);
	}


	return;
}

void testFactorRefinementDUQP(){
	cout << "------test FactorRefinement<DUQP>---------" << endl;
	vector<DUQP> rf1 = randomRF(3,4,0.5,3,'x');
	vector<DUQP> rf2 = randomRF(2,6,0.5,3,'x');

	vector<Factor<DUQP>> v1;
	v1.push_back(make_pair(rf1[0],1));
	v1.push_back(make_pair(rf1[1],1));

	vector<Factor<DUQP>> v2;

	DUQP one;
	one.one();
	v2.push_back(make_pair(rf1[0],1));
	v2.push_back(make_pair(rf1[1],1));
	//v2.push_back(make_pair(one,1));

	cout << "---test PolyRefine<DUQP>--------- " << endl;

	cout << "**input** " << endl;
	cout << rf1[0] << endl;
	cout << convertToDomain(v2) << endl;

	cout << "**output**" << endl;
	DUQP l;
	int m;
	vector<Factor<DUQP>> Q,S;
	FactorRefinement::MergeRefinePolySeq<DUQP>(rf1[0],1,v2,&l,&m,&Q,&S);
	cout <<"l :" <<l << endl;
	cout <<"m :" <<m << endl;
	cout <<"Q :" << convertToDomain(Q) << endl;
	cout <<"S :" << convertToDomain(S)  << endl;

	cout << "---end of test PolyRefine<DUQP>--------- " << endl;


	for (int i = 0; i < v1.size(); ++i)
	{
		cout <<"("<< v1[i].first << ")" ;
	}

	cout << endl;

	for (int i = 0; i < v2.size(); ++i)
	{
		cout <<"("<< v2[i].first << ")";
	}

	cout << endl;
	cout<< "---"<<endl;

	vector<Factor<DUQP>> ret1;
	vector<Factor<DUQP>> ret2;
	vector<Factor<DUQP>> ret3;

	FactorRefinement::MergeRefineTwoSeq<DUQP>(v1,v2,&ret1,&ret2,&ret3);

cout << "output:" <<endl;

cout << "size:" << ret1.size() << " " << ret2.size() << " " << ret3.size() << endl;
	cout << "ret1" << endl;
	for (int i = 0; i < ret1.size(); ++i)
	{
		cout << "(";
		ret1[i].first.print(std::cout);
		cout << ")";
		cout << "^";
		cout << ret1[i].second << endl;
	}

	cout << "ret2" << endl;
	for (int i = 0; i < ret2.size(); ++i)
	{
		cout << "(";
		ret2[i].first.print(std::cout);
		cout << ")";
		cout << "^";
		cout<< ret2[i].second << endl;
	}
	cout << "ret3" << endl;
		for (int i = 0; i < ret3.size(); ++i)
	{
		cout << "(";
		ret3[i].first.print(std::cout);
		cout << ")";
		cout << "^";
		cout<< ret3[i].second << endl;

	}

cout<< endl;




	return;
}


void testSmartFractionDUQP(){

	cout << "------test SmartFraction<DUQP>------" << endl;

	vector<DUQP> rf = randomRF(9,10,0.5,7,'x');
	vector<DUQP> rf2 = randomRF(6,7,0.5,4,'x');

	SmartFraction<DUQP> f1(rf[0],rf[1]);

	SmartFraction<DUQP> f2(rf2[0],rf2[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> r1(rf[0],rf[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> r2(rf2[0],rf2[1]);


	SmartFraction<DUQP> f3;
	if(f3.isZero()==true){
		cout << "SmartFraction<DUQP> isZero test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> isZero stest:\tFAIL" << endl;
		exit(1);
	}
	f3 = f1;
	f3.zero();
	if(f3.isZero()==true){
		cout << "SmartFraction<DUQP> zero test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> zero stest:\tFAIL" << endl;
		exit(1);
	}

	if(f3.isOne()==false){
		cout << "SmartFraction<DUQP> isOne test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> isOne test:\tFAIL" << endl;
		exit(1);
	}
	f3.one();


	if(f3.isOne()==true){
		cout << "SmartFraction<DUQP> one test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> one test:\tFAIL" << endl;
		exit(1);
	}

	SmartFraction<DUQP> mul_f = f1*f2;
	UnivariateRationalFunction<DUQP,RationalNumber> mul_r = r1*r2;

	if(mul_f.numerator() == mul_r.numerator()&&mul_f.denominator()==mul_r.denominator()){
		cout << "SmartFraction<DUQP> operator* test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> operator* test:\tFAIL" << endl;
		exit(1);
	}

	SmartFraction<DUQP> sum_f = f1+f2;

	UnivariateRationalFunction<DUQP,RationalNumber> sum_r = r1+r2;

	if(sum_f.numerator() == sum_r.numerator()&&sum_f.denominator()==sum_r.denominator()){
		cout << "SmartFraction<DUQP> operator+ test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> operator+ test:\tFAIL" << endl;
		exit(1);
	}
	SmartFraction<DUQP> sub_f = f1-f2;
	UnivariateRationalFunction<DUQP,RationalNumber> sub_r = r1-r2;

	if(sub_f.numerator() == sub_r.numerator()&&sub_f.denominator()==sub_r.denominator()){
		cout << "SmartFraction<DUQP> operator- test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> operator- test:\tFAIL" << endl;
		exit(1);
	}


	SmartFraction<DUQP> div_f = f1/f2;
	UnivariateRationalFunction<DUQP,RationalNumber> div_r = r1/r2;
	if(div_f.numerator() == div_r.numerator()&&div_f.denominator()==div_r.denominator()){
		cout << "SmartFraction<DUQP> operator/ test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> operator/ test:\tFAIL" << endl;
		exit(1);
	}

	SmartFraction<DUQP> inv_f = f1.inverse();
	UnivariateRationalFunction<DUQP,RationalNumber> inv_r = r1.inverse();

	//cout << "inv_f";
	//inv_f.print(cout);
	//cout << "inv_r" << inv_r << endl;

	if(inv_f.numerator() == inv_r.numerator()&&inv_f.denominator()==inv_r.denominator()){
		cout << "SmartFraction<DUQP> inverse test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> inverse test:\tFAIL" << endl;
		exit(1);
	}


	SmartFraction<DUQP> neg_f = -f1;
	UnivariateRationalFunction<DUQP,RationalNumber> neg_r = -r1;
	if(neg_f.numerator() == neg_r.numerator()&&neg_f.denominator()==neg_r.denominator()){
		cout << "SmartFraction<DUQP> operator-(neg) test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> operator-(neg) test:\tFAIL" << endl;
		exit(1);
	}
	SmartFraction<DUQP> f_copy = f1;

	if(f1 == f_copy&&!(f1==f2)){
		cout << "SmartFraction<DUQP> operator== test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> operator== test:\tFAIL" << endl;
		exit(1);
	}

	if(!(f1!=f_copy)&&(f1!=f2)){
		cout << "SmartFraction<DUQP> operator!= test:\tPASS" << endl;
	}
	else{
		cout << "SmartFraction<DUQP> operator!= test:\tFAIL" << endl;
		exit(1);
	}
return;
}

vector<DUQP> generateA_i(int i){
	vector<DUQP> ret;
	DUQP p1(11),p2(6),p3(21),p4(4);
	p1.setVariableName(Symbol('x'));
	p2.setVariableName(Symbol('x'));
	p3.setVariableName(Symbol('x'));
	p4.setVariableName(Symbol('x'));

	RN i_(i);
	RN one(1);
	RN i_1 = i_ + one;
	i_1 = i_1^10;

	p1.setCoefficient(0,i_1);

	i_1 = i_ + one;
	p1.setCoefficient(10,RN(1));
	p2.setCoefficient(5,RN(1));
	p2.setCoefficient(0,i_1);
	RN i_20(i^20);
	p3.setCoefficient(0,i_20);
	p3.setCoefficient(20,RN(1));

	p4.setCoefficient(3,RN(1));
	p4.setCoefficient(0,i_);

	p4 = p4*p4;

	p1*=p2;
	p3*=p4;

	ret.push_back(p1);
	ret.push_back(p3);

	return ret;
}

vector<DUQP> generateB_n(int n){
	vector<DUQP> A = generateA_i(1);
	DUQP first = A[0];
	DUQP second = A[1];
	for (int i = 2; i <= n; ++i)
	{
		vector<DUQP> tempA = generateA_i(i);
		first *= tempA[0];
		second *= tempA[1];
	}

	A[0] = first;
	A[1] = second;
	return A;
}



void experimentSmartFraction(){

	int input = 7;
	vector<DUQP> v1 = generateB_n(input);
	vector<DUQP> v2 = generateB_n(input);

	clock_t beginf1 = clock();
	SmartFraction<DUQP> f1(v1[0],v1[1]);
	SmartFraction<DUQP> f2(v2[0],v2[1]);
	SmartFraction<DUQP> f3 = f1+f2;
	clock_t endf1 = clock();

	clock_t beginr1 = clock();
	UnivariateRationalFunction<DUQP,RationalNumber> r1(v1[0],v1[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> r2(v2[0],v2[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> r3 = r1+r2;
  	clock_t endr1 = clock();

  	double elapsed_secs_SmartFraction = double(endf1 - beginf1) / CLOCKS_PER_SEC;
  	double elapsed_secs_RationNumber = double(endr1 - beginr1) / CLOCKS_PER_SEC;
	cout << "elapsed_secs_SmartFraction " << elapsed_secs_SmartFraction << endl;
	cout << "elapsed_secs_RationNumber " << elapsed_secs_RationNumber << endl;

	return;
}

vector<DUQP> generate_A2(int i){
	vector<DUQP> ret;
	DUQP p1(3),p2(2);

	p1.setVariableName(Symbol('x'));
	p2.setVariableName(Symbol('x'));
	RN i_(i);
	RN one(1);
	RN isqure2 = -(i_^2);
	RN i_plus1 = i_ + one;

	p1.setCoefficient(2,RN(1));
	p1.setCoefficient(0,isqure2);
	p2.setCoefficient(1,RN(1));
	p2.setCoefficient(0,i_plus1);
	p2 = p2^3;
	ret.push_back(p1);
	ret.push_back(p2);

	return ret;
}

vector<DUQP> generate_B2(int i){

	vector<DUQP> ret;
	DUQP p1(4),p2(1);

	p1.setVariableName(Symbol('x'));
	p2.setVariableName(Symbol('x'));

	RN i_(i);
	RN one(1);

	RN isqure3 = i_^3;

	p1.setCoefficient(3,one);
	p1.setCoefficient(0,isqure3);

	p2.setCoefficient(0,RN(1));
	ret.push_back(p1);
	ret.push_back(p2);

	return ret;
}


vector<DUQP> generate_C2(int i){

	vector<DUQP> ret;
	DUQP p1(2),p2(1);

	p1.setVariableName(Symbol('x'));
	p2.setVariableName(Symbol('x'));

	RN i_(i);
	RN one(1);

	RN iminus1 = i_-one;

	p1.setCoefficient(1,one);
	p1.setCoefficient(0,iminus1);

	p1 = p1^2;

	p2.setCoefficient(0,RN(1));
	ret.push_back(p1);
	ret.push_back(p2);
	return ret;

}


vector<DUQP> generate_F2(int n){

	vector<DUQP> ret;

	vector<DUQP> An = generate_A2(n);
	vector<DUQP> Bn = generate_B2(n);
	vector<DUQP> Cn = generate_C2(n);

	UnivariateRationalFunction<DUQP,RationalNumber> Anr(An[0],An[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> Bnr(Bn[0],Bn[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> Cnr(Cn[0],Cn[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> E;

	E = (Anr + Bnr)/Cnr;
	int i = n-1;
	while(i>0){
		An = generate_A2(i);
		Bn = generate_B2(i);

		UnivariateRationalFunction<DUQP,RationalNumber> Ai(An[0],An[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> Bi(Bn[0],Bn[1]);

		E = (Ai + Bi) /E; 
		i--;
	}

	ret.push_back(E.numerator());
	ret.push_back(E.denominator());


	return ret;


}





void experimentSmartFraction2(){

	for (int i = 0; i < 10; ++i)
	{
		cout << "n: " << i << endl;

		int input = i*10;
		vector<DUQP> v1 = generate_F2(input);
		vector<DUQP> v2 = generate_F2(input);

		clock_t beginf1 = clock();
		SmartFraction<DUQP> f1(v1[0],v1[1]);
		SmartFraction<DUQP> f2(v2[0],v2[1]);
		SmartFraction<DUQP> f3 = f1+f2;
		clock_t endf1 = clock();

		clock_t beginr1 = clock();
		UnivariateRationalFunction<DUQP,RationalNumber> r1(v1[0],v1[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> r2(v2[0],v2[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> r3 = r1+r2;
	  	clock_t endr1 = clock();

	  	double elapsed_secs_Smart_Fraction = double(endf1 - beginf1) / CLOCKS_PER_SEC;
	  	double elapsed_secs_Ration_Number = double(endr1 - beginr1) / CLOCKS_PER_SEC;
		cout << "elapsed_secs_SmartFraction " << elapsed_secs_Smart_Fraction << endl;
		cout << "elapsed_secs_RationNumber " << elapsed_secs_Ration_Number << endl;
	}

	return;

}



vector<DUQP> generate_A3(int i){
	vector<DUQP> ret;
	DUQP p1(3),p2(2);

	p1.setVariableName('x');
	p2.setVariableName('x');
	RN i_(i);
	RN one(1);
	RN isqure2 = -(i_^2);
	RN i_plus1 = i_ + one;

	p1.setCoefficient(2,RN(1));
	p1.setCoefficient(0,isqure2);

	p1 = p1^3;
	p2.setCoefficient(1,RN(1));
	p2.setCoefficient(0,i_plus1);
	p2 = p2^2;
	ret.push_back(p1);
	ret.push_back(p2);

	return ret;
}

vector<DUQP> generate_B3(int i){

	vector<DUQP> ret;
	DUQP p1(4),p2(1);

	p1.setVariableName('x');
	p2.setVariableName('x');

	RN i_(i);
	RN one(1);

	RN isqure3 = i_^3;

	p1.setCoefficient(3,one);
	p1.setCoefficient(0,isqure3);

	p2.setCoefficient(0,RN(1));
	ret.push_back(p1);
	ret.push_back(p2);

	return ret;
}


vector<DUQP> generate_C3(int i){

	vector<DUQP> ret;
	DUQP p1(2),p2(1);

	p1.setVariableName('x');
	p2.setVariableName('x');

	RN i_(i);
	RN one(1);

	RN iminus1 = i_-one;

	p1.setCoefficient(1,one);
	p1.setCoefficient(0,iminus1);

	p1 = p1^2;

	p2.setCoefficient(0,RN(1));
	ret.push_back(p1);
	ret.push_back(p2);
	return ret;

}


vector<DUQP> generate_F3(int n){

	vector<DUQP> ret;

	vector<DUQP> An = generate_A3(n);
	vector<DUQP> Bn = generate_B3(n);
	vector<DUQP> Cn = generate_C3(n);

	UnivariateRationalFunction<DUQP,RationalNumber> Anr(An[0],An[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> Bnr(Bn[0],Bn[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> Cnr(Cn[0],Cn[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> E;

	E = (Anr + Bnr)/Cnr;
	int i = n-1;
	while(i>0){
		An = generate_A3(i);
		Bn = generate_B3(i);

		UnivariateRationalFunction<DUQP,RationalNumber> Ai(An[0],An[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> Bi(Bn[0],Bn[1]);

		E = (Ai + Bi) /E; 
		i--;
	}

	ret.push_back(E.numerator());
	ret.push_back(E.denominator());


	return ret;


}



void experimentSmartFraction3(){
	for (int i = 0; i < 20; ++i)
	{
		cout << "n: " << i << endl;

		int input = i;
		vector<DUQP> v1 = generate_F3(input);
		vector<DUQP> v2 = generate_F3(input);

		clock_t beginf1 = clock();
		SmartFraction<DUQP> f1(v1[0],v1[1]);
		SmartFraction<DUQP> f2(v2[0],v2[1]);
		SmartFraction<DUQP> f3 = f1+f2;
		clock_t endf1 = clock();

		clock_t beginr1 = clock();
		
		UnivariateRationalFunction<DUQP,RationalNumber> r1(v1[0],v1[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> r2(v2[0],v2[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> r3 = r1+r2;
	  	clock_t endr1 = clock();

	  	double elapsed_secs_Smart_Fraction = double(endf1 - beginf1) / CLOCKS_PER_SEC;
	  	double elapsed_secs_Ration_Number = double(endr1 - beginr1) / CLOCKS_PER_SEC;
		cout << "elapsed_secs_SmartFraction " << elapsed_secs_Smart_Fraction << endl;
		cout << "elapsed_secs_RationNumber " << elapsed_secs_Ration_Number << endl;
	}


	return;
}


vector<DUQP> generate_A4(int i){
	vector<DUQP> ret;
	DUQP p1(3),p2(2);

	p1.setVariableName('x');
	p2.setVariableName('x');
	RN i_(i);
	RN one(1);
	RN isqure2 = -(i_^2);
	RN i_plus1 = i_ + one;

	p1.setCoefficient(2,RN(1));
	p1.setCoefficient(0,isqure2);

	p1 = p1^3;
	p2.setCoefficient(1,RN(1));
	p2.setCoefficient(0,i_plus1);
	p2 = p2^2;
	ret.push_back(p1);
	ret.push_back(p2);

	return ret;
}

vector<DUQP> generate_B4(int i){

	vector<DUQP> ret;
	DUQP p1(1),p2(1);

	p1.setVariableName('x');
	p2.setVariableName('x');

	RN i_(i);
	RN one(1);
	p1.setCoefficient(0,RN(1));
	p2.setCoefficient(0,RN(1));
	ret.push_back(p1);
	ret.push_back(p2);

	return ret;
}


vector<DUQP> generate_C4(int i){

	vector<DUQP> ret;
	DUQP p1(2),p2(1);

	p1.setVariableName('x');
	p2.setVariableName('x');

	RN i_(i);
	RN one(1);

	RN minus1 = -one;

	p1.setCoefficient(1,one);
	p1.setCoefficient(0,minus1);

	p1 = p1^2;

	p2.setCoefficient(0,RN(1));
	ret.push_back(p1);
	ret.push_back(p2);
	return ret;

}


vector<DUQP> generate_F4(int n){

	vector<DUQP> ret;

	vector<DUQP> An = generate_A4(n);
	vector<DUQP> Bn = generate_B4(n);
	vector<DUQP> Cn = generate_C4(n);

	UnivariateRationalFunction<DUQP,RationalNumber> Anr(An[0],An[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> Bnr(Bn[0],Bn[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> Cnr(Cn[0],Cn[1]);
	UnivariateRationalFunction<DUQP,RationalNumber> E;

	E = (Anr + Bnr)/Cnr;
	int i = n-1;
	while(i>0){
		An = generate_A4(i);
		Bn = generate_B4(i);

		UnivariateRationalFunction<DUQP,RationalNumber> Ai(An[0],An[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> Bi(Bn[0],Bn[1]);

		E = (Ai + Bi) /E; 
		i--;
	}

	ret.push_back(E.numerator());
	ret.push_back(E.denominator());


	return ret;


}



void experimentSmartFraction4(){


for (int i = 0; i < 20; ++i)
	{
		cout << "n: " << i << endl;

		int input = i;
		vector<DUQP> v1 = generate_F4(input);
		vector<DUQP> v2 = generate_F4(input);

		clock_t beginf1 = clock();
		SmartFraction<DUQP> f1(v1[0],v1[1]);
		SmartFraction<DUQP> f2(v2[0],v2[1]);
		SmartFraction<DUQP> f3 = f1+f2;
		clock_t endf1 = clock();

		clock_t beginr1 = clock();
		
		UnivariateRationalFunction<DUQP,RationalNumber> r1(v1[0],v1[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> r2(v2[0],v2[1]);
		UnivariateRationalFunction<DUQP,RationalNumber> r3 = r1+r2;
	  	clock_t endr1 = clock();

	  	double elapsed_secs_Smart_Fraction = double(endf1 - beginf1) / CLOCKS_PER_SEC;
	  	double elapsed_secs_Ration_Number = double(endr1 - beginr1) / CLOCKS_PER_SEC;
		cout << "elapsed_secs_SmartFraction " << elapsed_secs_Smart_Fraction << endl;
		cout << "elapsed_secs_RationNumber " << elapsed_secs_Ration_Number << endl;
	}

	return;
}



