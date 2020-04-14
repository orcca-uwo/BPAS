
#include "Ring/Integer.hpp"
#include "Ring/RationalNumber.hpp"
#include "Ring/ComplexRationalNumber.hpp"
#include "FiniteFields/SmallPrimeField.hpp"
#include "FiniteFields/BigPrimeField.hpp"
#include "FiniteFields/GeneralizedFermatPrimeField.hpp"
#include "IntegerPolynomial/uzpolynomial.h"
#include "RationalNumberPolynomial/urpolynomial.h"
#include "RingPolynomial/upolynomial.h"

#include <iostream>

using std::endl;
using std::cout;

//mpz_class GeneralizedFermatPrimeField::prime ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
mpz_class GeneralizedFermatPrimeField::characteristic ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);


// bool GeneralizedFermatPrimeField::isPrimeField = 1;
// bool GeneralizedFermatPrimeField::isSmallPrimeField = 0;
// bool GeneralizedFermatPrimeField::isComplexField = 0;
mpz_class& GeneralizedFermatPrimeField::prime = GeneralizedFermatPrimeField::characteristic;
unsigned long long int GeneralizedFermatPrimeField::r = 9223372054034644992ULL;
int GeneralizedFermatPrimeField::k = 8;

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField () {
    x = new usfixn64 [k]();
    int i;
    for (i = 0; i < k; i++) {
      x[i] = 0;
    }
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (mpz_class a) {
	setX(a);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (int a) {
	mpz_class a_mpz (a);
	setX (a_mpz);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const GeneralizedFermatPrimeField& c) {
	x = new usfixn64 [k]();
	// for (int i = 0; i < k; i++) {
	// 	x[i] = c.x[i];
	// }
	memcpy(x,c.x,(k) * sizeof(usfixn64));
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const Integer& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const RationalNumber& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const SmallPrimeField& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const BigPrimeField& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const DenseUnivariateIntegerPolynomial& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const DenseUnivariateRationalPolynomial& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<Integer>& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<RationalNumber>& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<ComplexRationalNumber>& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

template <class Ring>
GeneralizedFermatPrimeField::GeneralizedFermatPrimeField (const SparseUnivariatePolynomial<Ring>& c) {
	std::cerr << "Cannot convert input to GeneralizedFermatPrimeField! " << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField* GeneralizedFermatPrimeField::GPFpointer(GeneralizedFermatPrimeField* a) {
	return a;
}

GeneralizedFermatPrimeField* GeneralizedFermatPrimeField::GPFpointer(RationalNumber* a) {
	std::cout << "BPAS error, try to cast pointer to Rational Number to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField* GeneralizedFermatPrimeField::GPFpointer(SmallPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to BigPrimeField to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

GeneralizedFermatPrimeField* GeneralizedFermatPrimeField::GPFpointer(BigPrimeField* a) {
	std::cout << "BPAS error, try to cast pointer to GeneralizedFermatPrimeField to pointer to SmallPrimeField" << std::endl;
	exit(1);
}

void GeneralizedFermatPrimeField::setX (mpz_class a){
	x = new usfixn64 [k]();
	if (a >= prime){
		a = a % prime;
	}
	if ( a < 0){
		a = prime + a;
	}

	if (a == ( prime - 1)){

		x[k-1] = r;
// cout << "Object is prime - 1"<< endl;
	}
	else if (a == 0){

	}
	else{
			mpz_t q, r1;
			mpz_init_set_ui(q, 0);
			mpz_init_set_ui(r1, 0);

			mpz_set(q, a.get_mpz_t());
			for (int i = 0; i < k; i++){
		//      mpz_fdiv_qr (q, r, bigint, base);
				mpz_fdiv_qr_ui(q, r1, q, r);
				x[i] = mpz_get_ui(r1);
		//      mpz_set (bigint, q);
			}

			mpz_clear(q);
			mpz_clear(r1);
	}
}

mpz_class GeneralizedFermatPrimeField::Prime() const {
	return prime;
}

mpz_class GeneralizedFermatPrimeField::number() const{
	GeneralizedFermatPrimeField c(*this);
	//	char cmd_str[1024];
	mpz_t bigint;
	mpz_init(bigint);
	mpz_init_set_ui(bigint, 0);
//Horner method for evaluation
	mpz_init_set_ui(bigint, x[k - 1]);
	for (int i = k - 1; i > 0; i--)
	{
		mpz_mul_ui(bigint, bigint, r);
		mpz_add_ui(bigint, bigint, x[i - 1]);
	}


	mpz_class result(bigint);
	mpz_clear(bigint);
	return result;
}

GeneralizedFermatPrimeField GeneralizedFermatPrimeField::unitCanonical(GeneralizedFermatPrimeField* u, GeneralizedFermatPrimeField* v) const {
	if (isZero()) {
		if (u != NULL) {
			*u = 1;
		}
		if (v != NULL) {
			*v = 1;
		}
		return GeneralizedFermatPrimeField(0);
	} else {
		if (u != NULL) {
			*u = *this;
		}
		if (v != NULL) {
			*v = this->inverse();
		}
		return GeneralizedFermatPrimeField(1);
	}
}

GeneralizedFermatPrimeField& GeneralizedFermatPrimeField::operator= (const GeneralizedFermatPrimeField& c) {
	// for (int i = 0; i < k; i++) {
	// 	x[i] = c.x[i];
	// }
	memcpy(x,c.x,(k) * sizeof(usfixn64));


	return *this;
}

GeneralizedFermatPrimeField& GeneralizedFermatPrimeField::operator= (const mpz_class& c) {
	GeneralizedFermatPrimeField b (c);
	*this = b;
	return *this;
}

GeneralizedFermatPrimeField& GeneralizedFermatPrimeField::operator= (int c) {
	GeneralizedFermatPrimeField b (c);
	*this = b;
	return *this;
}

// computer zi = xi + yi for i in 0...k-1
// let zk = 0
// for i in 0...k-1, zi/r = qi*a + si
// zi+1 = zi+1 + qi
// if zk ==0 return
// if zk == 1 return (r,0,0...0)
GeneralizedFermatPrimeField& GeneralizedFermatPrimeField::operator+= (const GeneralizedFermatPrimeField& y) {

	short c = 0;
	short post = 0;
	usfixn64 sum = 0;
	int i = 0;

	for (i=0;i < k;i++) {
		sum = x[i] + y.x[i] + c;

		if (sum < x[i] || sum < y.x[i]) {
			c = 1;
			x[i] = ULMAX -r +1 + sum;
		}
		else if (sum >= r ) {
			c = 1;
			x[i] = sum - r;
		}
		else {
			x[i] = sum;
			c = 0;
		}
	}

	if (c > 0){
		post = -1;

		for (i = 0; i < k; i++) {
			if (x[i] != 0){
				post = i;
				break;
			}
		}

		if (post >= 0){
			for (i = 0; i < post; i++) {
				x[i] = r - 1;
			}
			x[post]--;
		}
		else {
			x[k-1] = r;

			for (i = 0;i < k-1; i++){
				x[i] = 0;
			}
		}
	}

	return *this;
}

// using similar algorithm as addition
GeneralizedFermatPrimeField& GeneralizedFermatPrimeField::operator-= (const GeneralizedFermatPrimeField& y) {

	int c = 0;
	int post = 0;
	usfixn64 sub = 0;
	int i = 0;

	for (i = 0; i < k; i++) {
		sub = y.x[i] + c;

		if (x[i] < sub) {
			c = 1;
			x[i] = r - sub + x[i];
		}
		else{
			c = 0;
			x[i] = x[i] - sub;
		}
	}

	if (c > 0){
		post = -1;
		for (i = 0;i<k;i++){
			if (x[i] < (r-1)) {
				post = i;
				break;
			}
		}

		if (post >= 0) {
			for (i = 0; i < post; i++) {
				x[i] = 0;
			}

			x[post] ++;
		}
		else{
			x[k-1] = r;

			for (i = 0; i < k - 1; i++){
				x[i] = 0;
			}
		}
	}
	return *this;
}

// part of multiplication
void GeneralizedFermatPrimeField::smallAdd2 (usfixn64 *xm, usfixn64* ym, short & c){
	c = 0;
	usfixn64 s = 0;

	s = xm[0] + ym[0];
	s < xm[0] || s < ym[0] ? c = 1 : c = 0;
	c > 0 ? s = s + RC : s = s;

	if (s >= r)
	{
		s = s - r;
		c = 1;
	}

	xm[0] = s;
  ym[1] = ym[1] + c;  //h1<r<2^64-1. This means no overflow
  s = xm[1] + ym[1];
  s < xm[1] || s < ym[1] ? c = 1 : c = 0;
  c > 0 ? s = s + RC : s = s;

  if (s >= r){
  	s = s - r;
  	c = 1;
  }
  xm[1] = s;
  xm[2] = xm[2] + ym[2] + c;
}

void GeneralizedFermatPrimeField::oneShiftRight (usfixn64 * xs){
	usfixn64 tmp;
	short i;
//	for (j = 0; j < shiftNo; j++)

	tmp = xs[k-1];
//#pragma unroll 7
	for (i = k - 1; i > 0; i--){
		xs[i] = xs[i - 1];
//		xs[0] = device_negate_plain(tmp);
	}
	xs[0] = tmp;
	if (xs[0] == r){
		xs[0]-=r;
		xs[1]++;
	}
}

// x*y = s0 + s1*r + s2 * r^2
void GeneralizedFermatPrimeField::mulLong_2 (usfixn64 x, usfixn64 y, usfixn64 &s0,usfixn64 &s1, usfixn64 &s2){
	std::stringstream str;
	str << x;
	mpz_class x_mpz(str.str());
	str.str( std::string() );
	str.clear();
	str << y;
	mpz_class y_mpz(str.str());
	str.str( std::string() );
	str.clear();
	str << r;
	mpz_class r_mpz(str.str());
	str.str( std::string() );
	str.clear();

	mpz_class xy;
	xy = x_mpz * y_mpz;
	str << (xy % r_mpz);
	str >> s0;
	str.str( std::string() );
	str.clear();
	xy = xy / r_mpz;
	str << (xy % r_mpz);
	str >> s1;
	str.str( std::string() );
	str.clear();
	xy = xy / r_mpz;
	str << (xy % r_mpz);
	str >> s2;
}

// special one for prime3
void GeneralizedFermatPrimeField::mulLong_3 (usfixn64 const &x, usfixn64 const &y, usfixn64 &s0,
	usfixn64 &s1, usfixn64 &s2)
{
	usfixn64 x1, y1;
	usfixn64 l, h, c, x0, y0, v2, v5, v9, v10, v11, v14, v15, v16, v17, q, t;
	usfixn64 a0, a1, b0, b1, c0, c1, c1prime, d0, d1, d2, e0, e1;

	if (x <= SQRTR && y <= SQRTR)
	{
		s0 = x * y;
		s1 = 0;
		s2 = 0;
		return;
	}


	x1 = (x >= r ? 1 : 0);
	x0 = (x1 > 0 ? x - r : x);
	y1 = (y >= r ? 1 : 0);
	y0 = (y1 > 0 ? y - r : y);

	v2 = x0 * y1; //[0,v2,0];
	v5 = x1 * y0; //[0,v5,0];
	v9 = x1 * y1; //[0,0,1];

	c = v9;
	l = 0;
	h = v5 + v2;
	h < v5 || h < v2 ? (c = c + 1) : (c = c);
	c > v9 ? (h = h + RC) : (h = h);

	if (x0 <= SQRTR && y0 <= SQRTR)
	{
		s0 = x0 * y0;
		s1 = h;
		s2 = c;
		return;
	}

//lhc
//x0*y0
	a1 = x0 >> 32;
	a0 = x0 - (a1 << 32);
	b1 = y0 >> 32;
	b0 = y0 - (b1 << 32);

	c0 = 0;
	c1 = a1 * b1;

	t = a0 * b1;
	q = t >> 32;
	t = (t - (q << 32)) << 32;
	c1 += q;
	c0 += t;  //safe

	t = a1 * b0;
	q = t >> 32;
	t = (t - (q << 32)) << 32;
	c1 += q;
	q = c0 + t;               //here, is not related to r.
	q < c0 || q < t ? (c1++) : (c1 = c1);  //c0=c0+t and carry, safe
	c0 = q;

	t = a0 * b0;
	q = c0 + t;
	q < c0 || q < t ? (c1++) : (c1 = c1);  //Now we finish [c0,c1]=x0*y0
	c0 = q;

	c1prime = c1 << 1;

	c0 >= r ? (v11 = 1) : (v11 = 0);
	v11 > 0 ? (v10 = c0 - r) : (v10 = c0);
//v12=0;

	q = l + v10;  //[l,h,c] + [v10,v11,0]
	q < l || q < v10 ? (v11 = v11 + 1) : (v11 = v11);
	q < l || q < v10 ? (l = q + RC) : (l = q);
	if (l >= r)
	{
		l = l - r;
		v11++;
	}
	q = h + v11;
	q < h || q < v11 ? (c = c + 1) : (c = c);
	q < h || q < v11 ? (h = q + RC) : (h = q);
	if (h >= r)
	{
		h = h - r;
		c++;
	}
//v13=0;
	c1prime >= r ? (v15 = 1) : (v15 = 0);
	v15 > 0 ? (v14 = c1prime - r) : (v14 = c1prime); //v13=0;

	q = h + v14;  //[l,h,c]+[0,v14,v15]
	q < h || q < v14 ? (c = c + v15 + 1) : (c = c + v15);
	q < h || q < v14 ? (h = q + RC) : (h = q);
	if (h >= r)
	{
		h = h - r;
		c++;
	}
//[l,h,c]

	d1 = c1prime >> 29;
	d0 = c1prime - (d1 << 29);
	if (d0 >= d1)
	{
		d2 = d0 - d1;
		e1 = d2 >> 29;
		e0 = d2 - (e1 << 29);
		e0 >= e1 ? (v16 = (e0 - e1) << 34) : (v16 = r - (e1 << 34) + (e0 << 34));
		e0 >= e1 ? (v17 = e1 + d1) : (v17 = e1 + d1 - 1);
		/*
		 if(e0>=e1)
		 {
		 v16=(e0-e1)<<34;
		 v17=e1+d1;
		 }
		 else
		 {
		 v17=e1+d1-1;
		 v16=R-(e1<<34)+(e0<<34);
		 }
		 */
		}
		else
		{
		//d1>d0
			d2 = d1 - d0;
			e1 = d2 >> 29;
			e0 = d2 - (e1 << 29);
			e0 >= e1 ? (v16 = r - ((e0 - e1) << 34)) : (v16 = (e1 - e0) << 34);
			e0 >= e1 ? (v17 = d1 - e1 - 1) : (v17 = d1 - e1);
		/*
		 if(e0>=e1)
		 {
		 v16=R-((e0-e1)<<34);
		 v17=d1-e1-1;
		 }
		 else
		 {
		 v16=(e1-e0)<<34;
		 v17=d1-e1;
		 }
		 */
		}
//[l,h,c]-[v16,v17,0]
//q
		q = 0;
		if (l >= v16)
		{
			l = l - v16;
		}
		else
		{
			l = r - v16 + l;
			q = 1;
		}
//t
		if (h < q + v17)
		{
			c = c - 1;
			h = r - q - v17 + h;
		}
		else
		{
			h = h - q - v17;
		}
		s0 = l;
		s1 = h;
		s2 = c;
	}

	void GeneralizedFermatPrimeField::multiplication (usfixn64* __restrict__ xs, const usfixn64* __restrict__ ys,
		usfixn64 permutationStride, usfixn64* lVector,
		usfixn64 *hVector, usfixn64* cVector,
		usfixn64* lVectorSub,
		usfixn64 *hVectorSub, usfixn64* cVectorSub ){
		short step = 0;

		short i = 0;
		usfixn64 tid = 0;
	//can change it to the following for a multithreaded implementation
	//usfixn64 tid = get_tid();

	//GeneralizedFermatPrimeField y;
		usfixn64 *y = (usfixn64*)malloc(k * sizeof(usfixn64));
		usfixn64 offset = tid;
		usfixn64 m0, m1, m2;


		usfixn64 *lhc = (usfixn64*)malloc(k * sizeof(usfixn64));
		usfixn64 *tmp = (usfixn64*)malloc(k * sizeof(usfixn64));
		usfixn64 *lhcSub = (usfixn64*)malloc(k * sizeof(usfixn64));
	//usfixn64 *lhc=new usfixn64 [k]();
	//usfixn64 *tmp=new usfixn64 [k]();
	//usfixn64 *lhcSub=new usfixn64 [k]();
		usfixn64 lArraySub;

		short k1,c;
	// for(tid=0;tid<n;tid++)
		c = 0;
		offset = tid;
		for (i = 0; i < k; i++){
			y[i] = ys[offset];
			offset += permutationStride;
		}
//	cout << "first for success"<<endl;
		offset = tid;
		k1 = 0;
//	cout << "begin mulLong_2"<<endl;
		for (k1 = 0; k1 < k; k1++){
			step = k - k1;
			offset = tid;
//		l = 0, h = 0, c = 0;
			lhc[0] = 0;
			lhc[1] = 0;
			lhc[2] = 0;

			lhcSub[0] = 0;
			lhcSub[1] = 0;
			lhcSub[2] = 0;

			if (k1 > 0){
				oneShiftRight (y);
			}

			offset = tid;
			c = 0;

			for (i = 0; i < step; i++){
				/*****************************
				* device_p3_mulLong_2_plain (xs[xIdx], ys[yIdx], m0, m1, m2)
				* multiplies digit xs_{xIdx} by digit ys_{yIdx}, then
				* rewrites the result in three machine words m0, m1, m2
				* such that
				*  xs_{xIdx} * digit ys_{yIdx} == m2*r^2 + m1*r + m0
				******************************/
				mulLong_3 (xs[offset], y[k -1 - i], m0, m1, m2);
//			cout << "m0:"<<m0<<endl;

				tmp[0] = m0;
				tmp[1] = m1;
				tmp[2] = m2;

				smallAdd2 (lhc, tmp, c);
				offset += permutationStride;
			}

			c = 0;
			for (i = step; i < k; i++){
	/*****************************
	* device_p3_mulLong_2_plain (xs[xIdx], ys[yIdx], m0, m1, m2)
	* multiplies digit xs_{xIdx} by digit ys_{yIdx}, then
	* rewrites the result in three machine words m0, m1, m2
	* such that
	*  xs_{xIdx} * digit ys_{yIdx} == m2*r^2 + m1*r + m0
	******************************/
	mulLong_3 (xs[offset], y[k -1 - i], m0, m1, m2);

	tmp[0] = m0;
	tmp[1] = m1;
	tmp[2] = m2;

	smallAdd2 (lhcSub, tmp, c);
	offset += permutationStride;
}
		// c=0;
		// device_p3_bigSubZero_3(lhcSub[0],lhcSub[1],lhcSub[2]);
		// device_p3_smallAdd2_plain (lhc,lhcSub,c);
		// printf( "=%llu \t =%llu \t =%llu \n",lhc[0],lhc[1],lhc[2]);

step = k1;
offset = tid + (permutationStride * (k -1 - k1));

lVector[offset] = lhc[0];
hVector[offset] = lhc[1];
cVector[offset] = lhc[2];

lVectorSub[offset] = lhcSub[0];
hVectorSub[offset] = lhcSub[1];
cVectorSub[offset] = lhcSub[2];
}
}

/**********************************************/

void GeneralizedFermatPrimeField::multiplication_step2 (usfixn64* __restrict__ xs, usfixn64 permutationStride,
	usfixn64* __restrict__ lVector,
	usfixn64 * __restrict__ hVector,
	usfixn64* __restrict__ cVector){
	short i = 0;
	usfixn64 offset;
	short c;
	GeneralizedFermatPrimeField v0;
	GeneralizedFermatPrimeField v1;
	GeneralizedFermatPrimeField v2;
	usfixn64 tid=0;
	offset = tid;
	i = 0;
	c = 0;
	memset (v0.x, 0x00, (k) * sizeof(usfixn64));
	memset (v1.x, 0x00, (k)* sizeof(usfixn64));
	memset (v2.x, 0x00, (k) * sizeof(usfixn64));
	//#################################################

	offset = tid;
	//move all values in lVEctor to v1 [0:7]
	for (i = 0; i < k; i++){
		v1.x[i] = lVector[tid + (i) * permutationStride];
		offset += permutationStride;
	}

	// v0[0]= last digit of hVector
//	memset (v0.x, 0x00, (k) * sizeof(usfixn64));
	v0.x[0] = hVector[tid + (k-1) * permutationStride];
	// v1[0:1]= last two digit of cVector
//	memset (v2.x, 0x00, (k) * sizeof(usfixn64));
	v2.x[0] = cVector[tid + (k-2) * permutationStride];
	v2.x[1] = cVector[tid + (k-1) * permutationStride];
	v2 = v2 + v0;
	v1 = v1 - v2;
  //  device_p3_bigPrimeAdd_plain (v2.x, v0.x);
	//	device_p3_bigSub_plain(v1.x,v2.x);

	offset = tid;
	memset (v0.x, 0x00, (k) * sizeof(usfixn64));
	//positive h's [1:7]
	for (i = 1; i < k; i++){
		v0.x[i] = hVector[tid + (i - 1) * permutationStride];
		offset += permutationStride;
	}
	v1 = v1 + v0;
	//device_p3_bigPrimeAdd_plain (v1, v0);

	offset = tid;
		//positive c's [2:7]
	memset (v0.x, 0x00, (k) * sizeof(usfixn64));
	for (i = 2; i < k; i++){
		v0.x[i] = cVector[tid + (i - 2) * permutationStride];
		offset += permutationStride;
	}
	v1 = v1 + v0;
	//device_p3_bigPrimeAdd_plain (v1, v0);
	//#################################################
	offset = tid;
	//#################################################
	for (i = 0; i < (k); i++){
		xs[offset] = v1.x[i];
		offset += permutationStride;
	}
}

// special multiplication for P3
//r = 9223372054034644992
// k = 8
GeneralizedFermatPrimeField GeneralizedFermatPrimeField::MultiP3 (GeneralizedFermatPrimeField ys) {
	GeneralizedFermatPrimeField xs (*this);
	if (x[k-1] == r){
		*this = -ys;
		return *this;
	}
	else if (ys.x[k-1] == r){
		*this = -xs;
		return *this;
	}
	else{
		int permutationStride = 1;
		usfixn64* lVector = new usfixn64[k]();
		usfixn64* hVector = new usfixn64[k]();
		usfixn64* cVector = new usfixn64[k]();

		usfixn64* lVectorSub = new usfixn64[k]();
		usfixn64* hVectorSub = new usfixn64[k]();
		usfixn64* cVectorSub = new usfixn64[k]();

		multiplication(&xs.x[0], &ys.x[0], permutationStride, &lVector[0],&hVector[0],&cVector[0],&lVectorSub[0],&hVectorSub[0],&cVectorSub[0]);
// //h:=h*r
// //c:=c*r^2
// //x:= l + h +c;
//xs: step 2 for positive values
		multiplication_step2(&xs.x[0], permutationStride,&lVector[0],&hVector[0],&cVector[0]);
//ys: step 2 for negative values
		multiplication_step2(&ys.x[0], permutationStride,&lVectorSub[0],&hVectorSub[0],&cVectorSub[0]);
//xs-ys: subtracting positive values from negative values
		xs = xs - ys;
//  device_p3_bigSub_plain(&xs[i*COEFFICIENT_SIZE],&ys[i*COEFFICIENT_SIZE]);
		*this = xs;
		return *this;
	}
}


//Multiple by some power of r
//using shift only
GeneralizedFermatPrimeField GeneralizedFermatPrimeField::MulPowR(int s){
	GeneralizedFermatPrimeField a(0);
	GeneralizedFermatPrimeField b(0);
	GeneralizedFermatPrimeField c(*this);
	s = s%(2 * k);
	if (s == 0){
		return c;
	}
	else if (s == k){
		return (-c);
	}
	else if ((s > k) && (s < (2 * k))){
		s = s - k;
		c = -c;
	}

	int i;
	for (i = 0; i < (k - s); i++){
		b.x[i + s] = c.x[i];
	}
	for (i = k - s; i < k; i++){
		a.x[i - (k - s)] = c.x[i];
	}
	if(c.x[k-1] == r){
		a.x[s-1] -=r;
		a.x[s] ++;
	}
	b = b - a;
	return b;
}


GeneralizedFermatPrimeField GeneralizedFermatPrimeField::euclideanDivision(const GeneralizedFermatPrimeField& b, GeneralizedFermatPrimeField* q) const {
	//TODO
	std::cerr << "GeneralizedFermatPrimeField::euclideanDivision not yet implemented" << std::endl;
	exit(1);
	return *this;
}

GeneralizedFermatPrimeField GeneralizedFermatPrimeField::extendedEuclidean(const GeneralizedFermatPrimeField& b, GeneralizedFermatPrimeField* s, GeneralizedFermatPrimeField* t) const {
	//TODO
	std::cerr << "GeneralizedFermatPrimeField::extendedEuclidean not yet implemented" << std::endl;
	exit(1);
	return *this;
}

GeneralizedFermatPrimeField GeneralizedFermatPrimeField::quotient(const GeneralizedFermatPrimeField& b) const {
	//TODO
	std::cerr << "GeneralizedFermatPrimeField::quotient not yet implemented" << std::endl;
	exit(1);
	return *this;
}

GeneralizedFermatPrimeField GeneralizedFermatPrimeField::remainder(const GeneralizedFermatPrimeField& b) const {
	//TODO

	std::cerr << "GeneralizedFermatPrimeField::remainder not yet implemented" << std::endl;
	exit(1);
	return *this;
}

void GeneralizedFermatPrimeField::egcd (const mpz_class& x, const mpz_class& y, mpz_class *ao, mpz_class *bo, mpz_class *vo, mpz_class P)
{
	mpz_class t, A, B, C, D, u, v, q;

	u = y;
	v = x;
	A = 1;
	B = 0;
	C = 0;
	D = 1;

	do
	{
		q = u / v;
		t = u;
		u = v;
		v = t - q * v;
		t = A;
		A = B;
		B = t - q * B;
		t = C;
		C = D;
		D = t - q * D;
	}
	while (v != 0);

	*ao = A;
	*bo = C;
	*vo = u;
}

// GMP inverse
GeneralizedFermatPrimeField GeneralizedFermatPrimeField::inverse() const {
	GeneralizedFermatPrimeField R(*this);
	mpz_t temp;
	mpz_init(temp);
	if (!mpz_invert(temp, R.number().get_mpz_t(), prime.get_mpz_t()))
		mpz_set_si(temp, 0);
	mpz_class temp_class(temp);
	mpz_clear(temp);
	R = temp_class;
	return R;
}
