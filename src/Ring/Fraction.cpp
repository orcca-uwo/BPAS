#include "Ring/Fraction.hpp"

template <class Domain>
void Fraction<Domain>::setNumerator(const Domain& b){
	num = b;
	canonicalize();
}

template <class Domain>
void Fraction<Domain>::setDenominator(const Domain& b){
	den = b;
	canonicalize();
}


template <class Domain>
void Fraction<Domain>::set(const Domain& a, const Domain& b){
	num = a;
	den = b;
	canonicalize();
}

template <class Domain>
Domain Fraction<Domain>::numerator() const{
	return num;
}

template <class Domain>
Domain Fraction<Domain>::denominator() const{
	return den;
}

template <class Domain>
bool Fraction<Domain>::operator!=(const Fraction<Domain> &b) const {
	return (!(num == b.num)||!(den == b.den));
}

template <class Domain>
bool Fraction<Domain>::operator==(const Fraction<Domain> &b) const {
	return ((num == b.num)&&(den == b.den));
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::operator+(const Fraction<Domain> &b) const{
	Fraction<Domain> r(*this);
	return (r+=b);
}

template <class Domain>
Fraction<Domain>& Fraction<Domain>::operator+=(const Fraction<Domain> &b){
	Domain g, r1(den), r2(b.den);
	g = den.gcd(b.den);
	r1/=g;
	r2/=g;

	den = r1*r2;
	r1*=b.num;
	r2*=num;
	r1 += r2; // r1 = ar2 + cr1;
	r2 = r1.gcd(g);
	r1 /= r2; // r1 = e;
	g /= r2; // g = g';
	den *= g;
	num = r1;
	
	canonicalize();

	normalize();

	return *this;

}

template <class Domain>
Fraction<Domain> Fraction<Domain>::operator-(const Fraction<Domain> &b) const{
	Fraction<Domain> r(*this);
	return (r-=b);
}

template <class Domain>
Fraction<Domain>& Fraction<Domain>::operator-=(const Fraction<Domain> &b){
	Domain g, r1(den), r2(b.den);
	g = den.gcd(b.den);
	r1/=g;
	r2/=g;
	den = r1*r2;
	r1*= -b.num;
	r2*=num;
	r1 += r2; // r1 = ar2 + cr1;
	r2 = r1.gcd(g);
	r1 /= r2; // r1 = e;
	g /= r2; // g = g';
	den *= g;
	num = r1;
	//canonicalize();
	normalize();

	return *this;
}






template <class Domain>
Fraction<Domain> Fraction<Domain>::operator*(const Fraction<Domain> &b) const{
	Fraction<Domain> r(*this);
	return (r*=b);
}

template <class Domain>
Fraction<Domain>& Fraction<Domain>::operator*=(const Fraction<Domain> &b){
	Domain g1,g2,r;
	g1 = num.gcd(b.den);
	g2 = den.gcd(b.num);
	num/=g1;
	r = b.num/g2;
	num *= r;
	den /= g2;
	r = b.den / g1;
	den *= r;
	//canonicalize();
	normalize();


	return *this;  
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::operator/(const Fraction<Domain> &b) const{
	Fraction<Domain> r(*this);
	return (r/=b);
}

template <class Domain>
Fraction<Domain>& Fraction<Domain>::operator/=(const Fraction<Domain> &b){
	if(b.isZero()){
		std::cout << "BPAS error: division by zero from Fraction<Domain> operator/=" << std::endl;
		exit(1);
	}

	Fraction<Domain> e(b.den,b.num);
	*this *=e;
	return *this;

}

template <class Domain>
Fraction<Domain> Fraction<Domain>::operator-() const{
	Fraction<Domain> r(-num,den);
	return r;
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::operator^(long long int e) const{
	Fraction<Domain> res;
	if(isZero()||isOne()||e == 1){
		res = *this;
	}
	else if(e==2){
		res = (*this) * (*this);
	}

	else if(e>2){
		Fraction<Domain> x (*this);
		res.one();
		while(e){
			if(e%2){res*=x;};
			x = x*x;
			e>>=1;
		}
	}
	else if(!e)
		res.one();
	else
	{
		
		res = *this;
	}
	return res;
}

template <class Domain>
Fraction<Domain>& Fraction<Domain>::operator^=(long long int e){
	*this = *this^e;
	return *this;
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::inverse() const{
	if (num.isZero()) {
	std::cout << "BPAS error: division by zero from Fraction<Domain> inverse()" << std::endl;
		exit(1);
	}
	Fraction<Domain> r(den, num);
	r.normalize();
	return r;
}

template <class Domain>
bool Fraction<Domain>::isZero() const{
	return num.isZero();
}

template <class Domain>
void Fraction<Domain>::zero(){
	num.zero();
	den.one();
	return;
}

template <class Domain>
bool Fraction<Domain>::isOne() const {
	if (num.isOne() && den.isOne())
	{
		return true;
	}
	else{
		return false;
	}
}

template <class Domain>
void Fraction<Domain>::one(){
	num.one();
	den.one();
	return;
}

template <class Domain>
bool Fraction<Domain>::isNegativeOne() const{
	return (num.isNegativeOne()&&den.isOne()) || (num.isOne()&&den.isNegativeOne());
}

template <class Domain>
void Fraction<Domain>::negativeOne(){
	num.negativeOne();
	den.one();
	return;
}

template <class Domain>
int Fraction<Domain>::isConstant() const{
	return num.isConstant()&&den.isConstant();
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::unitCanonical(Fraction<Domain>* u , Fraction<Domain>* v ) const{
	Domain du,dc,dv,g,temp;
	g = num.gcd(den);
	temp = den/g;
	dc = temp.unitCanonical(&du,&dv);
	temp = du*(num/g);
	Fraction<Domain> ret;
	ret.num = temp;
	ret.den = dc;
	if (u != NULL || v!= NULL) {
		Fraction<Domain> temp2;
		temp2.one();
		if (u != NULL) {
			*u = temp2;
		}
		if (v != NULL) {
			*v = temp2;
		}
	}
	return ret;
}

template <class Domain>
Fraction<Domain>& Fraction<Domain>::operator=(const Fraction<Domain>& b){
	if(this!=&b){
		num = b.num;
		den = b.den;
		characteristic = b.characteristic;
	}
	return *this;
}


template <class Domain>
void Fraction<Domain>::print(std::ostream& ostream) const {
	ostream << "(" << num << ")/(" << den << ")"; 
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::quotient(const Fraction<Domain>& b) const{
	return (*this / b);
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::remainder(const Fraction<Domain>& b) const{
	return Fraction<Domain>(0,1);
}

template <class Domain>
Fraction<Domain> Fraction<Domain>::operator%(const Fraction<Domain>& b) const{
	Fraction<Domain> ret = remainder(b);
	return ret;
}

template <class Domain>
Fraction<Domain>& Fraction<Domain>::operator%=(const Fraction<Domain>& b){
	*this = remainder(b);
	return *this;
}

template <class Domain>
void Fraction<Domain>::canonicalize(){
	Domain temp;
	temp = num.gcd(den);
	num/=temp;
	den/=temp;
	return;
}


template <class Domain>
void Fraction<Domain>::normalize(){

	Domain u, v; 
	Domain C = den.unitCanonical(&u,&v);
	num *= v;
	den = C;


}

template <class Domain>
void Fraction<Domain>::differentiate(){
	/*Domain D(den);
	Domain dD(den);
	Domain temp;
	dD.differentiate(1);
	temp = D.gcd(dD);
	D /= temp;
	dD /= temp;	
	temp = -num;
	temp *= dD;
	dD = num;
	dD.differentiate(1);
	dD *= D;
	temp += dD;
	D *= den; 
	num = temp;
	den = D;
	canonicalize();
	return;*/
}





template class Fraction<Integer>;
template class Fraction<DenseUnivariateRationalPolynomial>;



