#include "Ring/SmartFraction.hpp"






template <class Domain>
void SmartFraction<Domain>::setNumerator(const std::vector<std::pair<Domain, int>>& b){
	num.clear();
	for (auto p : b) {
		num.push_back(p);
	}
	//canonicalize();
}
template <class Domain>
void SmartFraction<Domain>::setDenominator(const std::vector<std::pair<Domain, int>>& b){
	den.clear();
	for (auto p : b) {
		den.push_back(p);
	}
	//canonicalize();
}


template <class Domain>
void SmartFraction<Domain>::set(const std::vector<std::pair<Domain, int>>& a, const std::vector<std::pair<Domain, int>>& b){
	this->setNumerator(a);
	this->setDenominator(a);
}

template <class Domain>
Domain SmartFraction<Domain>::numerator() const{
	Domain b =  convertToDomain(num);
	return b;
}

template <class Domain>
Domain SmartFraction<Domain>::denominator() const{
	Domain b = convertToDomain(den);
	return b;
}

template <class Domain>
bool SmartFraction<Domain>::operator!=(const SmartFraction<Domain> &b) const{
	return (!(num==b.num)||!(den==b.den));
}

template <class Domain>
bool SmartFraction<Domain>::operator==(const SmartFraction<Domain> &b) const{
	std::vector<Factor<Domain>> x,m1,y,z,m2,w;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(num,b.num,&x,&m1,&y);
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(den,b.den,&z,&m2,&w);
	if(smart_isOne(x)&&smart_isOne(y)&&smart_isOne(z)&&smart_isOne(w)){
		return true;
	}else{
		return false;
	}

	//return ((num==b.num)&&(den==b.den));
}

template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::operator*(const SmartFraction<Domain>& b) const{
	SmartFraction<Domain> r(*this);
	return (r*=b);

}


template <class Domain>
SmartFraction<Domain>& SmartFraction<Domain>::operator*=(const SmartFraction<Domain>& b){
		std::vector<Factor<Domain>> x,m1,y,z,m2,w, num1,num2,num3,den1,den2,den3;
		FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(num,b.den,&x,&m1,&y);

		FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(b.num,den,&z,&m2,&w);
		FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(x,z,&num1,&num2,&num3);
		FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(y,w,&den1,&den2,&den3);


		num1.insert(num1.end(),num2.begin(),num2.end());
		num1.insert(num1.end(),num3.begin(),num3.end());

		den1.insert(den1.end(),den2.begin(),den2.end());
		den1.insert(den1.end(),den3.begin(),den3.end());

		num = num1;
		den = den1;
		normalize();
		return *this;
}


template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::unitCanonical(SmartFraction<Domain> *u , SmartFraction<Domain> *v) const{
	std::vector<Factor<Domain>> g = smart_gcd(num,den);
	std::vector<Factor<Domain>> temp = smart_div(den,g);
	std::vector<Factor<Domain>> du,dc,dv;
	smart_unitCanonical(temp,&du,&dc,&dv);
	temp = smart_mul(du,smart_div(num,g));
	SmartFraction<Domain> ret;
	ret.num = temp;
	ret.den = dc;
	if(u!=NULL || v!=NULL){
		std::vector<Factor<Domain>> temp2;
		smart_one(&temp2);
		if(u!=NULL){
			u->den =  smart_one();
			u->num =  temp2;

		}
		if(v!=NULL){
			v->den = smart_one();
			v->num = temp2;
		}
	}
	return ret;
}

template <class Domain>
void SmartFraction<Domain>::canonicalize(){
	std::vector<Factor<Domain>> temp;
	temp = smart_gcd(num,den);
	num = smart_div(num,temp);
	den = smart_div(den,temp);
	return;
}

template <class Domain>
void SmartFraction<Domain>::normalize(){
	std::vector<Factor<Domain>>	u,C,v;
	smart_unitCanonical(den,&u,&C,&v);
	num = smart_mul(num,v);
	den = C;

}


template <class Domain>
void SmartFraction<Domain>::print(std::ostream& ostream) const{


	Domain n = convertToDomain(num);
	Domain d = convertToDomain(den);

	ostream << "(" << n << ")/(" << d << ")" << endl;


}


template <class Domain>
bool SmartFraction<Domain>::isZero() const{
	if(num.size() == 0){
		return true;
	}
	else{
		return false;
	}
}

template <class Domain>
bool SmartFraction<Domain>::isOne() const{
	if(num.size()!=1||den.size()!=1){
		return false;
	}
	else{
		if(num[0].first.isOne()&&den[0].first.isOne()){
			return true;
		}
		else{
			return false;
		}
	}
}


template <class Domain>
void SmartFraction<Domain>::zero(){
	num.clear();
	den.clear();
	return;
}


template <class Domain>
void SmartFraction<Domain>::one(){
	std::vector<Factor<Domain>> one;
	one = smart_one();
	den = one;
	num = one;
	return;
}



template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::operator+(const SmartFraction<Domain>& b) const{
	SmartFraction<Domain> r(*this);
	return (r += b);
}


template <class Domain>
SmartFraction<Domain>& SmartFraction<Domain>::operator+=(const SmartFraction<Domain>& b){


	std::vector<Factor<Domain>> da_,gd,db_;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(den,b.den,&da_,&gd,&db_);
	std::vector<Factor<Domain>> na_,gn,nb_;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(num,b.num,&na_,&gn,&nb_);
	std::vector<Factor<Domain>> m1,m2,s1;

	m1 = smart_mul(da_,nb_);

	m2 = smart_mul(db_,na_);

	s1 = smart_add(m1,m2);

	for (int i = 0; i < gn.size(); ++i)
	{
		gn[i].second = (gn[i].second/2);
	}

	for (int i = 0; i < gd.size(); ++i)
	{
		gd[i].second = (gd[i].second/2);
	}

	std::vector<Factor<Domain>> n_,n2,n3;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(gn,s1,&n_,&n2,&n3);
	n_.insert(n_.end(),n2.begin(),n2.end());
	n_.insert(n_.end(),n3.begin(),n3.end());


	std::vector<Factor<Domain>> num1,num2,df;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(n_,gd,&num1,&num2,&df);

	num.clear();
	num.insert(num.end(),num1.begin(),num1.end());

	std::vector<Factor<Domain>> d1,d2,d3;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(da_,db_,&d1,&d2,&d3);
	d1.insert(d1.end(),d2.begin(),d2.end());
	d1.insert(d1.end(),d3.begin(),d3.end());

	std::vector<Factor<Domain>> d4,d5,d6;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(d1,df,&d4,&d5,&d6);
	d4.insert(d4.end(),d5.begin(),d5.end());
	d4.insert(d4.end(),d6.begin(),d6.end());

	den.clear();
	den.insert(den.end(),d4.begin(),d4.end());

	normalize();

	return *this;

}

template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::operator-(const SmartFraction<Domain>& b) const{
	SmartFraction<Domain> r(*this);
	return (r -= b);

	return *this;


}


template <class Domain>
SmartFraction<Domain>& SmartFraction<Domain>::operator-=(const SmartFraction<Domain>& b){
	/*Domain dnum = convertToDomain(num);
	Domain dden = convertToDomain(den);
	Domain dbnum = convertToDomain(b.num);
	Domain dbden = convertToDomain(b.den);
	Domain g = dden.gcd(dbden);
	Domain r1 = dden;
	Domain r2 = dbden;
	r1/=g;
	r2/=g;
	dden = r1*r2;
	r1*= -dbnum;
	r2*=dnum;
	r1+=r2;
	r2 = r1.gcd(g);
	r1/=r2;
	dden *=g;
	dnum = r1;
	num = extractFactors(dnum);
	den = extractFactors(dden);
	normalize();*/

	std::vector<Factor<Domain>> da_,gd,db_;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(den,b.den,&da_,&gd,&db_);
	std::vector<Factor<Domain>> na_,gn,nb_;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(num,b.num,&na_,&gn,&nb_);
	std::vector<Factor<Domain>> m1,m2,s1;
	m1 = smart_mul(da_,nb_);
	m2 = smart_mul(db_,na_);
	s1 = smart_sub(m2,m1);
	for (int i = 0; i < gn.size(); ++i)
	{
		gn[i].second = (gn[i].second/2);
	}

	for (int i = 0; i < gd.size(); ++i)
	{
		gd[i].second = (gd[i].second/2);
	}

	std::vector<Factor<Domain>> n_,n2,n3;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(gn,s1,&n_,&n2,&n3);
	n_.insert(n_.end(),n2.begin(),n2.end());
	n_.insert(n_.end(),n3.begin(),n3.end());


	std::vector<Factor<Domain>> num1,num2,df;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(n_,gd,&num1,&num2,&df);

	num.clear();
	num.insert(num.end(),num1.begin(),num1.end());

	std::vector<Factor<Domain>> d1,d2,d3;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(da_,db_,&d1,&d2,&d3);
	d1.insert(d1.end(),d2.begin(),d2.end());
	d1.insert(d1.end(),d3.begin(),d3.end());

	std::vector<Factor<Domain>> d4,d5,d6;
	FactorRefinement::MergeRefineTwoSeqEmptyToIdentity<Domain>(d1,df,&d4,&d5,&d6);
	d4.insert(d4.end(),d5.begin(),d5.end());
	d4.insert(d4.end(),d6.begin(),d6.end());

	den.clear();
	den.insert(den.end(),d4.begin(),d4.end());

	normalize();

	return *this;
}

template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::inverse() const{
	if(num.size() == 0){
		std::cout << "BPAS error: division by zero from SmartFraction<Domain> inverse()" << std::endl;
		exit(1);
	}
	SmartFraction<Domain> r(den,num);
	r.normalize();
	return r;
}


template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::operator/(const SmartFraction<Domain>& b) const{
	SmartFraction<Domain> r(*this);
	return (r/=b);
}


template <class Domain>
SmartFraction<Domain>& SmartFraction<Domain>::operator/=(const SmartFraction<Domain>& b){

	if (b.isZero()) {
		std::cout << "BPAS error: division by zero from SmartFraction<Domain> operator/=" << std::endl;
		exit(1);
	}

	SmartFraction<Domain> e (b.den, b.num);
		*this *= e;
		return *this;
}

template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::operator^(long long int e) const{
	SmartFraction<Domain> res;
	if(isZero()||isOne()||e == 1){
		res = *this;
	}
	else if(e==2){
		res = (*this) * (*this);
	}

	else if(e>2){
		SmartFraction<Domain> x (*this);
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
SmartFraction<Domain>& SmartFraction<Domain>::operator^=(long long int e){
    *this = *this ^ e;
    return *this;

}


template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::operator-() const{
	Domain temp_num = convertToDomain(num);
	temp_num = -temp_num;
	std::vector<Factor<Domain>> f_num = extractFactors(temp_num);
	SmartFraction<Domain> r(f_num,den);
	return r;
}


template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::gcd(const SmartFraction<Domain>& b ) const{
	    	SmartFraction<Domain> ret;
    		if(isZero() == true &&b.isZero() == true){
    			ret.zero();

    		}
    		//otherwise is one
    		else{
    			ret.one();
    		}

    		return ret;
}


template <class Domain>
Integer SmartFraction<Domain>::euclideanSize() const{
	    	if(isZero() == true){
    			std::cerr << "in Fraction<Domain>, zero does not have a euclidean size" << std::endl;
    			exit(1);
    		}
    		else{
    			return Integer(1);
    		}
}


template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::euclideanDivision(const SmartFraction<Domain>&b, SmartFraction<Domain>* q ) const{
    	SmartFraction<Domain> ret;
    	ret.zero();
    	if(q!=NULL){
    		SmartFraction<Domain> curr(num,den);
    		*q = curr/b;
    	}
    	return ret;
}


template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::quotient(const SmartFraction<Domain>& b) const{
	return (*this / b);
}


template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::remainder(const SmartFraction<Domain>&b) const{
	SmartFraction<Domain> r;
	return r;
}



template <class Domain>
SmartFraction<Domain> SmartFraction<Domain>::operator%(const SmartFraction<Domain>&b) const{
	SmartFraction<Domain> ret = remainder(b);
	return ret;
}

template <class Domain>
SmartFraction<Domain>& SmartFraction<Domain>::operator%=(const SmartFraction<Domain> &b){
	*this = remainder(b);
	return *this;
}




template class SmartFraction<DenseUnivariateRationalPolynomial>;


