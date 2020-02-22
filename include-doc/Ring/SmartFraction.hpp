#ifndef _SMART_FRACTION_H_
#define _SMART_FRACTION_H_

#include <iostream>
#include "BPASFieldOfFractions.hpp"
#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"
#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include <algorithm>
#include <vector>
#include <utility>
#include "FactorRefinement.hpp"


/**
 * A field of fractions templated by an arbitrary BPASGCDDomain making
 * use of factor refinement.
 */
template <class Domain>
class SmartFraction: public BPASFieldOfFractions<Domain> {

private:
	std::vector<Factor<Domain>> num;
	std::vector<Factor<Domain>> den;
	void smart_unitCanonical(const std::vector<Factor<Domain>>&A, std::vector<Factor<Domain>>*u, std::vector<Factor<Domain>>* C, std::vector<Factor<Domain>>* v) const{

		for (int i = 0; i < A.size(); ++i)
		{
			Domain curr_Domain = A[i].first;
			int curr_exp = A[i].second;
			Domain temp_u;
			Domain temp_v;
			Domain temp_C = curr_Domain.unitCanonical(&temp_u,&temp_v);
			u->push_back(std::make_pair(temp_u,curr_exp));
			C->push_back(std::make_pair(temp_C,curr_exp));
			v->push_back(std::make_pair(temp_v,curr_exp));
		}
		return;
	}
 	std::vector<Factor<Domain>> smart_gcd(const std::vector<Factor<Domain>>&A,const std::vector<Factor<Domain>>&B) const{
 		std::vector<Factor<Domain>> left, mid, right;
 		FactorRefinement::MergeRefineTwoSeq<Domain>(A,B,&left,&mid,&right);
 		return mid;
 	}

 	static std::vector<Factor<Domain>> smart_add(const std::vector<Factor<Domain>>& A,const std::vector<Factor<Domain>>& B) {
 		Domain a;
 		a = convertToDomain(A);
 		Domain b;
 		b = convertToDomain(B);
 		a = a + b;
 		std::vector<Factor<Domain>> v;
 		v = extractFactors(a);
 		return v;
 	}

 	static std::vector<Factor<Domain>> smart_sub(const std::vector<Factor<Domain>>&A,const std::vector<Factor<Domain>>&B) {
 		Domain a;
 		a = convertToDomain(A);
 		Domain b;
 		b = convertToDomain(B);
 		a = a - b;
 		std::vector<Factor<Domain>> v;
 		v = extractFactors(a);
 		return v;
 	}
 	static std::vector<Factor<Domain>> smart_mul(const std::vector<Factor<Domain>>&A,const std::vector<Factor<Domain>>&B) {

		std::vector<Factor<Domain>> v1,v2,v3;
		FactorRefinement::MergeRefineTwoSeq<Domain>(A,B,&v1,&v2,&v3);
		//catVectors(&v1,v2,v3);

		v1.insert(v1.end(),v2.begin(),v2.end());
		v1.insert(v1.end(),v3.begin(),v3.end());
		return v1;
 	}
 	std::vector<Factor<Domain>> smart_div(const std::vector<Factor<Domain>>&A,const std::vector<Factor<Domain>>&B) const{
		std::vector<Factor<Domain>> v1,v2,v3;
		FactorRefinement::MergeRefineTwoSeq<Domain>(A,B,&v1,&v2,&v3);
		if(v3.size()!= 0){
			cout << "BPAS error: in SmartFraction<Domain>, not exact division" <<endl;
			exit(1);
		}
		return v1;
	}





 	static void smart_one(std::vector<Factor<Domain>>* b){
 		b->clear();
 		Domain one;
 		one.one();
 		b->push_back(std::make_pair(one,1));
 	}

 	static std::vector<Factor<Domain>> smart_one(){
 		std::vector<Factor<Domain>> b;
 		Domain one;
 		one.one();
 		b.push_back(std::make_pair(one,1));
 		return b;

 	}

 	void smart_zero(std::vector<Factor<Domain>>& b){
 		b.clear();
 	}

 	static bool smart_isOne(std::vector<Factor<Domain>>& b){
 		if(b.size() == 1 && b[0].first.isOne()){
 			return true;
 		}
 		else{
 			return false;
 		}

 	}





public:
	mpz_class characteristic;
	/**
	* Construct the zero fraction function
	*
	* @param
	**/ 
	SmartFraction<Domain>(){
		num;
		den;
		Domain c;
		characteristic = c.characteristic;
	}

	/**
	* Copy constructor
	*
	* @param b: A rational function
	**/ 
	SmartFraction<Domain>(const SmartFraction<Domain> &b){
		num = b.num;
		den = b.den;
		Domain c;
		characteristic = c.characteristic;
	}

	/**
	* constructor with two parameter
	*
	* @param a: the numerator
	* @param b: the denominator
	**/
	SmartFraction<Domain>(std::vector<Factor<Domain>> a, std::vector<Factor<Domain>> b){
		num = a;
		den = b;
		Domain c;
		characteristic = c.characteristic;
	}

	SmartFraction<Domain>(Domain a, Domain b){
		//cout << "print from constructor"<< endl;
		//a.print(cout);
		//cout << "num end" << endl;
		//b.print(cout);
		//cout << "den end" << endl;

		num = extractFactors(a);
		den = extractFactors(b);

		//cout << "const : num size " << num.size() << endl;	
		//cout << "const : den size " << den.size() << endl;

		//cout << "num contain:" << endl;
		//num[0].first.print(std::cout);
		//cout << endl;
		//normalize();


		Domain c;
		characteristic = c.characteristic;
	}

	~SmartFraction<Domain>(){};

	void setNumerator(const std::vector<std::pair<Domain, int>>& b);
	void setDenominator(const std::vector<std::pair<Domain, int>>& b);

	void set(const std::vector<std::pair<Domain, int>>& a, const std::vector<std::pair<Domain, int>>& b);

	Domain numerator() const;
	Domain denominator() const;

	bool operator!=(const SmartFraction<Domain> &b) const;
	bool operator==(const SmartFraction<Domain> &b) const;

	SmartFraction<Domain> operator*(const SmartFraction<Domain>& b) const;
	SmartFraction<Domain>& operator*=(const SmartFraction<Domain>& b);


	SmartFraction<Domain> unitCanonical(SmartFraction<Domain> *u = NULL, SmartFraction<Domain> *v = NULL) const;

	void canonicalize();
	void normalize();

	bool isZero() const;
	bool isOne() const;
	void zero();
	void one();
	SmartFraction<Domain> operator+(const SmartFraction<Domain>& b) const;
	SmartFraction<Domain>& operator+=(const SmartFraction<Domain>& b);
	SmartFraction<Domain> operator-(const SmartFraction<Domain>& b) const;
	SmartFraction<Domain>& operator-=(const SmartFraction<Domain>& b);
	SmartFraction<Domain> operator/(const SmartFraction<Domain>& b) const;
	SmartFraction<Domain>& operator/=(const SmartFraction<Domain>& b);

	SmartFraction<Domain> operator-() const;

	SmartFraction<Domain> inverse() const;

	SmartFraction<Domain> operator^(long long int e) const;

	SmartFraction<Domain>& operator^=(long long int e);

	Factors<SmartFraction<Domain>> squareFree() const{
		std::cerr << "SmartFraction<Domain>::squareFree NOT YET IMPLEMENTED" << std::endl;
		//TODO
		return Factors<SmartFraction<Domain>>(*this); 
	}
	ExpressionTree convertToExpressionTree() const{
			std::cerr << "SmartFraction<Domain>::convertToExpressionTree NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return ExpressionTree();

	}

	void print(std::ostream& ostream) const;
	SmartFraction<Domain> gcd(const SmartFraction<Domain>& b ) const;
	SmartFraction<Domain> euclideanSize() const;
	SmartFraction<Domain> euclideanDivision(const SmartFraction<Domain>&b, SmartFraction<Domain>* q = NULL) const;
	SmartFraction<Domain> quotient(const SmartFraction<Domain>& b) const;
	SmartFraction<Domain> remainder(const SmartFraction<Domain>&b) const;
	SmartFraction<Domain> extendedEuclidean(const SmartFraction<Domain> &b, SmartFraction<Domain>* s = NULL, SmartFraction<Domain>* t = NULL) const {}
	SmartFraction<Domain> operator%(const SmartFraction<Domain>&b) const;
	SmartFraction<Domain>& operator%=(const SmartFraction<Domain> &b);




};

	template <class Domain>
 	std::vector<Factor<Domain>> extractFactors(const Domain &p) {
 		std::vector<Factor<Domain>> ret;
 		Factors<Domain> v = p.squareFree();
 		if(v.size() == 1){
 			ret.emplace_back(v[0].first * v.ringElement(), 1);
 			return ret;
 		}

 		Domain v0 = v.ringElement();
 		if(!v0.isOne()){
 			ret.emplace_back(v[0].first * v0, v[0].second);
 		} else {
 			ret.emplace_back(v[0].first, v[0].second);
 		}
 		for (int i = 1; i < v.size(); ++i){
 			if(!v[i].first.isOne()){
 				ret.emplace_back(v[i].first, v[i].second);
 			}
 		}

 		return ret;
 	}




	template <class Domain>
 	Domain convertToDomain(const std::vector<Factor<Domain>>& b){
 		Domain ret,temp;
 		long long int a;
 		ret.one();

 		for (int i = 0; i < b.size(); ++i)
 		{
 			temp = b[i].first;
 			a = b[i].second;
 			//std::cout << "temp = " << temp << ", a = " << a << std::endl;
 			ret *= temp^a;
 		}

 		return ret;
 	}



	


#endif
