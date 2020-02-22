#ifndef _FACTOR_REFINEMENT_H_
#define _FACTOR_REFINEMENT_H_

#include <iostream>
#include "BPASFieldOfFractions.hpp"
#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"
#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include <algorithm>
#include <vector>
#include <utility>
#include "../Utils/TemplateHelpers.hpp"
#include "../DataStructures/Factors.hpp"

namespace FactorRefinement {
	/*
	*PolyRefine: it will do the factor refinement for two object a^e and b^f
	*input a, e, b, f, which are the input for FactorRefinement algorithm
	*input *ret_l,*ret_e,*G,*ret_r,*ret_f, should be empty pointers, the result of the algorithm will store in those pointers
	*ret_l^ret_e is the refinement result for a^e
	*ret_r^ret_f is the refinement result for b^f
	*G is the gcd parts for FactorRefinement algorithm 
	*/
	template <class Domain>
	static void PolyRefine(const Domain& a, int e, const Domain& b, int f, Domain* ret_l, int* ret_e, std::vector<Factor<Domain>>* ret_G, Domain* ret_r, int* ret_f) {
		Domain g = a.gcd(b);
		Domain a_ = a/g;
		Domain b_ = b/g;
		if(g.isOne() == true){
			*ret_l = a;
			//cout << "ret_l " << *ret_l << endl;
			*ret_e = e;
			//cout << "ret_e " << *ret_e << endl;
			//Domain one;
			//one.one();
			//ret_G->push_back(std::make_pair(one,1));
			ret_G->clear();
			*ret_r = b;
			*ret_f = f;
			return;

		}
		else if(a==b){
			Domain one;
			one.one();
			ret_G->clear();
			ret_G->push_back(std::make_pair(a,e+f));
			*ret_l = one;
			*ret_e = 1;
			*ret_r = one;
			*ret_f = 1;
			return;
		}
		else{
			Domain l1,l2;
			int e1,e2;
			Domain r1,r2;
			int f1,f2;	
			std::vector<Factor<Domain>> G1,G2;
			PolyRefine(a_,e,g,e+f,&l1,&e1,&G1,&r1,&f1);
			PolyRefine(r1,f1,b_,f,&l2,&e2,&G2,&r2,&f2);
			if(l2.isOne()==false){
				G2.push_back(std::make_pair(l2,e2));
			}
			*ret_l = l1;
			*ret_e = e1;
			ret_G->clear();
			ret_G->insert(ret_G->end(),G1.begin(),G1.end());
			ret_G->insert(ret_G->end(),G2.begin(),G2.end());
			*ret_r = r2;
			*ret_f = f2;
			return;

		}

	}

	/*
	*PolyRefine: it will do the factor refinement for one object a^e and one list B
	*input a, e, B which are the input for FactorRefinement algorithm
	*input *ret_l,*ret_m,*Q,*S, should be empty pointers, the result of the algorithm will store in those pointers
	*ret_l^ret_m is the refinement result for a^e
	*S is the refinement result for B
	*Q is the gcd parts for FactorRefinement algorithm 
	*/
	template <class Domain>
	static void MergeRefinePolySeq(const Domain& a, int e,const std::vector<Factor<Domain>>& B, Domain* ret_l, int* ret_m, std::vector<Factor<Domain>>* ret_Q,std::vector<Factor<Domain>>* ret_S){
		Domain l0 = a;
		int m0 = e;

		//cout << "a in 2: " << a << endl;
		//cout << "e in 2: " << e << endl;
		std::vector<Factor<Domain>> Q;
		std::vector<Factor<Domain>> S;


		//cout << "B.size(): " << B.size() << endl;

		for (int i = 0; i < B.size(); ++i)
		{

			Domain li;
			int mi;
			std::vector<Factor<Domain>> Gi;
			Domain di;
			int wi;


			PolyRefine(l0,m0,B[i].first,B[i].second,&li,&mi,&Gi,&di,&wi);

			//cout << "li in 2: " << li << endl;
			//cout << "di in 2: " << di << endl;

			Q.insert(Q.end(),Gi.begin(),Gi.end());
			if(di.isOne()==false){
				//cout << "di in 2: " << di << endl;
				//cout << "wi in 2: " << wi << endl; 
				S.push_back(std::make_pair(di,wi));
			}

			l0 = li;
			m0 = mi;

		}



		*ret_l = l0;
		*ret_m = m0;


		//cout << "ret_l in 2: " << *ret_l << endl;
		//cout << "ret_m in 2: " << *ret_m << endl;
		ret_Q->insert(ret_Q->end(),Q.begin(),Q.end());
		ret_S->insert(ret_S->end(),S.begin(),S.end());
		return;
	}
	/*
	*PolyRefine: it will do the factor refinement for two lists A and B
	*input A and B which are the input for FactorRefinement algorithm
	*input *ret_L,*ret_Q,*ret_S should be empty pointers, the result of the algorithm will store in those pointers
	*ret_L is the refinement result for A
	*ret_S is the refinement result for B
	*ret_Q is the gcd parts for FactorRefinement algorithm 
	*/
	template <class Domain>
	static void MergeRefineTwoSeq(const std::vector<Factor<Domain>>& A, const std::vector<Factor<Domain>>& B,std::vector<Factor<Domain>>* ret_L, std::vector<Factor<Domain>>* ret_Q,std::vector<Factor<Domain>>* ret_S){
		std::vector<Factor<Domain>> L,Q,S;
		S = B;
		for (int i = 0; i < A.size(); ++i)
		{
			Domain l;
			int m;
			std::vector<Factor<Domain>> tempQ,tempS;
			MergeRefinePolySeq(A[i].first,A[i].second,S,&l,&m,&tempQ,&tempS);
			S = tempS;
			Q.insert(Q.end(),tempQ.begin(),tempQ.end());
			if(!l.isOne()){
				L.push_back(std::make_pair(l,m));
				//cout << "l: " << l << endl;
				//cout << "m: " << m << endl;

			}
		}

		ret_L->insert(ret_L->end(),L.begin(),L.end());
		ret_Q->insert(ret_Q->end(),Q.begin(),Q.end());
		ret_S->insert(ret_S->end(),S.begin(),S.end());

		//(3^1)*(2^2)
		return;
	}

	/*
	*This function is similar with MergeRefineTwoSeq, it will call MergeRefineTwoSeq first
	*then it check the result of MergeRefineTwoSeq, if the result is an empty vector, it will change the vector to one
	*
	*PolyRefine: it will do the factor refinement for two lists A and B
	*input A and B which are the input for FactorRefinement algorithm
	*input *ret_L,*ret_Q,*ret_S should be empty pointers, the result of the algorithm will store in those pointers
	*ret_L is the refinement result for A
	*ret_S is the refinement result for B
	*ret_Q is the gcd parts for FactorRefinement algorithm 
	*/
	template <class Domain>
	static void MergeRefineTwoSeqEmptyToIdentity(const std::vector<Factor<Domain>>& A, const std::vector<Factor<Domain>>& B,
		std::vector<Factor<Domain>>* ret_L, std::vector<Factor<Domain>>* ret_Q,std::vector<Factor<Domain>>* ret_S){

		MergeRefineTwoSeq(A,B,&*ret_L,&*ret_Q,&*ret_S);

		Domain one;
		one.one();
		if(ret_L->size() == 0){
			ret_L->push_back(std::make_pair(one,1));
		}
		if(ret_Q->size() == 0){
			ret_Q->push_back(std::make_pair(one,1));
		}	
		if(ret_S->size() == 0){
			ret_S->push_back(std::make_pair(one,1));
		}

		return;
	}
};



//template <class Domain>
//void PolyRefine(const Domain& a, int e, const Domain& b, int f, Domain* ret_l, int* ret_e, std::vector<Factor<Domain>>* ret_G, Domain* ret_r, int* ret_f);










#endif