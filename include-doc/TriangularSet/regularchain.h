#ifndef _REGULARCHAIN_H_
#define _REGULARCHAIN_H_

#include "triangularset.h"

using std::endl;
using std::cout;

template <class Ring>
class RegularChain : public BPASRegularChain, public TriangularSet<Ring> {
	protected:
		using TriangularSet<Ring>::n;
		using TriangularSet<Ring>::set;
		using TriangularSet<Ring>::names;
		

		using TriangularSet<Ring>::algVariables;
		using TriangularSet<Ring>::stronglyNormalized;
		using TriangularSet<Ring>::isEmptyRC;

		//std::vector<std::string> algVariables;
		//bool isEmptyRC = true;
		bool isKnownToBePrime= true;
		bool isKnownToBeSquareFree = false;
		//bool stronglyNormalized = true;
		bool optLazard=false;


		

	public:

		/**
		 * Construct with number of polynomials and variable names
		 *
		 * @param s: The number of polynomials
		 * @param xs; Variable names
		 **/
		/*	RegularChain<Ring> (TriangularSet<Ring> tri) {
		  
		  n=tri.numberOfVariables();
		  
		  set = new SparseMultivariatePolynomial<Ring>[n];
		  //algVariables = new int[n];
		  std::vector<std::string> xs = tri.variables();
		    names = new std::string[n];
		    for (int i = 0; i < n; ++i){
		      names[i] = xs[n-i-1];
		      //set[i] = tri.select(names[i]);
		      //algVariables[i]=0;
		      *this+=tri.select(names[i]);
		    }
		    
		    
		 exit(1);
		  
		 }*/

		/**
		 * Default constructor
		 *
		 * @param
		 **/
		RegularChain<Ring> () { n = 0; }
		/**
		 * Construct with number of polynomials and variable names
		 *
		 * @param s: The number of polynomials
		 * @param xs; Variable names
		 **/
		RegularChain<Ring> (int s, std::vector<std::string> xs) {
                        if (s < 1) { n = 0; }
                        else {
                                n = s;
                                set = new SparseMultivariatePolynomial<Ring>[n];
				//algVariables = new int[n];
                                if (xs.size() == s) {
                                        names = new std::string[n];
                                        for (int i = 0; i < n; ++i){
                                                names[i] = xs[n-i-1];
						//algVariables[i]=0;
					}
					

                                }
                                else {
                                        std::cout << "BPAS: error, RegularChain<Ring> must be defined with " << s << " variables." << std::endl;
                                        exit(1);
                                }
                        }
		}
		
		/**
		 * Copy constructor
		 *
		 * @param a: A regular chain
		 **/
		RegularChain<Ring> (const RegularChain<Ring>& a) {
		  if (a.n) {
		    n = a.n;
		    set = new SparseMultivariatePolynomial<Ring>[n];
		    names = new std::string[n];
		    std::copy(a.set, a.set+n, set);
		    std::copy(a.names, a.names+n, names);
		    
		    
		    
		    
		    
		    algVariables = a.algVariables;
		    isEmptyRC = a.isEmptyRC;
		    stronglyNormalized = a.stronglyNormalized;
		    isKnownToBePrime = a.isKnownToBePrime;
		    isKnownToBeSquareFree = a.isKnownToBeSquareFree;
		    optLazard = a.optLazard;
		  }
		}
		
		/**
		 * Deconstructor
		 *
		 * @param
		 **/
		//~RegularChain<Ring>() {
		//	if (n) {
		//			delete [] set;
		//		delete [] names;
		//		delete [] algVariables;
		//		}
		//}

		/**
		 * Overload operator =
		 *
		 * Using TriangularSet<Ring>'s
		 *
		 * @param
		 **/
		//using TriangularSet<Ring>::operator=;
		inline RegularChain<Ring>& operator= (RegularChain<Ring> a) {
		  if (this != &a) {
		    if (n) {
		      delete [] set;
		      delete [] names;
		      
		    }
		    n = a.n;
		    set = new SparseMultivariatePolynomial<Ring>[n];
		    names = new std::string[n];
		    std::copy(a.set, a.set+n, set);
		    std::copy(a.names, a.names+n, names);
		    algVariables = a.algVariables;
		    isEmptyRC = a.isEmptyRC;
		    stronglyNormalized = a.stronglyNormalized;
		    isKnownToBePrime = a.isKnownToBePrime;
		    isKnownToBeSquareFree = a.isKnownToBeSquareFree;
		    optLazard = a.optLazard;
		  }
		  return *this;
		}





		inline RegularChain<Ring>& operator= (std::string &str) {
		  
		  std::vector<std::string > vars;
		  if(n>0){
		    
		    
		    for(int i=n-1 ; i>=0 ; i--){
		      //std::cout<< names[i] << std::endl;
		      vars.push_back(names[i]);
		    }
		    
		  }else{
		    std::cout<<"BPAS:: error: a list of variables has to be specified"<<std::endl;
		  }
		  
		  
		  
		  std::vector<SparseMultivariatePolynomial<Ring> > polys = this->parser(str);
		  for(int i =0; i<polys.size(); i++)
		    (*this) += polys[i];
		  
		  
		  return *this;
		  
		  
		}





		/**
		 * Overload operator +
		 *
		 * @param smp: A sparse multivariate polynomial
		 **/
		inline RegularChain<Ring> operator+ (SparseMultivariatePolynomial<Ring>& smp) {
			RegularChain<Ring> r (*this);
			return (r += smp);
		}
		/**
		 * Overload operator +=
		 *
		 * @param smp: A sparse multivariate polynomial
		 **/
		inline RegularChain<Ring>& operator+= (SparseMultivariatePolynomial<Ring> smp) {
		 
                        if (smp.isConstant()) {
                                std::cout << "BPAS: error, cannot adding a constant to RegularChain<Ring>" << std::endl;
                                exit(1);
                        }
			int v = smp.numberOfVariables();
                        if (v > n || !v) {
                                std::cout << "BPAS: error, cannot adding a polynomisl of " << v << " variables to a RegularChain<Ring> of " << n << " variables." << std::endl;
                                exit(1);
                        }
			// the array pos finds the position of variables of smp polynomials in the array names
			int* pos = new int[v];
			std::vector<std::string> xs = smp.variables();
			for (int i = 0; i < v; ++i) {
				pos[i] = -1;
				for (int j = n-1; j > -1; --j) {
					if (xs[i] == names[j]) {
						pos[i] = j;
						break;
					}
				}
			}
			bool isIt = 0;
			int k = pos[0];
			if (k < 0) { isIt = 1; }
			for (int i = 1; i < v; ++i) {
				if (pos[i] >= 0 && pos[i] < k) { k = pos[i]; }
				else { isIt = 1; break; }
			}
			if (isIt) {
				std::cout << "BPAS: error, RegularChain<Ring>[";
				for (int i = n-1; i > -1; --i) {
					std::cout << names[i];
					if (i) { std::cout << ", "; }
				}
				std::cout << "] cannot add a polynomial in Ring[";
				for (int i = 0; i < xs.size(); ++i) {
					std::cout << xs[i];
					if (i < xs.size()-1) { std::cout << ", "; }
				}
				std::cout << "]." << std::endl;
				exit(1);
			}
			k = pos[0];
			delete [] pos;

			// here we check whether or not there is a polynomial in rc with the same  main variable as of smp. 
                        if (!set[k].isZero()) {
				std::cout << "BPAS: error, a polynomial with the leading variable " << xs[0] << " already exists in RegularChain<Ring>[";
				for (int i = n-1; i > -1; --i) {
					std::cout << names[i];
					if (i) { std::cout << ", "; }
				}
				std::cout << "]." << std::endl;
				exit(1);
                        }
			// smp is added to polynomials of rc
			if(!smp.isZero()){
			  updateRegularChainsStates(smp,k);
			  //The following line was moved into updateRegularChainsStates function 
			  //set[k] += smp;
			  
			  // the string of main variables will become updated here
			  algVariables.push_back(names[k]);
			  //std::cout<<"AlgVars"<< algVariables.back()<<std::endl;
			  // check wether or not rc is strongly normalized
			  

			}
			// here we check whether or not rc union smp forms a regular chain or not. If not return an error.
			if (k < 2 && n > 1 && !set[0].isZero() && !set[1].isZero()) {
				SparseMultivariatePolynomial<Ring> l = set[1].leadingCoefficientInVariable(names[1]);

				if (!l.isConstant()) {
					SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> > p = l.convertToSUP(names[0]);
					SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> > q = set[0].convertToSUP(names[0]);

					SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> > r = q.resultant(p);

					if (r.isZero()) {
						std::cout << "BPAS: error, it is not a regular chain after adding a polynomial in RegularChain<Ring>." << std::endl;
						exit(1);
					}
				}
			}

			return *this;
		} 

		/**
		 * Get the number of variables
		 *
		 * @param
		 **/
		inline int numberOfVariables() {
			return n;
		}
		/**
		 * Get variable names
		 *
		 * @param
		 **/
		inline std::vector<std::string> variables() {
                        std::vector<std::string> xs;
                        for (int i = n-1; i > -1; --i)
                                xs.push_back(names[i]);
                        return xs;
		}
		/**
		 * Get the polynomials
		 *
		 * @param
		 **/
		
		inline std::vector<SparseMultivariatePolynomial<Ring> > polynomials(){
		  std::vector<SparseMultivariatePolynomial<Ring> > polys;
		  for(int i=0; i<n; i++){
		    if(!set[i].isZero()){
		      
		      polys.push_back(set[i]);
		    }
		  }
		  return polys;
		}



		
		


		/**
		 * Select a polynomial given the leading variable
		 *
		 * Using TriangularSet<Ring>'s
		 *
		 * @param 
		 **/
		using TriangularSet<Ring>::select;
		/**
		 * The regular chain under a variable
		 *
		 * @param x: The variable name
		 **/
		inline RegularChain<Ring> under(std::string x) {
                        int k = -1;
                        for (int i = 0; i < n; ++i) {
                                if (x == names[i]) {
                                        k = i;
                                        break;
                                }
                        }
                        if (k <= 0)
                                return RegularChain<Ring>();
                        else {

			 
			  RegularChain<Ring> r;
			  r.n = k ;
			  r.set = new SparseMultivariatePolynomial<Ring>[r.n];
			  r.names = new std::string[r.n];
			  std::copy(names, names+k, r.names);
			  
			  for(int i=0;i<k;i++){
			    if(!this->set[i].isZero())
			      r+=this->set[i];
			   
			  }
			   
			  return r;
                        }
		}	





/*




		inline RegularChain<Ring> under (std::string x) {
		  /*  TriangularSet<Ring> tri = this->under(x);
		  return *this;
		  exit(0);
		  //TriangularSet<Ring> tri = this->under(x);
		
		RegularChain<Ring> r(tri);
		return r;
		//return tri;
		*/
  /*int k = -1;
                        for (int i = 0; i < n; ++i) {
                                if (x == names[i]) {
                                        k = i;
                                        break;
                                }
                        }
			std::cout<< k<<std::endl;
                        if (k <= 0)
                                return RegularChain<Ring>();
                        else {
			  //std::cout<< "I am here"<< endl;
			  RegularChain<Ring> r;
			  r.n = k-1;
			  r.set = new SparseMultivariatePolynomial<Ring>[r.n];
			  r.names = new std::string[r.n];
			  std::copy(set, set+k-1, r.set);
			  //for(int i=0;i<k-1;i++){
			  //  r+=this->set[i];
			  //}
			  std::cout<< "the poly"<< set[k-1]<<std::endl;
			  std::copy(names, names+k-1, r.names);
			  return r;
			}
		}*/
		
		/**
		 * Get whether or not regular chain is empty
		 *
		 * @param
		 **/
		inline bool isEmpty() {
		  return isEmptyRC;
		  
		  // old code;
		  //for(int i=0;i<n;i++){
		  //  if(!set[i].isZero())
		  //    return false;
		  //}
		  //return true;
		}

		/**
		 * Get the value of isKnownToBePrime
		 *
		 * @param
		 **/
		inline bool isPrimeRegularChain() {
		  
		  return isKnownToBePrime;
		}


		/**
		 * The triangular set upper a variable 
		 *
		 * @param x: The variable name
		 **/
		inline RegularChain<Ring> upper(std::string x) {

		 
		  int k = -1;
		  for (int i = 0; i < n; ++i) {
		    if (x == names[i]) {
		      k = i; 
		      break;
		    }
		  }
		  if (k < 0 || k == n-1)
		    return RegularChain<Ring>();
		  else {
		    
		    std::vector<std::string> vars;
		    for(int i=k+1; i<n; ++i){
		      std::vector<std::string> varsTemp = this->set[i].variables();
		      vars =  setUnion(vars, varsTemp);
		      
		    }
		   
		    std::vector<std::string> rcVars(names, names+n);
		    std::vector<std::string> namesVec;
		    
		    for(int i=n-1; i>=0; --i){
		      std::vector<std::string>::iterator it = std::find(vars.begin(), vars.end(), rcVars[i]);
		     
		      if(it!=vars.end()){
			namesVec.push_back(names[i]);
			
		      }
		    }
		    
		    
		    RegularChain<Ring> r(namesVec.size(), namesVec);
		    
		    for(int i=k+1 ; i<n ; i++){
		      if(!this->set[i].isZero()){
			r+=this->set[i];
		      }
		    }
		    
		    return r;
		  }
		}
		
		
		/**
		 * Pseudo division
		 * Return the remainder
		 *
		 * @param p: An input polynomial
		 * @param quo: The quotients
		 * @param c: The constant multiplied to the input polynomial
		 **/
		using TriangularSet<Ring>::pseudoDivide;
		/**
		 * Monic division
		 * Return the remainder
		 *
		 * @param p: An input polynomial
		 * @param quo: The quotients
		 **/
		using TriangularSet<Ring>::monicDivide;
		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param a: A regular chain
		 **/
		/*inline friend std::ostream& operator<< (std::ostream& out, RegularChain<Ring>& a) {
			bool isFirst = 0;
			out << "[";
                        for (int i = 0; i < a.n; ++i) {
				if (isFirst) { out << ", "; }
                                out << a.set[i];
				isFirst = 1;
                        }
			out << "]";
                        return out;
			}*/
		
		/**
		 * Overload stream operator ==
		 *
		 *
		 * @param a: A triangular set
		 **/
		inline bool operator== (RegularChain<Ring>& a) {
		  if(n==a.n){
		    int i;
		    for(i=0; i<n;i++){
		      if(set[i]!=a.set[i])
			break;
		    }
		    if(i==n){
		      return 1;
		    }else{
		      
		      for(int j=0; j<n; j++){
			if(!this->pseudoDivide(a.set[j]).isZero())
			  return 0;
		      }

		      for(int j=0; j<n; j++){
			if(!a.pseudoDivide(this->set[j]).isZero())
			  return 0;
		      }
		      return 1;
		    }
		  }else
		    return 0;
		}
		





		/**
		 *
		 * isStronglyNormalized 
		 * Get the value of strongly normalized
		 *
		 **/
		inline bool isStronglyNormalized(){
		  return stronglyNormalized;
		}

		inline bool  isSquareFree(){
		  return isKnownToBeSquareFree;
		}
		/**
		 * Dimension of a regular chain
		 *
		 **/
		inline int dimension(){
		  
		  return (n - algVariables.size());
		  
		  
		}
		inline void updateRegularChainsStates(SparseMultivariatePolynomial<Ring> smp, int k){
		  
		  
		  //Check whether or not by adding smp to rc, rc remains strongly normalized
		  SparseMultivariatePolynomial<Ring> init = smp.initial();
		  //std::cout<<"Print smp:  "<<smp<< "  The leading variable"   <<   smp.leadingVariable()<<std::endl;
		  //std::cout<< "Print the initial smp:  " << init <<std::endl;
		  std::vector<std::string> xs = init.variables();
		  //std::cout<<"The size of init vars:  "<<xs.size()<<std::endl;
		  
		  //if(xs.size()!=0){
		  //  std::cout<<"The var in init vars:  "<<xs.at(0)<<std::endl;

		  //}
		  if(stronglyNormalized){
		    //if(xs.size()==0){
		      
		    //}	    
		    std::sort(xs.begin(),xs.end());
		    std::sort(algVariables.begin(),algVariables.end());
		    
		    std::vector<std::string> It;
		    std::set_intersection(xs.begin(),xs.end(),  algVariables.begin(),algVariables.end(), back_inserter(It));
		    
		    if(It.size()!=0 && It[0]!="9" && It[0]!="1"){
		      
		      stronglyNormalized = false;
		      //std::cout<<"strongly normalized state:  "<< stronglyNormalized<< std::endl;
		    }
		    
		    for (int i=0;i<n;++i){
		      std::vector<std::string> xxs = (this->set[i].initial()).variables();
		      std::sort(xxs.begin(),xxs.end());
		      std::set_intersection(xs.begin(),xs.end(),  algVariables.begin(),algVariables.end(), back_inserter(It));
		      
		      
		      
		      if(It.size()!=0&& It[0]!="9" && It[0]!="1"){
			stronglyNormalized = false;
			//std::cout<<" here 1"<<std::endl;
		      }
		    }
		  }
		  
		  if(isKnownToBePrime){

		    if(smp.leadingVariableDegree()!=1){
		      isKnownToBePrime = false;
		    }
		    //for(int i=0; i<=k; i++){
		    //if(!set[i].isZero()){
		    //	if(smp.leadingVariableDegree()!=1){
		    //	  isKnownToBePrime = false;
		    //	}
		    // }
		    //}
		  }


		  if(isKnownToBeSquareFree){
		    
		    SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > supsmp = smp.convertToSUP(names[k]);
		    
		    SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> >  supsmpDiff = supsmp;
		    supsmpDiff.differentiate(1);

		    //supsmp.monomialBasisSubresultantChain(supsmpDiff);
		    //regularGcd_by_src(supsmp,supsmpDiff, src );
		    // THE FOLLOWING LINE HAS TO BE CHANGED
		    set[k] += smp;
		    if(isEmptyRC)
		      isEmptyRC = false;
		    
		  }else{
		    
		    set[k] += smp;
		    if(isEmptyRC)
		      isEmptyRC = false;
		  }
		  
		  
		}
		




		

		

		
		/**
		 * principleSubresultantCoefficientOfIndex
		 *
		 * @param q: The other sparse univariate polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring> principleSubresultantCoefficientOfIndex(int k, std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src){
		  if(!k){ 
		    
		    return src[0].coefficient(0);
		  }else if(src[k].degree()<k){
		    SparseMultivariatePolynomial<Ring >poly;
		    poly.zero();
		    return poly;
		  }else if(!src[k].isZero()){
		    return  src[k].coefficient(src[k].degree()); 
		  }else{  
		    SparseMultivariatePolynomial<Ring >poly;
		    poly.zero();
		    return poly;
		  }
		}

		/**
		 *  subresultantOfIndex
		 *
		 * @param k: an index
		 * @param src: subresultant chains
		 **/
		inline SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > >  subresultantOfIndex(int k, std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src)
		{
		  return src[k];
		}
		
		
		inline int regularGcd(SparseMultivariatePolynomial<Ring> p, SparseMultivariatePolynomial<Ring> q ,std::vector<RegularChain<Ring> > &rc_out, std::vector<SparseMultivariatePolynomial<Ring> > &polys , std::vector<RegularChain<Ring> > &resultDrop){


		  cout<< "p  "<<p<<endl;
		  cout<< "q  "<<q<<endl;
		  this->Display();
		  
		  SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > upoly1 = p.convertToSUP(p.leadingVariable());
		  
		  SparseUnivariatePolynomial<SparseMultivariatePolynomial<RationalNumber> > upoly2 = q.convertToSUP(q.leadingVariable());
		  
		  
		  
		  if(p.variables().size()==1 && q.variables().size()==1){
		    SparseMultivariatePolynomial<Ring >  g = p.gcd(q);
		    
		    rc_out.push_back(*this);
		    polys.push_back(g);
		    return 1;
		  }
		  std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src = upoly2.monomialBasisSubresultantChain(upoly1);
		  SparseMultivariatePolynomial<Ring> r = this->principleSubresultantCoefficientOfIndex(0,src);
		  
		  cout<< "r  "<<r<<endl; 
		  // sted::vector<RegularChain<Ring> > resultDrop;
		  
		  if(this->isEmpty()){
		    if(!r.isZero()){
		      
		      rc_out.push_back(*this);
		      polys.push_back(r);
		      return 1;
		    }
		    
		  }
		  std::vector<RegularChain<Ring> > newRcOut;
		  std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
		  this->regularizeDim0(r, newRcOut, newPolys);

		  cout<< "in regulargcd size of out"  << newRcOut.size();
		  
		  for(int i=0; i<newPolys.size();++i){
		    
		   SparseMultivariatePolynomial<Ring> f = newPolys.at(i);
		   RegularChain<Ring> C =  newRcOut.at(i);
		   
		   cout<< "fpoly   "<< f<<endl;
		   cout<< "C  "<<endl;
		   newRcOut.at(i).Display();

		   //first ifstatment is commneted out for the moment
		   if(C.dimension()<this->dimension())
		     resultDrop.push_back(C);
		   else if(f.isZero()){
		     std::vector<RegularChain<Ring> > rc_outTemp;
		     std::vector<SparseMultivariatePolynomial<Ring> > polysTemp;
		     std::vector<RegularChain<Ring> > resultDropTemp;
		     
		     C.regularGcd_by_src(p,q,src, rc_outTemp, polysTemp, resultDropTemp);
		     rc_out.insert(rc_out.end(),rc_outTemp.begin(), rc_outTemp.end());
		     polys.insert(polys.end(),polysTemp.begin(), polysTemp.end());
		     
		     resultDrop.insert(resultDrop.end(), resultDropTemp.begin(),resultDropTemp.end() );

		     for(int j=0; j<polysTemp.size(); ++j){
		       cout<< "polyTemp: "<< polysTemp.at(j) <<endl;
		       cout<<"rcTemp  "<<endl;
		       rc_outTemp.at(i).Display();
		       }
		     cout<< "safe"<<endl;
		     cout<< "newPolys.size()" << newPolys.size()<<endl;
		   }else{  
		     rc_out.push_back(C);
		     polys.push_back(f);
		     cout<< "The end"<< endl;
		   }

		    
		  }
		  cout<< "The size of Output polys  "<< polys.size()<<endl;
		  cout<< "The size of Output rcs  "<< rc_out.size()<<endl;
		  cout<< "exit regulargcd"<<endl;
		  return 1;
		  

		}




		/**
		 *  regularGcd_by_src
		 *
		 * @param p: a multivariate polynomial
		 * @param q: a multivariate polynomial
		 * @param src: subresultant chains
		 * @param rc_out: a vector of regular chains rc_1,...rc_s which form a triangular decomposition for the current object (regular chain)
		 * @param polys: a vector of polynomials g_1,...,g_s, where each g_i is either zero or equal to regular gcd of polynomials p and q modulo rc_i  
		 **/
		inline int regularGcd_by_src(SparseMultivariatePolynomial<Ring> p, SparseMultivariatePolynomial<Ring> q, std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src, std::vector<RegularChain<Ring> > &rc_out, std::vector<SparseMultivariatePolynomial<Ring> > &polys , std::vector<RegularChain<Ring> > &resultDrop){

		  cout<< "entered regularGcd_by_src "<< endl;
		  

		  cout<< "p  "<<p<<endl;
		  cout<< "q  "<<q<<endl;
		  this->Display();


		  std::vector<RegularChain<Ring> > tasks_rc;
		  tasks_rc.push_back(*this);
		  std::vector<int > tasks_indice;
		  tasks_indice.push_back(1);
		  
		  int d = src.at(src.size()-2).degree();
		  
		  while(tasks_rc.size()!=0){
		    RegularChain<Ring> rc = tasks_rc[0];
		    tasks_rc.erase(tasks_rc.begin());
		    
		    int indice = tasks_indice[0];
		    tasks_indice.erase(tasks_indice.begin());
		    SparseMultivariatePolynomial<Ring> f;
		    if(indice>d || indice == q.leadingVariableDegree()){
		      polys.push_back(q);
		      rc_out.push_back(rc);
		      continue;
		    }else if(indice<d){
		      f = principleSubresultantCoefficientOfIndex(indice,src);
		      
		    }else{
		      f = p.initial();
		    }
		    cout<< "f"<< f<<endl;
		    std::vector<RegularChain<Ring> > newRcOut;
		    std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
		    rc.regularizeDim0(f, newRcOut, newPolys);
		   

		    cout<< "result size()" << newRcOut.size()<<endl;
  		    for(int i=0; i<newRcOut.size(); i++){

		      cout<< "g" <<  newPolys.at(i)<<endl;
		      cout<< "rc_C: "<<endl;
		      newRcOut.at(i).Display();
		      if(newRcOut.at(i).dimension()<rc.dimension())
			resultDrop.push_back(newRcOut.at(i));
		      else if(newPolys[i].isZero()){  
			
			tasks_rc.push_back(newRcOut[i]);
			tasks_indice.push_back(indice+1);
			
		      }else if(indice==d){
			rc_out.push_back(newRcOut[i]);
			polys.push_back(p);
			
		      }else{ 
			rc_out.push_back(newRcOut[i]);
			SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> >  tempUiniPoly = this->subresultantOfIndex(indice,src);
			cout<< "indice"<<indice<<endl;

			cout<<" tempUiniPoly   "  <<  tempUiniPoly  <<endl;
			SparseMultivariatePolynomial<Ring> polyTemp(tempUiniPoly);
			polys.push_back(polyTemp); 
	cout<< "here"<<endl;
		      }
		      
		      
		    }
		    
		    
		    }
		  cout<< "Exit regularGcd_by_src "<< endl;
		  return(1);
		  
		}
		inline int regularizeDim0(SparseMultivariatePolynomial<Ring> poly, std::vector<RegularChain<Ring> > &newRcOut, std::vector<SparseMultivariatePolynomial<Ring> > &newPolys ){
		  cout<< "enter RegularizeDim0:  "<<endl;
		  
		  cout<< "poly  "<< poly<<endl;
		  this->Display();

		  SparseMultivariatePolynomial<Ring> p = this->removeZero(poly);
		 
		  cout<< "poly  after remove zero"<< p<<endl;
		  /************************/
		  //The following lines are a temporary solution for the existing problem in monicdevide and PseudoDivide functions in which the order of variables would get changed...
		  std::stringstream ss;
		  ss << p;
		  std::string s = ss.str();
		  p = s;
		  /************************/


		  
		  if(p.isZero()){
		    SparseMultivariatePolynomial<Ring> zeroPoly;
		    zeroPoly.zero();
		    newPolys.push_back(zeroPoly);
		    		    
		    newRcOut.push_back(*this);
		    cout<< "Exit RegularizeDim0: "<< endl;
		    return(1);
		    
		  }else if(p.isConstant()|| this->isEmpty()){
		    
		    newRcOut.push_back(*this);
		    newPolys.push_back(poly);
		    cout<< "Exit RegularizeDim0: "<< endl;
		    return(1);
		  }else if(this->isPrimeRegularChain()){
		    if((this->pseudoDivide(p)).isZero()){
		      SparseMultivariatePolynomial<Ring> zeroPoly;
		      zeroPoly.zero();
		      newPolys.push_back(zeroPoly);
		      newPolys.push_back(zeroPoly);
		      newRcOut.push_back(*this);
		      cout<< "Exit RegularizeDim0: "<< endl;
		      return(1);  
		      
		    }else{
		      newRcOut.push_back(*this);
		      newPolys.push_back(poly);
		      cout<< "Exit RegularizeDim0: "<< endl;
		      return(1);
		    }

		  }

		  if((this->pseudoDivide(p)).isZero()){
		    SparseMultivariatePolynomial<Ring> zeroPoly;
		    zeroPoly.zero();
		    newPolys.push_back(zeroPoly);
		    newPolys.push_back(zeroPoly);
		    newRcOut.push_back(*this);
		    cout<< "Exit RegularizeDim0: "<< endl;
		    return(1);  
		    
		  }

		  
		  std::string mvar = p.leadingVariable();
		  
		  SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > >  up= p.convertToSUP(mvar);
		  
		 
		  cout<< "passed"<<endl;

		  RegularChain<Ring> under_rc = this->under(mvar);
		  
		  std::string * xx;
		  xx= std::find(names, names+n,mvar);
		  int indice = std::distance(names,xx);

		  SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > >  uq= set[indice].convertToSUP(mvar);
		  cout<< "uq"<< uq<<endl;

		  
		  std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src = up.monomialBasisSubresultantChain(uq);

		  SparseMultivariatePolynomial<Ring >  res = principleSubresultantCoefficientOfIndex(0,src);



		  cout<< "res:  "<<res<<endl;
		  
		  std::vector<RegularChain<Ring> > newRcOutTemp;
		  std::vector<SparseMultivariatePolynomial<Ring> > newPolysTemp;
		  under_rc.Display();
		  
		  under_rc.regularizeDim0(res, newRcOutTemp, newPolysTemp);
		 
		  for(int i=0; i<newPolysTemp.size();i++){
		    
		  RegularChain<Ring> rc = newRcOutTemp[i];
		  SparseMultivariatePolynomial<Ring> pNew= newPolysTemp[i];
		 

		  rc.Display();
		  cout<<"p" << p<<endl;

		  if(!pNew.isZero()){	    
		    newRcOut.push_back(*this);
		    newPolys.push_back(p);
		   
		  }
		  else{
		    
		    std::vector<RegularChain<Ring> > rcOutTemp;
		    std::vector<SparseMultivariatePolynomial<Ring> > polysTemp;
		    std::vector<RegularChain<Ring> > resultDrop;
		    rc.regularGcd_by_src(p, set[indice], src , rcOutTemp, polysTemp,resultDrop);
		    cout<< "PASS"<<endl;
		    for(int j=0;j<polysTemp.size();j++){
		      SparseMultivariatePolynomial<Ring> g = polysTemp[j];
		      cout<< "g : "<< g<<endl;
		      RegularChain<Ring> E = rcOutTemp[j];
		      cout<<"E  :"<<endl;
		      E.Display();
		      if(g.leadingVariableDegree()==set[indice].leadingVariableDegree()){
			
			/*   Build a new regular chain   by E C_v and rc_upper*/
			std::vector<std::string> xs(names, names+n);
			std::reverse(xs.begin(),xs.end());
			RegularChain<Ring> rcNew(n,xs);
			std::vector<SparseMultivariatePolynomial<Ring> > polySet = E.polynomials();
			for(int ll=0; ll<polySet.size();++ll ){
			  rcNew+=polySet[ll];
			}	
			rcNew+=g;
			RegularChain<Ring> rc_upper = this->upper(mvar);
			for(int ll=0; ll<rc_upper.n;++ll){
			  if(!rc_upper.set[ll].isZero())
			    rcNew+=rc_upper.set[ll];
			}

			newRcOut.push_back( rcNew  );
			SparseMultivariatePolynomial<Ring> zeroPoly;
			zeroPoly.zero();
			newPolys.push_back(zeroPoly);
			continue;
		      }
			cout<< "uq "<<uq<<endl;		 
		      SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > ug = g.convertToSUP(mvar);
		      

		      SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > q = uq/ug;
		      //SparseMultivariatePolynomial<Ring > q = this->set[indice]/g;  
		      
		      cout<< "q  "<<q<<endl;
		      /*   Build a new regular chain   by E g and rc_upper*/
		      std::vector<std::string> xs(names, names+n);
		      std::reverse(xs.begin(),xs.end());
		      

		      RegularChain<Ring> rcNew(n,xs);
		      
		      std::vector<SparseMultivariatePolynomial<Ring> > polySet = E.polynomials();
		      

		      RegularChain<Ring> rc_upper = this->upper(mvar);

		      for(int ll=0; ll<indice;ll++ ){
			rcNew+=polySet[ll];
		      }
		      rcNew+= g;
		     
		      for(int ll=0; ll<rc_upper.n;ll++){
			rcNew+=rc_upper.set[ll];
		      }
		
		      newRcOut.push_back(  rcNew  );
		      SparseMultivariatePolynomial<Ring> zeroPoly;
		      zeroPoly.zero();
		      newPolys.push_back(zeroPoly);
		      cout<< "poly: "<< zeroPoly <<endl;
		      cout<< "rc   "<<endl;
		      rcNew.Display();
			
		      if(q.degree()>0){
			std::vector<RegularChain<Ring> > newRcOutTemp2;
			std::vector<SparseMultivariatePolynomial<Ring> > newPolysTemp2;
			  
			  /*   Build a new regular chain   by E g and rc_upper*/
			std::vector<std::string> xs(names, names+n);
			std::reverse(xs.begin(),xs.end());
			RegularChain<Ring> rcNew2(n,xs);
			std::vector<SparseMultivariatePolynomial<Ring> > polySet = E.polynomials();
			for(int ll=0; ll<indice;ll++ ){
			  rcNew2+=polySet[ll];
			}

			SparseMultivariatePolynomial<Ring > tempPoly(q);
			rcNew2+= tempPoly;
			//rcNew+= q;
			for(int ll=0; ll<rc_upper.n;ll++){
			  rcNew2+=rc_upper.set[ll];
			}
		
			

			if(rcNew2.isSquareFree()){
			  newRcOut.push_back(rcNew2);
			  newPolys.push_back(p);

			}else{
			  cout<< "p added over q"<< p<<endl;
			  rcNew2.regularizeDim0(poly, newRcOutTemp2, newPolysTemp2 );
			    
			  for (int k=0; k<newPolysTemp2.size(); k++ ){
			    cout<< "poly  "<< newPolysTemp2.at(k);
			    cout<< "rc  "<< endl;
			    newRcOutTemp2[k].Display();
			    newRcOut.push_back(newRcOutTemp2[k]);
			    newPolys.push_back(newPolysTemp2[k]);
			  }
			}
			cout<< "Exit Regularize: "<< endl;
			return(1);  
			}
			
			
		    }
		      
		  }
		  }
		  cout<< "Exit RegularizeDim0: "<< endl;
		  return(0);
		  
		}	
		
	

		inline std::vector<RegularChain<Ring> > innerTriangularizeKalkbrener(std::vector< SparseMultivariatePolynomial<Ring> > lp){
		 
		  std::vector<RegularChain<Ring> > results;
		  std::vector<SparseMultivariatePolynomial<Ring> >  cleanSet;
		  bool inconsistent = this->cleanSet(lp,cleanSet);

		  if(inconsistent){
		    return results;
		  }
		  if(cleanSet.size()==0){
		    results.push_back(*this);
		  }
		  
		  std::vector< SparseMultivariatePolynomial<Ring> > tasksLp;
		  std::vector<RegularChain<Ring> > tasksRC;
		  tasksLp.push_back(cleanSet);
		  tasksRC.push_back(*this);
		  
		  while(tasksRC.size()!=0){
		    RegularChain<Ring> rc =  tasksRC.at(0);
		    tasksRC.erase(tasksRC.begin());
		    std::vector< SparseMultivariatePolynomial<Ring> > lpoly = tasksLp.at(0);

		    if(lpoly.size()==0){
		      results.push_back(rc);
		    }else{
		      int k = this->minRittSet(lpoly); 
		      SparseMultivariatePolynomial<Ring> poly = lpoly.at(k);
		      lpoly.erase(k);
			
		      std::vector<RegularChain<Ring> > intersectRes = rc.intersect(poly);
			for(int j=0; j<intersectRes.size(); j++){
			  std::vector<SparseMultivariatePolynomial<Ring> > cleanSetNew;
			  bool inconsistent2 = intersectRes.at(j).cleaSet(lpoly, cleanSetNew);
			  if(inconsistent2)
			    continue;
			  tasksLp.push_back(cleanSetNew);
			  tasksRC.push_back(intersectRes.at(j));
			}
		      
		    }
		    
		  }
		  
		  return results;
		  
		}
		
		inline std::vector<RegularChain<Ring> > intersect(SparseMultivariatePolynomial<Ring> poly){
		  
		  std::vector<RegularChain<Ring> > results;
		  if(poly.isZero()){
		    results.push_back(*this);
		    return results;
		  }else if(poly.isConstant()){
		    return results;
		  }else if(this->pseudoDivide(poly)==0){
		    results.push_back(*this);
		    return results;
		  }
		  
		  std::vector< SparseMultivariatePolynomial<Ring> > inequations = this->listOfInitials();
		  
		  inequations = mergeGcdFreeFactorization(inequations);
		  std::vector< SparseMultivariatePolynomial<Ring> > tasksIneqs;
		  std::vector< SparseMultivariatePolynomial<Ring> >  tasksEqs;
		  this->intersectProjection(poly, inequations, tasksIneqs, tasksEqs);
		  results = this->intersectExtension(tasksIneqs, tasksEqs);
		  return results;
		  
		}

		/**
		 * intersectProjection
		 * @param poly is a primitive and squarefree polynomial
		 * @param inequations list of polynomials
		 **/
		inline void intersectProjection(SparseMultivariatePolynomial<Ring> poly, std::vector< SparseMultivariatePolynomial<Ring> > inequations, std::vector< SparseMultivariatePolynomial<Ring> > &tasksIneqs, std::vector< SparseMultivariatePolynomial<Ring> >  &tasksEqs ){
		  
		  std::vector< SparseMultivariatePolynomial<Ring> >  Facs = poly.factors();

		  std::vector<std::vector< SparseMultivariatePolynomial<Ring> > > tasksFacs;
		  std::vector<std::vector< SparseMultivariatePolynomial<Ring> > > tasksSets;
		  for(int i=0; i< Facs.size(); i++){
		    std::vector< SparseMultivariatePolynomial<Ring> > eqs;
		    tasksSets.push_back(eqs);
		    eqs.push_back(Facs.at(i));
		    tasksFacs.push_back(eqs);
		  }
		  

		  while(tasksFacs.size()!=0){
		    std::vector< SparseMultivariatePolynomial<Ring> > lp = tasksFacs.at(0);
		    tasksFacs.erase(tasksFacs.begin());
		    std::vector< SparseMultivariatePolynomial<Ring> > s = tasksSets.at(0);
		    tasksFacs.erase(tasksFacs.begin());

		    SparseMultivariatePolynomial<Ring> r = lp.back();
		    lp.erase(lp.end()-1);
		      //r = remove_gcd_poly_wrt_polylist;
		    if(r.isConstant())
		      continue;
		    else
		      lp.push_back(r);
		    
		    std::string mvar = r.leadingVariable();
		    if(!this->isAlgebraic(mvar)){
		      tasksIneqs.push_back(lp);
		      tasksEqs.push_back(s);
		      continue;
		    }

		    std::vector< SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > > src = r.convertToSUP(mvar).subresultantChain(this->select(mvar).convertToSUP(mvar));
		    s.insert(s.end(), src.begin(),src.end());
		    r = principleSubresultantCoefficientOfIndex(0);
		    if(r.isZero()){
		      tasksIneqs.push_back(lp);
		      tasksEqs.push_back(s);
		      continue;
		    }
		      
		    if(r.isConstant())
		      continue;

		    std::vector< SparseMultivariatePolynomial<Ring> > rFacs = r.factors(); 
		    for(int j=0; j< rFacs.size(); j++){
		      std::vector< SparseMultivariatePolynomial<Ring> > eqs;
		      tasksSets.push_back(s);
		      eqs.push_back(rFacs.at(j));
		      tasksFacs.push_back(eqs);
		    }
		    
		  }
		  
		}
		
		inline std::vector<RegularChain<Ring> > intersectExtension(std::vector< SparseMultivariatePolynomial<Ring> > tasksIneqs, std::vector< SparseMultivariatePolynomial<Ring> >  tasksEqs){
		  
		  std::vector<RegularChain<Ring> > result;
		  for(int j=0; j<tasksIneqs.size(); ++j){
		    std::vector< SparseMultivariatePolynomial<Ring> > lp = tasksIneqs.at(j);
		    std::vector< SparseMultivariatePolynomial<Ring> > s = tasksEqs.at(j);
		    
		    std::vector<RegularChain<Ring> > dec;
		    RegularChain<Ring> emptyRC;
		    dec.push_back(emptyRC);
		    std::vector<RegularChain<Ring> > newDec;

		    
		    std::vector< std::string> lpMvars;
		    for(int i=0;i<=lp.size();i++){
		      lpMvars.push_back(lp.at(i).leadingVariable());

		    }
		    int i=0;
		    while(i<=n){
		      std::string var = names[i];
		      //A=??
		      for(int l=0; l<dec.size(); ++l){
			RegularChain<Ring> C = dec.at(l);
			if(!this->isAlgebraic(var) && std::find(lpMvars.begin(), lpMvars.end(), var) == lpMvars.end()){

			  std::vector<RegularChain<Ring> > tempResult = this->cleanChain(C, i-1);

			  newDec.insert(newDec.end(),tempResult.begin(),tempResult.end());
			}else if(find(lpMvars.begin(), lpMvars.end(), var) == lpMvars.end()){

			  std::vector<RegularChain<Ring> > tempRC = C.constructRegularChain(this->select(var));
			  for(int k=0; k<tempRC.size(); ++k){
			    RegularChain<Ring> D = tempRC.at(k);
			    std::vector<RegularChain<Ring> > tempResult = this->cleanChain(D, i+1);
			    newDec.insert(newDec.end(), tempResult.begin(),tempResult.end());
			  }
			  
			}else if(!this->isAlgebraic(var)){
			  SparseMultivariatePolynomial<Ring> pV;
			  for(int  k =0; k < lp.size();++k){
			    if(lp.at(k).leadingVariable == var){
			      pV = lp.at(k);
			      break;
			    }
			  }
			  std::vector<RegularChain<Ring> > tempResult = C.intersectFree(pV);
			  for(int k=0; k<tempResult.size();++k){
			   std::vector<RegularChain<Ring> > tempRC =  this->cleanChain(tempResult.at(k), i+1);
			   newDec.insert(newDec.end(),tempRC.begin(), tempRC.end());
			  }

			}else{
			  SparseMultivariatePolynomial<Ring> pV;
			  for(int  k =0; k < lp.size();++k){
			    if(lp.at(k).leadingVariable == var){
			      pV = lp.at(k);
			      break;
			    }
			  }
			  //THE FOLLOWING LINE IS NOT COMPLETE...
			  std::vector<RegularChain<Ring> > tempResult = C.intersectAlgebraic(pV, i+1);
			  for(int k=0; k<tempResult.size();++k){
			   std::vector<RegularChain<Ring> > tempRC =  this->cleanChain(tempResult.at(k), i+1);
			   newDec.insert(newDec.end(),tempRC.begin(), tempRC.end());
			  }

			}
			
		      }
		      

		      dec = newDec;
		      newDec.erase(newDec.begin(),newDec.end());
		      ++i;


		    }

		    result.insert(result.end(), dec.begin(),dec.end());
		    
		  }

		  return result;
		}


		inline std::vector<RegularChain<Ring> >	intersectAlgebraic(SparseMultivariatePolynomial<Ring> p,int indice ,std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src, RegularChain<Ring> C, int A){
		  

		  std::vector<RegularChain<Ring> > rc_out;
		  std::vector<SparseMultivariatePolynomial<Ring> > polys;
		  std::vector<RegularChain<Ring> > rcDrop;

		  C.regular_gcd_by_src_posDim(p, this->select(p.leadingVariable()),src,rc_out, polys, rcDrop);
		  std::vector<RegularChain<Ring> > finalResult;
		  for(int i=0; i<rcDrop.size();++i){
		    RegularChain<Ring> E = rcDrop.at(i);
		    std::vector<RegularChain<Ring> >result = this->cleanChain(E, indice);
		    for(int j=0; j<result.size();++j){
		      std::vector<RegularChain<Ring> >rcOut = this->intersectAlgebraic(p, indice, src, result.at(j),A-1);
		      finalResult.insert(finalResult.end(), rcOut.begin(), rcOut.end()); 
		      
		    }
		    
		  }

		  for(int i=0; i<polys.size();++i ){
		    //if( rc_out.at(i).dimension()==this->under(p.leadingVariable()).dimension()){
		    //replace the following line with constructRegularChain function...

		    RegularChain<Ring> E = rc_out.at(i);
		    E += polys.at(i);
		    finalResult.push_back(E);
		    
		    
		    
		    if(rc_out.at(i).height()<A-1){
		      if(rc_out.at(i).canComputeDim0(polys.initial()))
			continue;
		    }
		    std::vector<RegularChain<Ring> > result = rc_out.at(i).intersect(polys.initial());
		    for(int j=0;j<result.size();++j){
		    
		      std::vector<RegularChain<Ring> >result2 =this->cleanChain(result.at(j));
		      for(int l=0; l<result2.size(); ++l){
			std::vector<RegularChain<Ring> > rcOut = this->intersectAlgebraic(p, indice, src,result2.at(l),A );
			finalResult.insert(finalResult.end(), rcOut.begin(), rcOut.end());
		      }

		    } 
		  }
		  
		}

		std::vector<RegularChain<Ring> > intersectFree(SparseMultivariatePolynomial<Ring> poly, int A){

		  std::vector<RegularChain<Ring> > result;
		  SparseMultivariatePolynomial<Ring> ip = poly.initial();
		  std::vector<RegularChain<Ring> > newRcOut;
		  std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
		  this->regularize(ip, newRcOut, newPolys);
		  for(int i=0; i< newPolys.size(); ++i ){
		    if(newPolys.at(i).isZero()){
		      std::vector<RegularChain<Ring> > rcOut = newRcOut.at(i).intersect(poly.tail());
		      result.insert(result.end(), rcOut.begin(), rcOut.end());
		    }else if(newRcOut.at(i).height()<A){

		      std::vector<RegularChain<Ring> > Temp = newRcOut.at(i).constructRegularChain(poly );
		      
		      result.insert(result.end(), Temp.begin(), Temp.end());
		      if(newRcOut.at(i).canComputeDim0(ip))
			continue;

		      std::vector<RegularChain<Ring> > rcOut = newRcOut.at(i).intersect(ip);

		      for(int j=0;j<rcOut.size();j++){

			std::vector<RegularChain<Ring> > rcOutTemp = rcOut.at(i).intersect(poly.tail());
		      
			result.insert(result.end(), rcOutTemp.begin(), rcOutTemp.end());
		      }
		     
		    }
		  }
		  return result;
		}
		






		std::vector<RegularChain<Ring> > constructRegularChain(SparseMultivariatePolynomial<Ring> poly){
		  
		}



		int height(){

		  return algVariables.size();
		}

		inline int regular_gcd_by_src_posDim(SparseMultivariatePolynomial<Ring> p, SparseMultivariatePolynomial<Ring> q, std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src, std::vector<RegularChain<Ring> > &rc_out, std::vector<SparseMultivariatePolynomial<Ring> > &polys ,  std::vector<RegularChain<Ring> > &rcDrop){

		  std::vector<RegularChain<Ring> > tasks_rc;
		  tasks_rc.push_back(*this);
		  std::vector<int > tasks_indice;
		  tasks_indice.push_back(1);
		  
		  int d = src.at(src.size()-2).degree();
		  
		  while(tasks_rc.size()!=0){
		    RegularChain<Ring> rc = tasks_rc[0];
		    tasks_rc.erase(tasks_rc.begin());
		    
		    int indice = tasks_indice[0];
		    tasks_indice.erase(tasks_indice.begin());
		    SparseMultivariatePolynomial<Ring> f;
		    if(indice>d || indice == q.leadingVariableDegree()){
		      polys.push_back(q);
		      rc_out.push_back(rc);
		      continue;
		    }else if(indice<d){
		      f = principleSubresultantCoefficientOfIndex(indice,src);
		      
		    }else{
		      f = p.initial();
		    }
		    
		    std::vector<RegularChain<Ring> > newRcOut;
		    std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
		    rc.regularizeDim0(f, newRcOut, newPolys);
		    
		    for(int i=0; i<newRcOut.size(); i++){

		      if(newRcOut[i].dimension()< rc.dimensio()){
			rcDrop.push_back(newRcOut[i]);
		      }else if(newPolys[i].isZero()){

			tasks_rc.push_back(newRcOut[i]);
			tasks_indice.push_back(indice+1);

		      }else if(indice==d){
			rc_out.push_back(newRcOut[i]);
			polys.push_back(p);
			
		      }else{
			rc_out.push_back(newRcOut[i]);
			SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> >  tempUiniPoly = this->subresultantOfIndex(indice,src);
			SparseMultivariatePolynomial<Ring> polyTemp(tempUiniPoly);
			polys.push_back(polyTemp); 
		      }
		      
		      
		    }
		   
		    return(1);
		    
		    }
		  
		  return(1);
		  





		}












		inline  std::vector< SparseMultivariatePolynomial<Ring> >  mergeGcdFreeFactorization( std::vector< SparseMultivariatePolynomial<Ring> > lp){
		  
		  std::vector< SparseMultivariatePolynomial<Ring> > result;
		  for(int i =0; i<lp.size(); i++ ){
		    std::vector< SparseMultivariatePolynomial<Ring> > facs = lp.at(i).factors();
		    result.insert(result.end(),facs.begin(),facs.end());
		    
		  }
		  return result;
		}

		std::vector<RegularChain<Ring> > cleanChain(RegularChain<Ring> C,int i){
		  std::vector<RegularChain<Ring> > result;
		  if(i<1 || C.isEmpty()){
		    result.push_back(C);
		    return result;
		  }
		  
		  std::string var = names[i];
		  if(!this->isAlgebraic(var)){
		    result.push_back(C);
		    return result;
		  }
		  if(C.dimension() == this->under(var).dimension()){
		    result.push_back(C);
		    return result;
		  }
		  
		  SparseMultivariatePolynomial<Ring> ip = this->select(var).initial();
		  std::vector<RegularChain<Ring> > newRcOut;
		  std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
		  C.regularize(ip, newRcOut, newPolys);
		  for(int i=0; i<newPolys.size();++i){
		    if(!newPolys.at(i).isZero())
		      result.push_back(newRcOut.at(i));

		  }
		  return result;

		}
		

		inline int regularizePosDim(SparseMultivariatePolynomial<Ring> p, std::vector<RegularChain<Ring> > &newRcOut, std::vector<SparseMultivariatePolynomial<Ring> > &newPolys ){
		  cout<< "entered regularizePosDim"<< endl;
		  std::vector<SparseMultivariatePolynomial<Ring> > regL;
		  std::vector< SparseMultivariatePolynomial<Ring> > lp;
		  lp.push_back(p);
		  std::vector<SparseMultivariatePolynomial<Ring> > iregL = mergeGcdFreeFactorization(lp);
		  std::vector<RegularChain<Ring> > regRCOut;
		  std::vector<RegularChain<Ring> > iregRCOut;
		  
		  this->regularizeList(regL,iregL, regRCOut, iregRCOut);
		  
		  newRcOut.insert(newRcOut.end(),iregRCOut.begin(),iregRCOut.end());
		  newRcOut.insert(newRcOut.end(),regRCOut.begin(),regRCOut.end());
		  SparseMultivariatePolynomial<Ring> zero;
		  for(int i=0; i<iregRCOut.size(); i++)
		    newPolys.insert(newPolys.end(),zero);
		  for(int i=0; i<regRCOut.size(); i++)
		    newPolys.insert(newPolys.end(),p);


		}
		inline int regularizeList(std::vector<SparseMultivariatePolynomial<Ring> > regL,std::vector<SparseMultivariatePolynomial<Ring>  > iregL, std::vector<RegularChain<Ring> > &regRCOut,std::vector<RegularChain<Ring> > &iregRCOut ){
		  
		  if(iregL.size()==0){
		    iregRCOut.push_back(*this);
		    return(1);
		  }
		  SparseMultivariatePolynomial<Ring> p = iregL.at(0);
		  iregL.erase(iregL.begin());
		  
		  std::vector<RegularChain<Ring> > newRcOut;
		  std::vector<SparseMultivariatePolynomial<Ring> > newPolys;

		  this->regularizeSingle(p,newRcOut,newPolys);
		  

		  

		  for(int j=0;j<newRcOut.size(); j++){
		    std::vector<SparseMultivariatePolynomial<Ring> > regC;
		    std::vector<SparseMultivariatePolynomial<Ring> > iregC;
		    
		    
		    RegularChain<Ring> C = newRcOut.at(j);
		    SparseMultivariatePolynomial<Ring> pC = newPolys.at(j);
		    
		    if(pC.isZero()){
		      iregRCOut.push_back(C);
		      continue;
		    }else if(optLazard && C.dimension() < this->dimension() ){
		      regC.push_back(p);
		      iregC.insert(iregC.end(),regL.begin(),regL.end());
		      iregC.insert(iregC.end(),iregL.begin(),iregL.end());
		    }else{
		      regC.push_back(p);
		      iregC =  iregL;
		    }
		    std::vector<RegularChain<Ring> > newRcReg;
		    std::vector<RegularChain<Ring> > newRciReg;
		    

		    
		    regularizeList(regC, iregC, newRcReg,  newRciReg);
		    
		    regRCOut.insert(regRCOut.end(), newRcReg.begin(), newRcReg.end());
		    iregRCOut.insert(iregRCOut.end(), newRciReg.begin(), newRciReg.end());
		    
		  }


		  return (1);
		}
		
		int regularizeSingle(SparseMultivariatePolynomial<Ring> p, std::vector<RegularChain<Ring> > &newRcOut, std::vector<SparseMultivariatePolynomial<Ring> > &newPolys){
		  
		  if(p.isZero()){
		    SparseMultivariatePolynomial<Ring> zeroPoly;
		    zeroPoly.zero();
		    newPolys.push_back(zeroPoly);
		    newPolys.push_back(zeroPoly);		    
		    newRcOut.push_back(*this);
		    return(1);
		    
		  }else if(p.isConstant()|| this->isEmpty()){
		    newRcOut.push_back(*this);
		    newPolys.push_back(p);
		    return(1);
		  }else{
		    newRcOut.push_back(*this);
		    newPolys.push_back(p);
		    return(1);
		  }
		  
		  
		  
		  std::string mvar = p.leadingVariable();
		if(!this->isAlgebraic(mvar)){
		  
		  std::vector<RegularChain<Ring> > newRcOut2;
		  std::vector<SparseMultivariatePolynomial<Ring> > newPolys2;
		  
		  this->regularize(p.initial(), newRcOut2, newPolys2);
		  for(int i=0;i<newRcOut2.size();i++){
		    
		    RegularChain<Ring> C = newRcOut2.at(i);
		    SparseMultivariatePolynomial<Ring> pC = newPolys2.at(i);

		    if(pC.isZero()){

		      std::vector<RegularChain<Ring> > newRcOut3;
		      std::vector<SparseMultivariatePolynomial<Ring> > newPolys3;
		      C.regularize(p.tail(), newRcOut3, newPolys3);
		      
		      newRcOut.insert(newRcOut.end(),newRcOut3.begin(),newRcOut3.end());
		      newPolys.insert(newPolys.end(), newPolys.begin(), newPolys.end());
		      
		    }else{
		      newRcOut.push_back(C);
		      newPolys.push_back(p);
		    }

		    
		  }

		  return 1;
		    
		}
		
		SparseMultivariatePolynomial<Ring> Cv = this->select(mvar);




		}
		
		inline int regularize(SparseMultivariatePolynomial<Ring> p, std::vector<RegularChain<Ring> > &newRcOut, std::vector<SparseMultivariatePolynomial<Ring> > &newPolys ){
		  cout<< "entered regularize"<<endl;

		  if(p.isZero()){
		    SparseMultivariatePolynomial<Ring> zeroPoly;
		    zeroPoly.zero();
		    newPolys.push_back(zeroPoly);
		    newPolys.push_back(zeroPoly);		    
		    newRcOut.push_back(*this);
		    cout<< "exit regularize"<<endl;
		    return(1);
		    
		  }else if(p.isConstant()|| this->isEmpty()){
		    newRcOut.push_back(*this);
		    newPolys.push_back(p);
		    cout<< "exit regularize"<<endl;
		    return(1);
		  }else if(this->isPrimeRegularChain()){
		    if((this->pseudoDivide(p)).isZero()){
		      SparseMultivariatePolynomial<Ring> zeroPoly;
		      zeroPoly.zero();
		      newPolys.push_back(zeroPoly);
		      newPolys.push_back(zeroPoly);
		      newRcOut.push_back(*this);
		      cout<< "exit regularize"<<endl;
		      return(1);  
		      
		    }else{
		      newRcOut.push_back(*this);
		      newPolys.push_back(p);
		      cout<< "exit regularize"<<endl;
		      return(1);
		    }

		  }
		  
		  p = this->removeZero(p);

		  cout<< "p  "<<p<<endl;
		  if(p.isConstant()){
		    newRcOut.push_back(*this);
		    newPolys.push_back(p);
		    cout<< "exit regularize"<<endl;
		    return(1);
		  }
		  
		  if(this->isRegular(p)){
		    newRcOut.push_back(*this);
		    newPolys.push_back(p);
		    cout<< "exit regularize"<<endl;
		    return(1);
		  }


		  if(this->canComputeDim0(p)){
		    cout<< "dimension 0"<<endl;
		    this->regularizeDim0(p,newRcOut,newPolys );
		    cout<< "exit regularize"<<endl;
		    return(1);
		  }else{ cout<< "dimension +"<<endl;
		    this->regularizePosDim(p,newRcOut,newPolys );
		    cout<< "exit regularize"<<endl;
		    return(1);
		  }


		  
		}


		SparseMultivariatePolynomial<Ring> removeZero(SparseMultivariatePolynomial<Ring> p){
		  
		  if(p.isConstant() || this->isEmpty()) 
		    return p;
		  
		  if(this->isZeroDimensionalMathematically() && this->isStronglyNormalized()){
		    return this->normalForm(p);
		  }
		  if(this->pseudoDivide(p).isZero()){
		    return SparseMultivariatePolynomial<Ring>();
		  }
		  SparseMultivariatePolynomial<Ring> q = p;
		  
		  do{
		    if(q.isConstant())
		      return q;
		    SparseMultivariatePolynomial<Ring> init = q.initial();
		    cout<< "init"<<init<<endl;
		    SparseMultivariatePolynomial<Ring> t = q.tail();

		    std::string mvar = q.leadingVariable();
		    SparseMultivariatePolynomial<Ring> r(1);
		    //
		    //r.setVariableNames(vars);
		    int d[1] = {q.leadingVariableDegree()};
		    r.setCoefficient(1,d,Ring(mpq_class(1)));
		    std::vector<std::string> vars;
		    vars.push_back(mvar);
		    r.setVariableNames(vars);
		    SparseMultivariatePolynomial<Ring> h = this->removeZero(init);

		    //cout<<"*************"<<endl<<endl;
		    //cout<< "q:   "<< q<< endl;
		    //cout<< "init q:   "<<init<<endl;
		    //cout<< "r :  "<< r <<endl;
		    //cout<< "h"  << h<<endl;
		    //cout<<"*************"<<endl<<endl;


		    if(h.isZero())
		      q = t;
		    else{
		      h -= init;
		      
		      h*=r;
		      
		      q += h;
		      
		      return q;
		    }
		  }while(true);

		  
		}
		

		bool isRegular(SparseMultivariatePolynomial<Ring> p){

		  if(p.isZero()){
		    return false;
		  }
		  if(this->isEmpty() || p.isConstant()){
		    return true;
		  }

		  std::vector< SparseMultivariatePolynomial<Ring> > lp;
		  lp.push_back(p);
		  std::vector< SparseMultivariatePolynomial<Ring> > pList = mergeGcdFreeFactorization(lp);

		  //not completely follows the algorithm of is_regular in Maple
		  // optimizations are not applied here...

		  SparseMultivariatePolynomial<Ring> pi;
		  for(int i=0; i<pList.size();i++){
		    pi = pList.at(i);
		    if(this->dimension() == 0){
		      pi = this->normalForm(pi);
		      if(!pi.isConstant()){
			SparseMultivariatePolynomial<Ring> contentPart; 
			SparseMultivariatePolynomial<Ring> primPart;
			pi.primitiveFactorization(contentPart, primPart);
			pi = primPart;
		      }
		      if(pi.isZero())
			return false;
		      
		      if(!pi.isConstant()){
			pi = this->pseudoDivide(pi);
			if(pi.isZero())
			  return false;
		      }
		    }
		   
		  }
		  
		  return true;
		  
		}
		








		
		inline std::vector< SparseMultivariatePolynomial<Ring> > listOfInitials(){
		  
		  std::vector< SparseMultivariatePolynomial<Ring> > results;
		  for(int i=0; i<n; i++){
		    if(!set[i].isZero())
		      results.push_back(set[i].initial());
		  }
		  return results;
		}
		inline int minRittSet(std::vector< SparseMultivariatePolynomial<Ring> > lp){
		  
		  SparseMultivariatePolynomial<Ring> pmin = lp.at(0);
		  int k = 0;
		  for(int i =0; i<lp.size();i++){
		    if(strictlyLessRitt(lp.at(i), pmin)){
		      pmin = lp.at(i);
		      k=i;
		    }
		    
		  }
		  return k;
		}
		
		
		

		inline bool strictlyLessRitt(SparseMultivariatePolynomial<Ring> p1, SparseMultivariatePolynomial<Ring> p2){
		  
		  if(p1.isConstant())
		    return true;
		  if(p2.isConstant())
		    return false;
		  
		  
		  std::string v1 = p1.leadingVariable();
		  std::string v2 = p1.leadingVariable();
		  
		  std::string * xx;
		  xx= std::find(names, names+n,v1);
		  int indice1 = std::distance(names,xx);

		  xx= std::find(names, names+n,v2);
		  int indice2 = std::distance(names,xx);

		  if(indice1 < indice2 )
		    return true;
		  if(indice1 > indice2 )
		    return false;
		  if(p1.leadingVariableDegree() < p2.leadingVariableDegree() )
		    return true;
		  if(p1.leadingVariableDegree() > p2.leadingVariableDegree() )
		    return false;
		  if(strictlyLessRitt(p1.initial(), p2.initial()))
		    return true;
		  if(strictlyLessRitt(p2.initial(), p1.initial()))
		    return false;
		  
		  return strictlyLessRitt(p1.tail(), p2.tail());
		  
		}
		
		inline bool cleanSet(std::vector< SparseMultivariatePolynomial<Ring> > lp, std::vector<SparseMultivariatePolynomial<Ring> > &cleanSet){

		  for(int i=0; i<lp.size();i++){
		    SparseMultivariatePolynomial<Ring> cleanP = this->normalForm(lp.at(i));
		    if(!cleanP.isZero() && cleanP.isConstant()){
		      return true;
		    }
		    cleanSet.push_back(cleanP);
		    
		  }


		  return false;

		}
		bool canComputeDim0(SparseMultivariatePolynomial<Ring> poly){

		  std::vector<std::string> s1= poly.variables();
		  std::vector<std::string> s2;
		  s2.push_back(poly.leadingVariable());
		  std::vector<std::string> result;
		  std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
				      std::inserter(result, result.end()));

		  std::vector<std::string> vars(this->names, this->names+this->n);
		  //
		    bool dim0 =  (poly.isConstant() || IsSubset( result, vars) );
		    bool outPut = (this->isZeroDimensionalMathematically() && dim0);
		    //bool outPut = (this->dimension()==0 && dim0);
		    return outPut;
		   
//if()  
		  
		}

		bool isZeroDimensionalMathematically(){
		  int iter=0;
		  std::vector<std::string> vars;
		  for(int j=0; j<n;++j){
		    if(!this->set[j].isZero()){
		      vars = setUnion(vars, this->set[j].variables());
		      ++iter;
		    }
		  }
		  
		  
		  if(vars.size()==iter)
		    return true;
		  else 
		    return false;
		  


		}

		std::vector<RegularChain<Ring> > extend(SparseMultivariatePolynomial<Ring> p, int A){
		  SparseMultivariatePolynomial<Ring> init = p.initial();
		  std::vector<RegularChain<Ring> > result;
		  //What is the difference between regular and regularize?
		  std::vector<RegularChain<Ring> > newRcOut;
		  std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
		  regularize(init ,newRcOut,  newPolys, A-1);
		  for(int i=0;i<newPolys.size();++i){
		    if(newPolys.at(i).isConstant())
		      continue;
		    std::vector<RegularChain<Ring> > outRC =  newRcOut.at(i).constructRegularChain(p);
		    result.insert( result.end(),outRC.begin(), outRC.end());
		    
		  }

		  return result;

		}

		std::vector<RegularChain<Ring> > listExtend(std::vector<SparseMultivariatePolynomial<Ring> > lp, int A){
		  
		  std::vector<RegularChain<Ring> > result;
		  
		  std::vector<std::pair <std::vector<SparseMultivariatePolynomial<Ring> >  , RegularChain<Ring>  > > pile;
		  std::pair<std::vector<SparseMultivariatePolynomial<Ring> >  , RegularChain<Ring>  > p1(lp, *this);
		  pile.push_back(p1);

		  while(pile.size()!=0){
		    p1 = pile.at(pile.begin());
		    pile.erase(pile.begin());
		    std::vector<SparseMultivariatePolynomial<Ring> > lpNew = p1.first;
		    RegularChain<Ring> rc = p1.second;
		    if(lp.size()==0){
		      result.push_back(rc);
		      continue;
		    }
		    SparseMultivariatePolynomial<Ring> p =  lpNew.at( lpNew.begin());
		    lpNew.erase(lpNew.begin());
		    //A=?
		    
		    std::vector<RegularChain<Ring> > tempResult = rc.extend(p,A);
		    
		    for(int i=0; i< tempResult.size(); ++i){
		      p1 = (lpNew, tempResult.at(i));
		      pile.push_back(p1);
		    }
		    return result;
		  }
		  
		  
		}
		// helper functions:

		template <typename T>
		  bool IsSubset(std::vector<T> A, std::vector<T> B)
		  {
		    std::sort(A.begin(), A.end());
		    std::sort(B.begin(), B.end());
		    return std::includes(A.begin(), A.end(), B.begin(), B.end());
		  }
		

		std::vector<std::string> instersection(std::vector<std::string> &v1, std::vector<std::string> &v2)
		  {

		    std::vector<std::string> v3;
		    
		    std::sort(v1.begin(), v1.end());
		    std::sort(v2.begin(), v2.end());
		    
		    std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
		    
		    return v3;
		  }

		std::vector<std::string> setUnion(std::vector<std::string> v1, std::vector<std::string> v2){

		    std::vector<std::string> v3;
		    
		    std::sort(v1.begin(), v1.end());
		    std::sort(v2.begin(), v2.end());
		    
		    std::set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
		    
		    return v3;
		  }
		

};

#endif
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


