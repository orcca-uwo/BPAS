#ifndef _ZERODIMENSIONREGULARCHAIN_H_
#define _ZERODIMENSIONREGULARCHAIN_H_


#include "regularchain.h"


template <class Ring>

class ZeroDimensionRegularChain :  public RegularChain<Ring> {
 

  using RegularChain<Ring>::n;
  using RegularChain<Ring>::names;
  using RegularChain<Ring>::set;

  using RegularChain<Ring>::algVariables;
  using RegularChain<Ring>::stronglyNormalized;
  using RegularChain<Ring>::isEmptyRC;
  using RegularChain<Ring>::isKnownToBePrime;
  using RegularChain<Ring>::isKnownToBeSquareFree;
  using RegularChain<Ring>::optLazard;
  
 public:
  /**
   * Default constructor
   *
   * @param
   **/
  ZeroDimensionRegularChain<Ring> () { n = 0; }
  /**
   * Construct with number of polynomials and variable names
   *
   * @param s: The number of polynomials
   * @param xs; Variable names
   **/
  ZeroDimensionRegularChain<Ring> (int s, std::vector<std::string> xs) {
    if (s < 1) { n = 0; }
    else {
      n = s;
      set = new SparseMultivariatePolynomial<Ring>[n];
      
      if (xs.size() == s) {
	names = new std::string[n];
	for (int i = 0; i < n; ++i){
	  names[i] = xs[n-i-1];
	  
	}
	
	
      }
      else {
	std::cout << "BPAS: error, ZeroDimensionRegularChain<Ring> must be defined with " << s << " variables." << std::endl;
	exit(1);
      }
    }
  }
  
 /**
   * Copy constructor
   *
   * @param a: A regular chain
   **/
  ZeroDimensionRegularChain<Ring> (const ZeroDimensionRegularChain<Ring>& a) {
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
   * Overload operator =
   *
   * ZeroDimensionRegularChain<Ring>'s
   *
   * @param
   **/
 
  inline ZeroDimensionRegularChain<Ring>& operator= (ZeroDimensionRegularChain<Ring> a) {
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
  
  

  inline ZeroDimensionRegularChain<Ring>& operator= (std::string &str) {
    
    std::vector<std::string > vars;
    if(n>0){
      
      
      for(int i=n-1 ; i>=0 ; i--){
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
  inline ZeroDimensionRegularChain<Ring> operator+ (SparseMultivariatePolynomial<Ring>& smp) {
    ZeroDimensionRegularChain<Ring> r (*this);
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
      this->updateRegularChainsStates(smp,k);
      //The following line was moved into updateRegularChainsStates function 
      //set[k] += smp;
      
      // the string of main variables will become updated here
      algVariables.push_back(names[k]);
      
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
  


  inline ZeroDimensionRegularChain<Ring> under(std::string x) {
    int k = -1;
    for (int i = 0; i < n; ++i) {
      if (x == names[i]) {
	k = i;
	break;
      }
    }
    if (k <= 0)
      return ZeroDimensionRegularChain<Ring>();
    else {
      
      
      ZeroDimensionRegularChain<Ring> r;
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
  
  /**
   * The triangular set upper a variable 
   *
   * @param x: The variable name
   **/
  inline ZeroDimensionRegularChain<Ring> upper(std::string x) {
    
    
    int k = -1;
    for (int i = 0; i < n; ++i) {
      if (x == names[i]) {
	k = i; 
	break;
      }
    }
    if (k < 0 || k == n-1)
      return ZeroDimensionRegularChain<Ring>();
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
		    
      
      ZeroDimensionRegularChain<Ring> r(namesVec.size(), namesVec);
      
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
  //using RegularChain<Ring>::pseudoDivide;
  /**
   * Monic division
   * Return the remainder
   *
   * @param p: An input polynomial
   * @param quo: The quotients
   **/
  //using RegularChain<Ring>::monicDivide;
  
  //this function should be relocated as it does not necessarily belongs to regularchain class
  using RegularChain<Ring>::principleSubresultantCoefficientOfIndex;

//this function should be relocated as it does not necessarily belongs to regularchain class
  using RegularChain<Ring>::subresultantOfIndex;





  inline int regularizeDim0(SparseMultivariatePolynomial<Ring> poly, std::vector<ZeroDimensionRegularChain<Ring> > &newRcOut, std::vector<SparseMultivariatePolynomial<Ring> > &newPolys ){
   
    
    SparseMultivariatePolynomial<Ring> p = this->removeZero(poly);
    
    std::stringstream ss;
    ss << p;
    std::string s = ss.str();
    p = s;
    
    
    if(p.isZero()){
      SparseMultivariatePolynomial<Ring> zeroPoly;
      zeroPoly.zero();
      newPolys.push_back(zeroPoly);
      
      newRcOut.push_back(*this);
     
      return(1);
      
    }else if(p.isConstant()|| this->isEmpty()){
      
      newRcOut.push_back(*this);
      newPolys.push_back(poly);
      
      return(1);
    }else if(this->isPrimeRegularChain()){
      if((this->pseudoDivide(p)).isZero()){
	SparseMultivariatePolynomial<Ring> zeroPoly;
	zeroPoly.zero();
	newPolys.push_back(zeroPoly);
	newPolys.push_back(zeroPoly);
	newRcOut.push_back(*this);
	
	return(1);  
	
      }else{
	newRcOut.push_back(*this);
	newPolys.push_back(p);
	
	return(1);
      }
      
    }
    
    if((this->pseudoDivide(p)).isZero()){
      SparseMultivariatePolynomial<Ring> zeroPoly;
      zeroPoly.zero();
      newPolys.push_back(zeroPoly);
      newPolys.push_back(zeroPoly);
      newRcOut.push_back(*this);
      
      return(1);  
      
    }
    std::string mvar = p.leadingVariable();
    
    SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > >  up= p.convertToSUP(mvar);
    
    ZeroDimensionRegularChain<Ring> under_rc = this->under(mvar);
    
    std::string * xx;
    xx= std::find(names, names+n,mvar);
    int indice = std::distance(names,xx);
    SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > >  uq= set[indice].convertToSUP(mvar);
    std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src = up.monomialBasisSubresultantChain(uq);
    
    SparseMultivariatePolynomial<Ring >  res = principleSubresultantCoefficientOfIndex(0,src);
    
    std::vector<ZeroDimensionRegularChain<Ring> > newRcOutTemp;
    std::vector<SparseMultivariatePolynomial<Ring> > newPolysTemp;
    under_rc.Display();
    
    under_rc.regularizeDim0(res, newRcOutTemp, newPolysTemp);
    
    for(int i=0; i<newPolysTemp.size();i++){
      
      ZeroDimensionRegularChain<Ring> rc = newRcOutTemp[i];
      SparseMultivariatePolynomial<Ring> pNew= newPolysTemp[i];
      
      if(!pNew.isZero()){	    
	newRcOut.push_back(*this);
	newPolys.push_back(p);
	
      }
      else{
	
	std::vector<ZeroDimensionRegularChain<Ring> > rcOutTemp;
	std::vector<SparseMultivariatePolynomial<Ring> > polysTemp;
	rc.regularGcd_by_src(p, set[indice], src , rcOutTemp, polysTemp);
	
	for(int j=0;j<polysTemp.size();j++){
	  SparseMultivariatePolynomial<Ring> g = polysTemp[j];
	  
	  ZeroDimensionRegularChain<Ring> E = rcOutTemp[j];
	  
	  if(g.leadingVariableDegree()==set[indice].leadingVariableDegree()){
	    
	    /*   Build a new regular chain   by E C_v and rc_upper*/
	    std::vector<std::string> xs(names, names+n);
	    std::reverse(xs.begin(),xs.end());
	    ZeroDimensionRegularChain<Ring> rcNew(n,xs);
	    std::vector<SparseMultivariatePolynomial<Ring> > polySet = E.polynomials();
	    for(int ll=0; ll<polySet.size();++ll ){
	      rcNew+=polySet[ll];
	    }	
	    rcNew+=g;
	    ZeroDimensionRegularChain<Ring> rc_upper = this->upper(mvar);
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
	  
	  SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > ug = g.convertToSUP(mvar);
	  SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > q = uq.pseudoDivide(ug,NULL);
	  
	  
	  /*   Build a new regular chain   by E g and rc_upper*/
	  std::vector<std::string> xs(names, names+n);
	  std::reverse(xs.begin(),xs.end());
	  
	  
	  ZeroDimensionRegularChain<Ring> rcNew(n,xs);
	  
	  std::vector<SparseMultivariatePolynomial<Ring> > polySet = E.polynomials();
	  
	  
	  ZeroDimensionRegularChain<Ring> rc_upper = this->upper(mvar);
	  
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
	  
	  
	  if(q.degree()>0){
	    std::vector<ZeroDimensionRegularChain<Ring> > newRcOutTemp2;
	    std::vector<SparseMultivariatePolynomial<Ring> > newPolysTemp2;
	    
	    /*   Build a new regular chain   by E g and rc_upper*/
	    std::vector<std::string> xs(names, names+n);
	    std::reverse(xs.begin(),xs.end());
	    ZeroDimensionRegularChain<Ring> rcNew(n,xs);
	    std::vector<SparseMultivariatePolynomial<Ring> > polySet = E.polynomials();
	    for(int ll=0; ll<indice;ll++ ){
	      rcNew+=polySet[ll];
	    }
	    
	    SparseMultivariatePolynomial<Ring > tempPoly(q);
	    rcNew+= tempPoly;
	    for(int ll=0; ll<rc_upper.n;ll++){
	      rcNew+=rc_upper.set[ll];
	    }
	    
	    
	    
	    if(rcNew.isSquareFree()){
	      newRcOut.push_back(rcNew);
	      newPolys.push_back(p);
	      
	    }else{
	      
	      rcNew.regularizeDim0(p, newRcOutTemp2, newPolysTemp2 );
	      
	      for (int k=0; k<newPolysTemp2.size(); k++ ){
		newRcOut.push_back(newRcOutTemp2[k]);
		newPolys.push_back(newPolysTemp2[k]);
	      }
	    }
	   
	    return(1);  
	  }
	  
	  
	}
	
      }
    }
    
    return(0);
    
  }	
  
  
  inline int regularGcd(SparseMultivariatePolynomial<Ring> p, SparseMultivariatePolynomial<Ring> q ,std::vector<ZeroDimensionRegularChain<Ring> > &rc_out, std::vector<SparseMultivariatePolynomial<Ring> > &polys ){
    
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
    
    
    
    
    if(this->isEmpty()){
      if(!r.isZero()){
	
	rc_out.push_back(*this);
	polys.push_back(r);
	return 1;
      }
      
    }
    std::vector<ZeroDimensionRegularChain<Ring> > newRcOut;
    std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
    this->regularizeDim0(r, newRcOut, newPolys);
    for(int i=0; i<newPolys.size();++i){
      
      SparseMultivariatePolynomial<Ring> f = newPolys.at(i);
      ZeroDimensionRegularChain<Ring> C =  newRcOut.at(i);
      
      //first ifstatment is commneted out for the moment
      if(f.isZero()){
	C.regularGcd_by_src(p,q,src, rc_out, polys);
      }else{
	rc_out.push_back(C);
	polys.push_back(f);
      }
      
      
    }
    
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
  inline int regularGcd_by_src(SparseMultivariatePolynomial<Ring> p, SparseMultivariatePolynomial<Ring> q, std::vector<SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring > > > src, std::vector<ZeroDimensionRegularChain<Ring> > &rc_out, std::vector<SparseMultivariatePolynomial<Ring> > &polys ){
    

    
    std::vector<ZeroDimensionRegularChain<Ring> > tasks_rc;
    tasks_rc.push_back(*this);
    std::vector<int > tasks_indice;
    tasks_indice.push_back(1);
    
    int d = src.at(src.size()-2).degree();
    
    while(tasks_rc.size()!=0){
      ZeroDimensionRegularChain<Ring> rc = tasks_rc[0];
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
     
      std::vector<ZeroDimensionRegularChain<Ring> > newRcOut;
      std::vector<SparseMultivariatePolynomial<Ring> > newPolys;
      rc.regularizeDim0(f, newRcOut, newPolys);
      
      
      
      for(int i=0; i<newRcOut.size(); i++){
	if(newPolys[i].isZero()){  
	  
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
      
      
    }
    
    return(1);
    
  }
  
  
  
  

  //helper functions....
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


