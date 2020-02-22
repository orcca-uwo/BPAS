#ifndef _TRIANGULARSET_H_
#define _TRIANGULARSET_H_

//#include "Polynomial/mpolynomial.h"

using std::endl;
using std::cout;

// "Classical" way of enforcing constraints on template parameters.
template<class T, class B> struct Derived_from {
    static void constraints(T* p) { B* pb = p; }
    Derived_from() { void(*p)(T*) = constraints; }
};

template <class Field, class RecursivePoly>
class TriangularSet : public Derived_from<Field,BPASField>, Derived_from<RecursivePoly,BPASRecursivelyViewedPolynomial<Field,RecursivePoly>, BPASTriangularSet<Field,RecursivePoly,TriangularSet<Field,RecursivePoly>> {
	protected:
		int n;
		RecursivePoly* set;
		Symbol* names;
		std::vector<Symbol> algVariables;

		bool isEmptyRC = true;
		bool stronglyNormalized = true;

		/*inline std::vector<RecursivePoly > parser(const std::string& str){
			std::vector<RecursivePoly > polys;
			std::string poly;
			std::size_t pos = 0;
			std::size_t found;
			do{
				found = str.find(",",pos);
				if(found!=std::string::npos){
					poly = str.substr(pos, found-pos);
					RecursivePoly poly2;
					poly2  = poly; 
					polys.push_back(poly2);
					pos = found+1;
				}
				else{
					std::string poly = str.substr(pos, found-pos);
					RecursivePoly poly2;
					poly2  = poly;
					//tri += poly2;
					polys.push_back(poly2);
					break;
				}
			}while(true);
			return polys;
		}*/


		inline void updateRegularChainsStates(const RecursivePoly& rp, int k) {
			//Check whether or not by adding rp to rc, rc remains strongly normalized
			RecursivePoly init = rp.initial();
			std::vector<Symbol> xs = init.variables();
			if(stronglyNormalized){
				std::sort(xs.begin(),xs.end());
				std::sort(algVariables.begin(),algVariables.end());
				std::vector<Symbol> It;
				std::set_intersection(xs.begin(),xs.end(),algVariables.begin(),algVariables.end(),back_inserter(It));
				
				if(It.size()!=0 && It.at(0)!="9" && It.at(0)!="1"){
					stronglyNormalized = false;
				}
				for (int i=0;i<n;++i){
					std::vector<Symbol> xxs = (this->set.at(i).initial()).variables();
					std::sort(xxs.begin(),xxs.end());
					std::set_intersection(xs.begin(),xs.end(),algVariables.begin(),algVariables.end(),back_inserter(It));
					if(It.size()!=0 && It.at(0)!="9" && It.at(0)!="1")
						stronglyNormalized = false;
				}
			}
			set[k] += rp;
		}


	public:
		/**
		 * Default constructor
		 *
		 * @param
		 **/
		TriangularSet<Field,RecursivePoly> () : n(0) {}
		/**
		 * Construct the empty triangular set in the s decreasingly ordered 
		 variables given by xs
		 *
		 * @param s: The number of polynomials
		 * @param xs: The variable names
		 **/
		TriangularSet<Field,RecursivePoly> (int s,const std::vector<Symbol>& xs) {
			if (s < 1) { n = 0; }
			else {
				n = s;
				set = new RecursivePoly[n];
				if (xs.size() == s) {
					names = new Symbol[n];
					for (int i = 0; i < n; ++i)
						names[i] = xs.at(n-i-1); // why order reversing?
				}
				else {
					std::cout << "BPAS: error, TriangularSet<Field,RecursivelyViewedPolynomial> must be defined with " << s << " variables." << std::endl;
					exit(1);
				}
			}	
		}
		/**
		 * Copy constructor
		 *
		 * @param a: A triangular set
		 **/
		TriangularSet<Field,RecursivePoly> (const TriangularSet<Field,RecursivePoly>& a) : n (a.n) {
			if (n) {
				set = new RecursivePoly[n];
				names = new std::string[n];
				std::copy(a.set, a.set+n, set);
				std::copy(a.names, a.names+n, names);
			}
		}
		/**
		 * Deconstructor
		 *
		 * @param
		 **/
		~TriangularSet<Field,RecursivePoly>() {
			if (n) {
				delete [] set;
				delete [] names;
			}
		}
		/**
		 * Overload operator =
		 *
		 * @param a: A triangular set
		 **/
		inline TriangularSet<Field,RecursivePoly> operator= (const TriangularSet<Field,RecursivePoly>& a) {
			if (this != &a) {
				if (n) {
					delete [] set;
					delete [] names;
				}
				n = a.n;
				set = new RecursivePoly[n];
				names = new Symbol[n];
				//std::copy(a.set, a.set+n, set);
				for (int i=0; i<n; i++){
					*this += a.set[i];
				}
				std::copy(a.names, a.names+n, names);
				algVariables = a.algVariables;
				isEmptyRC = a.isEmptyRC;
				stronglyNormalized = a.stronglyNormalized;
			}
			return *this;
		}
		/**
		 * Parsing a string and assigning the result to
		 * a triangular set
		 * @param
		 **/
		/*virtual inline TriangularSet<Field,RecursivePoly>& operator= (const std::string& str) {
			std::vector<std::string > vars;
			if (n>0) {
				for(int i=n-1 ; i>=0 ; i--){
					vars.push_back(names[i]);
				}
			}
			else {
				std::cout<<"BPAS:: error: a list of variables has to be specified"<<std::endl;
			}
			TriangularSet<Field,RecursivePoly> tri(n,vars);
			std::vector<RecursivePoly> polys = this->parser(str);
			for(int i =0; i<polys.size(); i++)
				tri += polys[i];
			*this = tri;
			return *this;
		}*/
		/**
		 * Overload operator +
		 * Adds a polynomial to a triangular set and returns a new triangular set
		 * @param rp: A sparse multivariate polynomial
		 **/
		inline TriangularSet<Field,RecursivePoly> operator+ (const RecursivePoly& rp) {
			TriangularSet<Field,RecursivePoly> r(*this);
			return (r += rp);
		}
		/**
		 * Overload operator +=
		 * Adds a polynomial to a triangular set in a destructive manner
		 * @param rp: A sparse multivariate polynomial
		 **/
		inline TriangularSet<Field,RecursivePoly> operator+= (const RecursivePoly& rp) {
			if (rp.isConstant()) {
				std::cout << "BPAS: error, cannot add a constant to TriangularSet<Field,RecursivelyViewedPolynomial>" << std::endl;
				exit(1);
			}
			int v = rp.numberOfVariables();
			if (v>n || !v) {
				std::cout << "BPAS: error, cannot add a polynomial of " << v << " variables to a TriangularSet<Field,RecursivelyViewedPolynomial> of " << n << " variables." << std::endl;
				exit(1);
			}
			int* pos = new int[v];
			std::vector<Symbol> xs = rp.variables();
			for (int i=0; i<v; ++i) {
				pos[i] = -1;
				for (int j=n-1; j>-1; --j) {
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
				std::cout << "BPAS: error, TriangularSet<Field,RecursivelyViewedPolynomial>[";
				for (int i = n-1; i > -1; --i) {
		
					if (i) { std::cout << ", "; }	
				}
				std::cout << "] cannot add a polynomial in Field[";
				for (int i = 0; i < xs.size(); ++i) {
					std::cout << xs[i];
					if (i < xs.size()-1) { std::cout << ", "; }
				}
				std::cout << "]." << std::endl;
				exit(1);
			}
			k = pos[0];
			delete [] pos;
			if (!set[k].isZero()) {
				std::cout << "BPAS: error, a polynomial with the leading variable " << xs[0] << " already exists in TriangularSet<Field,RecursivePoly>[";
				for (int i = n-1; i > -1; --i) {
					if (i)
						std::cout << ", ";
				}
				std::cout << "]." << std::endl;
				exit(1);
			}
			updateRegularChainsStates(rp,k);
			//set[k] += rp;
			return *this;
		}
		/**
		 * Get the number of variables
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return n;
		}
		/**
		 * Get the variable names in decreasing order
		 *
		 * @param
		 **/
		inline std::vector<Symbol> variables() const {
			std::vector<Symbol> xs;
			for (int i=n-1; i>-1; --i)
				xs.push_back(names[i]);
			return xs;
		}
		/**
		 * Select a polynomial given the leading variable;
		 * if no such polynomial, 0 is returned
		 * @param x: The leading variable name
		 **/
		inline RecursivePoly select(const Symbol& x) const {
			int k = -1;
			for (int i=0; i<n; ++i) {
				if (x == names[i]) {
					k = i;
					break;
				}
			}
			if (k < 0)
				return RecursivePoly();
			else
				return set[k];
		}
		/**
		 * Returns the ts consisting of polynomials with 
		 * main variable strictly less than s
		 *
		 * @param s: Symbol of the main variable of specified element of the ts
		 **/
		inline TriangularSet<Field,RecursivePoly> under(const Symbol& s) const {
			int k = -1;
			for (int i=0; i<n; ++i) {
				if (s == names[i]) {
					k = i;
					break;
				}
			}
			if (k <= 0)
				return TriangularSet<Field,RecursivePoly>();
			else {
				TriangularSet<Field,RecursivePoly> r;
				r.n = k;
				r.set = new RecursivePoly[r.n];
				r.names = new Symbol[r.n];
				std::copy(names, names+k, r.names);
				//std::copy(set, set+k, r.set);
				for(int i=0; i<k; i++) {
					if(!this->set[i].isZero())
						r += this->set[i];
				}
				return r;
			}
		}
		/**
		 * Returns the ts consisting of polynomials with
		 * main variable strictly greater than s
		 *
		 * @param s: Symbol of the main variable of specified element of the ts
		 **/
		inline TriangularSet<Field,RecursivePoly> upper(const Symbol& s) const {
			int k = -1;
			for (int i=0; i<n; ++i) {
				if (s == names[i]) {
					k = i; 
					break;
				}
			}
			if (k<0 || k == n-1)
				return TriangularSet<Field,RecursivePoly>();
			else {
				std::vector<Symbol> vars;
					for (int i=k+1; i<n; ++i) {
					std::vector<Symbol> varsTemp = this->set[i].variables();
					vars =  setUnion(vars, varsTemp);
				}
				std::vector<Symbol> rcVars(names, names+n); // why not call this->variables()?
				std::vector<Symbol> namesVec;
				for(int i=n-1; i>=0; --i){
					std::vector<Symbol>::iterator it = std::find(vars.begin(), vars.end(), rcVars[i]);
					std::cout << rcVars[i] << endl;
					if(it!=vars.end())
						namesVec.push_back(names[i]);
				}
				TriangularSet<Field,RecursivePoly> r(namesVec.size(), namesVec);
				for(int i=k+1 ; i<n ; i++){
					if (!this->set[i].isZero())
						r += this->set[i];
				}

				return r;
			}
		}
		/**
		 * Pseudo division
		 * Return the pseudo-remainder, the pseudo-quotients and
		 * c such that c*p = ∑(q_i T_i) + r 
		 * @param p: An input polynomial
		 * @param quo: The array of quotients
		 * @param c: The constant multiplied to the input polynomial
		 **/
		inline RecursivePoly pseudoDivide (const RecursivePoly& p, std::vector<RecursivePoly>* quo=NULL, RecursivePoly* c=NULL) const {
			bool qflag = false;
			bool cflag = false;
			if (quo == NULL) {
				quo = new std::vector<RecursivePoly>;
				qflag = true;
			}
			if (c == NULL) {
				c = new RecursivePoly;
				cflag = true;
			}
			c->one();
			RecursivePoly r(p);
			// my code will assume a non-destructive pseudodivision routine for now.
			for (int i = n-1; i > -1 && !r.isZero(); --i) {
				if (!set[i].isZero()) {
					//SparseUnivariatePolynomial<RecursivePoly> s = set[i].convertToSUP(names[i]);
					//SparseUnivariatePolynomial<RecursivePoly> x = r.convertToSUP(names[i]);
					RecursivePoly d;
					RecursivePoly q;
					r = p.lazyPseudoDivide(set[i],names[i],&q,&d);// might need massaging to ensure that the main variables agree
					//SparseUnivariatePolynomial<RecursivePoly> q = x.lazyPseudoDivide(s, &d);
					//r = RecursivePoly(x);
					*c *= d;
					quo->push_back(q);
				}
			}
			if (qflag)
				delete quo;
			if (cflag)
				delete c;
			return r;
		}
		/**
		 * Monic division
		 * Returns the remainder; requires the triangular set be monic
		 *
		 * @param p: An input polynomial
		 * @param quo: The quotients
		 **/
		inline RecursivePoly monicDivide (const RecursivePoly& p, std::vector<RecursivePoly>* quo=NULL) const {
			// we need case checking here to ensure that the ts is monic 
			bool qflag = false;
			if (quo == NULL) {
				quo = new std::vector<RecursivePoly>;
				qflag = true;
			}
			RecursivePoly r(p);
			// my code will assume a non-destructive monic division routine for now.
			for (int i=n-1; i>-1 && !r.isZero(); --i) {
				if (!set[i].isZero()) {
					//SparseUnivariatePolynomial< RecursivePoly > s = set[i].convertToSUP(names[i]);
					//SparseUnivariatePolynomial< RecursivePoly > x = r.convertToSUP(names[i]);
					//SparseUnivariatePolynomial< RecursivePoly > q = x.monicDivide(s);
					//r = RecursivePoly(x);
					//cout<< "rpoly: "<< r << " lv " << r.leadingVariable()<<endl;
					RecursivePoly q;
					r = p.monicDivide(set[i],names[i],&q); // might need massaging to ensure that the main variables agree
					quo->push_back(q);
				}
			}
			if (qflag)
				delete quo;
			return r;
		}
		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param a: A triangular set
		 **/
		inline friend std::ostream& operator<< (std::ostream& out, TriangularSet<Field,RecursivePoly>& a) {
			bool isFirst = 0;
			out << "[";
			for (int i = 0; i < a.n; ++i) {
				if (!a.set[i].isZero()) {
					if (isFirst) { out << ", "; }
					out << a.set[i];
					isFirst = 1;
				}
			}
			out << "]";
			return out;
		}
		/**
		 * Overload stream operator ==
		 *
		 *
		 * @param a: A triangular set
		 **/
		inline bool operator== (TriangularSet<Field,RecursivePoly>& a) {
			if(n == a.n){
				for(int i=0; i<n; i++){
					if (set[i] != a.set[i])
						return 0;
				}
				return 1;
			}
			return 0;
		}
		/**
		 * Display the triangular set
		 *
		 **/
		inline Display() {
			if(this->n==0) {
				std::cout<<" /" << std::endl;
				std::cout << " | "<<std::endl;
				std::cout<<"<  "<<std::endl;
				std::cout << " | "<<std::endl;
				std::cout<<" \\" << std::endl;
			}
			else if(this->n==1) {
				std::cout<<" /" << std::endl;
				std::cout << " | "<<std::endl;
				std::cout<<"<  "<< this->set[this->n-1] << " = 0"<<std::endl;
				std::cout << " | "<<std::endl;
				std::cout<<" \\" << std::endl;
			}
			else {
				int half = this->n/2;
		    	int rem = this->n%2;
				std::cout << " /" << std::endl;
		    	//std::cout << " |" << std::endl;
		    	for (int i = 0; i < half; ++i) {
					if (!this->set[i].isZero()) {
						//if (isFirst) { std::cout << " | " << std::endl; }
						std::cout << " | "<<this->set[i] << " = 0"<< std::endl;
					}
				}
				if(rem==0) {
					std::cout<<"< "<<std::endl;
					for (int i = half; i < this->n; ++i) {
						if (!this->set[i].isZero()) {
							//if (isFirst) {  }
							std::cout << " | " << this->set[i] <<" = 0"<< std::endl;
						}
					}
				}
				else {
					std::cout<<"<  "<< this->set[half] << " = 0" <<std::endl;
					for (int i = half+1; i < this->n; ++i) {
						if (!this->set[i].isZero()) {
							std::cout << " | " << this->set[i] << " = 0"<<std::endl;
						}
					}
				}
				//std::cout << " | " << std::endl;
				std::cout << " \\ " << std::endl;
			}
		}
		/**
		 * normalForm in the sense of Groebner basis
		 *
		 * @param rp: A sparse multivaiate polynomial over a Ring  
		 * 
		 **/
		inline RecursivePoly normalForm(RecursivePoly rp) {
			if(!stronglyNormalized){
				std::cout<< "BPAS: error, The regular chain should be strongly normalized"<<std::endl;
				exit(0); 
			}
			// Compute the normal form in the sense of Groebner Basis
			for(int i=0; i<n; ++i){
				if(!this->set[i].isZero()) {
					std::vector<std::string> vars = this->set[i].initial().variables();
					if(vars.size()!=0){
						cout<<"salam"<<endl;
						std::vector<RecursivePoly > *quo;
						RecursivePoly c;
						RecursivePoly r = pseudoDivide(rp, quo, &c);
						r/=c;
						return r;
					}
				}
			}
			RecursivePoly r;
			this->makeMonic();
			r= this->monicDivide(rp);
			return r;
		}
		/**
		 * Make monic
		 *
		 * @param
		 **/
		 // Not sure this should be here?
		inline void makeMonic() {
			for(int i=0;i<n;i++) {
				if(!this->set[i].isZero())
					this->set[i] = this->set[i]/this->set[i].leadingCoefficient();
			}
		}
		/**
		 * Get algebraic variables
		 *
		 * @param
		 **/
		inline std::vector<Symbol> mainVariables() const {
			return algVariables;
		}
		
		inline bool isAlgebraic(std::string mvar) const {
			std::vector<std::string>::iterator it = std::find(algVariables.begin(),algVariables.end(), mvar);
			if (it!=algVariables.end())
				return true;
			else
				return false;
		}



		/*
		std::vector<std::string> setUnion(std::vector<std::string> v1, std::vector<std::string> v2){

		    std::vector<std::string> v3;
		    
		    std::sort(v1.begin(), v1.end());
		    std::sort(v2.begin(), v2.end());
		    
		    std::set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
		    
		    return v3;
		  }
		*/

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


