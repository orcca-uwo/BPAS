#include <bpas.h>
#include <string>
#include <sstream>

using namespace std;

int str2int (const string &str) {
  stringstream ss(str);
  int num;
  if((ss >> num).fail())
  { 
      cerr << "Cannot convert alphanumeric string to int" << endl;
  }
  return num;
}

void parseCoefficient(string coefficient, RN *coef){
	unsigned uint = coefficient.find("/");
	if (uint > coefficient.length()) {
		int icoef = str2int(coefficient);
		*coef = RN(mpq_class(icoef));
	}
	else {
		int num,den;
		string subcoef,subcoef2;
		unsigned uint2 = coefficient.find("(");
		if (uint2 > coefficient.length()) {
			subcoef = coefficient.substr(0,uint);
			num = str2int(subcoef);
			subcoef2 = coefficient.substr(uint+1,coefficient.length()-uint);
		}
		else {
			string subcoef = coefficient.substr(1,uint-1);
			num = str2int(subcoef);
			subcoef2 = coefficient.substr(uint+1,coefficient.length()-uint-2);
		}
		den = str2int(subcoef2);
		RN a(num,den);
		*coef = a;
	}
}

void parseTerm(int sign, string term, string variable, SparseUnivariatePolynomial<RN> *poly){
	unsigned uint;
	string coefficient;
	RN coef;
	uint = term.find(variable);
	if (uint > term.length()) {
		coefficient = term.substr(0,uint);
		parseCoefficient(coefficient,&coef);
		coef *= sign;
		poly->setCoefficient(0,coef);
	}
	else {
		if (term.length() == 1)
			poly->setCoefficient(1,RN(mpq_class(sign)));
		else {
			int exp;
			uint = term.find("*");
			if (uint > term.length())
				coef = 1;
			else {
				coefficient = term.substr(0,uint);
				parseCoefficient(coefficient,&coef);
			}
			uint = term.find("^");
			if (uint > term.length()){
				coef *= sign;
				poly->setCoefficient(1,coef);
			}
			else {
				string exponent = term.substr(uint+1,term.length());
				exp = str2int(exponent);
				coef *= sign;
				poly->setCoefficient(exp,coef);
			}
		}
	}
}

/*void parseTerm(int sign, string term, string variable, SparseUnivariatePolynomial<RN> *poly){
	unsigned uint;
	string coefficient;
	RN coef;
	uint = term.find(variable);
	if (uint > term.length()) {
		coefficient = term.substr(0,uint);
		coef = (RN)str2int(coefficient);
		coef *= sign;
		poly->setCoefficient(0,coef);
	}
	else {
		if (term.length() == 1)
			poly->setCoefficient(1,(RN)sign);
		else {
			int exp;
			uint = term.find("*");
			if (uint > term.length())
				coef = 1;
			else {
				coefficient = term.substr(0,uint);
				coef = (RN)str2int(coefficient);
			}
			uint = term.find("^");
			string exponent = term.substr(uint+1,term.length());
			exp = str2int(exponent);
			coef *= sign;
			poly->setCoefficient(exp,coef);
		}
	}
}*/

SparseUnivariatePolynomial<RN> polyParse(string s){
	SparseUnivariatePolynomial<RN> output;
	
	// find variable //
	string alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	unsigned uint = s.find_first_of(alpha);
	int sgn = 1;
	if (uint < s.length()){
		stringstream ss;
		string variable;
		ss << s.at(uint);
		ss >> variable;
		output.setVariableName(Symbol(variable));
	
		// remove spaces //
		s.erase(std::remove(s.begin(),s.end(),' '),s.end());
	
		// get and parse first term //
		string signs = "+-";
		string term;
		uint = s.find_first_of(signs);
		if (uint != 0) {
			term = s.substr(0,uint);
		}
		else {
			uint = s.find_first_of(signs,uint+1);
			term = s.substr(1,uint-1);
			if (s.at(0) == '-')
				sgn = -1;
			else
				sgn = 1;
		}
		parseTerm(sgn,term,variable,&output);
		
		// get and parse remaining terms, if any //
		unsigned uint2;
		while (uint < s.length()){
			uint2 = s.find_first_of(signs,uint+1);
			if (s.at(uint) == '-')
				sgn = -1;
			else
				sgn = 1;
			term = s.substr(uint+1,uint2-uint-1);
			parseTerm(sgn,term,variable,&output);
			uint = uint2;
		}
	}
	else {
		string signs = "+-";
		string coefficient;
		RN coef;
		
		uint = s.find_first_of(signs);
		if (uint == 0){
			if (s.at(0) == '-')
				sgn = -1;
			coefficient = s.substr(1,s.length()-1);
			parseCoefficient(coefficient,&coef);
			coef *= sgn;
			output.setCoefficient(0,coef);
		}
		else {	
			parseCoefficient(s,&coef);
			output.setCoefficient(0,coef);
		}
	}	
	return output;
}
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


