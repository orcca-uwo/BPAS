#ifndef _MPOLYNOMIAL_H_
#define _MPOLYNOMIAL_H_

/**
 * Data Structure for multivariate polynomial
 * stored in a sparse case
 **/

#include "../polynomial.h"
#include "../global.h"
#include "../DyadicRationalNumber/globals.h"
#include "upolynomial.h"
// #include "../RationalNumberPolynomial/globals.h" NO LONGER NEEDED
#include "../Utils/TemplateHelpers.hpp"

template <class Ring>
class SparseMultivariatePolynomial : public BPASMultivariatePolynomial<Ring,SparseMultivariatePolynomial<Ring>>,
									 private Derived_from<Ring, BPASRing<Ring>> {
	private:
		std::string mPoly;
		int var;		// Number of variables
		Symbol* names;	// Variable names
					// names[0] = "1" or "9" (user specified)
					// names[1] < names[2] < ...
		std::vector< MultivariateTerm<Ring> > terms;	// Each term of the polynomial

		/* Return the first non-zero variate in the ascending order */
		inline int firstNonZero(int* d) const {
			for (int i = 0; i < var; i++) {
				 if (d[i])
					return i;
			}
			return -1;
		}
		/* Return the first non-equal variate in the descending order */
		inline int checkDegs(int* d, int* b) const {
			for (int i = var-1; i > -1; --i) {
				if (d[i] != b[i])
					return i;
			}
			return -1;
		}
		inline bool isOrderedRing(const SparseMultivariatePolynomial<Ring>& b, std::vector<int>& xs) const {
			if (names[0] != b.names[0]) { return 0; }
			if (names[0] == "1") {
				int k = 1;
				for (; k <= var && k <= b.var; ++k) {
					xs.push_back(k);
					xs.push_back(k);
				}
				for (int i = k; i <= var; ++i) {
					xs.push_back(i);
					xs.push_back(0);
				}
				for (int i = k; i <= b.var; ++i) {
					xs.push_back(0);
					xs.push_back(i);
				}
				return 1;
			}
			if (!var) {
				for (int i = 1; i <= b.var; ++i) {
					xs.push_back(0);
					xs.push_back(i);
				}
				return 1;
			}
			if (!b.var) {
				for (int i = 1; i <= var; ++i) {
					xs.push_back(i);
					xs.push_back(0);
				}
				return 1;
			}

			bool isFound = 0;
			int* pos = new int[var];
			for (int i = 1; i <= var; ++i) {
				isFound = 0;
				for (int j = 1; j <= b.var; ++j) {
					if (names[i] == b.names[j]) {
						pos[i-1] = j;
						isFound = 1;
						break;
					}
				}
				if (!isFound) { pos[i-1] = 0; }
			}

			isFound = 0;
			int ak = 1, bk = pos[0];
			if (pos[0]) {
				for(int j = 1; j < pos[0]; ++j) {
					xs.push_back(0);
					xs.push_back(j);
				}
				xs.push_back(1);
				xs.push_back(pos[0]);
			}
			for (int i = 1; i < var; ++i) {
				if (pos[i] > bk) {
					while (pos[i] - bk > 1 && (i > ak || !bk)) {
						if (bk) { ak++; }
						bk++;
						if (names[ak] < b.names[bk]) {
							xs.push_back(ak);
							xs.push_back(0);
							xs.push_back(0);
							xs.push_back(bk);
						}
						else if (names[ak] > b.names[bk]) {
							xs.push_back(0);
							xs.push_back(bk);
							xs.push_back(ak);
							xs.push_back(0);
						}
					}
					for(int j = bk+1; j < pos[i]; ++j) {
						xs.push_back(0);
						xs.push_back(j);
					}
					for(int j = (!bk)? ak : ak+1; j <= i; ++j) {
						xs.push_back(j);
						xs.push_back(0);
					}
					xs.push_back(i+1);
					xs.push_back(pos[i]);
					bk = pos[i];
					ak = i + 1;
				}
				else if (pos[i] && pos[i] < bk) {
					isFound = 1;
					break;
				}
			}
			if (!isFound) {
				for (int i = (bk)? ak+1 : ak, j = bk+1; i <= var && j <= b.var; ++i, ++j) {
					if (names[i] < b.names[j]) {
						xs.push_back(i);
						xs.push_back(0);
						xs.push_back(0);
						xs.push_back(j);
					}
					else if (names[i] > b.names[j]) {
						xs.push_back(0);
						xs.push_back(j);
						xs.push_back(i);
						xs.push_back(0);
					}
					ak = i, bk = j;
				}
				for (int i = ak+1; i <= var; ++i) {
					xs.push_back(i);
					xs.push_back(0);
				}
				for (int j = bk+1; j <= b.var; ++j) {
					xs.push_back(0);
					xs.push_back(j);
				}
			}

			delete [] pos;

			if (isFound) { xs.clear(); return 0; }
			else { return 1; }
		}

		/**
		 * Compare two term degrees in the same Ring
		 *
		 * Return:
		 * 1: a > b
		 * 0: a = b
		 * -1: a < b
		 **/
		inline int compareTermDegs(const MultivariateTerm<Ring>& a, const MultivariateTerm<Ring>& b) const {
			for (int i = a.v-1; i > -1; --i) {
				if (a.degs[i] < b.degs[i])
					return -1;
				else if (a.degs[i] > b.degs[i])
					return 1;
			}
			return 0;
		}
		inline int compareTermDegs(const MultivariateTerm<Ring>& a, const MultivariateTerm<Ring>& b, std::vector<int> xs) const {
			int n = xs.size() / 2;
			for (int i = n-1; i > -1; --i) {
				if (xs[2*i] && xs[2*i+1]) {
					if (a.degs[xs[2*i]-1] < b.degs[xs[2*i+1]-1])
						return -1;
					else if (a.degs[xs[2*i]-1] > b.degs[xs[2*i+1]-1])
						return 1;
				}
				else if (!xs[2*i+1] && xs[2*i] && a.degs[xs[2*i]-1])
					return 1;
				else if (!xs[2*i] && xs[2*i+1] && b.degs[xs[2*i+1]-1])
					return -1;
			}
			return 0;
		}

		inline void basicOp(Ring c, bool isOp) {
			if (!c.isZero()) {
				if (terms.size()) {
					bool isIt = 1;
					for (int i = 0; i < var; ++i) {
						if (terms[0].degs[i]) {
							isIt = 0;
							break;
						}
					}
					if (isIt) {
						if (isOp)
							terms[0].coef -= c;
						else
							terms[0].coef += c;
						if (terms[0].coef == 0)
							terms.erase(terms.begin());
					}
					else {
						MultivariateTerm<Ring> a;
						if (isOp)
							a.coef = -c;
						else
							a.coef = c;
						a.v = var;
						a.degs = new int[var];
						for (int i = 0; i < var; ++i)
							a.degs[i] = 0;
						terms.insert(terms.begin(), a);
					}
				}
				else {
					MultivariateTerm<Ring> a;
					if (isOp)
						a.coef = -c;
					else
						a.coef = c;
					a.v = var;
					a.degs = new int[var];
					for (int i = 0; i < var; ++i)
						a.degs[i] = 0;
					terms.push_back(a);
				}
			}
		}

		inline void pomopo(const MultivariateTerm<Ring>& t, const SparseMultivariatePolynomial<Ring>& b, const std::vector<int>& xs) {
			for (int i = 0; i < b.terms.size(); ++i) {
				MultivariateTerm<Ring> a;
				a.coef = t.coef * b.terms[i].coef;
				a.v = var;
				a.degs = new int[var];
				for (int j = 0; j < var; ++j) {
					a.degs[j] = 0;
					if (xs[2*j])
						a.degs[j] += t.degs[xs[2*j]-1];
					if (xs[2*j+1])
						a.degs[j] += b.terms[i].degs[xs[2*j+1]-1];
				}
				int k = 0;
				while (k <= terms.size()) {
					if (k == terms.size()) {
						terms.push_back(a);
						break;
					}

					int s = compareTermDegs(terms[k], a);
					if (!s) {
						terms[k].coef += a.coef;
						if (terms[k].coef.isZero())
							terms.erase(terms.begin()+k);
						break;
					}
					else if (s > 0) {
						terms.insert(terms.begin()+k, a);
						break;
					}
					else { k++; }
				}
			}
		}

		/* Is equal to another polynomial */
		inline bool isEqual(const SparseMultivariatePolynomial<Ring>& b) const {
			if (terms.size() != b.terms.size()) { return 0; }

			std::vector<int> xs;
			bool isOrdered = isOrderedRing(b, xs);
			if (!isOrdered) { return 0; }

			int v = xs.size() / 2;
			for (int i = 0; i < terms.size(); ++i) {
				if (terms[i].coef != b.terms[i].coef)
					return 0;
				for (int j = 0; j < v; ++j) {
					if (xs[2*j] && xs[2*j+1] && (terms[i].degs[xs[2*j]-1] != b.terms[i].degs[xs[2*j+1]-1]))
						return 0;
					else if (!xs[2*j] && xs[2*j+1] && (b.terms[i].degs[xs[2*j+1]-1] != 0))
						return 0;
					else if (!xs[2*j+1] && xs[2*j] && (terms[i].degs[xs[2*j]-1] != 0))
						return 0;
				}
			}

			return 1;
		}

		/**
		 * private constructor to build the initial and tail polynomial
		 *
		 * @param
		 **/
		SparseMultivariatePolynomial<Ring>(int v, Symbol* xs, const std::vector<MultivariateTerm<Ring>>& newTerms) {

		  var = v;
		  names = new Symbol[var+1];
		  std::copy(xs, xs+var+1, names);


		  for (int i=0; i<newTerms.size(); i++){
		    terms.push_back(newTerms[i]);

		  }


		  Ring e;
		  characteristic = e.characteristic;

		}

		/**
		 * Second private constructor to build the initial and tail polynomial
		 *
		 * @param
		 **/
		inline void setTerms(const std::vector<MultivariateTerm<Ring>>& newTerms) {
		  for (int i=0; i<newTerms.size(); i++){
		    terms.push_back(newTerms[i]);

		  }
		  Ring e;
		  characteristic = e.characteristic;

		}


				//***********************************************amha parser source **********************************************************************************************


		/**
		 * This helper function check if a digit is numeric. Compatable with c++98 compiler.
		 *
		 * @param str : string value of number or ... .
		 * @return True if number else false.
		 */
		inline bool  is_digit_cpp98(const std::string &str)
		{
		    unsigned int i =0;
		    std::locale loc;
		    for(i=0; i<str.length(); i++)
		    {
			if(i==0 && str.length()>1 && str[i]=='-')
			{
			    continue;
			}
			if(!std::isdigit(str[i], loc))
			    return false;
		    }
		    return true;
		}

		inline bool isFloat( std::string myString ) {
			std::istringstream iss(myString);
			float f;
			iss >> std::noskipws >> f;
			return iss.eof() && !iss.fail();
		}

		/**
		 * This helper function compare c++ strings.
		 * @param str : string to compare
		 * @param str1 : second string to compare
		 * @return 0 if match else non zero value.
		 */
		inline bool  isMinus(const std::string &str, const std::string &str1)
		{
		    return !str.compare(str1);
		}

		/**
		 * This helper function modifies a string by removing all the space.
		 *
		 * @param str string to modifiy.
		 */
		inline void  remove_all_space(std::string &str)
		{
		    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
		}

		/**
		 * This helper function return the exact position of a character in a string.
		 *
		 * @param v String to search the character.
		 * @param s String format of the character to search.
		 * @return numeric position of the character in the string.
		 */
		inline size_t  findCharPosition(std::string &v, std::string s)
		{
		    size_t it = v.find(s);
		    if(it != std::string::npos)
			return it;
		    else
			return -1;
		}

		/**
		 * This helper function return a substring from a string.
		 *
		 * @param s :The original string where the substring is extracted from.
		 * @param spos : Starting position.
		 * @param epos : ending position.
		 * @return  Substring of the original string.
		 */
		inline std::string  getSubString(std::string &s, int spos, int epos)
		{
		    return (s.substr(spos, epos));
		}

		/**
		 * The helper function rearrange a polynomial term by putting the coefficient
		 * to the front and followed by the rest of the term components.
		 * Example: -x11^4*(56*x12^3) return -56*x11^4*x12^3.
		 *
		 * @param str string polynomial term.
		 * @return rearrange string polynomial term.
		 */
		inline std::string  TermReorder(std::string &str)
		{
		    //term contains bracket , split it
		    std::vector<std::string> ret = SimpleMultiDelimStringTokenizer(str, "()");
		    int tempcoeff = 1;              // hold the coeff of the term by collecting all the numbers within the bracket
		    std::string finalString = "";   // after all the reordering done everything concatinated to this string
		    std::locale loc;                // used in std::isdigit

		    std::vector<std::string>::iterator it;      //iterator to loop through the splitted tokens
		    // remove any "*" at the beginning , ending or if exists by itself.
		    for(it= ret.begin(); it<ret.end(); it++)
		    {
			std::string stemp = *it;

			if(stemp=="*")
			{
			    ret.erase(it);  // by itself
			}
			else if(findCharPosition(stemp, "*")==0)
			{
			    *it = getSubString(stemp, 1, stemp.length());   // at beginning
			}
			else if(findCharPosition(stemp, "*")==(stemp.length()-1))
			{
			    *it = getSubString(stemp, 0, stemp.length()-1); // at the end of the term
			}
		    }
		    for(it= ret.begin(); it<ret.end(); it++)
		    {
			std::string stemp = *it;
			// not one char, 1st char is '-' and 1st char is not digit, coeff is -1 and
			// '-' sign is removed from the token
			if(stemp.length()>1 && stemp[0]=='-' && !is_digit_cpp98(std::string(1, stemp[1])))
			{
			    tempcoeff *= -1;
			    *it = getSubString(stemp, 1, stemp.length());
			}
			// token start with digit strings, extract those digit
			// convert them and multply them to the tempcoeff temporary storage variable
			// then adjust the token  so it only contain the substring after the "*" sign.
			else if(stemp[0]!='-' && std::isdigit(stemp[0], loc))
			{
			    int tempint;
			    int pos = findCharPosition(stemp, "*");
			    std::string tempstring = getSubString(stemp, 0, pos);
			    std::stringstream ss(tempstring);
			    ss >> tempint;
			    tempcoeff = tempcoeff * tempint;
			    *it = getSubString(stemp, pos+1, stemp.length());
			}
			else if(stemp[0]=='-' && std::isdigit(stemp[1], loc))
			{
			    int tempint;
			    int pos = findCharPosition(stemp, "*");
			    std::string tempstring = getSubString(stemp, 0, pos);
			    std::stringstream ss(tempstring);
			    ss >> tempint;
			    tempcoeff = tempcoeff * tempint;
			    *it = getSubString(stemp, pos+1, stemp.length());
			}
			else if(stemp.size()==1 && stemp[0]=='-')
			{
			    tempcoeff = tempcoeff * -1;
			    //ret.erase(it);
			    //*it=getSubString(stemp, 3, stemp.length()-1);
			}
			else if(findCharPosition(stemp, "*")>1 && is_digit_cpp98(getSubString(stemp, findCharPosition(stemp, "*")+1, stemp.length())))
			{
			    int tempint;
			    int pos = findCharPosition(stemp, "*");
			    std::string tempstring = getSubString(stemp, findCharPosition(stemp, "*")+1, stemp.length());
			    std::stringstream ss(tempstring);
			    ss >> tempint;
			    tempcoeff = tempcoeff * tempint;
			    *it = getSubString(stemp, 0, pos);
			}
		    }
		    std::ostringstream oss;
		    oss << tempcoeff;
		    finalString += oss.str();
		    for(it= ret.begin(); it<ret.end(); it++)
		    {
			std::string stemp = *it;
			finalString+="*";
			finalString+=stemp;
		    }
		    return finalString;
		}


		inline std::vector<std::string>  SplitPolyToTerms(const std::string poly, std::string delim)
		{
		     //PrintToConsole(poly);
		    std::vector<std::string> tokens;
		    std::vector<std::string>::iterator tokenit;
		    std::size_t prev_pos = 0;
		    std::size_t cur_pos;
		    std::size_t prev_pos_temp = 0;

		    while((cur_pos = poly.find_first_of(delim, prev_pos)) != std::string::npos)
		    {
			if(poly[cur_pos-1]!='(')
			{
			    if((cur_pos > prev_pos_temp) )
			    {
				if(prev_pos != 0  && strncmp(&poly.at(prev_pos_temp-1), MINUS, 1) == 0)
				    tokens.push_back(poly.substr(prev_pos_temp-1, cur_pos-prev_pos_temp+1));
				else
				    tokens.push_back(poly.substr(prev_pos_temp, cur_pos-prev_pos_temp));
			    }
			    prev_pos_temp = cur_pos+1;
			}

			prev_pos = cur_pos+1;
		    }

		    if(prev_pos_temp < poly.length())
		    {
			if(prev_pos_temp!=0 && strncmp(&poly.at(prev_pos_temp-1), MINUS, 1) == 0 )
			    tokens.push_back(poly.substr(prev_pos_temp-1, std::string::npos));
			else
			    tokens.push_back(poly.substr(prev_pos_temp, std::string::npos));
		    }

		    for(tokenit=tokens.begin(); tokenit<tokens.end(); tokenit++)
		    {
			std::string tempstring = *tokenit;
			if(findCharPosition(tempstring, "(")!=-1 && findCharPosition(tempstring, ")")!=-1)
			{
			    *tokenit = TermReorder(tempstring);
			}
		    }

		    return tokens;
		}

		inline std::vector<std::string>  SimpleMultiDelimStringTokenizer(const std::string poly, std::string delim)
		{
		    std::vector<std::string> tokens;
		    std::size_t prev_pos = 0;
		    std::size_t cur_pos;

		    while((cur_pos = poly.find_first_of(delim, prev_pos)) != std::string::npos)
		    {
			if((cur_pos > prev_pos) )
			{
			    std::string stemp = poly.substr(prev_pos, cur_pos-prev_pos);
			    remove_all_space(stemp);
			    tokens.push_back(stemp);
			}
			prev_pos = cur_pos+1;
		    }

		    if(prev_pos < poly.length())
		    {
			std::string stemp = poly.substr(prev_pos, std::string::npos);
			remove_all_space(stemp);
			tokens.push_back(stemp);
		    }

		    return tokens;
		}

		inline std::vector<std::string>  AllPolyVaraibles(std::vector<std::string>& allTerms)
		{
		    std::vector<std::string> vars(0);

		    std::vector<std::string>::iterator iterm(0);
		    for( iterm=allTerms.begin(); iterm<allTerms.end(); iterm++)
		    {

			std::vector<std::string> r = SimpleMultiDelimStringTokenizer(*iterm, "-*^().");

			std::vector<std::string>::iterator it(0);
			for( it=r.begin(); it<r.end(); it++)
			{
			    std::string temp = *it;
			    remove_all_space(temp); // clean unnecessary space eg.-yyy when tokenized by "-" return " yyy" and !="yyy"
			    if(r.size() == 1)   //term split result is one, meaning no "*^", its is either eg.5x, x or -7y
			    {
				    if( (!is_digit_cpp98(temp)) && (std::find(vars.begin(), vars.end(), temp) == vars.end()) )
				    {
				        vars.push_back(temp/*ctos*/);
				    }
			    }
			    else if( (!is_digit_cpp98(temp)) && (std::find(vars.begin(), vars.end(), temp) == vars.end()) )
			    {
				vars.push_back(temp);
			    }
			}
		    }
		    std::sort(vars.begin(), vars.end());

		    return vars;
		}


		inline std::vector<std::vector<int> >   AllPolyTermsExponents(std::vector<std::string>& allvariables, std::vector<std::string>& polyterms)
		{
		    std::vector<std::vector<int> > exponent(polyterms.size(), std::vector<int>(0));

		    for(unsigned int i=0; i<polyterms.size(); i++)
		    {
			std::vector<std::string> termtokenized = SimpleMultiDelimStringTokenizer(polyterms[i], "^*");
			std::vector<std::string>::iterator ivars;
			for( ivars=allvariables.begin(); ivars<allvariables.end(); ivars++)
			{

			    unsigned int position;
			    //check if length of string is one, same for the second else if
			    //if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && std::isalpha(termtokenized[0][1]))
			    if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && !is_digit_cpp98(termtokenized[0].substr(1, termtokenized[0].length())))
			    {
				//termtokenized[0] = termtokenized[0][1];
				termtokenized[0] = termtokenized[0].substr(1, termtokenized[0].length());
				position = std::find(termtokenized.begin(), termtokenized.end(), *ivars) - termtokenized.begin();
			    }
			    //else if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && std::isdigit(termtokenized[0][1]) )
			    else if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && is_digit_cpp98(termtokenized[0].substr(1, termtokenized[0].length())))
			    {
				//termtokenized[0] = termtokenized[0][2]; // remove the '-' find position of variable in the token
				termtokenized[0] = termtokenized[0].substr(2, termtokenized[0].length());
				position = std::find(termtokenized.begin(), termtokenized.end(), *ivars) - termtokenized.begin();
			    }
			    else
			    {
				position = std::find(termtokenized.begin(), termtokenized.end(), *ivars) - termtokenized.begin();
			    }
			    if(position < termtokenized.size())
			    {
				if((position == termtokenized.size()-1) || !is_digit_cpp98(termtokenized[position+1])/*(std::isalpha(*termtokenized.at(position+1).c_str()))*/)
				{
				    exponent[i].push_back(1);
				}
				else
				{
				    exponent[i].push_back(std::atoi(termtokenized.at(position+1).c_str()));
				}
			    }
			    else if(position == termtokenized.size())
			    {
				exponent[i].push_back(0);
			    }
			}
		    }

		    return exponent;
		}
		/*
		inline std::vector<std::vector<int> >  AllPolyTermsCoeffs(std::vector<std::string>& allvariables, std::vector<std::string>& polyterms)
		{
		    // contain all the coeffs
		    std::vector<std::vector<int> > coeff(polyterms.size(), std::vector<int>(0));

		    // start with the first poly term that is provided in the argument
		    for(unsigned int i=0; i<polyterms.size(); i++)
		    {
			// split the term based the given delimeter and store them in a container
			std::vector<std::string> termtokenized = SimpleMultiDelimStringTokenizer(polyterms[i], "^*");
			    // first string is a minus sign by itself and second string a string made of alphabets
			    if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && !is_digit_cpp98(termtokenized[0].substr(1, termtokenized[0].length())))
			    {
				coeff[i].push_back(-1);
			    }
			    // first string a string made of alphabets
			    else if((termtokenized.size()==1) && !is_digit_cpp98(termtokenized[0]))
			    {
				coeff[i].push_back(1);
			    }
			    // termtokenized is not size=1 and first char is alpha
			    else if(!is_digit_cpp98(termtokenized[0]))
			    {
				coeff[i].push_back(1);
			    }
			    else if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && !is_digit_cpp98(termtokenized[0].substr(1, termtokenized[0].length())))
			    {
				coeff[i].push_back(-1);
			    }
			    // or termtokenizer is in the form -8*x^8 or 9*y^7 and so....
			    else
			    {
				int tempint;
				std::string tempstring;
				remove_all_space(termtokenized[0]);
				tempstring = termtokenized[0];
				std::stringstream(tempstring) >> tempint;
				coeff[i].push_back(tempint);
			    }
		    }

		    return coeff;
		}
*/
		inline std::vector<std::vector<mpq_class> > AllPolyTermsCoeffs(std::vector<std::string>& allvariables, std::vector<std::string>& polyterms)
		{
			// contain all the coeffs
			std::vector<std::vector<mpq_class> > coeff(polyterms.size(), std::vector<mpq_class>(0));

			// start with the first poly term that is provided in the argument
			for(unsigned int i=0; i<polyterms.size(); i++)
			{
				// split the term based the given delimeter and store them in a container
				std::vector<std::string> termtokenized = SimpleMultiDelimStringTokenizer(polyterms[i], "^*");
				//std::cout << "tok : " << termtokenized[0]<< std::endl;
				//std::cout << is_digit_cpp98(termtokenized[0]) << std::endl;
				    // first string is a minus sign by itself and second string a string made of alphabets
				    if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && !is_digit_cpp98(termtokenized[0].substr(1, termtokenized[0].length()))&&!isFloat(termtokenized[0]))
				    {
				        coeff[i].push_back(-1);
				    }
				    // first string a string made of alphabets
				    else if((termtokenized.size()==1) && !is_digit_cpp98(termtokenized[0]) &&!isFloat(termtokenized[0]))
				    {
				        coeff[i].push_back(1);
				    }
				    // termtokenized is not size=1 and first char is alpha
				    else if(!is_digit_cpp98(termtokenized[0]) && !isFloat(termtokenized[0]))
				    {
				        coeff[i].push_back(1);
				    }
				    else if(isMinus(std::string(1, termtokenized[0][0]), MINUS) && !is_digit_cpp98(termtokenized[0].substr(1, termtokenized[0].length())) && !isFloat(termtokenized[0]))
				    {
				        coeff[i].push_back(-1);
				    }
				    // or termtokenizer is in the form -8*x^8 or 9*y^7 and so....
				    else
				    {
				        //std::cout << "here : " << termtokenized[0] << std::endl;
				        mpq_class tempint;
				        std::string tempstring;
				        remove_all_space(termtokenized[0]);
				        tempstring = termtokenized[0];
				        std::istringstream iss(tempstring);
				        iss >> tempint;
				        coeff[i].push_back(tempint);
				    }
			}

			return coeff;
		}
		inline PolyAndExp  SplitPolyInBracket()
		{
		    /*potental bug in this function, if for example
		      f = (x^5+y-(x^6*5)+z^3)^4 is given, the function
		      will split it into 3 part based on the given delimeters, "()".
		      so not recommanded to nest bracket if a polynomial
		      has an outer bracket.

		    */
		  PolyAndExp retVal;
		  std::vector<std::string> vals(0);
		  if(mPoly.at(0) == '(')
		  {
		    vals = SimpleMultiDelimStringTokenizer(mPoly, "()");
		    vals[1]= getSubString(vals[1], 1, vals[1].length()); //remove the carrot sign
		  }
		  else
		  {
		    vals.push_back(mPoly);
		  }
		    // the returned token contain bracket and carrot sign or exponent
		    if((ValidPoly()) && (mPoly.find("(") != std::string::npos) && (mPoly.find(")") != std::string::npos) && (vals.size() == 2))
		    {

			retVal.exp = std::atoi(vals.at(1).c_str());
			retVal.poly = vals.at(0);

			return retVal;
		    }
		    // the return token only contain brakets so set the exponent to 1
		    else if((ValidPoly()) && (mPoly.find("(") != std::string::npos) && (mPoly.find(")") != std::string::npos))
		    {
			retVal.exp = 1;
			retVal.poly = vals.at(0);

			return retVal;
		    }
		    // no bracket but valid poly
		    else if(ValidPoly())
		    {
			retVal.exp = 1;
			retVal.poly = mPoly;

			return retVal;
		    }
		    // invalid polynomial
		    else
		    {
			std::cout << "Invalid Polynomial Format ! " << std::endl;
			retVal.exp = -1;
			retVal.poly = " ";
			return retVal;
		    }
		}


		inline bool  ValidPoly()
		{
		    if(mPoly[0] == '*')
			return false;
		    else if(mPoly[0] == '+')
			return false;
		    else if(mPoly.find("((") != std::string::npos)
			return false;
		    else if(mPoly.find("))") != std::string::npos)
			return false;
		    else if(mPoly.find("++") != std::string::npos)
			return false;
		    else if(mPoly.find("--") != std::string::npos)
			return false;
		    else if(mPoly.find("**") != std::string::npos)
			return false;
		    else if(mPoly[mPoly.length()-1] == '^')
			return false;
		    else
			return true;
		}



		inline int  PolyExponent()
		{
		    PolyAndExp exp = SplitPolyInBracket();
		    return exp.exp;
		}



	public:
		mpz_class characteristic;
		RingProperties properties;
		// static bool isPrimeField;
		// static bool isSmallPrimeField;
        // static bool isComplexField;
		/**
		 * Construct a multivariate polynomial
		 *
		 * @param
		 **/
		SparseMultivariatePolynomial<Ring>() : var(0) {
			names = new Symbol[1];
			names[0] = "1";
			Ring e;
			characteristic = e.characteristic;
		}
		/**
		 * Construct a multivariate polynomial with number of terms and variates
		 *
		 * @param v: Number of variables
		 **/
		SparseMultivariatePolynomial<Ring>(int v) {
			var = v;
			names = new Symbol[var+1];
			names[0] = "1";
			for (int i = 1; i <= var; ++i) {
				std::ostringstream convert;
				convert << var - i + 1;
				names[i] = "_";
				names[i] += convert.str();
			}
			Ring e;
			characteristic = e.characteristic;
		}
		/**
		 * Construct with a variable name
		 * such that f(x) = x
		 *
		 * @param x: The variable name
		 **/
		SparseMultivariatePolynomial<Ring> (const Symbol& x) {
			var = 1;
			names = new Symbol[2];
			names[0] = "9";
			names[1] = x;
			MultivariateTerm<Ring> t;
			t.coef.one();
			t.v = 1;
			t.degs = new int[1];
			t.degs[0] = 1;
			terms.push_back(t);
			Ring e;
			characteristic = e.characteristic;
		}

		/**
		 * Construct the constant polynomial.
		 */
		SparseMultivariatePolynomial<Ring>(const Ring& r) : var(0) {
			names = new Symbol[1];
			names[0] = "1";
			Ring e;
			characteristic = e.characteristic;
			MultivariateTerm<Ring> a;
			a.v = var;
			a.coef = r;
			a.degs = NULL;
			terms.push_back(a);
		}

		/**
		 * Copy Constructor
		 *
		 * @param b: A sparse multivariate polynomial
		 **/
		SparseMultivariatePolynomial<Ring>(const SparseMultivariatePolynomial<Ring>& b) : var(b.var), terms(b.terms) {
			names = new Symbol[var+1];
			std::copy(b.names, b.names+var+1, names);
			Ring e;
			characteristic = e.characteristic;
		}
		/**
		 * Construct from a SUP<SMQP> polynomial
		 *
		 * @param s: The SUP<SMQP> polynomial
		 **/
		SparseMultivariatePolynomial<Ring> (SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> >& s) {
			names = new Symbol[1];
			int d = s.degree().get_si();
			if(d==0){
				*this = s.coefficient(0);
			}else{

			for (int k = 0; k <= d; ++k) {
				SparseMultivariatePolynomial<Ring> t, c = s.coefficient(k);
				t.var = c.var + 1;
				t.names = new Symbol[t.var+1];
				if (c.var)
					std::copy(c.names, c.names+t.var, t.names);
				else
					t.names[0] = "9";
				t.names[t.var] = s.variable();
				for (int i = 0; i < c.terms.size(); ++i) {
					MultivariateTerm<Ring> a;
					a.coef = c.terms[i].coef;
					a.v = t.var;
					a.degs = new int[t.var];
					for (int j = 0; j < c.var; ++j)
						a.degs[j] = c.terms[i].degs[j];
					a.degs[c.var] = k;
					t.terms.push_back(a);
				}
				if (k) { *this = *this + t; }
				else { *this = t; }
			}
			Ring e;
			characteristic = e.characteristic;
		}
		}
		/**
		 * Destroy the polynomial
		 *
		 * @param
		 **/
		~SparseMultivariatePolynomial<Ring>() {
			terms.clear();
			delete [] names;
		}

		/**
		 * Get number of variables
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return variables().size();
		}

		/**
		 * Get the number of variables in this polynomial ring.
		 */
		inline int numberOfRingVariables() const {
			return ringVariables().size();
		}

		/**
		 * Get number of non-zero terms
		 *
		 * @param
		 **/
		inline Integer numberOfTerms() const{
			return terms.size();
		}
		/**
		 * Get the degree of a variable
		 *
		 * @param x: The variable name
		 **/
		inline Integer degree(const Symbol& x) const {
			int k = -1;
			for (int i = 1; i <= var; ++i) {
				if (names[i] == x) {
					k = i-1;
					break;
				}
			}
			if (k < 0) { return 0; }
			int d = 0;
			for (int i = terms.size()-1; i > -1; --i) {
				if (terms[i].degs[k] > d)
					d = terms[i].degs[k];
			}
			return d;
		}

		/**
		 * Get the total degree.
		 */
		inline Integer degree() const {
			int totalMax = 0, total;
			for (int i = 0; i < terms.size(); ++i) {
				total = 0;
				for (int k = 0; k < var; ++k) {
					total += terms[i].degs[k];
				}
				if (total > totalMax) {
					totalMax = total;
				}
			}
			return totalMax;
		}

		/**
		 * Get the leading variable
		 *
		 * @param
		 **/
		inline Symbol leadingVariable() const {
			return names[var];
		}
		/**
		 *  Get the leading term's leading variable's degree
		 *
		 *  @param
		 **/
		inline int leadingVariableDegree() const {
			int n = terms.size();
			if (n)
				return terms[n-1].degs[var-1];
			else
				return 0;
		}
		/**
		 * Get the leading coefficient
		 *
		 * @param
		 **/
		inline Ring leadingCoefficient() const {
			int n = terms.size();


			if (n){
			  return terms[n-1].coef;
			}
			else{
			  Ring r;
			  return r;
			}
		}

		inline Ring trailingCoefficient() const {
			int n = terms.size();
			if (n) {
				return terms[0].coef;
			}
			return Ring(0);
		}

		inline bool isConstantTermZero() const {
			int n = terms.size();
			if (n) {
				for (int i = 0; i < var; ++i) {
					if (terms[0].degs[i] != 0) {
						return true;
					}
				}

				//at this point we know terms[0] is a const term.
				return (terms[0].coef == 0);
			}
			return true;
		}

		/**
		 * Get the initial
		 *
		 * @param
		 **/
		inline SparseMultivariatePolynomial<Ring>  initial() const {


		  /*
SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > uPoly = this->convertToSUP(this->leadingVariable());
return uPoly.leadingCoefficient();*/


			  int n = terms.size();

			  if (n){
				    std::vector<MultivariateTerm<Ring> > initTerms;
				    Symbol name = this->leadingVariable();
				    int position=-1;
				    for(int j=1; j<=var ;j++){
				      if(name == names[j]){
				      	position = j;
				      	break;
				      }
				    }

				    for (int i=0 ; i< n; i++){
				      if(this->leadingVariableDegree() == terms[i].degs[position-1]){
						MultivariateTerm<Ring> term;
						term = terms[i];
						term.degs[position-1] = 0;
						initTerms.push_back(term);
				      }
				    }
					int newVar = var-1;
					std::vector<Symbol> newNames;
					for(int i=0; i<=var; i++){
						if(names[i] != this->leadingVariable() && names[i]!="9" && names[i]!="1"){
							
							newNames.push_back(names[i]);
						}
					}

					//std::cout << "Number of var: " << var << std::endl;
				    //SparseMultivariatePolynomial init(var, names, initTerms);
					SparseMultivariatePolynomial init(var, names, initTerms);

					//std::cout << "newNames size : " << newNames.size() << "  newvar size:"  << var << std::endl;
					//init.setVariableNames(newNames);
					//init.setTerms(initTerms);

				    std::stringstream ss;
				    ss << init;
				    std::string s1= ss.str();

				    SparseMultivariatePolynomial init2;
				    init2 = s1;
				    return init2;
				    //return init;
			  }else{
			    	SparseMultivariatePolynomial init(0);
			   		return init;
			  }

		}


		inline SparseMultivariatePolynomial<Ring>  tail() const {
			/*
			SparseMultivariatePolynomial<Ring> temp = this->initial();
			SparseMultivariatePolynomial<Ring> temp2 = *this;

			std::stringstream ss;
			ss << this->leadingVariableDegree ();
			std::string val = this->leadingVariable() + "^" + ss.str();
			SparseMultivariatePolynomial<Ring> temp3(0);
			temp3 = val;

			SparseMultivariatePolynomial<Ring> res = temp2 - (temp*temp3);
			std::cout << "inside: " << res <<std::endl;
			return res;
			*/

		  int n = terms.size();
		  if (n){
		    std::vector<MultivariateTerm<Ring> > tailTerms;
		    Symbol name = this->leadingVariable();
		    int position=-1;
		    for(int j=1; j<=var ;j++){
		      if(name == names[j]){
			position = j;
			break;
		      }
		    }
		    for (int i=0 ; i<n; i++){
		      if(this->leadingVariableDegree() != terms[i].degs[position-1]){
			tailTerms.push_back(terms[i]);
		      }
		    }
		    SparseMultivariatePolynomial tailResult(var, names,tailTerms);
			std::stringstream tostream;
			tostream << tailResult;
			std::string tostring = tostream.str();
			SparseMultivariatePolynomial finaltailResult(0);
			finaltailResult = tostring;
		    return finaltailResult;
		  }
		  else{
		    SparseMultivariatePolynomial tailResult(0);
		    return tailResult;
		  }

		}


		/**
		 * Get the leading coefficient over a variable
		 *
		 * @param x: The name of the variable
		 * @param e: The leading exponent of the variable
		 **/
		inline SparseMultivariatePolynomial<Ring> leadingCoefficientInVariable (const Symbol& x, int* e=NULL) {
			if (e == NULL)
				e = new int;
			*e = 0;

			if (isConstant()) {
				SparseMultivariatePolynomial<Ring> r;
				r.var = var;
				r.names = new Symbol[var+1];
				std::copy(names, names+var+1, r.names);
				return r;
			}

			int k = 0;
			for (int i = 1; i <= var; ++i) {
				if (names[i] == x) {
					k = i;
					break;
				}
			}
			int v = var - 1;
			SparseMultivariatePolynomial<Ring> r(v);
			for (int i = 0; i < var; ++i) {
				if (i < k)
					r.names[i] = names[i];
				else
					r.names[i] = names[i+1];
			}
			if (!k) { return r; }
			k--;
			for (int i = 0; i < terms.size(); ++i) {
				if (terms[i].degs[k] >= *e) {
					if (terms[i].degs[k] > *e) {
						*e = terms[i].degs[k];
						r.terms.clear();
					}
					MultivariateTerm<Ring> t;
					t.coef = terms[i].coef;
					t.v = v;
					t.degs = new int[v];
					for (int j = 0; j < v; ++j) {
						if (j < k)
							t.degs[j] = terms[i].degs[j];
						else
							t.degs[j] = terms[i].degs[j+1];
					}
					r.terms.push_back(t);
				}
			}
			return r;
		}

		/**
		 * Is zero polynomial
		 *
		 * @param
		 **/
		inline bool isZero() const {
			return !terms.size();
		}
		/**
		 * Zero polynomial
		 *
		 * @param
		 **/
		inline void zero() {
			terms.clear();
		}
		/**
		 * Is polynomial a constant 1
		 *
		 * @param
		 **/
		inline bool isOne() const {
			if (terms.size() == 1) {
				for (int i = 0; i < var; ++i) {
					if (terms[0].degs[i] != 0)
						return 0;
				}
				if (terms[0].coef.isOne())
					return 1;
			}
			return 0;
		}

		/**
		 * Set polynomial to 1
		 *
		 * @param
		 **/
		inline void one() {
			terms.clear();
			MultivariateTerm<Ring> t;
			t.coef.one();
			t.v = var;
			t.degs = new int[var];
			for (int i = 0; i < var; ++i)
				t.degs[i] = 0;
			terms.push_back(t);
		}
		
		/**
		 * Is polynomial a constant -1
		 *
		 * @param
		 **/
		inline bool isNegativeOne() const {
			if (terms.size() == 1) {
				for (int i = 0; i < var; ++i) {
					if (terms[0].degs[i] != 0)
						return 0;
				}
				if (terms[0].coef.isNegativeOne())
					return 1;
			}
			return 0;
		}

		/**
		 * Set polynomial to -1
		 *
		 * @param
		 **/
		inline void negativeOne() {
			terms.clear();
			MultivariateTerm<Ring> t;
			t.coef.negativeOne();
			t.v = var;
			t.degs = new int[var];
			for (int i = 0; i < var; ++i)
				t.degs[i] = 0;
			terms.push_back(t);
		}
		
		/**
		 * Is a constant
		 *
		 * @param
		 **/
		inline int isConstant() const {
			if (!terms.size())
				return 1;
			else if (terms.size() == 1) {
				for (int i = 0; i < var; ++i) {
					if (terms[0].degs[i] != 0)
						return 0;
				}
				if (terms[0].coef.isConstant() >= 0) { return 1; }
				else { return -1; }
			}
			return 0;
		}

		inline SparseMultivariatePolynomial<Ring> unitCanonical(SparseMultivariatePolynomial<Ring>* u = NULL, SparseMultivariatePolynomial<Ring>* v = NULL) const {
			Ring lead = leadingCoefficient();
			Ring ru, rv;
			Ring canon = lead.unitCanonical(&ru, &rv);

			SparseMultivariatePolynomial<Ring> ret = *this * ru;
			if (u != NULL) {
				*u = SparseMultivariatePolynomial<Ring>(ru);
			}
			if (v != NULL) {
				*v = SparseMultivariatePolynomial<Ring>(rv);
			}
			return ret;
		}
		
		/**
		 * Overload operator =
		 *
		 * @param b: A sparse multivariate polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator= (const SparseMultivariatePolynomial<Ring>& b) {
			if (this != &b) {
				terms.clear();
				delete [] names;

				var = b.numberOfVariables();
				names = new Symbol[var+1];
				std::copy(b.names, b.names+var+1, names);
				terms = b.terms;
				Ring e;
				characteristic = e.characteristic;
			}
			return *this;
		}

		inline SparseMultivariatePolynomial<Ring>& operator= (const Ring& r) {
			*this = SparseMultivariatePolynomial<Ring>(r);
			return *this;
		}

		SparseMultivariatePolynomial<Ring>& operator= (std::string &b) {

			mPoly = b;
			std::string valPoly;
			PolyAndExp pe = SplitPolyInBracket();
			valPoly = pe.poly;
			int e = pe.exp;

			std::vector<std::string> res = SplitPolyToTerms(valPoly, "+-");
			std::vector<std::string> allres = AllPolyVaraibles(res);
			std::vector<std::vector<int> > exp = AllPolyTermsExponents(allres, res);
			std::vector<std::vector<mpq_class> > coef = AllPolyTermsCoeffs(allres, res);

			SparseMultivariatePolynomial<Ring> t(allres.size());
			t.setRingVariables(allres);

			for (int i = 0; i < res.size(); ++i)
			{
				mpq_class elem = coef[i][0];
				int* d = new int[allres.size()];
				for (int j = 0; j < allres.size(); ++j)
				{
				        d[j] = exp[i][j];
				}
				t.setCoefficient(allres.size(), d, Ring(elem));
			}
			if(e==2)
			{
				*this *= t;
				return *this;
			}
			if(e >2)
			{
				*this = t;
				int ee = e -1 ;
				while(ee)
				{
					if(ee % 2) { *this *= t; }
					t *= t;
					ee >>=1;
				}
				return *this;
			}
			*this = t;
			return *this;
		}


		/**
		 * Overload operator ==
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline bool operator== (const SparseMultivariatePolynomial<Ring>& b) const {
			return isEqual(b);
		}
		
		/**
		 * Overload operator !=
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline bool operator!= (const SparseMultivariatePolynomial<Ring>& b) const {
			return !isEqual(b);
		}

		/**
		 * Overload operator ^
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
		inline SparseMultivariatePolynomial<Ring> operator^ (long long int e) const {
			SparseMultivariatePolynomial<Ring> res;
			res.var = var;
			res.names = new Symbol[var+1];
			std::copy(names, names+var+1, res.names);
			res.one();
			unsigned long int q = e / 2, r = e % 2;
			SparseMultivariatePolynomial<Ring> p2 = *this * *this;
			for (int i = 0; i < q; ++i)
				res *= p2;
			if (r) { res *= *this; }
			return res;
		}

		/**
		 * Overload operator ^=
		 * replace xor operation by exponentiation
		 *
		 * @param e: The exponentiation, e > 0
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator^= (long long int e) {
			*this = *this ^ e;
			return *this;
		}

		/**
		 * Overload operator +
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring> operator+ (const SparseMultivariatePolynomial<Ring>& b) const {
			if (!terms.size()) { return b; }
            if (!b.terms.size()) { return *this; }
            if (isConstant()) { return (b + terms[0].coef); }
            if (b.isConstant()) { return (*this + b.terms[0].coef); }
			std::vector<int> xs;
			bool isOrdered = isOrderedRing(b, xs);

			if (!isOrdered) {
				std::cout << "BPAS: error, trying to add between Ring[";
				for (int i = 1; i <= var; ++i) {
					std::cout << names[i];
					if (i < var) { std::cout << ", "; }
				}
				std::cout << "] and Ring[";
				for (int i = 1; i <= b.var; ++i) {
					std::cout << b.names[i];
					if (i < b.var) { std::cout << ", "; }
				}
				std::cout << "] from SparseMultivariatePolynomial<Ring>." << std::endl;
				exit(1);
			}

			int v = xs.size() / 2;
			SparseMultivariatePolynomial<Ring> r (v);
			r.names[0] = names[0];
			if (names[0] == "9") {
				for (int i = 1; i <= v; ++i) {
					if (xs[2*i-2])
						r.names[i] = names[xs[2*i-2]];
					else
						r.names[i] = b.names[xs[2*i-1]];
				}
			}

			int i = 0, j = 0;
			while (i < terms.size() && j < b.terms.size()) {
				MultivariateTerm<Ring> t;
				t.v = v;
				t.degs = new int[v];
				int c = compareTermDegs(terms[i], b.terms[j], xs);
				if (!c) {
					t.coef = terms[i].coef + b.terms[j].coef;
					if (!t.coef.isZero()) {
						for (int k = 0; k < v; ++k) {
							if (xs[2*k])
								t.degs[k] = terms[i].degs[xs[2*k]-1];
							else
								t.degs[k] = b.terms[j].degs[xs[2*k+1]-1];
						}
						r.terms.push_back(t);
					}
					i++, j++;
				}
				else if (c < 0) {
					t.coef = terms[i].coef;
					for (int k = 0; k < v; ++k) {
						if (xs[2*k])
							t.degs[k] = terms[i].degs[xs[2*k]-1];
						else
							t.degs[k] = 0;
					}
					r.terms.push_back(t);
					i++;
				}
				else {
					t.coef = b.terms[j].coef;
					for (int k = 0; k < v; ++k) {
						if (xs[2*k+1])
							t.degs[k] = b.terms[j].degs[xs[2*k+1]-1];
						else
							t.degs[k] = 0;
					}
					r.terms.push_back(t);
					j++;
				}
			}
			for (; i < terms.size(); ++i) {
				MultivariateTerm<Ring> t;
				t.v = v;
				t.coef = terms[i].coef;
				t.degs = new int[v];
				for (int k = 0; k < v; ++k) {
					if (xs[2*k])
						t.degs[k] = terms[i].degs[xs[2*k]-1];
					else
						t.degs[k] = 0;
				}
				r.terms.push_back(t);
			}
			for (; j < b.terms.size(); ++j) {
				MultivariateTerm<Ring> t;
				t.v = v;
				t.coef = b.terms[j].coef;
				t.degs = new int[v];
				for (int k = 0; k < v; ++k) {
					if (xs[2*k+1])
						t.degs[k] = b.terms[j].degs[xs[2*k+1]-1];
					else
						t.degs[k] = 0;
				}
				r.terms.push_back(t);
			}
			return r;
		}

		/**
		 * Overload operator +=
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator+= (const SparseMultivariatePolynomial<Ring>& b) {
			*this = *this + b;
			return *this;
		}
		/**
		 * Overload operator +
		 *
		 * @param c: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring> operator+ (const Ring& c) const {
			SparseMultivariatePolynomial<Ring> r (*this);
			return (r += c);
		}

		/**
		 * Overload operator +=
		 *
		 * @param c: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator+= (const Ring& c) {
			basicOp(c, 0);
			return *this;
		}

		inline friend SparseMultivariatePolynomial<Ring> operator+ (const Ring& c, const SparseMultivariatePolynomial<Ring>& p) {
			return (p + c);
		}

		/**
		 * Overload operator -
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring> operator- (const SparseMultivariatePolynomial<Ring>& b) const {
		 	if (!terms.size()) { return -b; }
		 	if (!b.terms.size()) { return *this; }
		 	if (isConstant()) { return (-b + terms[0].coef); }
		 	if (b.isConstant()) { return (*this - b.terms[0].coef); }

		 	std::vector<int> xs;
		 	bool isOrdered = isOrderedRing(b, xs);

		 	if (!isOrdered) {
		 		std::cout << "BPAS: error, trying to subtract between Ring[";
		 		for (int i = 1; i <= var; ++i) {
		 			std::cout << names[i];
		 			if (i < var) { std::cout << ", "; }
		 		}
		 		std::cout << "] and Ring[";
		 		for (int i = 1; i <= b.var; ++i) {
		 			std::cout << b.names[i];
		 			if (i < b.var) { std::cout << ", "; }
		 		}
		 		std::cout << "] from SparseMultivariatePolynomial<Ring>." << std::endl;
		 		exit(1);
		 	}

		 	int v = xs.size() / 2;
		 	SparseMultivariatePolynomial<Ring> r (v);
		 	r.names[0] = names[0];
		 	if (names[0] == "9") {
		 		for (int i = 1; i <= v; ++i) {
		 			if (xs[2*i-2])
		 				r.names[i] = names[xs[2*i-2]];
		 			else
		 				r.names[i] = b.names[xs[2*i-1]];
		 		}
		 	}

		 	int i = 0, j = 0;
		 	while (i < terms.size() && j < b.terms.size()) {
		 		MultivariateTerm<Ring> t;
		 		t.v = v;
		 		t.degs = new int[v];
		 		int c = compareTermDegs(terms[i], b.terms[j], xs);
		 		if (!c) {
		 			t.coef = terms[i].coef - b.terms[j].coef;
		 			if (!t.coef.isZero()) {
		 				for (int k = 0; k < v; ++k) {
		 					if (xs[2*k])
		 						t.degs[k] = terms[i].degs[xs[2*k]-1];
		 					else
		 						t.degs[k] = b.terms[j].degs[xs[2*k+1]-1];
		 				}
		 				r.terms.push_back(t);
		 			}
		 			i++, j++;
		 		}
		 		else if (c < 0) {
		 			t.coef = terms[i].coef;
		 			for (int k = 0; k < v; ++k) {
		 				if (xs[2*k])
		 					t.degs[k] = terms[i].degs[xs[2*k]-1];
		 				else
		 					t.degs[k] = 0;
		 			}
		 			r.terms.push_back(t);
		 			i++;
		 		}
		 		else {
		 			t.coef = -b.terms[j].coef;
		 			for (int k = 0; k < v; ++k) {
		 				if (xs[2*k+1])
		 					t.degs[k] = b.terms[j].degs[xs[2*k+1]-1];
		 				else
		 					t.degs[k] = 0;
		 			}
		 			r.terms.push_back(t);
		 			j++;
		 		}
		 	}
		 	for (; i < terms.size(); ++i) {
		 		MultivariateTerm<Ring> t;
		 		t.v = v;
		 		t.coef = terms[i].coef;
		 		t.degs = new int[v];
		 		for (int k = 0; k < v; ++k) {
		 			if (xs[2*k])
		 				t.degs[k] = terms[i].degs[xs[2*k]-1];
		 			else
		 				t.degs[k] = 0;
		 		}
		 		r.terms.push_back(t);
		 	}
		 	for (; j < b.terms.size(); ++j) {
		 		MultivariateTerm<Ring> t;
		 		t.v = v;
		 		t.coef = -b.terms[j].coef;
		 		t.degs = new int[v];
		 		for (int k = 0; k < v; ++k) {
		 			if (xs[2*k+1])
		 				t.degs[k] = b.terms[j].degs[xs[2*k+1]-1];
		 			else
		 				t.degs[k] = 0;
		 		}
		 		r.terms.push_back(t);
		 	}
		 	return r;
		}
		
		/**
		 * Overload operator -=
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator-= (const SparseMultivariatePolynomial<Ring>& b) {
			*this = *this - b;
			return *this;
		}
		
		/**
		 * Overload operator -
		 *
		 * @param c: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring> operator- (const Ring& c) const {
			SparseMultivariatePolynomial<Ring> r (*this);
			return (r -= c);
		}

		/**
		 * Overload operator -=
		 *
		 * @param c: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator-= (const Ring& c) {
			basicOp(c, 1);
			return *this;
		}

		inline friend SparseMultivariatePolynomial<Ring> operator- (const Ring& c, const SparseMultivariatePolynomial<Ring>& p) {
            return (-p + c);
        }

		/**
		 * Overload operator - (negate)
		 *
		 * @param
		 **/
		inline SparseMultivariatePolynomial<Ring> operator- () const {
			SparseMultivariatePolynomial<Ring> res(var);
			std::copy(names, names+var+1, res.names);
			for (int i = 0; i < terms.size(); ++i) {
				MultivariateTerm<Ring> t;
				t.coef = -terms[i].coef;
				t.v = var;
				t.degs = new int[var];
				for (int j = 0; j < var; ++j)
					t.degs[j] = terms[i].degs[j];
				res.terms.push_back(t);
			}
			return res;
		}

		/**
		 * Overload operator *
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring> operator* (const SparseMultivariatePolynomial<Ring>& b) const {
			if (!terms.size() || !b.terms.size()) {
				SparseMultivariatePolynomial<Ring> r;
				r.var = var;
				r.names = new Symbol[var+1];
				std::copy(names, names+var+1, r.names);
				return r;
			}
                        if (isConstant()) { return (b * terms[0].coef); }
                        if (b.isConstant()) { return (*this * b.terms[0].coef); }
                        std::vector<int> xs;
                        bool isOrdered = isOrderedRing(b, xs);

                        if (!isOrdered) {
                                std::cout << "BPAS: error, trying to multiply between Ring[";
                                for (int i = 1; i <= var; ++i) {
                                        std::cout << names[i];
                                        if (i < var) { std::cout << ", "; }
                                }
                                std::cout << "] and Ring[";
                                for (int i = 1; i <= b.var; ++i) {
                                        std::cout << b.names[i];
                                        if (i < b.var) { std::cout << ", "; }
                                }
                                std::cout << "] from SparseMultivariatePolynomial<Ring>." << std::endl;
                                exit(1);
                        }

                        int v = xs.size() / 2;
                        SparseMultivariatePolynomial<Ring> r (v);
                        r.names[0] = names[0];
                        if (names[0] == "9") {
                                for (int i = 1; i <= v; ++i) {
                                        if (xs[2*i-2])
                                                r.names[i] = names[xs[2*i-2]];
                                        else
                                                r.names[i] = b.names[xs[2*i-1]];
                                }
                        }
			for (int i = 0; i < terms.size(); ++i)
				r.pomopo(terms[i], b, xs);
			return r;
		}
		
		/**
		 * Overload operator *=
		 *
		 * @param b: A multivariate rational polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator*= (const SparseMultivariatePolynomial<Ring>& b) {
			*this = *this * b;
			return *this;
		}
		
		/**
		 * Overload operator *
		 *
		 * @param e: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring> operator* (const Ring& e) const {
			SparseMultivariatePolynomial<Ring> r (*this);
			return (r *= e);
		}

		// inline SparseMultivariatePolynomial<Ring> operator* (const sfixn& e) const {
            // SparseMultivariatePolynomial<Ring> r (*this);
            // return (r *= e);
        // }

		/**
		 * Overload operator *=
		 *
		 * @param e: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator*= (const Ring& e) {
			if (!e.isZero() && !e.isOne()) {
				for (int i = 0; i < terms.size(); ++i)
					terms[i].coef *= e;
			}
			else if (e.isZero()) { terms.clear(); }
			return *this;
		}

		inline friend SparseMultivariatePolynomial<Ring> operator* (const Ring& c, const SparseMultivariatePolynomial<Ring>& p) {
            return (p * c);
        }

		/**
		 * Overload operator *=
		 *
		 * @param e: A machine-word constant
		 **/
		// inline SparseMultivariatePolynomial<Ring>& operator*= (const sfixn& e) {
		// 	if (e != 0 && e != 1) {
		// 		for (int i = 0; i < terms.size(); ++i)
		// 			terms[i].coef *= e;
		// 	}
		// 	else if (e == 0) { terms.clear(); }
		// 	return *this;
		// }

		// inline friend SparseMultivariatePolynomial<Ring> operator* (const sfixn& c, const SparseMultivariatePolynomial<Ring>& p) {
  //           return (p * c);
  //       }

		/**
		 * Overload operator /
		 *
		 * @param b: A multivariate polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring> operator/ (const SparseMultivariatePolynomial<Ring>& b) const {
			SparseMultivariatePolynomial<Ring> r (*this);
			return (r /= b);
		}

		/**
		 * Overload operator /=
		 *
		 * @param b: A multivariate polynomial
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator/= (const SparseMultivariatePolynomial<Ring>& b) {
			if (!b.terms.size()) {
				std::cout << "BPAS: error, dividend is zero from SparseMultivariatePolynomial<Ring>." << std::endl;
				exit(1);
			}
			if (b.isConstant()) { return (*this /= b.terms[0].coef); }
			if (isConstant()) {
				terms.clear();
				return *this;
			}

			bool isIt = 0;
			int s = 0;
			std::vector<int> xs;
			for (int i = 1; i <= b.var; ++i) {
				bool isFound = 0;
				for (int j = i; j <= var; ++j)
					if (b.names[i] == names[j]) {
						if (j > i + s) {
							for (int k = i; k < j; ++k) {
								xs.push_back(k);
								xs.push_back(0);
							}
							s = j - i;
						}
						xs.push_back(j);
						xs.push_back(i);
						isFound = 1;
						break;
					}
				if (!isFound) {
					isIt = 1;
					break;
				}
			}
			for (int i = xs.size()/2; i <= var; ++i) {
				xs.push_back(i);
				xs.push_back(0);
			}

			if (isIt) {
				std::cout << "BPAS: error, trying to exact divide between Ring[";
				for (int i = 1; i <= var; ++i) {
					std::cout << names[i];
					if (i < var)
						std::cout << ", ";
				}
				std::cout << "] and Ring[";
				for (int i = 1; i <= b.var; ++i) {
					std::cout << b.names[i];
					if (i < var)
						std::cout << ", ";
				}
				std::cout << "] from SparseMultivariatePolynomial<Ring>." << std::endl;
				exit(1);
			}

			SparseMultivariatePolynomial<Ring> rem(*this);
			terms.clear();

			MultivariateTerm<Ring> bt = b.terms[b.terms.size()-1];
			s = compareTermDegs(rem.terms[rem.terms.size()-1], bt, xs);
			while (s >= 0) {
				MultivariateTerm<Ring> at = rem.terms[rem.terms.size()-1];
				MultivariateTerm<Ring> lc, nlc;
				lc.coef = at.coef / bt.coef;
				nlc.coef = - lc.coef;
				lc.v = var;
				nlc.v = var;
				lc.degs = new int[var];
				nlc.degs = new int[var];
				for (int i = 0; i < var; ++i) {
					nlc.degs[i] = lc.degs[i] = at.degs[xs[2*i]-1] - bt.degs[xs[2*i+1]-1];
					if (lc.degs[i] < 0) {
						std::cout << "BPAS: error, not exact division from SparseMultivariatePolynomial<Ring>." << std::endl;
						exit(1);
					}
				}
				rem.pomopo(nlc, b, xs);
				terms.insert(terms.begin(), lc);
				if (rem.terms.size())
					s = compareTermDegs(rem.terms[rem.terms.size()-1], bt, xs);
				else { s = -1; }
			}
			if (!rem.isZero()) {
				std::cout << "BPAS: error, not exact division from SparseMultivariatePolynomial<Ring>." << std::endl;
				exit(1);
			}

			return *this;
		}

		/**
		 * Overload operator /
		 *
		 * @param e: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring> operator/ (const Ring& e) const {
			SparseMultivariatePolynomial<Ring> r (*this);
			return (r /= e);
		}

		/**
		 * Overload operator /=
		 *
		 * @param e: A constant
		 **/
		inline SparseMultivariatePolynomial<Ring>& operator/= (const Ring& e) {
			if (e.isZero()) {
				std::cout << "BPAS: error, dividend is zero." << std::endl;
				exit(1);
			}
			else if (!e.isOne()) {
				for (int i = 0; i < terms.size(); ++i)
					terms[i].coef /= e;
			}
			return *this;
		}

		inline friend SparseMultivariatePolynomial<Ring> operator/ (const Ring& c, const SparseMultivariatePolynomial<Ring>& p) {
			if (p.isZero()) {
				std::cout << "BPAS: error, dividend is zero from SparseMultivariatePolynomial<Ring>." << std::endl;
				exit(1);
			}

			SparseMultivariatePolynomial<Ring> q;
			q.var = p.var;
			q.names = new Symbol[q.var+1];
			std::copy(p.names, p.names+p.var+1, q.names);

			if (p.isConstant()) {
				q += c;
				return (q /= p.terms[0].coef);
			}
            else { return q; }
        }

		/**
		 * Set variable names
		 *
		 * @param xs: Variable names
		 **/
		inline void setRingVariables (const std::vector<Symbol>& xs) {
			int ns = xs.size();
			if (ns != var) {
				std::cerr << "BPAS ERROR: SMP<Ring> shrinking and expanding polynomial ring NOT YET IMPLEMENTED" << std::endl;
				return;
			}
			if (names[0] == "1") {
				names[0] = "9";
				for (int i = var, j = 0; i > 0 && j < ns; --i, ++j)
					names[i] = xs[j];
			}
			else {
				bool isReorder = 0;
				int* pos = new int[var];
				for (int i = 0; i < ns; ++i) {
					bool isIt = 1;
					for (int j = 0; j < var; ++j) {
						if (names[j+1] == xs[ns-i-1]) {
							pos[i] = j;
							isIt = 0;
							break;
						}
					}
					if (isIt) { pos[i] = ns - i - 1; }
					else { isReorder = 1; }
				}

				if (isReorder) {
					SparseMultivariatePolynomial<Ring> r;
					r.var = var;
					r.names = new Symbol[var+1];
					r.names[0] = "9";
					for (int i = var, j = 0; i > 0 && j < ns; --i, ++j)
						r.names[i] = xs[j];

					int* d = new int[var];
					for (int i = 0; i < terms.size(); ++i) {
						for (int j = 0; j < var; ++j)
							d[var-j-1] = terms[i].degs[pos[j]];
						r.setCoefficient(var, d, terms[i].coef);
					}
					delete [] d;

					*this = r;
				}
				else {
					for (int i = var, j = 0; i > 0 && j < ns; --i, ++j)
						names[i] = xs[j];
				}

				delete []  pos;
			}
		}
		
		/**
		 * Get variable names
		 *
		 * @param
		 **/
		inline std::vector<Symbol> ringVariables() const {
			std::vector<Symbol> xs;
			for (int i = var; i > 0; --i)
				xs.push_back(names[i]);
			return xs;
		}

		inline std::vector<Symbol> variables() const {
			std::cerr << "BPAS ERROR: SDMP::variables() NOT YET IMPLEMENTED" << std::endl;
			return ringVariables();
		}
		
		/**
		 * Get the coefficient of one term
		 *
		 * @param d: The exponet of each variable
		 * @param v: Number of variables
		 **/
		inline Ring coefficient(int v, const int* d) const {
			int n = terms.size();
			for (int i = 0; i < n; ++i) {
				bool isIt = 1;
				for (int j = var-1, k = 0; j > -1 && k < v; --j, ++k) {
					if (terms[i].degs[j] != d[k]) {
						isIt = 0;
						break;
					}
				}
				for (int j = v-var-1; j > -1; --j) {
					if (d[j] != 0) {
						isIt = 0;
						break;
					}
				}
				if (isIt)
					return terms[i].coef;
			}
			return Ring();
		}

		inline Ring coefficient(const std::vector<int>& v) const {
			return coefficient(v.size(), v.data());
		}

		/**
		 * Set the coefficients of one term
		 *
		 * @param v: Number of variables
		 * @param d: Its exponent of each variable
		 * @param val: Coefficient
		 **/
		inline void setCoefficient(int v, const int* d, const Ring& val) {
			if (v != var) {
				std::cout << "BPAS: error, SMQP(" << var << "), but trying to setCoefficient with " << v << " variables." << std::endl;
				exit(1);
			}
			MultivariateTerm<Ring> a;
			a.coef = val;
			a.v = var;
			a.degs = new int[var];
			for (int i = 0; i < var; ++i)
				a.degs[i] = d[v-i-1];

			int i = 0;
			for (; i < terms.size(); ++i) {
				int k = compareTermDegs(a, terms[i]);
				if (!k) {
					if (val.isZero())
						terms.erase(terms.begin()+i);
					else
						terms[i].coef = val;
					break;
				}
				else if (k < 0 && !val.isZero()) {
					terms.insert(terms.begin()+i, a);
					break;
				}
			}
			if (i == terms.size() && !val.isZero())
				terms.push_back(a);
		}

		inline void setCoefficient(const std::vector<int>& v, const Ring& val) {
			setCoefficient(v.size(), v.data(), val);
		}

		//amha's function
		inline void addCoefficient(MultivariateTerm<Ring> a, int i) {
		    int k = 1;
		    for (; i < terms.size(); ++i) {
			k = compareTermDegs(a, terms[i]);
			if (!k) {
			    terms[i].coef += a.coef;
			    if (terms[i].coef == 0)
				terms.erase(terms.begin()+i);
			    break;
			}
			else if (k < 0 && a.coef != 0) {
			    terms.insert(terms.begin()+i, a);
			    break;
			}
		    }
		    if (k && i == terms.size() && a.coef != 0)
			terms.push_back(a);
		}

		/**
		 * Negate all the coefficients
		 *
		 * @param
		 **/
		inline void negate() {
			for (int i = 0; i < terms.size(); ++i) {
				terms[i].coef = -terms[i].coef;
			}
		}

		/**
		 * Convert to SUP<SMQP>
		 *
		 * @param x: The variable name
		 **/
		inline SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> > convertToSUP(const Symbol& x) const {
			SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> > r;
			r.setVariableName(x);

			int k = 0;
			for (int i = 1; i <= var; ++i)
				if (names[i] == x) { k = i; break; }

			if (k) {
				int d = 0;
				for (int i = terms.size()-1; i > -1; --i)
					if (terms[i].degs[k-1] > d) { d = terms[i].degs[k-1]; }
				int v = var - 1;
				SparseMultivariatePolynomial<Ring>* t = new SparseMultivariatePolynomial<Ring>[d+1];
				for (int i = 0; i <= d; ++i) {
					t[i].var = v;
					delete [] t[i].names;
					t[i].names = new Symbol[v+1];
					for (int j = 0; j < var; ++j) {
						if (j < k) { t[i].names[j] = names[j]; }
						else { t[i].names[j] = names[j+1]; }
					}
				}
				k--;
				for (int i = 0; i < terms.size(); ++i) {
					MultivariateTerm<Ring> a;
					a.v = v;
					a.coef = terms[i].coef;
					a.degs = new int[v];
					for (int j = 0; j < v; ++j) {
						if (j < k) { a.degs[j] = terms[i].degs[j]; }
						else { a.degs[j] = terms[i].degs[j+1]; }
					}
					t[terms[i].degs[k]].terms.push_back(a);
				}
				for (int i = 0; i <= d; ++i)
					r.setCoefficient(i, t[i]);
				delete [] t;
			}
			else { r.setCoefficient(0, *this); }

			return r;
		}

		/**
 		 * Overload stream operator <<
 		 *
 		 * @param out: Stream object
 		 * @param b: A multivariate rational polynoial
 		 **/
 		inline void print(std::ostream& out) const {
 			int n = this->numberOfTerms().get_si();
			int var = this->numberOfVariables();
			if (!n) { out << "0"; }
			for (int i = 0; i < n; i++) {
				if (i && this->terms[i].coef.isConstant() >= 0)
					out << "+";
				else if (this->terms[i].coef.isNegativeOne())
					out << "-";
				if (!this->terms[i].coef.isOne() && !this->terms[i].coef.isNegativeOne())
					out << this->terms[i].coef;
				bool isIt = 1;
				int* d = this->terms[i].degs;
				for (int j = 0; j < this->terms[i].v; ++j) {
					if (d[j]) {
						if ((!this->terms[i].coef.isOne() && !this->terms[i].coef.isNegativeOne() && isIt) || !isIt)
							out << "*";
						out << this->names[j+1];
						if (d[j] > 1)
							out << "^" << d[j];
						isIt = 0;
					}
				}
				if (isIt && (this->terms[i].coef.isOne() || this->terms[i].coef.isNegativeOne()))
					out << "1";
			}
		}

		inline ExpressionTree convertToExpressionTree() const {
			std::cerr << "SMP::convertToExpressionTree() NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return ExpressionTree();
		}



//***************************************************************multiplication section added by amha**************************************************

		inline std::vector<int> gammaMatrixByIndex(SparseMultivariatePolynomial<Ring> &b, int i, int j)
		{
			if(i>=terms.size() || j>=b.terms.size() || i<0 || j<0)
			{
				std::cout << "Gamma Matrix dimenstion out of bound."<< std::endl;
				exit(EXIT_FAILURE);
			}

			std::vector<int> tempGammaResult(0);
			for(unsigned int k=0; k<terms[0].v; k++)
			{
				tempGammaResult.push_back(terms[i].degs[k]+b.terms[j].degs[k]);
			}
			return tempGammaResult;
		}

		//potential bug, l will not function if it is greater than the number of term in the polynomial
		inline int JZeroI(int i, int na, int nb, int l)
		{
			return 1+((i/(na/l)) % 2)*std::floor(nb/(2*l));
		}

		inline int NStar(int l)
		{
		    return std::pow((l+1), 2);
		}

		//used by sort2DVector (if c++11 and above this function can be embbeded in sort algorithm as a lamda )
		static bool sortCompareHelper(const std::vector<int> &a, const std::vector<int> &b)
		{
			for(int i=a.size()-1; i>=0; i--)
			{
				if(a[i]!=b[i])
					return a[i]<b[i];
			}
			return false;
		}

		void sort2DVector(std::vector<std::vector<int> > &v)
		{
			std::sort(v.begin(), v.end(), sortCompareHelper);
		}

		inline std::vector<std::vector<int> > SList(SparseMultivariatePolynomial<Ring> &b, int l)
		{
			std::vector<std::vector<int> > S;
			std::vector<int> INF(0);
			//int infinity = std::numeric_limits<int>::max();
			int i, j;

			i=1;
			j=JZeroI(i, terms.size(), b.terms.size(), l);
			while(i<terms.size()-1 && j<b.terms.size()-1)
			{
				S.push_back(gammaMatrixByIndex(b, i, j));
				i=i+std::floor((terms.size())/l);
				j=j+std::floor((b.terms.size())/l);
			}

			j=1;
			for(; j<b.terms.size(); j+=std::floor((b.terms.size())/l))
			{
				S.push_back(gammaMatrixByIndex(b, terms.size()-1, j));
			}

			i=1;
			for(; i<terms.size(); i+=std::floor((terms.size())/l))
			{
				S.push_back(gammaMatrixByIndex(b, i, b.terms.size()-1));
			}

			S.push_back(gammaMatrixByIndex(b, 0, 0));
			sort2DVector(S);
			S.erase(std::unique(S.begin(), S.end()), S.end());
			for(unsigned int i=0; i<S[0].size(); i++)
			{
				INF.push_back(std::numeric_limits<int>::max());
			}
			S.push_back(INF);

			return S;
		}

		inline int compareList(std::vector<int> S, std::vector<int> gammaij)
		{
			std::vector<std::vector<int> > tempList;
			tempList.push_back(S);
			tempList.push_back(gammaij);
			sort2DVector(tempList);

			if(S==gammaij)
				return 0;
			else if(S==tempList[0])
				return -1;
			else
				return 1;
		}

		inline LminAndLmax findEdge(SparseMultivariatePolynomial<Ring> &b, std::vector<int> Sk, std::vector<int> Sk1)
		{
			LminAndLmax ret(terms.size());
			ret.Lmax.resize(terms.size());
			ret.Lmin.resize(terms.size());
			std::vector<int> gammaij;
			int i=0, j, k;

			for(;i<terms.size(); i++)
			{
				j=b.terms.size()-1;
				while(j>=0)
				{
					gammaij=gammaMatrixByIndex(b, i, j);
					if((compareList(gammaij, Sk1)==-1) && ((compareList(Sk, gammaij)==-1) || (compareList(Sk, gammaij)==0)))
					{
						ret.Lmax[i] = std::make_pair(i, j);

						while(j>=0)
						{
						    gammaij = gammaMatrixByIndex(b, i, j);
						    if((compareList(Sk, gammaij)==-1) || (compareList(Sk, gammaij)==0))
						    {
							ret.Lmin[i] = std::make_pair(i, j);
						    }
						    j--;
						}
						break;
					}
					else
					{
						j--;
					}
				}
			}
			gammaij.clear(); //new
			return ret;
		}
		//remove duplicate results with the same coefficients, by adding the coefficient.
		inline void uniqueExpAndCoef(std::vector<std::vector<int> > &D, std::vector<lfixq> &C)
		{
			for(int i=0; i<D.size(); i++)
			{
			    for(int j=i+1; j<D.size(); j++)
			    {
				if(compareList(D[i], D[j])==0)
				{
				    C[i]=C[i]+C[j];
				    C.erase(C.begin()+(j));
				    D.erase(D.begin()+(j));
				}
			    }
			    if(C[i]==0)
			    {
				C.erase(C.begin()+(i));
				D.erase(D.begin()+(i));
			    }
			}
		}
		/*
		inline SparseMultivariatePolynomial<Ring>& operator= (std::string &b) {

			mPoly = b;
			std::string valPoly;
			PolyAndExp pe = SplitPolyInBracket();
			valPoly = pe.poly;
			int e = pe.exp;

			std::vector<std::string> res = SplitPolyToTerms(valPoly, "+-");
			std::vector<std::string> allres = AllPolyVaraibles(res);
			std::vector<std::vector<int> > exp = AllPolyTermsExponents(allres, res);
			std::vector<std::vector<int> > coef = AllPolyTermsCoeffs(allres, res);

			SparseMultivariatePolynomial<Ring> t(allres.size());
			t.setVariableNames(allres);

			for (int i = 0; i < res.size(); ++i)
			{
				mpq_class elem = coef[i][0];
				int* d = new int[allres.size()];
				for (int j = 0; j < allres.size(); ++j)
				{
				        d[j] = exp[i][j];
				}
				t.setCoefficient(allres.size(), d, elem);
			}
			if(e==2)
			{
				*this *= t;
				return *this;
			}
			if(e >2)
			{
				*this = t;
				int ee = e -1 ;
				while(ee)
				{
					if(ee % 2) { *this *= t; }
					t *= t;
					ee >>=1;
				}
				return *this;
			}
			*this = t;
			return *this;
		}
		*/
		//parallel multiplication cilk_for
		inline SparseMultivariatePolynomial<Ring> sparsePolyMultiplication(const SparseMultivariatePolynomial<Ring>& b, int l) const
		{
		//	int l = (terms.size() < b.terms.size())? terms.size() : b.terms.size();
		//	if(l>12)
		//		l = 12;


			std::vector< std::vector<int> > S = SList(b, l);
			int nsize = S.size() - 1;
			SparseMultivariatePolynomial<Ring> *forcilk_res = new SparseMultivariatePolynomial<Ring>[nsize];


			cilk_for(int k = 0; k < nsize; k++)
			{
				LminAndLmax minMax = findEdge(b, S[k], S[k+1]);
				// merge all calculated exponents and coefficients to D and C repectively.

				for(int i=0; i<minMax.Lmin.size(); i++)
				{
					int pos = 0;
					for(int j=minMax.Lmin[i].second; j<=minMax.Lmax[i].second; j++)
					{
						if((minMax.Lmin[i].first != -1) && (minMax.Lmax[i].first != -1))
						{
							std::vector<int> tmp = gammaMatrixByIndex(b, i, j);
							MultivariateTerm<Ring> a;
							a.coef = terms[i].coef*b.terms[j].coef;
							a.v = var;
							a.degs = new int[var];
							std::copy(tmp.begin(), tmp.end(), a.degs);
							forcilk_res[k].addCoefficient(a, pos++);
						}
					}
				}
			}

			//cleaning all the duplicate result, this section is major performance bottle neck.
			SparseMultivariatePolynomial<Ring> res(var);
			std::copy(names, names+var+1, res.names);

			for(int k = 0; k < nsize; ++k)
			{
				std::copy(forcilk_res[k].terms.begin(), forcilk_res[k].terms.end(), std::back_inserter(res.terms));
			}

			delete [] forcilk_res;
			return res;
		}


		inline SparseMultivariatePolynomial<Ring> sparsePolyMultiplicationSerial(const SparseMultivariatePolynomial<Ring>& b, int l) const
		{
			std::vector< std::vector<int> > S = SList(b, l);
			int nsize = S.size() - 1;
			SparseMultivariatePolynomial<Ring> *forcilk_res = new SparseMultivariatePolynomial<Ring>[nsize];


			for(int k = 0; k < nsize; k++)
			{
				LminAndLmax minMax = findEdge(b, S[k], S[k+1]);
				// merge all calculated exponents and coefficients to D and C repectively.

				for(int i=0; i<minMax.Lmin.size(); i++)
				{
				        int pos = 0;
				        for(int j=minMax.Lmin[i].second; j<=minMax.Lmax[i].second; j++)
				        {
				                if((minMax.Lmin[i].first != -1) && (minMax.Lmax[i].first != -1))
				                {
				                        std::vector<int> tmp = gammaMatrixByIndex(b, i, j);
				                        MultivariateTerm<Ring> a;
				                        a.coef = terms[i].coef*b.terms[j].coef;
				                        a.v = var;
				                        a.degs = new int[var];
				                        std::copy(tmp.begin(), tmp.end(), a.degs);
				                        forcilk_res[k].addCoefficient(a, pos++);
				                }
				        }
				}
			}

			//cleaning all the duplicate result, this section is major performance bottle neck.
			SparseMultivariatePolynomial<Ring> res(var);
			std::copy(names, names+var+1, res.names);

			for(int k = 0; k < nsize; ++k)
			{
				std::copy(forcilk_res[k].terms.begin(), forcilk_res[k].terms.end(), std::back_inserter(res.terms));
			}

			delete [] forcilk_res;
			return res;
		}

		inline std::vector< SparseMultivariatePolynomial<Ring> > factors(){
		  std::vector< SparseMultivariatePolynomial<Ring> > facs;
		  //TODO 
		  std::cerr << "SparseMultivariatePolynomial<Ring>::factors NOT YET IMPLEMENTED" << std::endl;
		  // to be implemented
		  return facs;
		}

		inline SparseMultivariatePolynomial<Ring> primitivePart() const {
			std::cerr << "SparseMultivariatePolynomial<Ring>::primitivePart NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		inline Ring content() const {
			std::cerr << "SparseMultivariatePolynomial<Ring>::content() NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return leadingCoefficient();
		}

		inline SparseMultivariatePolynomial<Ring> content(const Symbol& sym) const {
			SparseMultivariatePolynomial<Ring> c;

		  //Ring c;
			if(this->isConstant()){
				if (terms.size() == 0) {
					return Ring(0);
				}
				return terms[0].coef;
			}
			SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> > sup = this->convertToSUP(sym);


		  //std::cout<< "In the content....." <<std::endl;
		  //std::cout<< "The Poly is :   "<<*this <<std::endl;

			int n = sup.degree().get_si();
			if (n>=0) {

				c = sup.coefficient(0);

				for (int i = 1; i <= n; ++i) {

					c = c.gcd(sup.coefficient(i));
					if (c.isOne())
						break;
				}
			}
			return c;
		}

		// TODO
		// inline void primitiveFactorization(SparseMultivariatePolynomial<Ring> &contentPart, SparseMultivariatePolynomial<Ring> &primPart){
		// 	primPart = *this;
		// 	contentPart =  this->content();
		// 	primPart/=contentPart;
		// }

		inline SparseMultivariatePolynomial<Ring> gcd(const SparseMultivariatePolynomial<Ring>& poly) const {

			SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > b = poly.convertToSUP(poly.leadingVariable());

			SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > a = this->convertToSUP(this->leadingVariable());

			if (this->isZero()) { return poly; }
			if (b.isZero()) { return *this; }

			if (a.variable() != b.variable()) {
				std::cout << "BPAS: error, trying to compute GCD between Ring[" << a.variable() << "] and Ring[" << b.variable() << "]." << std::endl;
				exit(1);
			}


			if (a.degree() == 0 || b.degree() == 0) {
				a.one();
				return a;
			}

			if(this->variables().size()<2){

				SparseUnivariatePolynomial<Ring> p1;
				for(int i=0; i<=a.degree();i++){
					p1.setCoefficient(i,a.coefficient(i).leadingCoefficient());
		      //std::cout<< "The required poly is: "<< p1<<std::endl;
				}
				SparseUnivariatePolynomial<Ring> p2;
				for(int i=0; i<=b.degree();i++){
					p2.setCoefficient(i,b.coefficient(i).leadingCoefficient());
		      //std::cout<< "The required poly is: "<< p2<<std::endl;
				}

				SparseUnivariatePolynomial<Ring> g = p1.gcd(p2);
		    //std::cout << "The gcd is: "<< g<<std::endl;
				if(g.degree()==0){
					SparseMultivariatePolynomial<Ring> temp(0);
					temp.setCoefficient(1,0, Ring(g.coefficient(0)));
				}else{
					SparseMultivariatePolynomial<Ring> temp(1);
					temp.setRingVariables(this->ringVariables());
					for(int j=0; j<=g.degree();j++){
						int d [1] = {j};
						temp.setCoefficient(1,d, Ring(g.coefficient(j)));
					}
					return temp;
				}
		    //SparseMultivariatePolynomial<Ring> temp(g);

			}




			SparseUnivariatePolynomial< SparseMultivariatePolynomial<Ring> > r;
			SparseMultivariatePolynomial<Ring> ca, cb, cr;
			ca = this->content(this->leadingVariable());
		    //std::cout<< "ca: "<<ca<<std::endl;
			a /=ca.convertToSUP(this->leadingVariable());
		    //std::cout<< "a: "<<a<<std::endl;
			cb = poly.content(poly.leadingVariable());
		    //std::cout<< "cb: "<<cb<<std::endl;
			b /= cb.convertToSUP(poly.leadingVariable());
		    //std::cout<< "b: "<<b<<std::endl;
			std::vector< SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > > R = b.subresultantChain(a);
		    //std::cout << "the size of src:  "<< R.size() << std::endl;

			r.setCoefficient(0, ca.gcd(cb));

		    //r *= cb;
			int n = R.size();
			bool isZero = 0;
			if (n) {

				isZero = 1;
				for (int i = 0; i < n; ++i) {
					if (!R[i].isZero()) {
						SparseMultivariatePolynomial<Ring> temp(R[i]);
						cr = temp.content(temp.leadingVariable());
			  //std::cout<<" the dividor is: " << temp<<std::endl;
			  //std::cout<<" the dividen is: " << cr<<std::endl;
						R[i]/=cr.convertToSUP(R[i].variable());

						SparseMultivariatePolynomial<Ring> temp2(R[i]);

						r *= temp2.convertToSUP(temp2.leadingVariable());
						isZero = 0;
						break;
					}
				}
			}
			if (isZero) {
				if (a.degree() <= b.degree()) { r *= a; }
				else { r *= b; }
			}
			SparseMultivariatePolynomial<Ring> temp2(r);
			return temp2;
		}

		inline Factors<SparseMultivariatePolynomial<Ring>> squareFree() const {
			// if(this->leadingVariableDegree()<=1){
			// 	return *this;
			// }

			// SparseUnivariatePolynomial<SparseMultivariatePolynomial<Ring> > diffTemp = this->convertToSUP(this->leadingVariable());
			// diffTemp.differentiate(1);

			// SparseMultivariatePolynomial<Ring> diff(diffTemp);
			// SparseMultivariatePolynomial<Ring> g = this->gcd(diff);
		 //  //std::cout<< "g:"<<g << std::endl;
			// SparseMultivariatePolynomial<Ring> sqf = *this;
			// sqf /=g;
			// return sqf;

			std::cerr << "SparseMultivariatePolynomial::squareFree NOT YET IMPLEMENTED" << std::endl;
			return Factors<SparseMultivariatePolynomial<Ring>>(*this);

		}

		inline void differentiate(const Symbol& s,int k) {
			std::cerr << "SparseMultivariatePolynomial::differentiate NOT YET IMPLEMENTED" << std::endl;
			//TODO
		}

		inline void differentiate(const Symbol& sym) {
			return this->differentiate(sym, 1);
		}

		inline SparseMultivariatePolynomial<Ring> derivative(const Symbol& sym, int k) const {
			SparseMultivariatePolynomial<Ring> ret(*this);
			ret.differentiate(sym, k);
			return ret;
		}

		inline SparseMultivariatePolynomial<Ring> derivative(const Symbol& sym) const {
			SparseMultivariatePolynomial<Ring> ret(*this);
			ret.differentiate(sym,1);
			return ret;
		}

		inline SparseMultivariatePolynomial<Ring> evaluate(const std::vector<Symbol>& syms, const std::vector<Ring>& vals) const {
			std::cerr << "SparseMultivariatePolynomial::evaluate NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		inline SparseMultivariatePolynomial<Ring> evaluate(int n, const Symbol* syms, const Ring* vals) const {
			std::vector<Symbol> vecSym;
			std::vector<Ring> vecVals;
			vecSym.resize(n);
			vecVals.resize(n);
			for (int i = 0; i < n; ++i) {
				vecSym[i] = syms[i];
				vecVals[i] = vals[i];
			}
			return evaluate(vecSym, vecVals);
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


