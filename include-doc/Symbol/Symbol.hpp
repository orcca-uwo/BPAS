#ifndef _SYMBOL_H_
#define _SYMBOL_H_

#include<locale>
#include<codecvt>
#include<string>
#include "../ExpressionTree/ExprTreeNode.hpp"
#include "../ExpressionTree/ExpressionTree.hpp"

union SymbolString {
	std::string* str;	// ASCII string
	std::wstring* wstr;	// Unicode string
};

typedef enum SymbolType{
	NULL_TYPE = 0x0,
	STRING_TYPE,
	WSTRING_TYPE
} SymbolType;

extern std::wstring string2wstring(const std::string& str);
extern std::string wstring2string(const std::wstring& wstr);

// Do we need a Unicode processing library (such as ICU) to normalize strings,
//   i.e., can multiple unicode strings be equivalent?
// This is currently mixing hpp and cpp code
/** 
 * An encapsulation of a mathematical symbol. 
 * This symbol can be an indeterminant, or more generally, any unicode string.
 * Provides comparison, concatenation, and printing.
 */
class Symbol {
	private:
		SymbolString s;
		SymbolType type;

	public:
		Symbol () : type(NULL_TYPE) {
			s.str = new std::string("");
		}
		Symbol (char c) : type(STRING_TYPE) {
			char localC[2];
			localC[0] = c;
			localC[1] = '\0';
			s.str = new std::string(localC);
		}
		explicit Symbol (const char* c) : type(STRING_TYPE) {
			s.str = new std::string(c);
		}
		explicit Symbol (const wchar_t* c) : type(WSTRING_TYPE) {
			s.wstr = new std::wstring(c);
		}
		explicit Symbol (const std::string& a) : type(STRING_TYPE) {
			s.str = new std::string(a);
		}
		explicit Symbol (const std::wstring& a) : type(WSTRING_TYPE) {
			s.wstr = new std::wstring(a);
		}
		Symbol (const Symbol& a) : type(a.type) {
			switch (type) {
				case NULL_TYPE :
				case STRING_TYPE : {
					s.str = new std::string(*a.s.str);
					break;
				}
				case WSTRING_TYPE: {
					s.wstr = new std::wstring(*a.s.wstr);
					break;
				}
			}
		}
		~Symbol () {
			switch (type) {
				case NULL_TYPE:
				case STRING_TYPE: {
					delete s.str;
					break;
				}
				case WSTRING_TYPE: {
					delete s.wstr;
					break;
				}
			}
		}
		
		std::string getType() {
			switch (type) {
				case NULL_TYPE: {
					return "NULL_TYPE";
				}
				case STRING_TYPE: {
					return "STRING_TYPE";
				}
				case WSTRING_TYPE: {
					return "WSTRING_TYPE";
				}
			}
		}
		inline Symbol& operator= (char c) {
			char localC[2];
			localC[0] = c;
			localC[1] = '\0';
			return (*this = std::string(localC));
		}

		inline Symbol& operator= (const char* c) {
			std::string s(c);
			return (*this = s);
		}

		inline Symbol& operator= (const std::string& a) {
			switch (type) {
				case NULL_TYPE: {
					if (a != "") {
						type = STRING_TYPE;
						*s.str = a;	
					}
					break;
				}
				case STRING_TYPE: {
					if (a == "") {
						type = NULL_TYPE;
					}
					*s.str = a;
					break;
				}
				case WSTRING_TYPE: {
					delete s.wstr;
					s.str = new std::string(a);
					type = STRING_TYPE;
					break;
				}
			}
			return *this;
		}
		inline Symbol& operator= (const Symbol& a) {
			switch (a.type) {
				case NULL_TYPE: {
					switch (type) {
						case NULL_TYPE: {
							break;
						}
						case STRING_TYPE: {
							*s.str = *a.s.str;
							type = NULL_TYPE;
							break;
						}
						case WSTRING_TYPE: {
							delete s.wstr;
							s.str = new std::string("");
							type = NULL_TYPE;
							break;
						}
					}
					break;
				}
				case STRING_TYPE: {
					switch (type) {
						case NULL_TYPE: {
							*s.str = *a.s.str;
							type = STRING_TYPE;
							break;
						}
						case STRING_TYPE: {
							*s.str = *a.s.str;
							break;
						}
						case WSTRING_TYPE: {
							delete s.wstr;
							s.str = new std::string(*a.s.str);
							type = STRING_TYPE;
							break;
						}
					}
					break;
				}
				case WSTRING_TYPE: {
					switch (type) {
						case NULL_TYPE:
						case STRING_TYPE: {
							delete s.str;
							s.wstr = new std::wstring(*a.s.wstr);
							type = WSTRING_TYPE;
							break;
						}
						case WSTRING_TYPE: {
							*s.wstr = *a.s.wstr;
							break;
						}
					}
				}
			}
			return *this;
		}
		inline Symbol& operator+= (const std::string& a) {
			if (a != "") {
				switch (type) {
					case NULL_TYPE: {
						type = STRING_TYPE;
					}
					case STRING_TYPE: {
						*s.str += a;
						break;
					}
					case WSTRING_TYPE: {
						std::wstring b;
						b = string2wstring(a);
						*s.wstr += b;
						break;
					}
				}
			}
			return *this;
		}
		inline Symbol& operator+= (const Symbol& a) {
			switch (type) {
				case NULL_TYPE: {
					switch (a.type) {
						case NULL_TYPE: {
							break;
						}
						case STRING_TYPE: {
							*s.str += *a.s.str;
							type = STRING_TYPE;
							break;
						}
						case WSTRING_TYPE: {
							delete s.str;
							s.wstr = new std::wstring(*a.s.wstr);
							type = WSTRING_TYPE;
						}
					}
					break;
				}
				case STRING_TYPE: {
					switch (a.type) {
						case NULL_TYPE: {
							break;
						}
						case STRING_TYPE: {
							*s.str += *a.s.str;
							break;
						}
						case WSTRING_TYPE: {
							std::wstring b;
							b = string2wstring(*s.str);
							delete s.str;
							s.wstr = new std::wstring(b);
							*s.wstr += *a.s.wstr;
							type = WSTRING_TYPE;
							break;
						}
					}
					break;
				}
				case WSTRING_TYPE: {
					switch (a.type) {
						case NULL_TYPE: {
							break;
						}
						case STRING_TYPE: {
							std::wstring b;
							b = string2wstring(*a.s.str);
							*s.wstr += b;
							break;
						}
						case WSTRING_TYPE: {
							*s.wstr += *a.s.wstr;
						}
					}
				}
			}
			return *this;
		}
		
		friend bool operator<(const Symbol& s1,const Symbol& s2);
		friend bool operator==(const Symbol& s1,const Symbol& s2);
		friend bool operator== (const std::string& s2, const Symbol& s1);		
		friend bool operator== (const Symbol& s1, const std::string& s2);		
		friend bool operator!= (const std::string& s2, const Symbol& s1);		
		friend bool operator!= (const Symbol& s1, const std::string& s2);		
		//not needed as we can convert form string to symbol
		// friend bool operator==(const Symbol& s1, const std::string& s2){
		// 	return s1 == Symbol(s2);
		// }
		// friend bool operator==(const std::string& s1, const Symbol& s2) {
			// return Symbol(s1) == s2;
		// }
		friend bool operator<=(const Symbol& s1,const Symbol& s2) {
			return (s1<s2 || s1==s2);
		}
		friend bool operator>(const Symbol& s1,const Symbol& s2) {
			return !(s1<=s2);
		}
		friend bool operator>=(const Symbol& s1,const Symbol& s2) {
			return !(s1<s2);
		}
		friend bool operator!=(const Symbol& s1,const Symbol& s2) {
			return !(s1==s2);
		}
		
		static Symbol randomElement();
		static std::vector<Symbol> randomElements(int n);
		
		// need to extend ExpressionTree to cover Symbols
		/*ExpressionTree toExpressionTree() const {
			ExprTreeNode etn(*this);
			ExpressionTree et(&etn);
			return et;
		}*/
		std::string toString() const {
			std::string out;
			switch (type) {
				case NULL_TYPE:
				case STRING_TYPE: {
					out = *(s.str);
					break;
				}
				case WSTRING_TYPE: {
					out = wstring2string(*s.wstr);
					break;
				}
			}
			return out;
		}

		inline friend std::ostream& operator<< (std::ostream &out, const Symbol& b) {
			out << b.toString();
			return out;
		}
		
		inline ExpressionTree convertToExpressionTree() const {
			ExprTreeNode* etn = new ExprTreeNode(*this);
			ExpressionTree et(etn);
			return et;
		}

};

#endif
