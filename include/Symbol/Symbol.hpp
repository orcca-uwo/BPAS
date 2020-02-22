#ifndef _SYMBOL_H_
#define _SYMBOL_H_

#include <locale>
#include <codecvt>
#include <string>
#include "../ExpressionTree/ExprTreeNode.hpp"
#include "../ExpressionTree/ExpressionTree.hpp"


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
	/**
	 * Simple union for different 
	 */
	union SymbolString {
		std::string* str;	// ASCII string
		std::wstring* wstr;	// Unicode string
	};
	
	/** 
	 * Enum to differentiate between symbol types.
	 */
	typedef enum SymbolType{
		NULL_TYPE = 0x0,
		STRING_TYPE,
		WSTRING_TYPE
	} SymbolType;

	SymbolString s;
	SymbolType type;

	public:
		/**
		 * Construct an empty Symbol.
		 */
		Symbol () : type(NULL_TYPE) {
			s.str = new std::string("");
		}

		/**
		 * Construct a Symbol from the given character.
		 * @param c the character
		 */
		Symbol (char c) : type(STRING_TYPE) {
			char localC[2];
			localC[0] = c;
			localC[1] = '\0';
			s.str = new std::string(localC);
		}

		/** 
		 * Construct a Symbol from a c-string.
		 * @param c the c-string.
		 */
		explicit Symbol (const char* c) : type(STRING_TYPE) {
			s.str = new std::string(c);
		}

		/**
		 * Construct a Symbol from a wide character c-string.
		 */
		explicit Symbol (const wchar_t* c) : type(WSTRING_TYPE) {
			s.wstr = new std::wstring(c);
		}

		/**
		 * Construct a Symbol from a std::string.
		 * @param a the string
		 */
		explicit Symbol (const std::string& a) : type(STRING_TYPE) {
			s.str = new std::string(a);
		}

		/**
		 * Construct a Symbol from a wide string, std::wstring.
		 * @param a the wide string
		 */
		explicit Symbol (const std::wstring& a) : type(WSTRING_TYPE) {
			s.wstr = new std::wstring(a);
		}

		/**
		 * Copy constructor.
		 */
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

		/**
		 * Destructor.
		 */
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
		
		/**
		 * Get the type of this Symbol.
		 */
		// std::string getType() {
		// 	switch (type) {
		// 		case NULL_TYPE: {
		// 			return "NULL_TYPE";
		// 		}
		// 		case STRING_TYPE: {
		// 			return "STRING_TYPE";
		// 		}
		// 		case WSTRING_TYPE: {
		// 			return "WSTRING_TYPE";
		// 		}
		// 	}
		// }

		/**
		 * Assignment operator from a character.
		 *
		 * @param c the character to assign from.
		 */
		inline Symbol& operator= (char c) {
			char localC[2];
			localC[0] = c;
			localC[1] = '\0';
			return (*this = std::string(localC));
		}

		/**
		 * Assignment operator from a c-string.
		 *
		 * @param c the c-string to assign from.
		 */
		inline Symbol& operator= (const char* c) {
			std::string s(c);
			return (*this = s);
		}

		/**
		 * Assignment operator from a std::string.
		 *
		 * @param a the string to assign from.
		 */
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

		/**
		 * Copy assignment.
		 * @param a the other Symbol to copy from.
		 */
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

		/**
		 * Concatenate a string to this Symbol.
		 * @param a the string to concatenate.
		 */
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

		/**
		 * Concatenate a Symbol to this Symbol.
		 * @param a the Symbol to concatenate.
		 */
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
		
		/**
		 * Less than comparison operator between two Symbol.
		 * @param s1 the first Symbol
		 * @param s2 the second Symbol
		 * @return true iff s1 < s2
		 */
		friend bool operator<(const Symbol& s1,const Symbol& s2);

		/**
		 * Equality comparison operator between two Symbol.
		 * @param s1 the first Symbol
		 * @param s2 the second Symbol
		 * @return true iff s1 = s2
		 */
		friend bool operator==(const Symbol& s1,const Symbol& s2);

		/**
		 * Equality comparison operator between a Symbol and a string.
		 * @param s1 a Symbol
		 * @param s2 a string
		 * @return true iff s1 = s2
		 */
		friend bool operator== (const std::string& s2, const Symbol& s1);		

		/**
		 * Equality comparison operator between a Symbol and a string.
		 * @param s1 a Symbol
		 * @param s2 a string
		 * @return true iff s1 = s2
		 */
		friend bool operator== (const Symbol& s1, const std::string& s2);	

		/**
		 * Inequality comparison operator between a Symbol and a string.
		 * @param s1 a Symbol
		 * @param s2 a string
		 * @return true iff s1 != s2
		 */
		friend bool operator!= (const std::string& s2, const Symbol& s1);		

		/**
		 * Inequality comparison operator between a Symbol and a string.
		 * @param s1 a Symbol
		 * @param s2 a string
		 * @return true iff s1 != s2
		 */
		friend bool operator!= (const Symbol& s1, const std::string& s2);		
		
		//not needed as we can convert form string to symbol
		// friend bool operator==(const Symbol& s1, const std::string& s2){
		// 	return s1 == Symbol(s2);
		// }
		// friend bool operator==(const std::string& s1, const Symbol& s2) {
			// return Symbol(s1) == s2;
		// }
		

		/**
		 * Less-than-equal comparison operator between two Symbol.
		 * @param s1 the first Symbol
		 * @param s2 the second Symbol
		 * @return true iff s1 <= s2
		 */
		friend bool operator<=(const Symbol& s1,const Symbol& s2) {
			return (s1<s2 || s1==s2);
		}

		/**
		 * Greater-than comparison operator between two Symbol.
		 * @param s1 the first Symbol
		 * @param s2 the second Symbol
		 * @return true iff s1 > s2
		 */
		friend bool operator>(const Symbol& s1,const Symbol& s2) {
			return !(s1<=s2);
		}

		/**
		 * Greater-than-equal comparison operator between two Symbol.
		 * @param s1 the first Symbol
		 * @param s2 the second Symbol
		 * @return true iff s1 >= s2
		 */
		friend bool operator>=(const Symbol& s1,const Symbol& s2) {
			return !(s1<s2);
		}

		/**
		 * Inequality comparison operator between two Symbol.
		 * @param s1 the first Symbol
		 * @param s2 the second Symbol
		 * @return true iff s1 != s2
		 */
		friend bool operator!=(const Symbol& s1,const Symbol& s2) {
			return !(s1==s2);
		}
		
		/**
		 * Get a random Symbol.
		 * @return the random Symbol.
		 */
		static Symbol randomElement();

		/**
		 * Get a vector of pair-wise different random Symbols.
		 * n: the size of the vector.
		 * @return the vector of Symbols
		 */
		static std::vector<Symbol> randomElements(int n);
		
		/**
		 * Conver a Symbol to a std::string.
		 * @return the Symbol as a string.
		 */
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

		/**
		 * Output operator. Outputs the Symbol to the output stream.
		 * out: the output stream.
		 * b: the Symbol to output.
		 */
		inline friend std::ostream& operator<< (std::ostream &out, const Symbol& b) {
			out << b.toString();
			return out;
		}
		
		/**
		 * Convert a Symbol to an ExpressionTree.
		 * @return the ExpressionTree encoding this Symbol.
		 */
		inline ExpressionTree convertToExpressionTree() const {
			ExprTreeNode* etn = new ExprTreeNode(*this);
			ExpressionTree et(etn);
			return et;
		}

};

#endif
