#include "Symbol/Symbol.hpp"
#include <locale>
#include <codecvt>
#include <string>
#include <iostream>
#include "ExpressionTree/ExprTreeNode.hpp"
#include "ExpressionTree/ExpressionTree.hpp"
#include "Utils/RandomHelpers.hpp"

// Convert std::string to std::wstring
std::wstring string2wstring(const std::string& str) {
    using convert_typeX = std::codecvt_utf8<wchar_t>;
    std::wstring_convert<convert_typeX, wchar_t> converterX;

    return converterX.from_bytes(str);
}

// Convert std::wstring to std::string
std::string wstring2string(const std::wstring& wstr) {
    using convert_typeX = std::codecvt_utf8<wchar_t>;
    std::wstring_convert<convert_typeX, wchar_t> converterX;

    return converterX.to_bytes(wstr);
}

bool operator< (const Symbol& s1,const Symbol& s2) {
	switch (s1.type) {
		case Symbol::WSTRING_TYPE: {
			switch (s2.type) {
				case Symbol::WSTRING_TYPE: {
					return *(s1.s.wstr)<*(s2.s.wstr);
				}
				case Symbol::NULL_TYPE:
				case Symbol::STRING_TYPE: {
					return *(s1.s.wstr)<string2wstring(*(s2.s.str));
				}
			}
		}
		case Symbol::NULL_TYPE:
		case Symbol::STRING_TYPE: {
			switch (s2.type) {
				case Symbol::WSTRING_TYPE: {
					return string2wstring(*(s1.s.str))<*(s2.s.wstr);
				}
				case Symbol::NULL_TYPE:
				case Symbol::STRING_TYPE: {
					return *(s1.s.str)<*(s2.s.str);
				}
			}
		}
	}
	return 0;
}

bool operator== (const Symbol& s1,const Symbol& s2) {
	switch (s1.type) {
		case Symbol::WSTRING_TYPE: {
			switch (s2.type) {
				case Symbol::WSTRING_TYPE: {
					return *(s1.s.wstr)==*(s2.s.wstr);
				}
				case Symbol::NULL_TYPE:
				case Symbol::STRING_TYPE: {
					return *(s1.s.wstr)==string2wstring(*(s2.s.str));
				}
			}
		}
		case Symbol::NULL_TYPE:
		case Symbol::STRING_TYPE: {
			switch (s2.type) {
				case Symbol::WSTRING_TYPE: {
					return string2wstring(*(s1.s.str))==*(s2.s.wstr);
				}
				case Symbol::NULL_TYPE:
				case Symbol::STRING_TYPE: {
					return *(s1.s.str)==*(s2.s.str);
				}
			}
		}
	}
	return 0;
}

bool operator!= (const std::string& s2, const Symbol& s1) {
	switch (s1.type) {
		case Symbol::WSTRING_TYPE: {
			return *(s1.s.wstr)!=string2wstring(s2);
		}
		case Symbol::NULL_TYPE:
		case Symbol::STRING_TYPE: {
			return *(s1.s.str)!=(s2);
		}
	}
	return 0;
}

bool operator!= (const Symbol& s1,const std::string& s2) {
	switch (s1.type) {
		case Symbol::WSTRING_TYPE: {
			return *(s1.s.wstr)!=string2wstring(s2);
		}
		case Symbol::NULL_TYPE:
		case Symbol::STRING_TYPE: {
			return *(s1.s.str)!=(s2);
		}
	}
	return 0;
}

bool operator== (const std::string& s2, const Symbol& s1) {
	switch (s1.type) {
		case Symbol::WSTRING_TYPE: {
			return *(s1.s.wstr)==string2wstring(s2);
		}
		case Symbol::NULL_TYPE:
		case Symbol::STRING_TYPE: {
			return *(s1.s.str)==(s2);
		}
	}
	return 0;
}

bool operator== (const Symbol& s1,const std::string& s2) {
	switch (s1.type) {
		case Symbol::WSTRING_TYPE: {
			return *(s1.s.wstr)==string2wstring(s2);
		}
		case Symbol::NULL_TYPE:
		case Symbol::STRING_TYPE: {
			return *(s1.s.str)==(s2);
		}
	}
	return 0;
}

Symbol Symbol::randomElement() {
	// Generate n unicode string of length 'len' whose characters are in range [start, end]
	std::vector<int> vals;
	int index;
	for (int i=65; i<=90; ++i) {
		if (i == 73) {
			//skip I as it often encodes mathematical i.
			continue; 
		}
		vals.push_back(i);
	}
	for (int i=97; i<=122; ++i)
		vals.push_back(i);
	size_t len(2);
	size_t start(0);
	size_t end(vals.size()-1);
//	wchar_t* ustr = new wchar_t[len+1];      // +1 for '\0'
//	size_t intervalLength = end - start + 1; // +1 for inclusive range

//	srand(time(NULL));
//	for (auto i = 0; i < len; i++) {
//		ustr[i] = (rand() % intervalLength) + start;
//		std::cout << "ustr[" << i << "] = " << ustr[i] << std::endl;
//	}
//	ustr[len] = L'\0';
	char* ustr = new char[len+1];      // +1 for '\0'
	size_t intervalLength = end - start + 1; // +1 for inclusive range

//	srand(time(NULL));
	for (auto i = 0; i < len; i++) {
		index = (rand() % intervalLength) + start;
		ustr[i] = vals[index];
	}
	std::string str(ustr,len);
	return Symbol(str);
}

std::vector<Symbol> Symbol::randomElements(int n) {
	// Generate n unicode characters in a locally generated list //
	std::vector<int> vals,indices;
	int index;
	for (int i=65; i<=90; ++i) {
		if (i == 73) {
			//skip I as it often encodes mathematical i.
			continue; 
		}
		vals.push_back(i);
	}
	for (int i=97; i<=122; ++i)
		vals.push_back(i);
		
	if (n > vals.size()/2) {
		std::cerr << "BPAS: error, n must have a value less than " << vals.size()/2 << std::endl;
	}
	std::vector<Symbol> out;
	size_t size(2);
	char str[size];
	str[size-1] = '\0';

//	ustr[size] = L'\0';
	indices = randValsInRange(0,vals.size()-1,n);
	for (auto i=0; i<n; ++i) {
		str[0] = vals[indices[i]];
		out.push_back(Symbol(str));
	}
	return out;
//	size_t len(2);
//	size_t start(0);
//	size_t end(vals.size()-1);
////	wchar_t* ustr = new wchar_t[len+1];      // +1 for '\0'
////	size_t intervalLength = end - start + 1; // +1 for inclusive range

////	srand(time(NULL));
////	for (auto i = 0; i < len; i++) {
////		ustr[i] = (rand() % intervalLength) + start;
////		std::cout << "ustr[" << i << "] = " << ustr[i] << std::endl;
////	}
////	ustr[len] = L'\0';
//	char* ustr = new char[len+1];      // +1 for '\0'
//	size_t intervalLength = end - start + 1; // +1 for inclusive range

////	srand(time(NULL));
//	for (auto i = 0; i < len; i++) {
//		index = (rand() % intervalLength) + start;
//		ustr[i] = vals[index];
//	}
//	std::string str(ustr,len);
//	return Symbol(str);
}
