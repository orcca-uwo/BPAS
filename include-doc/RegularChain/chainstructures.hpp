#ifndef _CHAINSTRUCTURES_H_
#define _CHAINSTRUCTURES_H_

#include <string>
#include "../ExpressionTree/ExprTreeNode.hpp"
#include "../ExpressionTree/ExpressionTree.hpp"


template <class RecursivePoly, class RegularChainType>
class PolyChainPair {
	public:
		RecursivePoly poly;
		RegularChainType chain;
		
		PolyChainPair () {
			poly = RecursivePoly();
			chain = RegularChainType();
		}
		
		PolyChainPair (const RecursivePoly& p, const RegularChainType& rc) {
			poly = RecursivePoly(p);
			chain = RegularChainType(rc);
		}
		
		PolyChainPair (const PolyChainPair<RecursivePoly,RegularChainType>& a) {
			poly = a.poly;
			chain = a.chain;
		}
		
		~PolyChainPair () {}
		
		inline PolyChainPair& operator= (const PolyChainPair<RecursivePoly,RegularChainType>& a) {
			poly = a.poly;
			chain = a.chain;
			return *this;
		}
		
//		friend bool operator==(const PolyChainPair& s1,const PolyChainPair& s2);
//		friend bool operator!=(const PolyChainPair& s1,const PolyChainPair& s2) {
//			return !(s1==s2);
//		}
		// need to extend ExpressionTree to cover PolyChainPairs
		/*ExpressionTree toExpressionTree() const {
			ExprTreeNode etn(*this);
			ExpressionTree et(&etn);
			return et;
		}*/
		std::string toString() const {
			std::stringstream ss;
			std::string out;
			ss << "<" << poly << "," << chain << ">";
			out = ss.str();
			return out;
		}

		inline friend std::ostream& operator<< (std::ostream &out, const PolyChainPair<RecursivePoly,RegularChainType>& a) {
			out << a.toString();
			return out;
		}

};

template <class RegularChainType>
class BoolChainPair {
	public:
		bool isTrue;
		RegularChainType chain;
		
		BoolChainPair () {
			isTrue = false;
			chain = RegularChainType();
		}
		
		BoolChainPair (bool b, const RegularChainType& rc) {
			isTrue = b;
			chain = RegularChainType(rc);
		}
		
		BoolChainPair (const BoolChainPair<RegularChainType>& a) {
			isTrue = a.isTrue;
			chain = a.chain;
		}
		
		~BoolChainPair () {}
		
		inline BoolChainPair& operator= (const BoolChainPair<RegularChainType>& a) {
			isTrue = a.isTrue;
			chain = a.chain;
			return *this;
		}
		
//		friend bool operator==(const BoolChainPair& s1,const BoolChainPair& s2);
//		friend bool operator!=(const BoolChainPair& s1,const BoolChainPair& s2) {
//			return !(s1==s2);
//		}
		// need to extend ExpressionTree to cover BoolChainPairs
		/*ExpressionTree toExpressionTree() const {
			ExprTreeNode etn(*this);
			ExpressionTree et(&etn);
			return et;
		}*/
		std::string toString() const {
			std::stringstream ss;
			std::string out;
			ss << "<";
			if (isTrue)
				ss << "true";
			else
				ss << "false";
			ss << "," << chain << ">";
			out = ss.str();
			return out;
		}

		inline friend std::ostream& operator<< (std::ostream &out, const BoolChainPair<RegularChainType>& a) {
			out << a.toString();
			return out;
		}

};

//template <class Field, class RecursivePoly>
//union DynamicChain {
//	RegularChain<Field,RecursivePoly>* rc;
//	ZeroDimensionalRegularChain<Field,RecursivePoly>* zdrc;
//};

//typedef enum ChainType{
//	RC_TYPE,
//	ZDRC_TYPE
//} ChainType;

//template <class Field, class RecursivePoly>
//class PolyChainPair {
//	private:
//		RecursivePoly poly;
//		DynamicChain<Field,RecursivePoly> chain;
//		ChainType type;

//	public:
//		PolyChainPair () : type(RC_TYPE) {
//			poly = RecursivePoly();
//			chain.rc = new RegularChain<Field,RecursivePoly>();
//		}
//		PolyChainPair (const RecursivePoly& p, const RegularChain<Field,RecursivePoly>& rc) : type(RC_TYPE) {
//			poly = RecursivePoly(p);
//			chain.rc = new RegularChain<Field,RecursivePoly>(rc);
//		}
//		PolyChainPair (const RecursivePoly& p, const ZeroDimensionalRegularChain<Field,RecursivePoly>& zdrc) : type(ZDRC_TYPE) {
//			poly = RecursivePoly(p);
//			chain.zdrc = new ZeroDimensionalRegularChain<Field,RecursivePoly>(zdrc);
//		}
//		PolyChainPair (const PolyChainPair& a) : type(a.type), poly(a.poly) {
//			switch (type) {
//				case RC_TYPE : {
//					chain.rc = new RegularChain<Field,RecursivePoly>(*a.chain.rc);
//					break;
//				}
//				case ZDRC_TYPE: {
//					chain.zdrc = new ZeroDimensionalRegularChain<Field,RecursivePoly>(*a.chain.zdrc);
//					break;
//				}
//			}
//		}
//		~PolyChainPair () {
//			switch (type) {
//				case RC_TYPE: {
//					delete chain.rc;
//					break;
//				}
//				case ZDRC_TYPE: {
//					delete chain.zdrc;
//					break;
//				}
//			}
//		}
//		
//		std::string getType() {
//			switch (type) {
//				case RC_TYPE: {
//					return "RC_TYPE";
//				}
//				case ZDRC_TYPE: {
//					return "ZDRC_TYPE";
//				}
//			}
//		}
//		inline PolyChainPair& operator= (const PolyChainPair<Field,RecursivePoly>& a) {
//			switch (a.type) {
//				case RC_TYPE: {
//					switch (type) {
//						case RC_TYPE: {
//							poly = a.poly;
//							*chain.rc = *a.chain.rc;
//							break;
//						}
//						case ZDRC_TYPE: {
//							delete chain.zdrc;
//							poly = a.poly;
//							chain.rc = new RegularChain<Field,RecursivePoly>(*a.chain.rc);
//							type = RC_TYPE;
//							break;
//						}
//					}
//					break;
//				}
//				case ZDRC_TYPE: {
//					switch (type) {
//						case RC_TYPE: {
//							delete chain.rc;
//							poly = a.poly;
//							chain.zdrc = new ZeroDimensionalRegularChain<Field,RecursivePoly>(*a.chain.zdrc);
//							type = ZDRC_TYPE;
//							break;
//						}
//						case ZDRC_TYPE: {
//							poly = a.poly;
//							*chain.zdrc = *a.chain.zdrc;
//							break;
//						}
//					}
//				}
//			}
//			return *this;
//		}
//		
////		friend bool operator==(const PolyChainPair& s1,const PolyChainPair& s2);
////		friend bool operator!=(const PolyChainPair& s1,const PolyChainPair& s2) {
////			return !(s1==s2);
////		}
//		// need to extend ExpressionTree to cover PolyChainPairs
//		/*ExpressionTree toExpressionTree() const {
//			ExprTreeNode etn(*this);
//			ExpressionTree et(&etn);
//			return et;
//		}*/
//		std::string toString() const {
//			std::stringstream ss;
//			std::string out;
//			switch (type) {
//				case RC_TYPE: {
//					ss << "<" << poly << "," << chain.rc << ">";
//					out = ss.str();
//					break;
//				}
//				case ZDRC_TYPE: {
//					ss << "<" << poly << "," << chain.zdrc << ">";
//					out = ss.str();
//					break;
//				}
//			}
//			return out;
//		}

//		inline friend std::ostream& operator<< (std::ostream &out, const PolyChainPair<Field,RecursivePoly>& b) {
//			out << b.toString();
//			return out;
//		}

//};

//template <class Field, class RecursivePoly>
//class BoolChainPair {
//	private:
//		bool isTrue;
//		DynamicChain chain;
//		ChainType type;

//	public:
//		BoolChainPair () : type(RC_TYPE) {
//			isTrue = false;
//			chain.rc = new RegularChain<Field,RecursivePoly>();
//		}
//		BoolChainPair (bool b, const RegularChain<Field,RecursivePoly>& rc) : type(RC_TYPE), isTrue(b) {
//			chain.rc = new RegularChain<Field,RecursivePoly>(rc);
//		}
//		BoolChainPair (bool b, const ZeroDimensionalRegularChain<Field,RecursivePoly>& zdrc) : type(ZDRC_TYPE), isTrue(b) {
//			chain.zdrc = new ZeroDimensionalRegularChain<Field,RecursivePoly>(zdrc);
//		}
//		BoolChainPair (const BoolChainPair& a) : type(a.type), isTrue(a.isTrue) {
//			switch (type) {
//				case RC_TYPE : {
//					chain.rc = new RegularChain<Field,RecursivePoly>(*a.chain.rc);
//					break;
//				}
//				case ZDRC_TYPE: {
//					chain.zdrc = new ZeroDimensionalRegularChain<Field,RecursivePoly>(*a.chain.zdrc);
//					break;
//				}
//			}
//		}
//		~BoolChainPair () {
//			switch (type) {
//				case RC_TYPE: {
//					delete chain.rc;
//					break;
//				}
//				case ZDRC_TYPE: {
//					delete chain.zdrc;
//					break;
//				}
//			}
//		}
//		
//		std::string getType() {
//			switch (type) {
//				case RC_TYPE: {
//					return "RC_TYPE";
//				}
//				case ZDRC_TYPE: {
//					return "ZDRC_TYPE";
//				}
//			}
//		}
//		inline BoolChainPair& operator= (const BoolChainPair<Field,RecursivePoly>& a) {
//			switch (a.type) {
//				case RC_TYPE: {
//					switch (type) {
//						case RC_TYPE: {
//							isTrue = a.isTrue;
//							*chain.rc = *a.chain.rc;
//							break;
//						}
//						case ZDRC_TYPE: {
//							delete chain.zdrc;
//							isTrue = a.isTrue;
//							chain.rc = new RegularChain<Field,RecursivePoly>(*a.chain.rc);
//							type = RC_TYPE;
//							break;
//						}
//					}
//					break;
//				}
//				case ZDRC_TYPE: {
//					switch (type) {
//						case RC_TYPE: {
//							delete chain.rc;
//							isTrue = a.isTrue;
//							chain.zdrc = new ZeroDimensionalRegularChain<Field,RecursivePoly>(*a.chain.zdrc);
//							type = ZDRC_TYPE;
//							break;
//						}
//						case ZDRC_TYPE: {
//							isTrue = a.isTrue;
//							*chain.zdrc = *a.chain.zdrc;
//							break;
//						}
//					}
//				}
//			}
//			return *this;
//		}
//		
////		friend bool operator==(const BoolChainPair& s1,const BoolChainPair& s2);
////		friend bool operator!=(const BoolChainPair& s1,const BoolChainPair& s2) {
////			return !(s1==s2);
////		}
//		// need to extend ExpressionTree to cover BoolChainPairs
//		/*ExpressionTree toExpressionTree() const {
//			ExprTreeNode etn(*this);
//			ExpressionTree et(&etn);
//			return et;
//		}*/
//		std::string toString() const {
//			std::stringstream ss;
//			std::string out;
//			switch (type) {
//				case RC_TYPE: {
//					ss << "<";
//					if (isTrue)
//						ss << "true";
//					else
//						ss << "false";
//					ss << "," << chain.rc << ">";
//					out = ss.str();
//					break;
//				}
//				case ZDRC_TYPE: {
//					ss << "<";
//					if (isTrue)
//						ss << "true";
//					else
//						ss << "false";
//					ss << "," << chain.zdrc << ">";
//					out = ss.str();
//					break;
//				}
//			}
//			return out;
//		}

//		inline friend std::ostream& operator<< (std::ostream &out, const BoolChainPair<Field,RecursivePoly>& b) {
//			out << b.toString();
//			return out;
//		}

//};

#endif
