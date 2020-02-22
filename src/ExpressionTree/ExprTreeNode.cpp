
#include "ExpressionTree/ExprTreeNode.hpp"
#include "Symbol/Symbol.hpp"
#include <sstream>
#include <iostream>


ExprTreeNode::ExprTreeNode() : left(NULL), right(NULL), type(EXPR_NONE) {
	val.genericData = NULL;
}

ExprTreeNode::ExprTreeNode(ExprTreeType tType, ExprTreeVal* vVal, ExprTreeNode* lNode, ExprTreeNode* rNode) :
	type(tType), left(lNode), right(rNode) {
	if (vVal != NULL) {
		this->val = *vVal;
	} else {
		this->val.genericData = NULL;
	}
} 


ExprTreeNode::ExprTreeNode(const std::string& s) : left(NULL), right(NULL), type(EXPR_VAR) {
	val.var = new std::string(s);
}

ExprTreeNode::ExprTreeNode(long int i) : left(NULL), right(NULL), type(EXPR_INT) {
	val.i = i;
}

ExprTreeNode::ExprTreeNode(const mpz_class& z) : left(NULL), right(NULL), type(EXPR_MPZ) {
	val.z = new mpz_class(z);
}

ExprTreeNode::ExprTreeNode(const mpq_class& q) : left(NULL), right(NULL), type(EXPR_MPQ) {
	val.q = new mpq_class(q);
}

ExprTreeNode::ExprTreeNode(const Symbol& s) : left(NULL), right(NULL), type(EXPR_SYM) {
	val.sym = new Symbol(s);
}

ExprTreeNode::~ExprTreeNode() {
	switch (type) {
		case EXPR_VAR: {
			delete val.var;
			break;
		}
		case EXPR_INT: {
			break;
		}
		case EXPR_MPZ: {
			delete val.z;
			break;
		}
		case EXPR_MPQ: {
			delete val.q;
			break;
		}
		case EXPR_SYM: {
			delete val.sym;
			break;
		}
	}

	delete left;
	delete right;
}


bool ExprTreeNode::findChildType(ExprTreeType searchType) {
	if (type == searchType) {
		return true;
	}
	bool found = 0;
	if (left != NULL) {
		found = left->findChildType(searchType);
	}
	if (!found && right != NULL) {
		found = right->findChildType(searchType);
	}
	return found;
}

std::string ExprTreeNode::toString() const {
	std::string s;
	switch (type) {
		case EXPR_NONE: {
			s = "";
			break;
		}
		case EXPR_VAR: {
			s = *(val.var);
			break;
		}
		case EXPR_INT: {
			s = std::to_string(val.i);
			break;
		}
		case EXPR_MPZ: {
			s = (val.z)->get_str();	
			break;
		}
		case EXPR_MPQ: {
			s = (val.q)->get_str();
			break;
		}
		case EXPR_ADD: {
			s = " + ";
			break;
		}
		case EXPR_SUB: {
			s = " - ";
			break;
		}
		case EXPR_MULT: {
			s = " * ";
			break;
		}
		case EXPR_DIV: {
			s = " / ";
			break;
		}
		case EXPR_EXP: {
			s = "^";
			break;
		}
		case EXPR_NEG: {
			s = " -";
			break;
		}
		case EXPR_SYM: {
			Symbol tempSym = *(val.sym);
			std::string tempS = tempSym.toString();
			s = (val.sym)->toString();
			break;
		}
		case EXPR_ARRAY: {
			break;
		}
		case EXPR_DATA: {
			std::stringstream ss;
			ss << val.genericData;
			s = ss.str();
			break;
		}
	}
	return s;
}

bool ExprTreeNode::isVar() const {
	if (type == EXPR_VAR) {
		return true;
	}
	if (type == EXPR_EXP 
			&& left != NULL && left->type == EXPR_VAR
			&& right != NULL && right->type == EXPR_INT) {
		return true;
	}

	return false;
}

bool ExprTreeNode::isConstant() const {
	return (type == EXPR_INT || type == EXPR_MPZ || type == EXPR_MPQ);
}

bool ExprTreeNode::isPolynomialTerm() const {
	if (type == EXPR_MULT) {
		if (left == NULL || right == NULL) {
			return false;
		}
		if (left->isConstant()) {
			return right->isPolynomialTerm();
		} else {
			return left->isPolynomialTerm() && right->isPolynomialTerm();
		}
	} else {
		return (isVar() || isConstant());
	}
}

ExprTreeNode* ExprTreeNode::deepCopy() const {
	ExprTreeNode* r = new ExprTreeNode();
	r->type = type;
	switch (r->type) {
		case EXPR_VAR: {
			r->val.var = new std::string(*(val.var));
			break;
		}
		case EXPR_INT: {
			r->val.i = val.i;
			break;
		}
		case EXPR_MPZ: {
			r->val.z = new mpz_class(*(val.z));
			break;
		}
		case EXPR_MPQ: {
			r->val.q = new mpq_class(*(val.q));
			break;
		}
		case EXPR_SYM: {
			r->val.sym = new Symbol(*(val.sym));
			break;
		}
		case EXPR_NONE:
		case EXPR_DATA: {
			r->val.genericData = val.genericData;
			break;
		}
	}

	if (left != NULL) {
		r->left = left->deepCopy();
	}
	if (right != NULL) {
		r->right = right->deepCopy();
	}

	return r;
}



///////////////////
// Static Helpers
///////////////////

ExprTreeNode* ExprTreeNode::combineExprTreeNodes(ExprTreeNode* lNode, ExprTreeNode* rNode, ExprTreeType tType, ExprTreeVal* vVal) {
	ExprTreeNode* ret = new ExprTreeNode();
	ret->left = lNode;
	ret->right = rNode;
	ret->type = tType;
	if (vVal != NULL) {
		ret->val = *vVal;
	} else {
		ret->val.genericData = NULL;
	}

	return ret;
}
