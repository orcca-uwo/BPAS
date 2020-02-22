
#ifndef _SMQP_CPP_SUPPORT_
#define _SMQP_CPP_SUPPORT_

#include <sstream>
#include <string>
#include <gmpxx.h>

#include "SMQP_Support.h"

typedef mpq_class ratNum_class;

/**
 * Convert a polynomial given its head node, node, to a string.
 * nvar: the number of variables in the polynomial
 * vars: the strings representing each variable. e.g. {"x1","x2","x3"}  
 */
static std::string polyToString(Node* node, int nvar, const char** vars) {
	
	if (node == NULL) {
		return "0";
	}

	std::stringstream ss;

	bool first = true;
	bool needsMult = false;
	bool isConst = true;
	ratNum_class coef;
	while (node != NULL) {
		coef = ratNum_class(node->coef);
		// coef = node->coef;
		isConst = true;
		if (coef < 0) {
			coef *= -1; 
			ss << " - ";
		} else if (!first) {
			ss << " + ";
		}

		if (coef != 1) {
			ss << coef;
			needsMult = true;
		}
		degrees_t degs = node->degs;
		for (int i = 0; i < nvar; ++i) {
			if (degs[i] == 0) {
				continue;
			}
			isConst = false;
			if (needsMult) {
				ss << "*";
			}
			ss << vars[i];
			if (degs[i] > 1) {
				ss << "^" << degs[i];
			}
			needsMult = true;
		}

		node = node->next;
		first = false;
		needsMult = false;
	}
	if (isConst && coef == 1) {
		ss << coef;
	}

	return ss.str();
}

static std::string polyToString(Node* node, int nvar, const std::string* vars) {
	const char* charvars[nvar];
	for (int i = 0; i < nvar; ++i) {
		charvars[i] = vars[i].c_str();
	}
	return polyToString(node, nvar, charvars);
}

#endif