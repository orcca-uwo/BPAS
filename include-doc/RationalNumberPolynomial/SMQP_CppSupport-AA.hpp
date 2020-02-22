
#ifndef _SMQP_CPP_SUPPORT_AA_
#define _SMQP_CPP_SUPPORT_AA_

#include <sstream>
#include <string>
#include <gmpxx.h>
#include <iomanip>
#include <ios>

#include "SMQP_Support-AA.h"

typedef mpq_class ratNum_class;

static bool degsToString(std::stringstream& ss, const char** vars, degrees_t degs, bool needsMult, int nvar) {
    int* sizes = getExpOffsetArray(nvar);
    unsigned long long int* masks = getExpMaskArray(nvar);

	bool ret = 1;
	for (int i = 0; i < nvar; ++i) {
		unsigned long int deg = (degs & masks[i]) >> sizes[i];
		if (deg > 0) {
			if (needsMult) {
				ss << "*";
			}
			ss << vars[i];
			if (deg > 1) {
				ss << "^" << deg;
			}
			needsMult = true;
		
			ret = 0;
		}
	}

	free(sizes);
	free(masks);
	return ret;
}

/**
 * Convert a polynomial given its head node, node, to a string.
 * nvar: the number of variables in the polynomial
 * vars: the strings representing each variable. e.g. {"x1","x2","x3"}  
 */
static std::string polyToString_AA(AltArr_t* aa, const char** vars) {
	if (aa == NULL || aa->size == 0) {		
		return "0";
	}

	std::stringstream ss;

	int nvar = aa->nvar;
	bool first = true;
	bool needsMult = false;
	bool isConst = true;
	ratNum_class coef;
	for (int i = 0; i < aa->size; ++i) {		
		coef = ratNum_class(aa->elems[i].coef);
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
		if (nvar > 0) {
			isConst = degsToString(ss, vars, aa->elems[i].degs, needsMult, nvar);
		} else {
			isConst = 1;
		}

		first = false;
		needsMult = false;
	}
	if (isConst && coef == 1) {
		ss << coef;
	}

	return ss.str();
}

static std::string polyToString_AA(AltArr_t* aa, const std::string* vars) {
	int nvar = aa->nvar;
	const char* charvars[nvar];
	for (int i = 0; i < nvar; ++i) {
		charvars[i] = vars[i].c_str();
	}
	return polyToString_AA(aa, charvars);
}

#endif