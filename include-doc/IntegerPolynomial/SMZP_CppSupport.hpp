
#ifndef _SMZP_CPP_SUPPORT_AA_
#define _SMZP_CPP_SUPPORT_AA_

#include <sstream>
#include <string>
#include <gmpxx.h>
#include <iomanip>
#include <ios>

#include "SMZP_Support.h"
#include "../RationalNumberPolynomial/SMQP_CppSupport-AA.hpp"

/**
 * Convert a polynomial given its head node, node, to a string.
 * nvar: the number of variables in the polynomial
 * vars: the strings representing each variable. e.g. {"x1","x2","x3"}  
 */
static std::string polyToString_AAZ(AltArrZ_t* aa, const char** vars) {
	if (aa == NULL || aa->size == 0) {		
		return "0";
	}

	std::stringstream ss;

	int nvar = aa->nvar;
	bool first = true;
	bool needsMult = false;
	bool isConst = true;
	mpz_class coef;
	for (int i = 0; i < aa->size; ++i) {		
		coef = mpz_class(aa->elems[i].coef);
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

static std::string polyToString_AAZ(AltArrZ_t* aa, const std::string* vars) {
	int nvar = aa->nvar;
	const char* charvars[nvar];
	for (int i = 0; i < nvar; ++i) {
		charvars[i] = vars[i].c_str();
	}
	return polyToString_AAZ(aa, charvars);
}

#endif