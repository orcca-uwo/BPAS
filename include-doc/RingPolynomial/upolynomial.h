#ifndef _UPOLYNOMIAL_H_
#define _UPOLYNOMIAL_H_

//#include <type_traits>
#include "../Polynomial/BPASUnivarPolynomial.hpp"
#include "../IntegerPolynomial/uzpolynomial.h"
#include "../RationalNumberPolynomial/urpolynomial.h"
#include "../Utils/TemplateHelpers.hpp"
#include "../Utils/RandomHelpers.hpp"


/**
 * A univariate polynomial over an arbitrary BPASRing represented sparsely.
 *
 *
 * Inheritance of proper base class, and exporting of proper functions,
 * is automatic when the Ring template parameter is specified at compile time
 * by means of std::conditional.
 */
template <class Ring>
class SparseUnivariateTempPoly : public virtual BPASUnivariatePolynomial<Ring>,
								   private Derived_from<Ring, BPASRing> {

									   // public std::conditional<std::is_base_of<BPASGCDDomain<Ring>, Ring>::value,
										  // SparseUnivariateDomainPolynomial<Ring, Derived>,
										  // SparseUnivariateTempPolyBase<Ring, Derived> >::type {

private:
	Symbol name;
	std::vector< UnivariateTerm<Ring> > terms;

	inline bool isEqual(const Derived& b) const {
		if (degree() > 0 && b.degree() > 0 && (name != b.name))
			return 0;
		int n = terms.size();
		int m = b.terms.size();
		if (n != m)
			return 0;
		for (int i = 0; i < n; ++i) {
			UnivariateTerm<Ring> t = b.terms[i];
			if (terms[i].coef != t.coef || terms[i].exp != t.exp)
				return 0;
		}
		return 1;
	}

	inline bool isEqual(const DenseUnivariateRationalPolynomial& b) const {
		if (name != b.variable())
			return 0;
		int n = degree(), m = b.degree();
		if (n != m)
			return 0;
		n = terms.size(), m++;
		int k = 0;
		for (int i = 0; i < m; ++i) {
			if (k < n && terms[k].exp == i) {
				if (terms[k].coef != b.coefficient(i))
					return 0;
				k++;
			}
			else if (b.coefficient(i) != 0)
				return 0;
		}
		return 1;
	}

	inline bool isEqual(const DenseUnivariateIntegerPolynomial& b) const {
		if (name != b.variable())
			return 0;
		int n = degree().get_si(), m = b.degree().get_si();
		if (n != m)
			return 0;
		n = terms.size(), m++;
		int k = 0;
		for (int i = 0; i < m; ++i) {
			if (k < n && terms[k].exp == i) {
				if (terms[k].coef != b.coefficient(i))
					return 0;
				k++;
			}
			else if (b.coefficient(i) != 0)
				return 0;
		}
		return 1;
	}

/* this += t * b */
	inline void pomopo(UnivariateTerm<Ring> t, const Derived& b) {
		int ai = 0, m = b.terms.size();
		int n = terms.size();
		int i;

		// //product
		// std::vector<UnivariateTerm<Ring>> prodTerms;
		// prodTerms.reserve(m);
		// Ring coef;
		// int exp;
		// for (i = 0; i < m; ++i) {
		// 	coef = t.coef + b.terms[i].coef;
		// 	exp = t.exp + b.terms[i].exp;
		// 	prodTerms.emplace_back(coef, exp);
		// }

		// //add
		// // UnivariateTerm<Ring>* mergeOut = (UnivariateTerm<Ring>*) malloc((terms.size() + m)*sizeof(UnivariateTerm<Ring>));
		// UnivariateTerm<Ring>* mergeOut = new UnivariateTerm<Ring>[n + m];
		// UnivariateTerm<Ring>* adata = terms.data(), *proddata = prodTerms.data();
		// int j, curIdx = 0;
		// int aExp = 0, bExp = 0;
		// for (i = 0, j = 0; i < n && j < m; ++curIdx) {
		// 	aExp = adata[i].exp;
		// 	bExp = proddata[j].exp;

		// 	if (aExp == bExp) {
		// 		mergeOut[curIdx] = adata[i];
		// 		mergeOut[curIdx].coef += proddata[j].coef;
		// 		++i;
		// 		++j;
		// 	}
		// 	else if (aExp < bExp) {
		// 		mergeOut[curIdx] = adata[i];
		// 		++i;
		// 	} else {
		// 		mergeOut[curIdx] = proddata[j];
		// 		++j;
		// 	}
		// }

		// for ( ; i < n; ++i, ++curIdx) {
		// 	mergeOut[curIdx] = adata[i];
		// }
		// for ( ; j < m; ++j, ++curIdx) {
		// 	mergeOut[curIdx] = proddata[j];
		// }


		// // std::vector<UnivariateTerm<Ring>> newVec;
		// // newVec.reserve(curIdx);
		// terms.clear();
		// terms.reserve(curIdx);
		// for (i = 0; i < curIdx; ++i) {
		// 	terms.push_back(mergeOut[i]);
		// }

		// // free(mergeOut);
		// delete[] mergeOut;

		// UnivariateTerm<Ring> product;
		const UnivariateTerm<Ring>* bterms = b.terms.data();
		Ring prodCoef, tcoef = t.coef;
		int prodExp, tExp = t.exp, aExp;
		for (int i = 0; i < m; ++i) {
			prodCoef = tcoef * bterms[i].coef;
			prodExp = tExp + bterms[i].exp;

			if (ai >= terms.size()) {
				terms.emplace_back(prodCoef, prodExp);
				ai++;
			} else {
				aExp = terms[ai].exp;
				while (aExp < prodExp) {
					++ai;
					if (ai < terms.size())
						aExp = terms[ai].exp;
					else {
						terms.emplace_back(prodCoef, prodExp);
						++ai;
						break;
					}
				}
				if (aExp == prodExp) {
					terms[ai].coef += prodCoef;
					if (!terms[ai].coef.isZero()) {
						++ai;
					} else {
						terms.erase(terms.begin()+ai);
					}
				}
				else if (aExp > prodExp) {
					terms.emplace(terms.begin()+ai, prodCoef, prodExp);
					++ai;
				}
			}
		}
	}

/* this = c*this + t*b */
	inline void pomopo(Ring c, UnivariateTerm<Ring> t, const Derived& b) {
		int ai = 0, m = b.terms.size();

		for (int i = 0; i < m; ++i) {
			UnivariateTerm<Ring> product;
			product.coef = t.coef * b.terms[i].coef;
			product.exp = t.exp + b.terms[i].exp;

			if (ai >= terms.size()) {
				terms.push_back(product);
				ai++;
			}
			else {
				int e1 = terms[ai].exp, e2 = product.exp;
				while (e1 < e2) {
					terms[ai].coef *= c;
					ai++;
					if (ai < terms.size())
						e1 = terms[ai].exp;
					else {
						terms.push_back(product);
						ai++;
						break;
					}
				}
				if (e1 == e2) {
					terms[ai].coef *= c;
					terms[ai].coef += product.coef;
					if (terms[ai].coef.isZero())
						terms.erase(terms.begin()+ai);
					else { ai++; }
				}
				else if (e1 > e2) {
					terms.insert(terms.begin()+ai, product);
					ai++;
				}
			}
		}
		for (int i = ai; i < terms.size(); ++i)
			terms[i].coef *= c;
	}

	// For subresultant
	// y = s;
	inline Derived LazardSe(Derived& sd, Derived& sdm, Ring& y) const {
		int n = (sd.degree() - sdm.degree()).get_si() - 1;
		if (!n) { return sdm; }
		Ring x = sdm.leadingCoefficient();
		// Ring y = sd.leadingCoefficient();
		int a = (int) pow(2, floor(log2(n)));
		Ring c = x;
		Derived se = sdm;

		n -= a;
		/* int orign = n; */
		while (a != 1) {
			a >>= 1;
			c *= c;
			c /= y;
			if (n >= a) {
				c *= x;
				c /= y;
				n -= a;
			}
		}
		se *= c;
		/* if (orign != 0) { */
			se /= y;
		/* } */
		return se;
	}

	inline Derived DucosSem(Derived& a, Derived& sdm, Derived& se, Ring sd) const {

	    Integer d = a.degree();
	    int e = sdm.degree().get_si();
		Derived res;
		res.name = name;
		Ring ec = se.leadingCoefficient();

		/* std::cout << "Sd := " << a << std::endl; */
		/* std::cout << "Sdm := " << sdm << std::endl; */
		/* std::cout << "Se := " << se << std::endl; */
		/* std::cout << "sd := " << sd << std::endl; */
		/* std::cout << "se := " << ec << std::endl; */

		if (d == e){
	// Streamlined resultant computation for regular case
			Ring dc = a.leadingCoefficient();
			Ring rTemp;
			Derived supTemp;
			supTemp.name = name;
	// compute first of two main terms
			res = a;
			res *= ec;
			rTemp = a.coefficient(e);
			supTemp = se;
			supTemp *= rTemp;
			res -= supTemp;
			res /= dc;
			res *= ec;
	// compute second of two main terms
			rTemp = se.coefficient(e-1);
			supTemp.zero();
			supTemp.setCoefficient(0,rTemp);
			rTemp = -ec;
			supTemp.setCoefficient(1,rTemp);
			supTemp *= se;
			res += supTemp;
	// divide out dc to obtain resultant
			res /= dc;
			return res;
		}
		else {
	// Ducos algorithm for defective case
			Ring mc = sdm.leadingCoefficient();
			Derived* P = new Derived[d.get_si()];
			for (int i = 0; i < e; ++i) {
				P[i].setVariableName(name);
				P[i].setCoefficient(i, ec);
			}
			for (int i = e+1; d > i; ++i)
				P[i].setVariableName(name);
			P[e].setVariableName(name);
			P[e].setCoefficient(e, ec);
			P[e] -= se;
			ec.one();
			res.setCoefficient(1, ec);
			UnivariateTerm<Ring> t (ec, 1);
			for(int i = e+1; i < d; ++i) {
				if (!P[i-1].coefficient(e-1).isZero()) {
					ec = -P[i-1].coefficient(e-1);
				}
				else {
					ec.zero(); // prevents poly name error
				}
				ec /= mc;
				P[i] = sdm;
				P[i] *= ec;
				//P[i] += res * P[i-1];
				P[i].pomopo(t, P[i-1]);
			}
			res *= P[d.get_si()-1];

			Derived D;
			for (int i = 0; i < d; ++i) {
				P[i] *= a.coefficient(i);
				D += P[i];
			}
			D /= a.leadingCoefficient();
			delete [] P;

			ec = -res.coefficient(e);
			//res += D;
			//res *= mc;
			t.coef = mc;
			t.exp = 0;
			res.pomopo(mc, t, D);
			sdm *= ec;
			res += sdm;
			res /= sd;

			if ((d - e + 1) % 2 == 1) {
				return -res;
			}
			return res;
		}
	}

	inline UnivariateTerm<Ring> leadingTerm() const {
		int n = terms.size();
		if (n) { return terms[n-1]; }
		else { return UnivariateTerm<Ring>(); }
	}

public:
		mpz_class characteristic;
	// static bool isPrimeField;
	// static bool isSmallPrimeField;
	// static bool isComplexField;
	/**
	* Construct a polynomial
	*
	* @param
	**/
	SparseUnivariateTempPoly<Ring, Derived> () : name("%"), terms() {
		Ring e;
		characteristic = e.characteristic;
	}

	/**
	* Copy constructor
	*
	* @param b: A sparse univariate polynomial
	**/
	SparseUnivariateTempPoly<Ring, Derived> (const SparseUnivariateTempPoly<Ring, Derived>& b) : name(b.name), terms(b.terms) {
		Ring e;
		characteristic = e.characteristic;
	}

	SparseUnivariateTempPoly<Ring, Derived> (int a) {
		UnivariateTerm<Ring> t;
		t.coef = Ring(a);
		terms.push_back(t);
	}

	SparseUnivariateTempPoly<Ring, Derived> (const Integer& c) {
		UnivariateTerm<Ring> t;
		t.coef = Ring(c);
		terms.push_back(t);
	}

	SparseUnivariateTempPoly<Ring, Derived> (const RationalNumber& c) {
		UnivariateTerm<Ring> t;
		t.coef = Ring(c);
		t.exp = 0;
		terms.push_back(t);
	}

	SparseUnivariateTempPoly<Ring, Derived> (const ComplexRationalNumber& c) {
		UnivariateTerm<Ring> t;
		t.coef = Ring(c);
		terms.push_back(t);
	}

	SparseUnivariateTempPoly<Ring, Derived> (const DenseUnivariateIntegerPolynomial& b) {
		name = b.variable();
		for (int i = 0; i <= b.degree().get_si(); ++i) {
			Ring e (Integer(b.coefficient(i)));
			if (!e.isZero()) {
				UnivariateTerm<Ring> t;
				t.coef = e;
				t.exp = i;
				terms.push_back(t);
			}
		}
	}

	SparseUnivariateTempPoly<Ring, Derived> (const DenseUnivariateRationalPolynomial& b) {
		name = b.variable();
		for (int i = 0; i <= b.degree().get_si(); ++i) {
			Ring e (RationalNumber(b.coefficient(i)));
			if (!e.isZero()) {
				UnivariateTerm<Ring> t;
				t.coef = e;
				t.exp = i;
				terms.push_back(t);
			}
		}
	}

	SparseUnivariateTempPoly<Ring, Derived> (Symbol sym) : name(sym), terms() {
		this->one();
		terms[0].exp = 1;
	}

	/**
	 * Destroy the polynomial
 	 *
     * @param
     **/
	~SparseUnivariateTempPoly<Ring, Derived> () {
		terms.clear();
	}

	/**
	 * Get the number of terms
	 *
	 * @param
	 **/
	inline Integer numberOfTerms() const {
		return terms.size();
	}

	/**
	 * Get the degree of the polynomial
	 *
	 * @param
	 **/
	inline Integer degree() const {
		int n = terms.size();
		if (n) { return terms[n-1].exp; }
		else { return 0; }
	}

	/**
	 * Get the leading coefficient
	 *
	 * @param
	 **/
	inline Ring leadingCoefficient() const {
		int n = terms.size();
		if (n) { return terms[n-1].coef; }
		else { return Ring(); }
	}

	inline Ring trailingCoefficient() const {
		int n = terms.size();
		if (n) {
			return terms[0].coef;
		} else {
			return Ring();
		}
	}

	/**
	 * Get a coefficient
	 *
	 * @param k: The exponent
	 **/
	inline Ring coefficient(int k) const {
		int n = terms.size();
		for (int i = 0; i < n; ++i) {
			if (k == terms[i].exp)
				return terms[i].coef;
			else if (k < terms[i].exp)
				break;
		}
		return Ring();
	}

	/**
	 * Set a coeffcient with its exponent
	 *
	 * @param e: The exponent
	 * @param c: The coefficient
	 **/
	inline void setCoefficient(int e, const Ring& c) {
		UnivariateTerm<Ring> b;
		b.coef = c;
		b.exp = e;


		int k = terms.size() - 1;
		if ((k < 0 || b.exp > terms[k].exp) && !c.isZero())
			terms.push_back(b);
		else {
			for (int i = 0; i <= k; ++i) {
				int e1 = terms[i].exp, e2 = b.exp;
				if (e1 == e2) {
					terms[i].coef = b.coef;
					if (terms[i].coef.isZero())
						terms.erase(terms.begin()+i);
					break;
				}
				else if (e1 > e2) {
					if (!c.isZero())
						terms.insert(terms.begin()+i, b);
					break;
				}
			}
		}
	}

	/**
	 * Get the variable name
	 *
	 * @param
	 **/
	inline Symbol variable() const {
		return name;
	}

	/**
	 * Set the variable name
	 *
	 * @param c: Variable's name
	 **/
	inline void setVariableName (const Symbol& c) {
		name = c;
	}

	inline Derived unitCanonical(Derived* u = NULL, Derived* v = NULL) const {
		Ring lead = leadingCoefficient();
		Ring ru, rv;
		Ring canon = lead.unitCanonical(&ru, &rv);
		Derived ret = *this * ru;
		if (u != NULL) {
			*u = ru;
		}
		if (v != NULL) {
			*v = rv;
		}

		return ret;
	}

	/**
	 * Overload operator =
	 *
	 * @param b: A sparse univariate polynomial
	 **/
	inline Derived& operator= (const Derived& b) {
		if (this != &b) {
			terms.clear();
			name = b.name;
			terms = b.terms;
			Ring e;
			characteristic = e.characteristic;
		}
		return dynamic_cast<Derived&>(*this);
	}

	/**
	 * Overload operator =
	 *
	 * @param r: A base ring element
	 */
	inline Derived& operator= (const Ring& r) {
		terms.clear();
		UnivariateTerm<Ring> t;
		t.coef = Ring(r);
		terms.push_back(t);
		return dynamic_cast<Derived&>(*this);
	}


	/**
	 * Overload operator !=
	 *
	 * @param b: A sparse univariate polynomial
	 **/
	inline bool operator!= (const Derived& b) const {
		return !(isEqual(b));
	}

	/**
	 * Overload operator ==
	 *
	 * @param b: A sparse univariate polynomial
	 **/
	inline bool operator== (const Derived& b) const {
		return isEqual(b);
	}

	inline bool operator== (const DenseUnivariateRationalPolynomial& b) const {
		return isEqual(b);
	}

	inline bool operator== (const DenseUnivariateIntegerPolynomial& b) const {
		return isEqual(b);
	}

	/**
	 * Is zero polynomial
	 *
	 * @param
	 **/
	inline bool isZero() const {
		return !(terms.size());
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
		if (terms.size() == 1 && !terms[0].exp)
			return terms[0].coef.isOne();
		return 0;
	}

	/**
	 * Set polynomial to 1
	 *
	 * @param
	 **/
	inline void one() {
		terms.clear();
		UnivariateTerm<Ring> t;
		t.coef.one();
		t.exp = 0;
		terms.push_back(t);
	}

	/**
	 * Is polynomial a constant -1
	 *
	 * @param
	 **/
	inline bool isNegativeOne() const {
		if (terms.size() == 1 && !terms[0].exp)
			return terms[0].coef.isNegativeOne();
		return 0;
	}

	/**
	 * Set polynomial to -1
	 *
	 * @param
	 **/
	inline void negativeOne() {
		terms.clear();
		UnivariateTerm<Ring> t;
		t.coef.negativeOne();
		t.exp = 0;
		terms.push_back(t);
	}

	/**
	 * Is a constant
	 *
	 * @param
	 **/
	inline int isConstant() const {
		if (terms.size() == 1 && !terms[0].exp)
			return terms[0].coef.isConstant();
		return 0;
	}

	/**
	 * Convert to a constant
	 *
	 * @param
	 **/
	inline Ring convertToConstant() {
		if (terms.size() && !terms[0].exp)
			return terms[0].coef;
		else {
			Ring e;
			e.zero();
			return e;
		}
	}

	/**
	 * Content of the polynomial
	 *
	 * @param
	 **/
	inline Ring content() const override {
		Ring c;
		int n = terms.size();
		if (n) {
			c = terms[0].coef;
			for (int i = 1; i < n; ++i) {
				c = c.gcd(terms[i].coef);
				if (c.isOne())
					break;
			}
		}
	//else { c.zero(); }
		return c;
	}

	inline Derived primitivePart() const {
		//TODO
		std::cerr << "BPAS ERROR: SUP<Ring>::primitivePart NOT YET IMPLEMENTED" << std::endl;
		return (*this);
	}

	/**
	 * Overload operator ^
	 * replace xor operation by exponentiation
	 *
	 * @param e: The exponentiation, e > 0
	 **/
	inline Derived operator^ (long long int e) const {
		Derived res;
		res.name = name;
	//res.one();
	//unsigned long int q = e / 2, r = e % 2;
	//Derived power2 = *this * *this;
	//for (int i = 0; i < q; ++i)
	//	res *= power2;
	//if (r) { res *= *this; }
		if (isZero() || isOne() || e == 1)
			res = *this;
		else if (e == 2)
			res = *this * *this;
		else if (e > 2) {
			Derived x (*this);
			res.one();

			while (e) {
				if (e % 2) { res *= x; }
				x = x * x;
				e >>= 1;
			}
		}
		else if (!e)
			res.one();
		else {
			res = *this;
		}
		return res;
	}

	/**
	 * Overload operator ^=
	 * replace xor operation by exponentiation
	 *
	 * @param e: The exponentiation, e > 0
	 **/
	inline Derived& operator^= (long long int e) {
		*this = *this ^ e;
		return dynamic_cast<Derived&>(*this);
	}

	/**
	 * Overload operator <<
	 * replace by muplitying x^k
	 *
	 * @param k: The exponent of variable, k > 0
	 **/
	inline Derived operator<< (int k) const {
		Derived r (*this);
		return (r <<= k);
	}

	/**
	 * Overload operator <<=
	 * replace by muplitying x^k
	 *
	 * @param k: The exponent of variable, k > 0
	 **/
	inline Derived& operator<<= (int k) {
		for (int i = 0; i < terms.size(); ++i)
			terms[i].exp += (unsigned long int) k;
		return dynamic_cast<Derived&>(*this);
	}

	/**
	 * Overload operator >>
	 * replace by dividing x^k, and
	 * return the quotient
	 *
	 * @param k: The exponent of variable, k > 0
	 **/
	inline Derived operator>> (int k) const {
		Derived r (*this);
		return (r >>= k);
	}

	/**
	 * Overload operator >>=
	 * replace by dividing x^k, and
	 * return the quotient
	 *
	 * @param k: The exponent of variable, k > 0
	 **/
	inline Derived& operator>>= (int k) {
		int i = 0;
		unsigned long int e = (unsigned long int) k;
		while (i < terms.size()) {
			if (terms[i].exp >= e) {
				terms[i].exp -= e;
				i++;
			}
			else { terms.erase(terms.begin()); }
		}
		return dynamic_cast<Derived&>(*this);
	}

	/**
	 * Overload operator +
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived operator+ (const Derived& b) const {
		Derived res(*this);
		return (res += b);
	}

	/**
	 * Overload operator+=
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived& operator+= (const Derived& b) {
		if (!terms.size()) { return (*this = b); }
		if (!b.terms.size()) { return dynamic_cast<Derived&>(*this); }
		if (isConstant()) { return (*this = b + terms[0].coef); }
		if (b.isConstant()) { return (*this += b.terms[0].coef); }
		if (name != b.name) {
			std::cout << "BPAS: error, trying to add between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}

		int ai = 0;
		for (int i = 0; i < b.terms.size(); ++i) {
			UnivariateTerm<Ring> bt = b.terms[i];
			if (ai >= terms.size()) {
				terms.push_back(bt);
				ai++;
			}
			else {
				int e1 = terms[ai].exp, e2 = bt.exp;
				while (e1 < e2) {
					ai++;
					if (ai < terms.size())
						e1 = terms[ai].exp;
					else {
						terms.push_back(bt);
						ai++;
						break;
					}
				}
				if (e1 == e2) {
					terms[ai].coef += bt.coef;
					if (terms[ai].coef.isZero())
						terms.erase(terms.begin()+ai);
					else { ai++; }
				}
				else if (e1 > e2)
					terms.insert(terms.begin()+ai, bt);
			}
		}
		return dynamic_cast<Derived&>(*this);
	}

	/**
	 * Overload operator +
	 *
	 * @param e: A coefficient constant
	 **/
	inline Derived operator+ (const Ring& e) const {
		Derived r (*this);
		return (r += e);
	}

	/**
	 * Overload operator +=
	 *
	 * @param e: A coefficient constant
	 **/
	inline Derived& operator+= (const Ring& e) {
		if (!e.isZero()) {
			if (terms.size()) {
				UnivariateTerm<Ring> a = terms[0];
				if (a.exp) {
					a.coef = e;
					a.exp = 0;
					terms.insert(terms.begin(), a);
				}
				else {
					terms[0].coef += e;
					if (terms[0].coef.isZero())
						terms.erase(terms.begin());
				}
			}
			else {
				UnivariateTerm<Ring> a;
				a.coef = e;
				a.exp = 0;
				terms.push_back(a);
			}
		}
		return dynamic_cast<Derived&>(*this);
	}

	inline friend Derived operator+ (const Ring& c, const Derived& p) {
		return (p + c);
	}

	/**
	 * Overload operator -, negate
	 *
	 * @param
	 **/
	inline Derived operator- () const {
		Derived res;
		res.name = name;
		for (int i = 0; i < terms.size(); ++i) {
			UnivariateTerm<Ring> t;
			t.coef = -terms[i].coef;
			t.exp = terms[i].exp;
			res.terms.push_back(t);
		}
		return res;
	}

	/**
	 * Subtract another polynomial
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived operator- (const Derived& b) const {
		Derived res(*this);
		return (res -= b);
	}

	/**
	 * Overload operator -=
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived& operator-= (const Derived& b) {
		if (!terms.size()) { return (*this = -b); }
		if (!b.terms.size()) { return dynamic_cast<Derived&>(*this); }
		if (isConstant()) { return (*this = -b + terms[0].coef); }
		if (b.isConstant()) { return (*this -= b.terms[0].coef); }
		if (name != b.name) {
			std::cout << "BPAS: error, trying to subtract between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}

		int ai = 0;
		for (int i = 0; i < b.terms.size(); ++i) {
			UnivariateTerm<Ring> t = b.terms[i];
			t.coef = -t.coef;

			if (ai >= terms.size()) {
				terms.push_back(t);
				ai++;
			}
			else {
				int e1 = terms[ai].exp, e2 = t.exp;
				while (e1 < e2) {
					ai++;
					if (ai < terms.size())
						e1 = terms[ai].exp;
					else {
						terms.push_back(t);
						ai++;
						break;
					}
				}
				if (e1 == e2) {
					terms[ai].coef += t.coef;
					if (terms[ai].coef.isZero())
						terms.erase(terms.begin()+ai);
					else { ai++; }
				}
				else if (e1 > e2)
					terms.insert(terms.begin()+ai, t);
			}
		}
		return dynamic_cast<Derived&>(*this);
	}

	/**
	 * Overload operator -
	 *
	 * @param e: A coefficient constant
	 **/
	inline Derived operator- (const Ring& e) const {
		Derived r (*this);
		return (r -= e);
	}

	/**
	 * Overload operator -=
	 *
	 * @param e: A coefficient constant
	 **/
	inline Derived& operator-= (const Ring& e) {
		if (!e.isZero()) {
			if (terms.size()) {
				UnivariateTerm<Ring> t = terms[0];
				if (t.exp) {
					t.coef = -e;
					t.exp = 0;
					terms.insert(terms.begin(), t);
				}
				else {
					terms[0].coef -= e;
					if (terms[0].coef.isZero())
						terms.erase(terms.begin());
				}
			}
			else {
				UnivariateTerm<Ring> t;
				t.coef = -e;
				t.exp = 0;
				terms.push_back(t);
			}
		}
		return dynamic_cast<Derived&>(*this);
	}

	inline friend Derived operator- (const Ring& c, const Derived& p) {
		return (-p + c);
	}

	/**
	 * Multiply another polynomial
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived operator* (const Derived& b) const {
		int n = terms.size(), m = b.terms.size();
		if (!n)
			return *this;
		if (!m)
			return b;

		Derived res;
		if (degree() == 0) {
			for (int i = 0; i < m; ++i) {
				UnivariateTerm<Ring> bt = b.terms[i];
				UnivariateTerm<Ring> t;
				t.coef = terms[0].coef * bt.coef;
				t.exp = bt.exp;
				res.terms.push_back(t);
			}
			res.name = b.name;
			return res;
		}
		res.name = name;
		if (b.degree() == 0) {
			UnivariateTerm<Ring> bt = b.terms[0];
			for (int i = 0; i < n; ++i) {
				UnivariateTerm<Ring> t;
				t.coef = terms[i].coef * bt.coef;
				t.exp = terms[i].exp;
				res.terms.push_back(t);
			}
			return res;
		}

		if (name != b.name) {
			std::cout << "BPAS: error, trying to multiply between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}

		if (n + m < 64) {
			if (n <= m) {
				for (int i = 0; i < n; ++i)
					res.pomopo(terms[i], b);
			}
			else {
				for (int i = 0; i < m; ++i)
					res.pomopo(b.terms[i], *this);
			}
		}
		else {
			int s = (m > n) ? m : n;
			s >>= 1;
			Derived f0, f1;
			f0.name = f1.name = name;
			for (int i = 0; i < n; ++i) {
				if (terms[i].exp < s)
					f0.terms.push_back(terms[i]);
				else {
					UnivariateTerm<Ring> t (terms[i].coef, terms[i].exp-s);
					f1.terms.push_back(t);
				}
			}
			Derived g0, g1;
			g0.name = g1.name = name;
			for (int i = 0; i < m; ++i) {
				if (b.terms[i].exp < s)
					g0.terms.push_back(b.terms[i]);
				else {
					UnivariateTerm<Ring> t (b.terms[i].coef, b.terms[i].exp-s);
					g1.terms.push_back(t);
				}
			}
			Derived t0, t1;
			t0.name = t1.name = name;
			n = f0.terms.size(), m = g0.terms.size();
			if (n <= m) {
				for (int i = 0; i < n; ++i)
					res.pomopo(f0.terms[i], g0);
			}
			else {
				for (int i = 0; i < m; ++i)
					res.pomopo(g0.terms[i], f0);
			}
			n = f1.terms.size(), m = g1.terms.size();
			if (n <= m) {
				for (int i = 0; i < n; ++i)
					t0.pomopo(f1.terms[i], g1);
			}
			else {
				for (int i = 0; i < m; ++i)
					t0.pomopo(g1.terms[i], f1);
			}
			f0 += f1, g0 += g1;
			n = f0.terms.size(), m = g0.terms.size();
			if (n <= m) {
				for (int i = 0; i < n; ++i)
					t1.pomopo(f0.terms[i], g0);
			}
			else {
				for (int i = 0; i < m; ++i)
					t1.pomopo(g0.terms[i], f0);
			}
			t1 -= res + t0;
			for (int i = 0; i < t1.terms.size(); ++i)
				t1.terms[i].exp += s;
			s <<= 1;
			for (int i = 0; i < t0.terms.size(); ++i)
				t0.terms[i].exp += s;
			res += t0 + t1;
		}
		return res;
	}

	/**
	 * Overload operator *=
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived& operator*= (const Derived& b) {
		*this = *this * b;
		return dynamic_cast<Derived&>(*this);
	}

	/**
	 * Overload operator *
	 *
	 * @param c: A coefficient
	 **/
	inline Derived operator* (const Ring& c) const {
		Derived r (*this);
		return (r *= c);
	}

	inline Derived operator* (const sfixn& e) const {
		Derived r (*this);
		return (r *= e);
	}

	/**
	 * Overload operator *=
	 *
	 * @param c: A coefficient
	 **/
	inline Derived& operator*= (const Ring& c) {
		if (!isZero()) {
			if (!c.isZero() && !c.isOne()) {
				for (int i = 0; i < terms.size(); ++i)
					terms[i].coef *= c;
			}
			else if (c.isZero())
				terms.clear();
		}
		return dynamic_cast<Derived&>(*this);
	}

	inline Derived& operator*= (const sfixn& e) {
		if (e != 0 && e != 1) {
			for (int i = 0; i < terms.size(); ++i)
				terms[i].coef *= e;
		}
		else if (e == 0) { zero(); }
		return dynamic_cast<Derived&>(*this);
	}

	inline friend Derived operator* (const Ring& e, const Derived& p) {
		return (p * e);
	}

	inline friend Derived operator* (const sfixn& e, const Derived& p) {
		return (p * e);
	}

	/**
	 * Overload operator /
	 * EdeDivision
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived operator/ (const Derived& b) const {
		Derived rem(*this);
		rem /= b;
		return rem;
	}

	/**
	 * Overload operator /=
	 * ExactDivision
	 *
	 * @param b: A univariate polynomial
	 **/
	inline Derived& operator/= (const Derived& b) {
		if (b.isZero()) {
			std::cout << "BPAS: error, dividend is zero from Derived." << std::endl;
			exit(1);
		}
		if (b.isConstant()) { return (*this /= b.terms[0].coef); }
		if (isConstant()) {
			zero();
			return dynamic_cast<Derived&>(*this);
		}
		if (name != b.name) {
			std::cout << "BPAS: error, trying to exact divide between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}

		Derived q;
		q.name = name;

		Integer db = b.degree();
		UnivariateTerm<Ring> bt = b.leadingTerm();

		if (db == 1 && bt.coef.isOne()) {
			if (b.terms.size() > 1) {
				int k = 0;
				Ring e, prev;
				UnivariateTerm<Ring> t, at = leadingTerm();
				if (!terms[0].exp) {
					prev = t.coef = terms[0].coef / b.terms[0].coef;
					t.exp = 0;
					q.terms.push_back(t);
					k++;
				}
				else { prev.zero(); }
				for (int i = 1; i < at.exp; ++i) {
					if (k < terms.size() && terms[k].exp == i) {
						e = terms[k].coef;
						k++;
					}
					else { e.zero(); }
					t.coef = (e - prev) / b.terms[0].coef;
					if (!t.coef.isZero()) {
						t.exp = i;
						q.terms.push_back(t);
					}
					prev = t.coef;
				}
				if (prev == at.coef)
					return (*this = q);
				else {
					std::cout << "BPAS: error, not exact division in Derived." << std::endl;
					exit(1);
				}
			}
			else {
				if (!terms[0].exp) {
					std::cout << "BPAS: error, not exact division in Derived." << std::endl;
					exit(1);
				}
				else {
					for (int i = 0; i < terms.size(); ++i)
						terms[i].exp--;
					return dynamic_cast<Derived&>(*this);
				}
			}
		}

		while (!isZero() && degree() >= db) {
			UnivariateTerm<Ring> at = leadingTerm();
			UnivariateTerm<Ring> lc, nlc;
			lc.coef = at.coef / bt.coef;
			lc.exp = at.exp - bt.exp;
			nlc.coef = -lc.coef;
			nlc.exp = lc.exp;
			pomopo(nlc, b);
			q.terms.insert(q.terms.begin(), lc);
		}
		if (!isZero()) {
			std::cout << "BPAS: error, not exact division in Derived." << std::endl;
			exit(1);
		}
		return (*this = q);
	}

	/**
	 * Overload operator /
	 *
	 * @param e: A coefficient constant
	 **/
	inline Derived operator/ (const Ring& e) const {
		Derived r (*this);
		r /= e;
		return r;
	}

	/**
	 * Overload operator /=
	 *
	 * @param e: A coefficient constant
	 **/
	inline Derived& operator/= (const Ring& e) {
		if (e.isZero()) {
			std::cout << "BPAS: error, dividend is zero from Derived." << std::endl;
			exit(1);
		}
		else if (!e.isOne()) {
			int i = 0;
			while (i < terms.size()) {
				terms[i].coef /= e;
				if (terms[i].coef.isZero())
					terms.erase(terms.begin()+i);
				else { ++i; }
			}
		}
		return dynamic_cast<Derived&>(*this);
	}

	inline friend Derived operator/ (const Ring& e, const Derived& p) {
		if (p.isZero()) {
			std::cout << "BPAS: error, dividend is zero from Derived." << std::endl;
			exit(1);
		}
		Derived q;
		q.name = p.name;

		if (p.isConstant()) {
			q += e;
			return (q /= p.terms[0].coef);
		}
		else { return q; }
	}

	/**
 	 * Negate the polynomial
	 *
	 * @param
	 **/
	inline void negate() {
		for (int i = 0; i < terms.size(); ++i)
			terms[i].coef = -terms[i].coef;
	}

	/**
	 * Monic division
	 * Assuming the leading coefficient of dividend is 1
	 * Return quotient and itself becomes remainder
	 *
	 * @param b: The dividend polynomial
	 **/
	inline Derived monicDivide(const Derived& b) {
		if (b.isZero()) {
			std::cout << "BPAS: error, dividend is zero from Derived." << std::endl;
			exit(1);
		}
		else if (!b.leadingCoefficient().isOne()) {
			std::cout << "BPAS: error, leading coefficient is not one in monicDivide() from Derived." << std::endl;
			exit(1);
		}

		if (b.isConstant()) {
			Derived r (*this);
			zero();
			return r;
		}
		if (isConstant()) {
			Derived r;
			r.zero();
			return r;
		}
		if (name != b.name) {
			std::cout << "BPAS: error, trying to monic divide between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}

		Derived quo;
		quo.name = name;
		UnivariateTerm<Ring> bt = b.leadingTerm();
		while (degree() >= b.degree()) {
			UnivariateTerm<Ring> at = leadingTerm();
			UnivariateTerm<Ring> nlc;
			nlc.coef = -at.coef;
			nlc.exp = at.exp - bt.exp;
			pomopo(nlc, b);
			at.exp = nlc.exp;
			quo.terms.insert(quo.terms.begin(), at);
		}
		return quo;
	}

	/**
	 * Monic division
	 * Assuming the leading coefficient of dividend is 1
	 * Return quotient
	 *
	 * @param b: The dividend polynomial
	 * @param rem: The remainder polynomial
	 **/
	inline Derived monicDivide(const Derived& b, Derived* rem) const {
		std::cout << "*this " <<  *this << std::endl;
		*rem = *this; std::cout<< "salam"<<std::endl;
		return rem->monicDivide(b);
	}

	/**
	 * Lazy pseudo dividsion
	 * Return the quotient and itself becomes remainder
	 * e is the exact number of division steps
	 *
	 * @param b: The dividend polynomial
	 * @param c: The leading coefficient of b to the power e
	 * @param d: That to the power deg(a) - deg(b) + 1 - e
	 **/
	inline Derived lazyPseudoDivide (const Derived& b, Ring* c, Ring* d=NULL) {
		if (d == NULL)
			d = new Ring;
		Integer da = degree(), db = b.degree();
		if (b.isZero() || db == 0) {
			std::cout << "BPAS: error, dividend is zero or constant from Derived." << std::endl;
			exit(1);
		}
		c->one(), d->one();
		if (isConstant()) {
			Derived r;
			r.zero();
			return r;
		}
		if (name != b.name) {
			std::cout << "BPAS: error, trying to pseudo divide between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}

		if (da < db) {
			Derived r;
			r.name = name;
			return r;
		}

		Derived quo;
		quo.name = name;
		Ring blc = b.leadingTerm().coef;

		int e = 0;
		Integer diff = da - db;
		while (degree() >= db) {
			UnivariateTerm<Ring> at = leadingTerm();
			UnivariateTerm<Ring> nlc;
			nlc.coef = -at.coef;
			nlc.exp = at.exp - db.get_si();

			*c *= blc;
			e++;
			pomopo(blc, nlc, b);
			at.exp = nlc.exp;
			quo.terms.insert(quo.terms.begin(), at);
		}
		for (int i = e; diff >= i; ++i)
			*d *= blc;
		return quo;
	}

	/**
	 * Lazy pseudo dividsion
	 * Return the quotient
	 * e is the exact number of division steps
	 *
	 * @param b: The divident polynomial
	 * @param rem: The remainder polynomial
	 * @param c: The leading coefficient of b to the power e
	 * @param d: That to the power deg(a) - deg(b) + 1 - e
	 **/
	inline Derived lazyPseudoDivide (const Derived& b, Derived* rem, Ring* c, Ring* d) const {
		*rem = *this;
		return rem->lazyPseudoDivide(b, c, d);
	}

	/**
	 * Pseudo dividsion
	 * Return the quotient and itself becomes remainder
	 *
	 * @param b: The divident polynomial
	 * @param d: The leading coefficient of b
	 * 	     to the power deg(a) - deg(b) + 1
	 **/
	inline Derived pseudoDivide (const Derived& b, Ring* d=NULL) {
		Ring c;
		if (d == NULL)
			d = new Ring;
		Derived quo = lazyPseudoDivide(b, &c, d);
		quo *= *d;
		*this *= *d;
		*d *= c;
		return quo;
	}

	/**
	 * Pseudo dividsion
	 * Return the quotient
	 *
	 * @param b: The divident polynomial
	 * @param rem: The remainder polynomial
	 * @param d: The leading coefficient of b
	 *           to the power deg(a) - deg(b) + 1
	 **/
	inline Derived pseudoDivide (const Derived& b, Derived* rem, Ring* d) const {
		Ring c;
		Derived quo = lazyPseudoDivide(b, rem, &c, d);
		quo *= *d;
		*rem *= *d;
		*d *= c;
		return quo;
	}

	/**
	 * Compute k-th derivative
	 *
	 * @param k: k-th derivative
	 **/
	inline void differentiate(int k) {
		if (k <= 0) { return; }
		int i = 0;
		while (i < terms.size()) {
			if (terms[i].exp >= k) {
				for (int j = 0; j < k; ++j)
					terms[i].coef *= Ring(terms[i].exp - j);
				terms[i].exp -= k;
				i++;
			}
			else
				terms.erase(terms.begin());
		}
	}

	/**
	 * Convert current object to its derivative
	 *
	 **/
	inline void differentiate() {
		this->differentiate(1);
	}

	/**
	 * Return k-th derivative
	 *
	 * @param k: k-th derivative, k > 0
	 **/
	inline Derived derivative(int k) const {
		Derived a(*this);
		a.differentiate(k);
		return a;
	}

	/**
	 * Compute derivative
	 *
	 **/
	inline Derived derivative() const {
		return this->derivative(0);
	}

	/**
	 * Compute integral with constant of integration set to 0
	 *
	 * @param
	 **/
	inline void integrate() {
		int i = terms.size()-1;
		while (i > -1) {
			terms[i].coef /= (terms[i].exp + 1);
			terms[i].exp += 1;
			i--;
		}
	}

	/**
	 * Compute integral with constant of integration 0
	 *
	 **/
	inline Derived integral() const {
		Derived a(*this);
		a.integrate();
		return a;
	}

	/**
	 * Is trailing coefficient zero
	 *
	 * @param
	 **/
	inline bool isConstantTermZero() const {
		if (isZero())
			return 1;
		int n = terms.size();
		if (n && terms[0].exp == 0 && terms[0].coef == Ring(0))
			return 1;
		return 0;
	}

	/**
	 * Evaluate f(x)
	 *
	 * @param x: Evaluation point
	 **/
	inline Ring evaluate(const Ring& x) const {
		int d = terms.size() - 1;
		if (d < 0) { return Ring(); }
		int e = terms[d].exp - 1;
		Ring px = terms[d].coef;
		d--;
		for (int i = e; i > -1; --i) {
			px *= x;
			if (i == terms[d].exp && d > -1) {
				px += terms[d].coef;
				d--;
			}
		}
		return px;
	}

	/**
	 * Evaluate f(x)
	 *
	 * @param x: Evaluation point in larger ring, i.e. a ring in which the Ring of SUP<Ring> can be embedded
	 **/
	template <class LargerRing>
	inline LargerRing evaluate(const LargerRing& x) const {
	// we might need a way of checking that this is always possible
		int d = terms.size() - 1;
		if (d < 0) { return LargerRing(); }
		int e = terms[d].exp - 1;
		LargerRing px = (LargerRing)terms[d].coef;
		LargerRing a;
		d--;
		for (int i = e; i > -1; --i) {
			px *= x;
			if (i == terms[d].exp && d > -1) {
				a = (LargerRing)terms[d].coef;
				px += a;
				d--;
			}
		}
		return px;
	}

	inline void fillChain (std::vector<Derived>& chain) const {
		Derived zero;
		zero.zero();
		int fullSize(chain[chain.size()-2].degree().get_ui()+2);
		int delta;
	//	std::cerr << "chain.size() = " << chain.size() << std::endl;
	//	std::cerr << "fullSize = " << fullSize << std::endl;
		if (chain.size() < fullSize) {
			chain.reserve(fullSize);
			for (int i=chain.size()-2; i>0; --i) {
				if (chain[i].degree() != chain[i-1].degree()+1) {
					delta = chain[i].degree().get_ui() - chain[i-1].degree().get_ui();
					if (i > 1) {
						i = i-1;
						for (int j=0; j<delta-2; ++j)
							chain.insert(chain.begin()+i,zero);
					}
					else {
						for (int j=0; j<delta-1; ++j)
							chain.insert(chain.begin()+i,zero);
					}
				}
			}
			if (chain[0].degree() != 0) {
					for (int j=0; j<chain[0].degree(); ++j)
						chain.insert(chain.begin(),zero);
			}
		}
	//	std::cerr << "chain.size() = " << chain.size() << std::endl;
	}

	/**
	 * Subresultant Chain
	 * Return the list of subresultants
	 *
	 * @param q: The other sparse univariate polynomial
	 **/
	inline std::vector< Derived > subresultantChain (const Derived& q, int filled=0) const {
		if (name != q.name) {
			std::cout << "BPAS: error, trying to compute subresultant chains between Ring[" << name << "] and Ring[" << q.name << "]." << std::endl;
			exit(1);
		}

		if (degree() == 0 || q.degree() == 0){
			std::cout << "BPAS: error, Input polynomials to subresultantChain must have positive degree." << std::endl;
			exit(1);
		}

		std::vector< Derived > S;
		Derived a, b;
		if (q.degree() > degree()) {
			a = q;
			b = *this;
		}
		else {
			a = *this;
			b = q;
		}

		int k = (a.degree() - b.degree()).get_si();
		Ring s = b.leadingCoefficient() ^ k;

		Derived A = b, B = a, C = -b;
		if (k > 1) {
			b *= b.leadingCoefficient()^(k-1); // converts b to S_{degree(b)}
		}
		S.push_back(b);
		S.push_back(a);
		B.pseudoDivide(C);
		Integer delta = 0;
		while (true) {
		    if (B.isZero())
				break;
			S.insert(S.begin(), B);
			delta = A.degree() - B.degree();
			if (delta > 1) {
			    C = LazardSe(S[1], S[0], s);
				S.insert(S.begin(), C);
			}
			else { C = B; }
			if (B.degree() == 0)
			    break;
			B = DucosSem(A, B, C, s);
			A = C;
			s = A.leadingCoefficient();
		}
		// if resultant is 0, add it to subresultantChain
		if (S.at(0).degree() > 0) {
			S.insert(S.begin(), B);
		}
		if (filled) {
			std::cerr << "filling chain..." << std::endl;
			fillChain(S);
		}
		return S;
	}

	/**
 	 * monomialBasisSubResultantChain
	 *
	 * @param q: The other sparse univariate polynomial
	 **/
	inline std::vector<Derived > monomialBasisSubresultantChain(const Derived& q) {
		std::vector< Derived > s = this->subresultantChain(q);
		Derived sup;
		int delta,n;
		for (int i=s.size()-2; i>0; --i) {
			delta = s.at(i).degree() - s.at(i-1).degree();
			if (delta > 1) {
				if (i == 1 && s.at(i-1).isZero())
					n = delta-1;
				else
					n = delta-2;
				for (int j=0; j<n; j++)
					s.insert(s.begin()+i-1,sup);
			}
		}
		return s;
		/*int delta;
		Derived pp;
		std::vector< Derived >    src;
		if(this->degree()>q.degree()){
		src.push_back(*this);
		src.push_back(q);
		delta = q.degree();
		pp = q;
		}
		else{
		src.push_back(q);
		src.push_back(*this);
		delta = this->degree();
		pp = *this;
		}


		int i = s.size()-2;



		Derived qq = s[s.size()-1];
		src.push_back(qq);
		delta= delta - qq.degree();
		while(true){

		if(delta>=2){

		for(int j=0;j<delta-2;j++){
		Derived poly;
		poly.zero();
		src.push_back(poly);
		}

		}

		if(i==0){
		src.push_back(src[0]);
		std::reverse(src.begin(),src.end());
		return src;

		}
		if(i<0){
		std::reverse(src.begin(),src.end());
		return src;

		}


		pp= s[i];
		src.push_back(pp);

		qq = s[i-1];
		src.push_back(qq);

		delta  = pp.degree() - qq.degree();
		i = i-2;




		}*/

	}




	/**
	 * Resultant
	 *
	 * @param q: The other sparse univariate polynomial
	 **/
	inline Derived resultant (const Derived& q) {
		std::vector< Derived > s = subresultantChain(q);
		return s[0];
	}

	/**
	 * GCD(p, q)
	 *
	 * @param q: The other polynomial
	 **/
	inline Derived gcd (const Derived& q) const {
		if (isZero()) { return q; }
		if (q.isZero()) { return *this; }
		if (name != q.name) {
			std::cout << "BPAS: error, trying to compute GCD between Ring[" << name << "] and Ring[" << q.name << "]." << std::endl;
			exit(1);
		}

		Derived a(*this), b(q);
		if (a.degree() == 0 || b.degree() == 0) {
			a.one();
			return a;
		}

		Derived r;
		r.name = name;

		//std::cout << "is_same: " << std::is_same<Ring, RationalNumber>::value << std::endl;
		Ring rng;
		//TODO properly implement the ring properties.
		// if (Ring::properties.has(PRIME_FIELD)) {
		// 	DenseUnivariateRationalPolynomial f = a.convertToDUQP();
		// 	DenseUnivariateRationalPolynomial g = b.convertToDUQP();
		// 	DenseUnivariateRationalPolynomial z = f.gcd(g);
		// 	r = Derived (z);
		// }
		// else {
			Ring ca, cb, cr;
			ca = a.content();
			a /= ca;
			cb = b.content();
			b /= cb;
			std::vector< Derived > R = a.subresultantChain(b);

			r.setCoefficient(0, ca.gcd(cb));
			//r *= cb;
			int n = R.size();
			bool isZero = 0;
			if (n) {
				isZero = 1;
				for (int i = 0; i < n; ++i) {
					if (!R[i].isZero()) {
						cr = R[i].content();
						R[i] /= cr;
						r *= R[i];
						isZero = 0;
						break;
					}
				}
			}
			if (isZero) {
				if (a.degree() <= b.degree()) { r *= a; }
				else { r *= b; }
			}
		// }
		return r;
	}

	/**
	 * Square free
	 *
	 * @param
	 **/
	inline Factors<Derived> squareFree() const {
		std::vector< Derived > sf;
		int d = terms.size()-1;
		if (!terms[d].exp)
			sf.push_back(*this);
		else if (terms[d].exp == 1) {
			Derived t;
			t.name = name;
			t += terms[d].coef;
			sf.push_back(t);
			t = *this / terms[d].coef;
			sf.push_back(t);
		}
		else {
			Derived a (*this), b(*this);
			b.differentiate(1);
			Derived g = a.gcd(b);
			g /= g.content();
			Derived x = a / g;
			Derived y = b / g;
			Derived z = -x;
			z.differentiate(1);
			z += y;

			while (!z.isZero()) {
				g = x.gcd(z);
				g /= g.content();
				sf.push_back(g);
				x /= g;
				y = z / g;
				z = -x;
				z.differentiate(1);
				z += y;
			}
			sf.push_back(x);

			Ring e;
			e.one();
			for (int i = 0; i < sf.size(); ++i) {
				e *= sf[i].leadingCoefficient();
				sf[i] /= sf[i].leadingCoefficient();
			}
			Derived t;
			t.name = name;
			t += e;
			sf.insert(sf.begin(), t);
		}

		Factors<Derived> f;
		f.setRingElement(sf[0]);
		for (int i = 1; i < sf.size(); ++i) {
			f.addFactor(sf[i], i);
		}

		return f;
	}

	/**
	 * Overload stream operator <<
	 *
	 * @param out: Stream object
	 * @param b: The univariate polynomial
	 **/
	inline void print (std::ostream &out) const {
		int n = terms.size();
		if (!n) { out << "0"; }
		for (int i = 0; i < n; ++i) {
			if (this->terms[i].exp) {
				if (this->terms[i].coef.isNegativeOne())
					out << "-";
				else if (i && this->terms[i].coef.isConstant() >= 0)
					out << "+";
				if (!this->terms[i].coef.isConstant())
					out << "(" << this->terms[i].coef << ")*";
				else if (!this->terms[i].coef.isOne() && !this->terms[i].coef.isNegativeOne())
					out << this->terms[i].coef << "*";
				out << this->name;
				if (this->terms[i].exp > 1)
					out << "^" << this->terms[i].exp;
			}
			else {
				if (this->terms[i].coef.isConstant()) { out << this->terms[i].coef; }
				else { out << "(" << this->terms[i].coef << ")"; }
			}
		}
	}

	inline ExpressionTree convertToExpressionTree() const {
		//TODO
		std::cerr << "BPAS ERROR: SMP<Ring>::convertToExpressionTree NOT YET IMPLEMENTED" << std::endl;
		return ExpressionTree();
	}

	inline DenseUnivariateRationalPolynomial convertToDUQP() {
		bool isDense = 1;
		int k = 0, n = terms.size(), d = 0;
		if (n) { d = terms[n-1].exp; }
		DenseUnivariateRationalPolynomial res(d+1);
		res.setVariableName(name);
		for (int i = 0; i <= d; ++i) {
			if (k < n) {
				if (!terms[k].coef.isConstant()) {
					isDense = 0;
					break;
				}
				else if (terms[k].exp == i) {
					res.setCoefficient(i, RationalNumber(terms[k].coef));
					k++;
				}
			}
		}
		if (!isDense) { res.zero(); }
		return res;
	}

	inline DenseUnivariateIntegerPolynomial convertToDUZP() {
		bool isDense = 1;
		int k = 0, n = terms.size(), d = 0;
		if (n) { d = terms[n-1].exp; }
		DenseUnivariateIntegerPolynomial res(d+1);
		res.setVariableName(name);
		for (int i = 0; i <= d; ++i) {
			if (k < n) {
				if (!terms[k].coef.isConstant()) {
					isDense = 0;
					break;
				}
				else if (terms[k].exp == i) {
					res.setCoefficient(i, Integer(terms[k].coef));
					k++;
				}
			}
		}
		if (!isDense) { res.zero(); }
		return res;
	}







};


template <class Field>
class SparseUnivariateTempFieldPoly : public virtual SparseUnivariateTempPoly<Field> {

	Integer euclideanSize() const {
		return this->degree();
	}

	Derived euclideanDivision(
		const Derived& b,
		Derived* q = NULL) const
	{
		Field lc = b.leadingCoefficient();
		Derived monicb = b * lc.inverse();
		if (q != NULL) {
			Derived rem;
			*q = this->monicDivide(monicb, &rem);
			*q *= lc;
			return rem;
		} else {
			Derived rem = *this;
			return rem.monicDivide(b);
		}

	}

	Derived extendedEuclidean(
		const Derived& b,
		Derived* s,
		Derived* t) const
	{
		//TODO NOT YET IMPLEMENTED

		return *this;
	}


	Derived quotient(const Derived& b) const {
		Derived q;
		this->euclideanDivision(b, &q);
		return q;
	}

	Derived remainder(const Derived& b) const {
		return this->euclideanDivision(b);
	}

	Derived operator%(const Derived& b) const {
		Derived ret = *this;
		ret %= b;
		return ret;
	}

	Derived& operator%=(const Derived& b) {
		*this = this->remainder(b);
		return dynamic_cast<Derived&>(*this);
	}
};

template <class Ring>
class SparseUnivariatePolynomial : public std::conditional<std::is_base_of<BPASField, Ring>::value, SparseUnivariateTempFieldPoly<Ring>, SparseUnivariateTempPoly<Ring>> >::type {

public:
	typedef SparseUnivariateTempPoly<Ring, SparseUnivariatePolynomial<Ring>> Base;

	SparseUnivariatePolynomial<Ring>() : Base() {}

	SparseUnivariatePolynomial<Ring>& operator=(const SparseUnivariatePolynomial<Ring>& other) {
		Base::operator=(other);
		return *this;
	}

	SparseUnivariatePolynomial<Ring>& operator=(SparseUnivariatePolynomial<Ring>&& other) {
		Base::operator=(other);
		return *this;
	}

	SparseUnivariatePolynomial<Ring>& operator=(const Ring& r) {
		Base::operator=(r);
		return *this;
	}

	SparseUnivariatePolynomial<Ring>(const SparseUnivariatePolynomial<Ring>& other) : Base(other) {}
	SparseUnivariatePolynomial<Ring>(SparseUnivariatePolynomial<Ring>&& other) : Base(other) {}

	SparseUnivariatePolynomial<Ring>(int a) : Base(a) {}


	SparseUnivariatePolynomial<Ring>(const Integer& a) : Base(a) {}

	SparseUnivariatePolynomial<Ring>(const RationalNumber& a) : Base(a) {}

	SparseUnivariatePolynomial<Ring>(const ComplexRationalNumber& a) : Base(a) {}

	SparseUnivariatePolynomial<Ring>(const DenseUnivariateIntegerPolynomial& b) :  Base(b) {}

	SparseUnivariatePolynomial<Ring> (const DenseUnivariateRationalPolynomial& b)  : Base(b) {}


	SparseUnivariatePolynomial<Ring> (Symbol sym) : Base(sym) {}

	SparseUnivariatePolynomial<Ring> (const Base& b) : Base(b) {}

	/**
	 * Destroy the polynomial
 	 *
     * @param
     **/
	~SparseUnivariatePolynomial<Ring> () {}


};




//TODO: Develop random element generator for all BPASRings so that this can be placed in SUP<Ring>
/**
 * Generate random polynomial
 *
 * @param n: degree of the random polynomial
 * @param sparsity: the proportion of non-zero elements
 * @param bits: maximum number of bits for each coefficient
 **/
static SparseUnivariatePolynomial<RationalNumber> randomSUPQPolynomial(int n, double sparsity, unsigned long int coefBound){
	int k;
	mpq_t randVal;
	rand_mpq_t(coefBound,1,randVal);
	mpq_class coef(randVal);
	SparseUnivariatePolynomial<RationalNumber> P;
	P.setCoefficient(n, coef);
	int nTerms = ceil(sparsity*n);
	int index;
	for(int i = 0; i < nTerms; i++) {
		// Set random coefficients with sparsity
		index = rand() % n;
		rand_mpq_t(coefBound,1,randVal);
		coef = mpq_class(randVal);
		P.setCoefficient(index, coef);
	}

	return P;
}




















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


