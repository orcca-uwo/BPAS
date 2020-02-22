

#ifndef _FACTORS_HPP_
#define _FACTORS_HPP_

#include <iostream>
#include <utility>
#include "../Ring/BPASRing.hpp"
#include "../ExpressionTree/ExpressionTree.hpp"
#include "../Utils/TemplateHelpers.hpp"

/**
 * A Factor is a pair of a BPASRing element and an integer exponent.
 * Template Ring should derive from BPASRing.
 */
template <class Ring>
class Factor : public std::pair<Ring, int>, public ExpressionTreeConvert {

public:
	Factor() : std::pair<Ring,int>() {}

	Factor(const Ring& r, int e) : std::pair<Ring,int>(r,e) {}

	Factor(Ring&& r, int e) : std::pair<Ring,int>(r,e) {}

	Factor(const std::pair<Ring, int>& p) : std::pair<Ring, int>(p) {}

	bool operator==(const Factor<Ring>& f) {
		return (this->first == f.first && this->second == f.second);
	}

	bool operator!=(const Factor<Ring>& f) {
		return (this->first != f.first || this->second != f.second);
	}

	/**
     * Output operator. 
     *
     * Defines a to string conversion. 
     */
    friend std::ostream& operator<<(std::ostream& out, const Factor<Ring>& f) {
		out << "[" << f.first << ", " << f.second << "]";
		return out;

	}

	ExpressionTree convertToExpressionTree() const {
		std::vector<ExpressionTree> trees;
		trees.push_back(this->first.convertToExpressionTree());
		trees.emplace_back(new ExprTreeNode(this->second));

		return ExpressionTree(trees);
	}
};

/**
 * A simple data structure for encapsulating a collection of Factor elements.
 * Factors are encoded as a list of Factor items. That is, a pair of a Ring 
 * element and an integer representing its exponent. This list of Factor elements
 * is augmented by an invertible element of the ring (if one such exists) as a constant
 * multiplicative factor.  
 */
template <class Ring>
class Factors : public ExpressionTreeConvert {
	
private: 
	Ring u;
	std::vector<Factor<Ring>> facts;

public:
	
	/**
	 * Default constructor. Constructs an empty list of factors.
	 */
	Factors();

	/**
	 * Construct a Factors with a single factor. 
	 */
	Factors(const Ring& r);

	/**
	 * Create a Factors from a list of ring elements. Each factor is assumed
	 * to have an exponent of 1. 
	 */
	Factors(const std::vector<Ring>& v);
	
	/**
	 * Create a Factors from a list of ring elements and a parallel array of integers
	 * as exponents for each ring element.
	 */
	Factors(const std::vector<Ring>& v, const std::vector<int>& e);

	/**
	 * Create a Factors from a list of ring elements, a parallel array of integers
	 * as exponents for each ring element, and a ring element representing the
	 * invertible unit of the Ring. 
	 */
	Factors(const std::vector<Ring>& v, const std::vector<int>& e, const Ring& u);

	/**
	 * Create a Factors from a vector of Factor elements.
	 */
	Factors(const std::vector<Factor<Ring>>& v);

	/**
	 * Create a Factors from a vector of Factor elements.
	 */
	Factors(const std::vector<Factor<Ring>>& v, const Ring& u);

	/**
	 * Copy consturctor.
	 */
	Factors(const Factors& f);
	
	/**
	 * Move consturctor.
	 */
	Factors(Factors&& f);

	/**
	 * Destructor;
	 */
	~Factors();

	inline void setRingElement(const Ring& r) {
		u = r;
	}

	inline Ring ringElement() const {
		return u;
	}

	inline void multiplyRingElement(const Ring& r) {
		u *= r;
	}

	inline std::vector<Factor<Ring>> factors() const {
		return facts;
	}
	
	inline Factor<Ring> factor(int i) const {
		if (i < facts.size())
			return facts[i];
		else {
			std::cerr << "BPAS: error, i exceeds array bounds of list of factors" << std::endl;
			exit(1);
		}
	}

	inline void setFactors(const std::vector<Factor<Ring>>& v) {
		facts = v;
	}

	inline void addFactor(const Factor<Ring>& f) {
		facts.push_back(f);
	}

	inline void addFactor(const Ring& r, int e) {
		facts.emplace_back(r, e);
	}

	inline void addFactors(const Factors<Ring> f) {
		facts.insert(facts.end(), f.facts.begin(), f.facts.end());
	}

	inline size_t size() const {
		return facts.size();
	}

	/**
	 * Copy assignment.
	 */
	Factors<Ring>& operator=(const Factors<Ring>& f);
	
	/**
	 * Move assignment.
	 */
	Factors<Ring>& operator=(Factors<Ring>&& f);

	/**
	 * Equality testing for Factors. Does NOT check if, when factors are distributed,
	 * that the underlying elements are the same. Only if the factoirzation is the same.
	 * returns true iff all factors are equal in both Factors.
	 */
	bool operator==(const Factors<Ring>& f) const;

	/**
	 * Inequality testing for Factors. Does NOT check if, when factors are distributed,
	 * that the underlying elements are the same. Only if the factoirzation is the same.
	 * returns false iff all factors are equal in both Factors.
	 */
	inline bool operator!=(const Factors<Ring>& f) const {
		return !(*this == f);
	}

	inline const Factor<Ring>& operator[](size_t idx) const {
		return facts[idx];
	}

	/**
     * Output operator. 
     *
     * Defines a to string conversion. 
     */
    friend std::ostream& operator<<(std::ostream& out, const Factors<Ring>& f) {
		out << "[" << f.u << ", [";


		for (int i = 0; i < f.facts.size(); ++i) {
			out << f.facts[i];
			if (i + 1 < f.facts.size()) {
				out << ", ";
			}
		}

		out << "]]";
		return out;
	}

	ExpressionTree convertToExpressionTree() const {
		std::vector<ExpressionTree> trees;
		ExpressionTree t;
		t.fromVector(facts);

		trees.push_back(u.convertToExpressionTree());
		trees.push_back(t);

		return ExpressionTree(trees);
	}

};	


//Include the implementation
#include "Factors_impl.hxx"

#endif
