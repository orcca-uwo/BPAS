

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
class Factor : public std::pair<Ring, int>, public ExpressionTreeConvert, private Derived_from<Ring, BPASRing<Ring>> {

public:
	/**
	 * Construct an empty factor.
	 */
	Factor() : std::pair<Ring,int>() {}

	/**
	 * Construct a factor from a Ring element and an exponent.
	 * @param r the Ring element
	 * @param e the exponent
	 */
	Factor(const Ring& r, int e) : std::pair<Ring,int>(r,e) {}

	/**
	 * Move-construct a factor from a Ring element and an exponent.
	 * @param r the Ring element
	 * @param e the exponent
	 */
	Factor(Ring&& r, int e) : std::pair<Ring,int>(r,e) {}

	/**
	 * Construct a Factor from a std::pair of a Ring element and integer exponent.
	 * @param p the pair
	 */
	Factor(const std::pair<Ring, int>& p) : std::pair<Ring, int>(p) {}

	/** 
	 * Equality comparison operator.
	 * @param f the Factor to compare
	 * @return true iff *this == f.
	 */
	bool operator==(const Factor<Ring>& f) {
		return (this->first == f.first && this->second == f.second);
	}

	/** 
	 * Inequality comparison operator.
	 * @param f the Factor to compare
	 * @return true iff *this != f.
	 */
	bool operator!=(const Factor<Ring>& f) {
		return (this->first != f.first || this->second != f.second);
	}

	/**
     * Output operator, defines a to string conversion. 
     * @param out the output stream
     * @param f the Factor
     */
    friend std::ostream& operator<<(std::ostream& out, const Factor<Ring>& f) {
		out << "[" << f.first << ", " << f.second << "]";
		return out;

	}

	/**
	 * Convert the Factor to an ExpressionTree.
	 * @return the ExpressionTree
	 */
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
class Factors : public ExpressionTreeConvert, private Derived_from<Ring, BPASRing<Ring>> {
	
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
	 * @param r the Ring element
	 */
	Factors(const Ring& r);

	/**
	 * Create a Factors from a list of Ring elements. Each factor is assumed
	 * to have an exponent of 1. 
	 * @param v the vector of Ring elements.
	 */
	Factors(const std::vector<Ring>& v);
	
	/**
	 * Create a Factors from a list of Ring elements and a parallel array of integers
	 * as exponents for each ring element.
	 * @param v the vector of Ring elements
	 * @param e the parallel vector of integer exponents
	 */
	Factors(const std::vector<Ring>& v, const std::vector<int>& e);

	/**
	 * Create a Factors from a list of ring elements, a parallel array of integers
	 * as exponents for each ring element, and a ring element representing the
	 * invertible unit of the Ring. 
	 * @param v the vector of Ring elements
	 * @param e the parallel vector of integer exponents
	 * @param u the multicative factor
	 */
	Factors(const std::vector<Ring>& v, const std::vector<int>& e, const Ring& u);

	/**
	 * Create a Factors from a vector of Factor elements.
	 * @param v the vector of Factor elements.
	 */
	Factors(const std::vector<Factor<Ring>>& v);

	/**
	 * Create a Factors from a vector of Factor elements and an invertible Ring element.
	 * @param v the vector of Factor elements.
	 * @param u the invertible Ring element.
	 */
	Factors(const std::vector<Factor<Ring>>& v, const Ring& u);

	/**
	 * Copy consturctor.
	 * @param f the Factors to copy.
	 */
	Factors(const Factors& f);
	
	/**
	 * Move consturctor.
	 * @param f the Factors to move.
	 */
	Factors(Factors&& f);

	/**
	 * Destructor.
	 */
	~Factors();

	/**
	 * Set the invertible Ring element.
	 * @param r the new invertible Ring element.
	 */
	inline void setRingElement(const Ring& r) {
		u = r;
	}

	/**
	 * Get the invertible Ring element.
	 * @return the invertible Ring element of this Factors
	 */
	inline Ring ringElement() const {
		return u;
	}

	/**
	 * Multiply this Factors invertible Ring element by another Ring element.
	 * @param r the other Ring element.
	 */
	inline void multiplyRingElement(const Ring& r) {
		u *= r;
	}

	/**
	 * Get the list of Factor elements.
	 * @return the vector of Factor elements.
	 */
	inline std::vector<Factor<Ring>> factors() const {
		return facts;
	}
	
	/**
	 * Get the i'th Factor.
	 * @param i the index
	 * @return the Factor at index i
	 */
	inline Factor<Ring> factor(int i) const {
		if (i < facts.size())
			return facts[i];
		else {
			std::cerr << "BPAS: error, i exceeds array bounds of list of factors" << std::endl;
			exit(1);
		}
	}

	/**
	 * Set the list of Factor elements.
	 * @param v the vector of Factor elements.
	 */
	inline void setFactors(const std::vector<Factor<Ring>>& v) {
		facts = v;
	}

	/**
	 * A Factor to the list of Factor elements.
	 * @param f the Factor to add.
	 */
	inline void addFactor(const Factor<Ring>& f) {
		if (f.first.isOne()) {
			return;
		}
		for (int i = 0; i < facts.size(); ++i) {
			if (f.first == facts[i].first) {
				facts[i].second += f.second;
				return;
			}
		}
		facts.push_back(f);
	}

	/**
	 * Add a Ring element and its corresponding exponent as a Factor.
	 * @param r the ring element.
	 * @param e the exponent
	 */
	inline void addFactor(const Ring& r, int e) {
		if (r.isOne()) {
			return;
		}
		for (int i = 0; i < facts.size(); ++i) {
			if (r == facts[i].first) {
				facts[i].second += e;
				return;
			}
		}
		facts.emplace_back(r, e);
	}


	/**
	 * Add a Ring element and its corresponding exponent as a Factor.
	 * @param r the ring element.
	 * @param e the exponent
	 */
	inline void addFactor(Ring&& r, int e) {
		if (r.isOne()) {
			return;
		}
		for (int i = 0; i < facts.size(); ++i) {
			if (r == facts[i].first) {
				facts[i].second += e;
				return;
			}
		}
		facts.emplace_back(r, e);
	}

	/**
	 * Combine another Factors' list of Factor elements with this Factors. This does
	 * not modify the invertible ring element.
	 * @param f the Factors object
	 */
	inline void addFactors(const Factors<Ring> f) {
		for (int i = 0; i < f.facts.size(); ++i) {
			this->addFactor(f.facts[i]);
		}
	}

	/**
	 * Get the number of Factor elements in this Factors object.
	 * @return the number of Factor elements.
	 */
	inline size_t size() const {
		return facts.size();
	}

	/**
	 * Copy assignment.
	 * @param f the Factors to copy from.
	 */
	Factors<Ring>& operator=(const Factors<Ring>& f);
	
	/**
	 * Move assignment.
	 * @param f the Factors to move from.
	 */
	Factors<Ring>& operator=(Factors<Ring>&& f);

	/**
	 * Equality testing for Factors. Does NOT check if, when factors are distributed,
	 * that the underlying elements are the same. Only if the factoirzation is the same.
	 * @param f the Factors to test equality with.
	 * @return true iff all factors are equal in both Factors.
	 */
	bool operator==(const Factors<Ring>& f) const;

	/**
	 * Inequality testing for Factors. Does NOT check if, when factors are distributed,
	 * that the underlying elements are the same. Only if the factoirzation is the same.
	 * @param f the Factors to test equality with.
	 * @return false iff all factors are equal in both Factors.
	 */
	inline bool operator!=(const Factors<Ring>& f) const {
		return !(*this == f);
	}

	/**
	 * Indexing operator.
	 * @param idx the index 
	 * @return the Factor at index idx.
	 */
	inline Factor<Ring>& operator[](size_t idx) {
		return facts[idx];
	}

	/**
     * Output operator. Defines a to string conversion. 
     * @param out the output stream
     * @param f the Factors to output
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

	/**
	 * Convert the Factors object to an ExpressionTree.
	 * @return the ExpressionTree encoding the Factors.
	 */ 
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
