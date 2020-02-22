
#ifndef _FACTORS_HXX_
#define _FACTORS_HXX_

/**
 * A simple data structure for encapsulating factors of a ring element.
 * Factors are encoded as a list of Factor items. That is, a pair of a Ring 
 * element and an integer representing its exponent. This list of Factor elemements
 * is augmented by an invertible element of the ring (if one such exists) as a constant
 * multiplicative factor.  
 */
	
/**
 * Default constructor. Constructs an empty list of factors.
 */
template <class Ring>
Factors<Ring>::Factors () {
	u.one();
}

/**
 * Construct a Factors with a single factor. 
 */
template <class Ring>
Factors<Ring>::Factors(const Ring& r) {
	u.one();
	facts.emplace_back(r, 1);
}

/**
 * Create a Factors from a list of ring elements. Each factor is assumed
 * to have an exponent of 1. 
 */
template <class Ring>
Factors<Ring>::Factors (const std::vector<Ring>& v) {
	u.one();
	for (int i = 0; i < v.size(); ++i) {
		facts.emplace_back(v[i], 1);
	}
}

/**
 * Create a Factors from a list of ring elements and a parallel array of integers
 * as exponents for each ring element.
 */
template <class Ring>
Factors<Ring>::Factors (const std::vector<Ring>& v, const std::vector<int>& e) {
	u.one();
	for (int i = 0; i < v.size() && i < e.size(); ++i) {
		facts.emplace_back(v[i], e[i]);
	}
}

/**
 * Create a Factors from a list of ring elements, a parallel array of integers
 * as exponents for each ring element, and a ring element representing the common
 * factor among all factors, usually a unit. 
 */
template <class Ring>
Factors<Ring>::Factors (const std::vector<Ring>& v, const std::vector<int>& e, const Ring& u) {
	this->u = u;
	for (int i = 0; i < v.size() && i < e.size(); ++i) {
		facts.emplace_back(v[i], e[i]);
	}
}

/**
 * Create a Factors from a vector of Factor elements.
 */
template <class Ring>
Factors<Ring>::Factors (const std::vector<Factor<Ring>>& v) {
	u.one();
	facts = v;
}

/**
 * Create a Factors from a vector of Factor elements and a common ring element factor.
 */
template <class Ring>
Factors<Ring>::Factors (const std::vector<Factor<Ring>>& v, const Ring& u) {
	this->u = u;
	facts = v;
}

/**
 * Copy consturctor.
 */
template <class Ring>
Factors<Ring>::Factors (const Factors<Ring>& f) {
	this->u = f.u;
	facts = f.facts;
}

/**
 * Move consturctor.
 */
template <class Ring>
Factors<Ring>::Factors (Factors&& f ) {
	this->u = f.u;
	f.u.one();

	this->facts = std::move(f.facts);
}

/**
 * Destructor;
 */
template <class Ring>
Factors<Ring>::~Factors() {
	facts.clear();
}

/**
 * Copy assignment.
 */
template <class Ring>
Factors<Ring>& Factors<Ring>::operator=(const Factors<Ring>& f) {
	if (this != &f) {
		this->u = f.u;
		this->facts = f.facts;
	}
	return *this;
}

/**
 * Move assignment.
 */
template <class Ring>
Factors<Ring>& Factors<Ring>::operator=(Factors<Ring>&& f) {
	if (this != &f) {
		this->u = f.u;	
		f.u.one();

		this->facts = std::move(f.facts);
	}
	return *this;
}

/**
 * Equality testing for Factors. Does NOT check if, when factors are distributed,
 * that the underlying elements are the same. Only if the factoirzation is the same.
 * returns true iff all factors are equal in both Factors.
 */
template <class Ring>
bool Factors<Ring>::operator==(const Factors<Ring>& f) const {
	if (u != f.u) {
		return false;
	}

	if (facts.size() != f.facts.size()) {
		return false;
	}

	for (int i = 0; i < facts.size(); ++i) {
		bool found = false;
		for (int j = 0; j < f.facts.size(); ++j) {
			if (facts[i].first == f.facts[j].first &&
					facts[i].second == f.facts[j].second) {
				found = true;
				break;
			}
		}
		if (found == false) {
			return false;
		}
	}

	return true;
}

#endif