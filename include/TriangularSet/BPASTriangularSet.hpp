

#ifndef _BPAS_TRISET_H_
#define _BPAS_TRISET_H_

#include "../Polynomial/BPASRecursivePolynomial.hpp"

/**
 * An abstract class defining the interface of a triangular set.
 * A BPASTriangularSet is templated by a BPASRecursivelyViewedPolynomial with coefficients
 * in a BPASField.
 */
template <class Field, class RecursiveFieldPoly> // TODO: add Derived
class BPASTriangularSet :
	private Derived_from<Field, BPASField<Field>>,
	private Derived_from<RecursiveFieldPoly, BPASRecursivelyViewedPolynomial<Field,RecursiveFieldPoly>>
{
	public:

		virtual BPASTriangularSet<Field,RecursiveFieldPoly>& operator= (const BPASTriangularSet<Field,RecursiveFieldPoly>&) = 0;
		virtual BPASTriangularSet<Field,RecursiveFieldPoly>& operator= (BPASTriangularSet<Field,RecursiveFieldPoly>&&) = 0;
//		virtual void copy(const BPASTriangularSet<Field,RecursiveFieldPoly>&) = 0;
		virtual int numberOfVariables() const = 0;
		virtual std::vector<Symbol> variables() const = 0;

		virtual RecursiveFieldPoly select(const Symbol&) const = 0;
		virtual void lower(const Symbol&, BPASTriangularSet<Field,RecursiveFieldPoly>&) const = 0;
		virtual void upper(const Symbol&, BPASTriangularSet<Field,RecursiveFieldPoly>&) const = 0;
		virtual RecursiveFieldPoly pseudoDivide (const RecursiveFieldPoly&, std::vector<RecursiveFieldPoly>*, RecursiveFieldPoly*) const = 0;
		virtual RecursiveFieldPoly normalForm (const RecursiveFieldPoly&, std::vector<RecursiveFieldPoly>*) const = 0;
};


#endif