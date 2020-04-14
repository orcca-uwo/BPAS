

#ifndef _BPAS_REGCHAIN_H_
#define _BPAS_REGCHAIN_H_

/**
 * An abstract class defining the interface of a regular chain.
 * See also BPASTriangularSet.
 */
template <class Field, class RecursiveFieldPoly>
class BPASRegularChain : public virtual BPASTriangularSet<Field,RecursiveFieldPoly>
{
	public:

		virtual BPASRegularChain<Field,RecursiveFieldPoly>& operator= (const BPASRegularChain<Field,RecursiveFieldPoly>&) = 0;
		virtual BPASRegularChain<Field,RecursiveFieldPoly>& operator= (BPASRegularChain<Field,RecursiveFieldPoly>&&) = 0;
//		virtual RecursiveFieldPoly normalize (const RecursiveFieldPoly&) = 0;
};

#endif