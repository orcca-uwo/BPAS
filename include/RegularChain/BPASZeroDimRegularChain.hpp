

#ifndef _BPAS_ZERODIM_REGCHAIN_H_
#define _BPAS_ZERODIM_REGCHAIN_H_

#include "BPASRegularChain.hpp"

/**
 * An abstract class defining the interface of a zero-dimensional regular chain.
 */
template <class Field, class RecursiveFieldPoly>
class BPASZeroDimensionalRegularChain : public virtual BPASRegularChain<Field,RecursiveFieldPoly>
{
	public:

		virtual BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>& operator= (const BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>&) = 0;
		virtual BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>& operator= (BPASZeroDimensionalRegularChain<Field,RecursiveFieldPoly>&&) = 0;
};

#endif