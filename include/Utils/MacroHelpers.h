
#ifndef _MACROHELPERS_H_
#define _MACROHELPERS_H_


/*
 * A macro which SAFELY casts the variable X
 * of type T** to T const*const*. 
 * If the conversion cannot be performed safely then
 * a compiler error is thrown.
 *
 * Note that this conversion is automatic in C++.
 *
 * How it works: 
 *   - (T const*){ X[0] } determines if the implicit conversion
 *     of *X to the type T const* is valid using a compound literal.
 *   - All of that inside of sizeof ensures that the 
 *     compound literal isn't actually created.
 *   - The comma operator is then used to evaluate the actual cast
 *     and return it. 
 * 
 */
#define CONSTCONSTCAST(T, X) (                 \
   (void)sizeof((T const*){ (X)[0] }), \
   (T const*const*)(X)                 \
)



#endif