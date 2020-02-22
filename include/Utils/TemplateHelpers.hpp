
#ifndef _TEMPLATE_HELPERS_H_
#define _TEMPLATE_HELPERS_H_

/**
 * Template to enforce that a template paramater T is a dervied class of 
 * (or the class itself) T.
 *
 * To use, make a template class inherit Derived_from and pass it the template
 * template paramater as T and the restricting class as B.
 */
template<class T, class B> struct Derived_from {
	static void constraints(T* p) { 
		B* pb = p;
	}
	Derived_from() {
		void(*p)(T*) = constraints;
	}
};

#endif