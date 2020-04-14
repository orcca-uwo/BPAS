#ifndef _COMPLEXDOUBLE
#define _COMPLEXDOUBLE
#include <complex>
#include "../ring.h"

std::complex<double> CRN_to_complex_double(ComplexRationalNumber &a);

void complexAdd(double a_rp, double a_ip, double b_rp, double b_ip, double *out_rp, double *out_ip);

void complexMultiply(double a_rp, double a_ip, double b_rp, double b_ip, double *out_rp, double *out_ip);

void complexDivide(double a_rp, double a_ip, double b_rp, double b_ip, double *out_rp, double *out_ip);

#endif
