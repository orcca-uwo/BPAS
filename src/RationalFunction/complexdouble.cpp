#include "../../include/RationalFunction/complexdouble.h"

std::complex<double> CRN_to_complex_double(ComplexRationalNumber &a) {
	double rp,ip;
	rp = a.realPart().get_mpq().get_d();
	ip = a.imaginaryPart().get_mpq().get_d();
	std::complex<double> out(rp,ip);
	return out;
}

void complexAdd(double a_rp, double a_ip, double b_rp, double b_ip, double *out_rp, double *out_ip) {
	*out_rp = a_rp + b_rp;
	*out_ip = a_ip + b_ip;
}

void complexMultiply(double a_rp, double a_ip, double b_rp, double b_ip, double *out_rp, double *out_ip) {
	double rp;
	double ip;
	rp = a_rp*b_rp - a_ip*b_ip;
	ip = a_rp*b_ip + a_ip*b_rp;
	*out_rp = rp;
	*out_ip = ip;
}

void complexDivide(double a_rp, double a_ip, double b_rp, double b_ip, double *out_rp, double *out_ip) {
	double rp;
	double ip;
	double modb;
	modb = b_rp*b_rp + b_ip*b_ip;
	rp = a_rp*b_rp + a_ip*b_ip;
	ip = -a_rp*b_ip + a_ip*b_rp;
	rp /= modb;
	ip /= modb;
	*out_rp = rp;
	*out_ip = ip;
}
