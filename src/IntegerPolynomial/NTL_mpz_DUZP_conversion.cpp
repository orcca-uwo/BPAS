#include "IntegerPolynomial/NTL_mpz_DUZP_conversion.hpp"
#include <sstream>

void NTLZZX_set_DUZP_t (NTL::ZZX& g, const DUZP_t* f) {

	if (isZero_DUZP(f)) {
		SetCoeff (g, 0, NTL::ZZ(0));
		return;
	}
	
	int deg = f->lt;
	
	g = NTL::ZZX (deg);
	
	char* cc;
	NTL::ZZ nc;
   
	for (int i = 0; i <= deg; i+=1) {
		
		cc = (char*) malloc(sizeof(char)*(mpz_sizeinbase(f->coefs[i], 10) + 2));
		cc = mpz_get_str(cc,10,f->coefs[i]);
		conv (nc, cc);
		SetCoeff (g, i, nc);	
		free(cc);
		
	}
}

void DUZP_t_set_NTLZZX (DUZP_t** f, const NTL::ZZX& g) {

	if (IsZero (g)) {
		freePolynomial_DUZP(*f);
		*f = makeConstPolynomial_DUZP(1,0);
		return;
	}
	
	int d = deg (g);

	DUZP_t* ff = makePolynomial_DUZP(d+1);
	for (int i=1; i<=d; ++i) // these three lines should be incorporated into a routine in DUZP_Support
		mpz_init(ff->coefs[i]);
	ff->lt = d;
	
	NTL::ZZ nc;
	std::stringstream ssc;
	mpz_t mpc;
	mpz_init (mpc);
	
	for (int i = 0; i <= d; i++) {
		
		GetCoeff (nc, g, i);
		ssc << nc;
		mpz_set_str(mpc, ssc.str().c_str(), 10); // base = 10
		
		mpz_set(ff->coefs[i],mpc); // f->coefs[i] = g[i]
		
		ssc.str("");
	}

	mpz_clear(mpc);
	freePolynomial_DUZP(*f);
	*f = ff;
}

void NTLZZ_set_mpz_class (NTL::ZZ& g, const mpz_class& f) {

	std::string sc;
	
	sc = f.get_str (10); // base = 10
	const char* cc = sc.c_str();
	conv (g, cc);
}

void NTLZZ_set_mpz_t (NTL::ZZ& g, const mpz_t f) {

	char* cc;
	cc = (char*) malloc(sizeof(char)*(mpz_sizeinbase(f, 10) + 2));
	cc = mpz_get_str(cc,10,f);
	conv (g, cc);
	free(cc);
}

void mpz_class_set_NTLZZ (mpz_class& f, const NTL::ZZ& g) {

	std::stringstream ssc;
	mpz_t mpc;
	
	ssc << g;
	mpz_init (mpc);
	mpz_set_str(mpc, ssc.str().c_str(), 10); // base = 10

	f = mpz_class (mpc);
	
	mpz_clear(mpc);
	ssc.str("");
}

void mpz_t_set_NTLZZ (mpz_t& f, const NTL::ZZ& g) {

	mpz_clear(f);
	std::stringstream ssc;
	
	ssc << g;
	mpz_init (f);
	mpz_set_str(f, ssc.str().c_str(), 10); // base = 10
	
	ssc.str("");
}
