#include "IntegerPolynomial/DUZP_Support_Factoring.hpp"
#include "IntegerPolynomial/NTL_mpz_DUZP_conversion.hpp"
#include "ModularPolynomial/DUSP_NTL_Support.h"
#include "NTL/ZZXFactoring.h"
#include "NTL/LLL.h"
#include "Utils/MacroHelpers.h"

//#include <bpas.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include "Utils/Unix_Timer.h"
#include "../MapleTestTool/MapleTestTool.hpp"

polysize_t maxDeg = 5;
int coefBoundBits = 5ul;
float sparsity = 0;
int includeNeg = 1;

prime_t pr = 11;
polysize_t nfactors = 5;

char xvar[1] = {'x'};

using namespace DUZP::Factoring;

extern void naive_search_test();
extern void recombine_test();
extern void vec_mpz_test();
extern void vec_DUZP_test();
extern void vec_duspoly_test();
extern void vec_DUZP_long_test();
extern void factor_test(int deg, int n, int monic);
extern void factor_test(int* sizes, int n, int monic);
extern void factor_timing(int deg, int n, int monic);

extern void ddf_test();
extern void edf_test();

int main() {

//	naive_search_test();
//	vec_vec_test();
//	vec_mpz_test();
//	vec_DUZP_test();
//	vec_duspoly_test();
//	vec_DUZP_long_test();

	factor_timing(25,10,0);

//	factor_test(2,10,1);	
//	int sizes[2] = {5,2};
//	factor_test(sizes,2,0);
//	int sizes[5] = {1,4,2,5,4};
//	factor_test(sizes,5,0);
//	int sizes[10] = {1,4,2,5,4,9,2,7,9,12};
//	factor_test(sizes,10,0);
//	int sizes[8] = {1,2,3,4,5,6,7,8};
//	factor_test(sizes,8,0);

//	Recombine_test(); // compare with NTL
//	recombine_test(); // pure BPAS

//	ddf_test();
//	edf_test();
	
	return 0;
}

void ddf_test() {

	DUZP_t* f = buildRandomPoly_DUZP(12,coefBoundBits,sparsity,includeNeg);
	mpz_set_ui(f->coefs[f->lt],1);
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(3);
	duspoly_t* ff = convertToDUSP_DUZP(f,Pptr);
	factsll_t* factors;
	distinctDegFactorizationInFormll1_spX(ff,&factors,Pptr);

	freePolynomial_DUZP(f);
	freePolynomial_spX(&ff);
	freeFactsll_spX(&factors);
	free(Pptr);

}

void edf_test() {

	DUZP_t* f = buildRandomPoly_DUZP(12,coefBoundBits,sparsity,includeNeg);
	mpz_set_ui(f->coefs[f->lt],1);
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(3);
	duspoly_t* ff = convertToDUSP_DUZP(f,Pptr);
	factsll_t* ddfactors;
	factsll_t* edfactors;
	
	distinctDegFactorizationInFormll1_spX(ff,&ddfactors,Pptr);
	
	equalDegFactorsInFormll1_spX(ddfactors,&edfactors,Pptr);
	
	freePolynomial_DUZP(f);
	freePolynomial_spX(&ff);
	freeFactsll_spX(&ddfactors);
	freeFactsll_spX(&edfactors);
	free(Pptr);

}

DUZP_t** makePolynomials(int deg, long r, int monic) {
	
	DUZP_t** F = (DUZP_t**) malloc(sizeof(DUZP_t*)*r);
	mpz_t one;
	mpz_init(one);
	mpz_set_ui(one,1);
	
	for (long k=0; k<r; ++k) {
		F[k] = buildRandomPoly_DUZP(deg,coefBoundBits,sparsity,includeNeg);
		if (monic) mpz_set(F[k]->coefs[F[k]->lt],one);
	}
	
	mpz_clear(one);
	return F;
}

DUZP_t** makePolynomials(int* sizes, long r, int monic) {
	
	DUZP_t** F = (DUZP_t**) malloc(sizeof(DUZP_t*)*r);
	mpz_t one;
	mpz_init(one);
	mpz_set_ui(one,1);
	
	for (long k=0; k<r; ++k) {
		F[k] = buildRandomPoly_DUZP(sizes[k],coefBoundBits,sparsity,includeNeg);
		if (monic) mpz_set(F[k]->coefs[F[k]->lt],one);
	}
	
	mpz_clear(one);
	return F;
}

DUZP_t* makeExamplePolynomial(int* sizes, long r) {
	
	DUZP_t* f = makeConstPolynomial_DUZP(1,1);
	DUZP_t* tmp;
	mpz_t one;
	mpz_init(one);
	mpz_set_ui(one,1);
	
	for (long k=0; k<r; ++k) {
		tmp = buildRandomPoly_DUZP(sizes[k],coefBoundBits,sparsity,includeNeg);
		mpz_set(tmp->coefs[tmp->lt],one);
		multiplyPolynomials_DUZP_inp(tmp,&f);
		freePolynomial_DUZP(tmp);
	}
	
	mpz_clear(one);
	return f;
}

void factor_test(int deg, int n, int monic) {

	DUZP_t* f;
	if (n > 1) {
		DUZP_t** F = makePolynomials(deg,n,monic);
		
		f = makeConstPolynomial_DUZP(1,1);
		
		for (int i=0; i<n; ++i) {
			multiplyPolynomials_DUZP_inp(F[i],&f);
			freePolynomial_DUZP(F[i]);
		}
		free(F);
			
	}
	else if (n == 1) {
		f = buildRandomPoly_DUZP(deg,coefBoundBits,sparsity,includeNeg);
		if (monic) {
			mpz_t one;
			mpz_init(one);
			mpz_set_ui(one,1);
			mpz_set(f->coefs[f->lt],one);
			mpz_clear(one);
		}
	}
	else {
		fprintf(stderr,"Invalid number of factors!");
		exit(1);
	}
	
	printPoly_DUZP(f,xvar);
	
	NTL::ZZX ff;
	NTLZZX_set_DUZP_t(ff,f);
	NTL::ZZ c;
	NTL::vec_pair_ZZX_long factors;
	
//	factor_data_t fd;
	factor(c, factors, ff, 1, 0);
//	factor(c, factors, ff, 1, 0, fd);
	
	DUZP_t* g = makeConstPolynomial_DUZP(1,0);
	for (int i=1; i<=factors.length(); ++i) {
		DUZP_t_set_NTLZZX(&g,factors(i).a);
		printPoly_DUZP(g,xvar);
	}
	
	freePolynomial_DUZP(g);
	
	mpz_t cc;
	mpz_init(cc);
	vec_DUZP_long_t factors2;
	vec_DUZP_long_t_init2(&factors2,0);
	
	long long unsigned int start(0);
	float elapsed;
	
	fprintf(stderr,"calling BPAS factor:\n");
	
	_startTimer(&start);
	
	factor(&cc,&factors2,f);
	
	_stopTimer(&start,&elapsed);
	fprintf(stderr,"time: %f\n",elapsed);
	
	g = makeConstPolynomial_DUZP(1,1);
	DUZP_t* h;
	for (int i=0; i<factors2.size; ++i) {
		h = makeConstPolynomial_DUZP(1,1);
		fprintf(stderr,"exp = %ld\n",factors2.pairs[i].b);
		for (int j=0; j<factors2.pairs[i].b; ++j)
			multiplyPolynomials_DUZP_inp(factors2.pairs[i].a,&h);
		multiplyPolynomials_DUZP_inp(h,&g);
		printPoly_DUZP(h,xvar);
		freePolynomial_DUZP(h);
	}
	multiplyByInteger_DUZP_inp(g,cc);
	
	fprintf(stderr,"g = ");
	printPoly_DUZP(g,xvar);
	fprintf(stderr,"f = ");
	printPoly_DUZP(f,xvar);
	g = subtractPolynomials_DUZP(f,g);
	fprintf(stderr,"f - g = ");
	printPoly_DUZP(g,xvar);
	
	
//	print_factor_data(fd);
	
	mpz_clear(cc);
	freePolynomial_DUZP(f);
	freePolynomial_DUZP(g);
	vec_DUZP_long_t_free(&factors2);

}

void factor_test(int* sizes, int n, int monic) {

	DUZP_t** F = makePolynomials(sizes,n,monic);
	
	DUZP_t* f = makeConstPolynomial_DUZP(1,1);
	
	for (int i=0; i<n; ++i)
		multiplyPolynomials_DUZP_inp(F[i],&f);

//	DUZP_t* f = buildRandomPoly_DUZP(16,coefBoundBits,sparsity,includeNeg);
	
	printPoly_DUZP(f,xvar);
	
	NTL::ZZX ff;
	NTLZZX_set_DUZP_t(ff,f);
	NTL::ZZ c;
	NTL::vec_pair_ZZX_long factors;
	
//	factor_data_t fd;
	factor(c, factors, ff, 1, 0);
	
	DUZP_t* g = makeConstPolynomial_DUZP(1,0);
	for (int i=1; i<=factors.length(); ++i) {
		DUZP_t_set_NTLZZX(&g,factors(i).a);
		printPoly_DUZP(g,xvar);
	}
	
	freePolynomial_DUZP(g);
	
	mpz_t cc;
	mpz_init(cc);
	vec_DUZP_long_t factors2;
	vec_DUZP_long_t_init2(&factors2,0);
	
	factor(&cc,&factors2,f);
	
	g = makeConstPolynomial_DUZP(1,1);
	DUZP_t* h;
	for (int i=0; i<factors2.size; ++i) {
		h = makeConstPolynomial_DUZP(1,1);
		fprintf(stderr,"exp = %ld\n",factors2.pairs[i].b);
		for (int j=0; j<factors2.pairs[i].b; ++j)
			multiplyPolynomials_DUZP_inp(factors2.pairs[i].a,&h);
		multiplyPolynomials_DUZP_inp(h,&g);
		printPoly_DUZP(h,xvar);
		freePolynomial_DUZP(h);
	}
	multiplyByInteger_DUZP_inp(g,cc);
	
	fprintf(stderr,"g = ");
	printPoly_DUZP(g,xvar);
	fprintf(stderr,"f = ");
	printPoly_DUZP(f,xvar);
	g = subtractPolynomials_DUZP(f,g);
	fprintf(stderr,"f - g = ");
	printPoly_DUZP(g,xvar);
	
//	print_factor_data(fd);
	
	freePolynomial_DUZP(f);
	freePolynomial_DUZP(g);

}



void factor_timing(int deg, int n, int monic) {

	int showOutput = 0;
	std::streambuf *oldcout,*oldcerr;
	std::fstream fs("/dev/null");
	int bakfs1, bakfs2, newfs1, newfs2, bakcout, newcout;
	
	long long unsigned int start(0);
	float elapsed;

	if (!showOutput) {
//		fprintf(stderr,"hiding output...\n");
		oldcout = std::cout.rdbuf(); // <-- save  
		oldcerr = std::cerr.rdbuf(); // <-- save  

		std::cout.rdbuf(fs.rdbuf());
		std::cerr.rdbuf(fs.rdbuf());
		
		fflush(stdout);
		bakfs1 = dup(1);
		newfs1 = open("/dev/null", O_WRONLY);
		dup2(newfs1, 1);
		fflush(stderr);
		bakfs2 = dup(2);
		newfs2 = open("/dev/null", O_WRONLY);
		dup2(newfs2, 2);
		
	}

	DUZP_t* f;
	if (n > 1) {
		DUZP_t** F = makePolynomials(deg,n,monic);
		
		f = makeConstPolynomial_DUZP(1,1);
		
		for (int i=0; i<n; ++i) {
			multiplyPolynomials_DUZP_inp(F[i],&f);
			freePolynomial_DUZP(F[i]);
		}
		free(F);
			
	}
	else if (n == 1) {
		f = buildRandomPoly_DUZP(deg,coefBoundBits,sparsity,includeNeg);
		if (monic) {
			mpz_t one;
			mpz_init(one);
			mpz_set_ui(one,1);
			mpz_set(f->coefs[f->lt],one);
			mpz_clear(one);
		}
	}
	else {
		fprintf(stderr,"Invalid number of factors!");
		exit(1);
	}
	
	printPoly_DUZP(f,xvar);
	
	/// NTL Test ///
	
	NTL::ZZX ff;
	NTLZZX_set_DUZP_t(ff,f);
	NTL::ZZ c;
	NTL::vec_pair_ZZX_long factors;
	
	_startTimer(&start);
	
//	factor_data_t fd;

	factor(c, factors, ff, 1, 0);
	
//	factor(c, factors, ff, 1, 0, fd);
	
	_stopTimer(&start,&elapsed);
	
	if (!showOutput) {
		std::cout.rdbuf(oldcout);   // <-- restore
		std::cerr.rdbuf(oldcerr);   // <-- restore
		fflush(stdout);
		dup2(bakfs1, 1);
		fflush(stderr);
		dup2(bakfs2, 2);
	}
	
	fprintf(stderr,"NTL time: %f\n",elapsed);
	
	if (!showOutput) {
		std::cout.rdbuf(fs.rdbuf());
		std::cerr.rdbuf(fs.rdbuf());
		fflush(stdout);
		dup2(newfs1, 1);
		fflush(stderr);
		dup2(newfs2, 2);
	}
	
	DUZP_t* g = makeConstPolynomial_DUZP(1,0);
	for (int i=1; i<=factors.length(); ++i) {
		DUZP_t_set_NTLZZX(&g,factors(i).a);
		printPoly_DUZP(g,xvar);
	}
	
	freePolynomial_DUZP(g);
	
	/// BPAS Test ///
	
	mpz_t cc;
	mpz_init(cc);
	vec_DUZP_long_t factors2;
	vec_DUZP_long_t_init2(&factors2,0);
	
	fprintf(stderr,"calling BPAS factor:\n");
	
	_startTimer(&start);
	
	factor(&cc,&factors2,f);
	
	_stopTimer(&start,&elapsed);
	
	if (!showOutput) {
		std::cout.rdbuf(oldcout);   // <-- restore
		std::cerr.rdbuf(oldcerr);   // <-- restore
		fflush(stdout);
		dup2(bakfs1, 1);
		fflush(stderr);
		dup2(bakfs2, 2);
//		close(bakfs);
	}
	fprintf(stderr,"BPAS time: %f\n",elapsed);
	
	if (!showOutput) {
		std::cout.rdbuf(fs.rdbuf());
		std::cerr.rdbuf(fs.rdbuf());
		fflush(stdout);
		dup2(newfs1, 1);
		fflush(stderr);
		dup2(newfs2, 2);
		close(newfs1);
		close(newfs2);
	}
	
	g = makeConstPolynomial_DUZP(1,1);
	DUZP_t* h;
	for (int i=0; i<factors2.size; ++i) {
		h = makeConstPolynomial_DUZP(1,1);
		fprintf(stderr,"exp = %ld\n",factors2.pairs[i].b);
		for (int j=0; j<factors2.pairs[i].b; ++j)
			multiplyPolynomials_DUZP_inp(factors2.pairs[i].a,&h);
		multiplyPolynomials_DUZP_inp(h,&g);
		printPoly_DUZP(h,xvar);
		freePolynomial_DUZP(h);
	}
	multiplyByInteger_DUZP_inp(g,cc);
	
	fprintf(stderr,"g = ");
	printPoly_DUZP(g,xvar);
	fprintf(stderr,"f = ");
	printPoly_DUZP(f,xvar);
	g = subtractPolynomials_DUZP(f,g);
	fprintf(stderr,"f - g = ");
	printPoly_DUZP(g,xvar);
	
	/// MAPLE Test ///
	
	
    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    MKernelVector kv = mapleTest->getMKernelVector();
    
    FILE *fp;
	fprintf(stderr,"This will display on the screen.\n");
	
	if (!showOutput) {
		std::cout.rdbuf(oldcout);   // <-- restore
		std::cerr.rdbuf(oldcerr);   // <-- restore
		fflush(stdout);
		dup2(bakfs1, 1);
		fflush(stderr);
		dup2(bakfs2, 2);
//		close(bakfs);
	}

	fflush(stderr);
	bakfs2 = dup(2);
	dup2(newfs2, 2);
	int fp2 = open("./duz_tst_tmp_OUT", O_RDWR | O_CREAT | O_TRUNC, 0600);
//	if((fp=freopen("/home/rmoir2/BPAS/Main.new/tests/IntegerPolynomial/OUT", "w" ,stderr))==NULL) {
//		printf("Cannot open file.\n");
//		exit(1);
//	}
	
	fflush(stderr);
	bakfs2 = dup(2);
	dup2(fp2, 2);

//	fprintf(stderr,"This will be written to the file OUT.");
	printPoly_DUZP(f,xvar);

	close(fp2);
	
	fflush(stderr);
	dup2(bakfs2, 2);
	
	if (!showOutput) {
		std::cout.rdbuf(fs.rdbuf());
		std::cerr.rdbuf(fs.rdbuf());
		fflush(stdout);
		dup2(newfs1, 1);
		fflush(stderr);
		dup2(newfs2, 2);
	}
	
	std::ifstream file("./duz_tst_tmp_OUT");
	if (!file) {
		std::cout << "unable to open file";
		exit(1);
	}
	
	std::string s;
	std::getline(file, s);
	file.close();
    	
	std::cout << "this is a test" << std::endl;
	
	s += ":";
	
	std::cout << s << std::endl;
	
	fprintf(stderr,"STDERR is working again...\n");
    
    char* cstr;
    std::string evalStr;
    evalStr = "kernelopts(numcpus=1);";
    cstr = new char[evalStr.length()+1];
    std::strcpy (cstr, evalStr.c_str());
    EvalMapleStatement(kv, cstr);
    delete[] cstr;
    
    cstr = new char[s.length()+1];
    std::strcpy(cstr, s.c_str());
    ALGEB res = EvalMapleStatement(kv, cstr);
    delete[] cstr;

	std::string procStr = "factorTest := proc (f::polynom) local a; a := factors(f); end proc:";
	
    cstr = new char[procStr.length()+1];
    std::strcpy (cstr, procStr.c_str());
    ALGEB testProc = EvalMapleStatement(kv, cstr);
    delete[] cstr;
	
	
	if (!showOutput) {
		std::cout.rdbuf(oldcout);   // <-- restore
		std::cerr.rdbuf(oldcerr);   // <-- restore
		fflush(stdout);
		dup2(bakfs1, 1);
		fflush(stderr);
		dup2(bakfs2, 2);
	}
	
	_startTimer(&start);
	
//	fprintf(stderr,"computing factorization in Maple...\n");
	
    ALGEB result = EvalMapleProc(kv, testProc, 1, res);
	
	_stopTimer(&start,&elapsed);
	
	fprintf(stderr,"Maple time: %f\n",elapsed);
	
	if (!showOutput) {
		std::cout.rdbuf(fs.rdbuf());
		std::cerr.rdbuf(fs.rdbuf());
		fflush(stdout);
		dup2(newfs1, 1);
		fflush(stderr);
		dup2(newfs2, 2);
		close(newfs1);
		close(newfs2);
	}
	
	
//	print_factor_data(fd);
	
	mpz_clear(cc);
	freePolynomial_DUZP(f);
	freePolynomial_DUZP(g);
	vec_DUZP_long_t_free(&factors2);
	
	if (!showOutput) {
		std::cout.rdbuf(oldcout);   // <-- restore
		std::cerr.rdbuf(oldcerr);   // <-- restore
		
		fflush(stdout);
		dup2(bakfs1, 1);
		fflush(stderr);
		dup2(bakfs2, 2);
//		close(bakfs1);
//		close(bakfs2);
	}

}

void make_Recombine_input(int* sizes, int r, long_t p, DUZP_t** f, vec_DUZP_t* padicFactors, mpz_t* P, int* e, long_t* bound) {


	mpz_t one;
	mpz_init(one);
	mpz_set_ui(one,1);
	DUZP_t* tmp;
	
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(pr);
	duspoly_t* modf;
	factsll_t* modularFactorsll = initFactsll_spX(1);
	factors_t* modularFactors;
	
	int allOne = 0;
	long n;
	
	while (!allOne) {
		*f = makeConstPolynomial_DUZP(1,1);
		fprintf(stderr,"start build f for loop\n");
		for (int i=0; i<r; ++i) {
			tmp = buildRandomPoly_DUZP(sizes[i],coefBoundBits,sparsity,includeNeg);
			mpz_set(tmp->coefs[tmp->lt],one);
			printPoly_DUZP(tmp,xvar);
			*f = multiplyPolynomials_DUZP(*f,tmp);
		}
		printPoly_DUZP(*f,xvar);
		
		modf = convertToDUSP_DUZP(*f,Pptr);
		modularFactorsll = initFactsll_spX(1);
		modFactorizationInFormll1_spX(modf, &modularFactorsll, Pptr);

		modularFactors = convertToFactorsList_spX(modularFactorsll,Pptr);
		n = modularFactors->alloc;
		allOne = 1;
		fprintf(stderr,"start print modular factors loop\n");
		for (int i=0; i<n; ++i) {
//			printPolynomialOutForm_spX(modularFactors->polys[i],Pptr);
			fprintf(stderr,"exp = %ld\n",modularFactors->exps[i]);
			if (modularFactors->exps[i] > 1) {
				allOne = 0;
				freePolynomial_DUZP(*f);
				break;
			}
		}
	}
	DUZP_t** liftedModularFactors = (DUZP_t**) malloc(sizeof(DUZP_t*)*n);
	
	fprintf(stderr,"check1\n");
//	*e = multiTermPadicLiftOpt_Iterative(*f, (const duspoly_t**) modularFactors->polys, liftedModularFactors, n, Pptr);
	*e = multiTermPadicLiftOpt_Iterative(*f, CONSTCONSTCAST(duspoly_t,modularFactors->polys), liftedModularFactors, n, Pptr);
	
	fprintf(stderr,"check2\n");
	for (int i=0; i<n; ++i) {
		printPoly_DUZP(liftedModularFactors[i],xvar);
	}
	
	long* idxs = (long*) malloc(sizeof(long)*n);
	for (int i=0; i<n; ++i) {
		idxs[i] = i;
	}
	
	mpz_t pp;
	mpz_init_set_ui(pp,p);
	mpz_pow_ui(*P,pp,*e);
	DUZP_t* g;
	ComputePotentialFactorRational(&g, liftedModularFactors, idxs, n, *P);
	
	printPoly_DUZP(g,xvar);
	
	*bound = LandauMignotteBoundInBits(*f,n-1);
	fprintf(stderr,"log2(LM_Bound) = %Ld\n",*bound);
	
	vec_DUZP_t_free(padicFactors);
	vec_DUZP_t_init_set(padicFactors,n);
	padicFactors->polys = liftedModularFactors;
	padicFactors->alloc = n;
	padicFactors->size = n;

	mpz_clear(one);
	mpz_clear(pp);	
	free(idxs);
//	free(liftedModularFactors);

}

//void recombine_test() {

//	mpz_t one;
//	mpz_init(one);
//	mpz_set_ui(one,1);
//	DUZP_t* f; // random polynomial to factor
////	DUZP_t* tmp;
//	
//	Prime_ptr* Pptr = smallprimefield_get_prime_constants(pr);
//	duspoly_t* modf;
////	factsll_t* modularFactorsll; // = initFactsll_spX(1);
//	factors_t* modularFactors;
//	
//	int r = 5; // minimum number of rational factors
//	int sizes[5] = {2,3,4,5,6}; // sizes of factors
//	
//	int allOne = 0;
//	long n;
//	
//	f = makeExamplePolynomial(sizes,r);
//	
////	while (!allOne) {
////		f = makeConstPolynomial_DUZP(1,1);
////		fprintf(stderr,"start build f for loop\n");
////		for (int i=0; i<r; ++i) {
////			tmp = buildRandomPoly_DUZP(sizes[i],coefBoundBits,sparsity,includeNeg);
////			mpz_set(tmp->coefs[tmp->lt],one);
////			printPoly_DUZP(tmp,xvar);
////			multiplyPolynomials_DUZP_inp(tmp,&f);
////			freePolynomial_DUZP(tmp);
////		}
//////		printPoly_DUZP(f,xvar);
//////		return;
////		
////		modf = convertToDUSP_DUZP(f,Pptr);
//////		modularFactorsll = initFactsll_spX(1);
////		modFactorizationInFormll1_spX(modf, &modularFactorsll, Pptr);
////		
////		modularFactors = convertToFactorsList_spX(modularFactorsll,Pptr);
////		freeFactsll_spX(&modularFactorsll);
////		n = modularFactors->alloc;
////		allOne = 1;
////		fprintf(stderr,"start print modular factors loop\n");
////		for (int i=0; i<n; ++i) {
//////			printPolynomialOutForm_spX(modularFactors->polys[i],Pptr);
////			fprintf(stderr,"exp = %ld\n",modularFactors->exps[i]);
////			if (modularFactors->exps[i] > 1) {
////				allOne = 0;
////				freePolynomial_DUZP(f);
////				freePolynomial_spX(&modf);
////				freeFactors_spX(&modularFactors);
////				break;
////			}
////		}
////	}
////	printPoly_DUZP(f,xvar);
////	return;

//	factoring_info_t fi;
//	factoring_info_t_init(&fi);
//	SmallPrimeFactorize(&modularFactors,&fi,f);
//	n = modularFactors->alloc;

//	DUZP_t** liftedModularFactors = (DUZP_t**) malloc(sizeof(DUZP_t*)*n);
//	duspoly_t ** sigmas;
//	long bnd = LandauMignotteBoundInBits(f,f->lt);
////	long bnd = 6;
//	mpz_t bound;
//	mpz_init_set_ui(bound,1);
//	mpz_mul_2exp(bound,bound,bnd);
//	gmp_fprintf(stderr,"bound = %Zd\n",bound);
//	
//	fprintf(stderr,"check1\n");
////	int e = multiTermPadicLiftOpt_Iterative(f, (const duspoly_t**) modularFactors->polys, liftedModularFactors, n, Pptr);
////	int e = multiTermPadicLiftOpt_Iterative(f, CONSTCONSTCAST(duspoly_t,modularFactors->polys), &liftedModularFactors, n, Pptr);

//	int e = multiTermPadicLiftStart(f, CONSTCONSTCAST(duspoly_t,modularFactors->polys), liftedModularFactors, &sigmas, n, bound, fi.Pptr);

//	fprintf(stderr,"check2\n");
//	for (int i=0; i<n; ++i) {
//		printPoly_DUZP(liftedModularFactors[i],xvar);
//	}
//	
//	long* idxs = (long*) malloc(sizeof(long)*n);
//	for (int i=0; i<n; ++i) {
//		idxs[i] = i;
//	}
////	
//	mpz_t P,p;
//	mpz_init(P);
//	mpz_init(p);
//	mpz_set_si(p,11);
//	mpz_pow_ui(P,p,e);
////	DUZP_t* g;
////	ComputePotentialFactorRational(&g, liftedModularFactors, idxs, n, P);
////	
////	printPoly_DUZP(g,xvar);
////	
//	vec_DUZP_t padicFactors;
//	vec_DUZP_t_init2(&padicFactors,n);
//	for (int i=0; i<n; ++i)
//		vec_DUZP_t_push(&padicFactors,liftedModularFactors[i]);
////	padicFactors.polys = liftedModularFactors;
////	padicFactors.alloc = n;
////	padicFactors.size = n;
////	
//	vec_DUZP_t rationalFactors;
//	vec_DUZP_t_init(&rationalFactors);
////	rationalFactors.polys = (DUZP_t**) malloc(sizeof(DUZP_t*));
////	rationalFactors.alloc = 1;
////	rationalFactors.size = 0;
////	
//	DUZP_t* sigmas2;
////	
//	Recombine(&rationalFactors,f,&padicFactors,&sigmas2,P,11,e,mpz_sizeinbase(P,2),modularFactors,sigmas);
////	
////	fprintf(stderr,"\n%ld rationalFactors:\n",rationalFactors.size);
////	for (int i=0; i<rationalFactors.size; ++i) {
////		printPoly_DUZP(rationalFactors.polys[i],xvar);
////	}
//////	fprintf(stderr,"alloc = %ld:\n",rationalFactors.alloc);
//////	fprintf(stderr,"size = %ld:\n",rationalFactors.size);
////	fprintf(stderr,"\n%ld padicFactors:\n",padicFactors.size);
////	for (int i=0; i<padicFactors.size; ++i) {
////		printPoly_DUZP(padicFactors.polys[i],xvar);
////	}

//	factoring_info_t_free(&fi);
//	vec_DUZP_t_free(&rationalFactors);
//	vec_DUZP_t_free(&padicFactors);
//	for (int i=0; i<n; ++i)
//		freePolynomial_DUZP(liftedModularFactors[i]);
//	free(liftedModularFactors);
//	for (int i=0; i<n; ++i) {
//		freePolynomial_spX(&sigmas[i]);
//	}
//	free(sigmas);
//	free(sigmas2);
//	freeFactors_spX(&modularFactors);
//	freePolynomial_DUZP(f);
//	freePolynomial_spX(&modf);
//	free(idxs);
//	mpz_clear(one);
//	mpz_clear(bound);
//	mpz_clear(P);
//	mpz_clear(p);
//	free(Pptr);

//}

void vec_mpz_test() {

	int size = 10;
	vec_mpz_t a;
	vec_mpz_t_init2(&a,size);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	mpz_t n;
	mpz_init(n);
	
	for (int i=0; i<size+1; ++i) {
		mpz_set_ui(n,i+1);
		vec_mpz_t_push(&a,n);
		fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	}
	
	print_vec_mpz(a);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	for (int i=0; i<size+1; ++i) {
		vec_mpz_t_pop(&n,&a);
		gmp_fprintf(stderr,"n = %Zd, alloc=%ld, size=%ld\n",n,a.alloc,a.size);
	}
	
	print_vec_mpz(a);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	mpz_clear(n);
	vec_mpz_t_free(&a);

}

void vec_DUZP_test() {

	int alloc = 10;
	vec_DUZP_t a;
	vec_DUZP_t_init2(&a,alloc);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	DUZP_t* f;
	
	for (int i=0; i<alloc+1; ++i) {
		f = buildRandomPoly_DUZP(5,coefBoundBits,sparsity,includeNeg);
		vec_DUZP_t_push(&a,f);
		freePolynomial_DUZP(f);
		fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	}
	
	print_vec_DUZP(a);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	for (int i=0; i<alloc+1; ++i) {
		f = vec_DUZP_t_pop(&a);
		printPoly_DUZP(f,xvar);
		freePolynomial_DUZP(f);
		fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	}
	
	print_vec_DUZP(a);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	vec_DUZP_t_free(&a);

}

void vec_duspoly_test() {

	int alloc = 11;
	vec_duspoly_t a;
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(11);
	vec_duspoly_t_init2(&a,Pptr,alloc);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
//	duspoly_t* f;
//	
//	for (int i=0; i<alloc+1; ++i) {
//		f = binPolynomialInForm_spX(i+1,a.ptr);
//		vec_duspoly_t_push(&a,f);
//		freePolynomial_spX(f);
//		fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
//	}
//	
//	print_vec_duspoly(a);
//	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
//	
//	for (int i=0; i<alloc+1; ++i) {
//		f = vec_duspoly_t_pop(&a);
//		printPolynomialOutForm_spX(f,a.ptr);
//		freePolynomial_spX(f);
//		fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
//	}
//	
//	print_vec_duspoly(a);
//	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	vec_duspoly_t_free(&a);

}

void vec_DUZP_long_test() {

//	DUZP_t* f = buildRandomPoly_DUZP(5,coefBoundBits,sparsity,includeNeg);
//	printPoly_DUZP(f,xvar);
//	differentiate_DUZP_inp(&f);
//	printPoly_DUZP(f,xvar);
//	DUZP_t* df = differentiate_DUZP(f);
//	printPoly_DUZP(df,xvar);
//	freePolynomial_DUZP(df);
//	freePolynomial_DUZP(f);

	int alloc = 10;
	vec_DUZP_long_t a;
	vec_DUZP_long_t_init2(&a,alloc);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	DUZP_t* f;
	DUZP_long_t tmp;
	
	for (int i=0; i<alloc+1; ++i) {
		f = buildRandomPoly_DUZP(5,coefBoundBits,sparsity,includeNeg);
		tmp.a = f;
		tmp.b = i+1;
		vec_DUZP_long_t_push(&a,tmp);
		freePolynomial_DUZP(f);
		fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	}
	
	print_vec_DUZP_long(a);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	DUZP_t* tmp2;
	
	for (int i=0; i<alloc-1; ++i) {
		tmp = vec_DUZP_long_t_pop(&a);
		printPoly_DUZP(tmp.a,xvar);
		differentiate_DUZP_inp(&tmp.a);
		printPoly_DUZP(tmp.a,xvar);
		tmp2 = differentiate_DUZP(tmp.a);
		printPoly_DUZP(tmp2,xvar);
		freePolynomial_DUZP(tmp.a);
		freePolynomial_DUZP(tmp2);
		fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	}
	
	print_vec_DUZP_long(a);
	fprintf(stderr,"alloc=%ld, size=%ld\n",a.alloc,a.size);
	
	vec_DUZP_long_t_free(&a);

}

void vec_vec_test() {
	
	int length = 10;
	int depth = 2;
	vec_vec_mpz_t vec;
	vec_vec_mpz_t_init(&vec,length);
	fprintf(stderr,"one\n");
	vec_vec_mpz_t_set(&vec,length,depth);
	fprintf(stderr,"two\n");
	
	mpz_t t;
	mpz_init(t);
	
	// filling the array
	for (int i=0; i<vec.length; ++i)
		for (int j=0; j<vec.depth; ++j) {
			mpz_add_ui(t,t,1);
//			fprintf(stderr,"i=%d, j=%d\n",i,j);
			mpz_set(vec.array[i][j],t);
		}
	print_vec_vec(&vec);
	
	// decrease the depth
	vec_vec_mpz_t_set(&vec,length,depth-1);
	print_vec_vec(&vec);
	
	// increase the depth
	vec_vec_mpz_t_set(&vec,length,depth+1);
	
	// filling the array
	for (int i=0; i<vec.length; ++i)
		for (int j=depth-1; j<vec.depth; ++j) {
			mpz_add_ui(t,t,1);
//			fprintf(stderr,"i=%d, j=%d\n",i,j);
			mpz_set(vec.array[i][j],t);
		}
	print_vec_vec(&vec);
	
	// decrease depth and length
	vec_vec_mpz_t_set(&vec,length-5,depth);
	print_vec_vec(&vec);
	
	// increase depth and length
	vec_vec_mpz_t_set(&vec,length,depth+2);
	
	// filling the array
	for (int i=0; i<length-5; ++i)
		for (int j=depth; j<vec.depth; ++j) {
			mpz_add_ui(t,t,1);
//			fprintf(stderr,"i=%d, j=%d\n",i,j);
			mpz_set(vec.array[i][j],t);
		}
	for (int i=length-5; i<vec.length; ++i)
		for (int j=0; j<vec.depth; ++j) {
			mpz_add_ui(t,t,1);
//			fprintf(stderr,"i=%d, j=%d\n",i,j);
			mpz_set(vec.array[i][j],t);
		}
	print_vec_vec(&vec);
	
	vec_vec_mpz_t_free(&vec);
	mpz_clear(t);
}

void naive_search_test() {

	mpz_t one;
	mpz_init(one);
	mpz_set_ui(one,1);
	DUZP_t* f; // random polynomial to factor
	DUZP_t* tmp;
	
	Prime_ptr* Pptr = smallprimefield_get_prime_constants(pr);
	duspoly_t* modf;
	factsll_t* modularFactorsll = initFactsll_spX(1);
	factors_t* modularFactors;
	
	int r = 3; // minimum number of rational factors
	int sizes[3] = {2,3,4}; // sizes of factors
	
	int allOne = 0;
	long n;
	
	while (!allOne) {
		f = makeConstPolynomial_DUZP(1,1);
		fprintf(stderr,"start build f for loop\n");
		for (int i=0; i<r; ++i) {
			tmp = buildRandomPoly_DUZP(sizes[i],coefBoundBits,sparsity,includeNeg);
			mpz_set(tmp->coefs[tmp->lt],one);
			printPoly_DUZP(tmp,xvar);
			f = multiplyPolynomials_DUZP(f,tmp);
		}
		printPoly_DUZP(f,xvar);
		
		modf = convertToDUSP_DUZP(f,Pptr);
		modularFactorsll = initFactsll_spX(1);
		modFactorizationInFormll1_spX(modf, &modularFactorsll, Pptr);
		

		modularFactors = convertToFactorsList_spX(modularFactorsll,Pptr);
		n = modularFactors->alloc;
		allOne = 1;
		fprintf(stderr,"start print modular factors loop\n");
		for (int i=0; i<n; ++i) {
//			printPolynomialOutForm_spX(modularFactors->polys[i],Pptr);
			fprintf(stderr,"exp = %ld\n",modularFactors->exps[i]);
			if (modularFactors->exps[i] > 1) {
				allOne = 0;
				freePolynomial_DUZP(f);
				break;
			}
		}
	}
	DUZP_t** liftedModularFactors = (DUZP_t**) malloc(sizeof(DUZP_t*)*n);
	
	fprintf(stderr,"check1\n");
//	int e = multiTermPadicLiftOpt_Iterative(f, (const duspoly_t**) modularFactors->polys, liftedModularFactors, n, Pptr);
	int e = multiTermPadicLiftOpt_Iterative(f, CONSTCONSTCAST(duspoly_t,modularFactors->polys), liftedModularFactors, n, Pptr);
	
	fprintf(stderr,"check2\n");
	for (int i=0; i<n; ++i) {
		printPoly_DUZP(liftedModularFactors[i],xvar);
	}
	
	long* idxs = (long*) malloc(sizeof(long)*n);
	for (int i=0; i<n; ++i) {
		idxs[i] = i;
	}
	
	mpz_t P,p;
	mpz_init(P);
	mpz_init(p);
	mpz_set_si(p,11);
	mpz_pow_ui(P,p,e);
	DUZP_t* g;
	ComputePotentialFactorRational(&g, liftedModularFactors, idxs, n, P);
	
	printPoly_DUZP(g,xvar);
	
	long_t bound = LandauMignotteBoundInBits(f,n-1);
	fprintf(stderr,"log2(LM_Bound) = %Ld\n",bound);
	
	vec_DUZP_t padicFactors;
	padicFactors.polys = liftedModularFactors;
	padicFactors.alloc = n;
	padicFactors.size = n;
	
	vec_DUZP_t rationalFactors;
	rationalFactors.polys = (DUZP_t**) malloc(sizeof(DUZP_t*));
	rationalFactors.alloc = 1;
	rationalFactors.size = 0;
	
	NaiveFactorSearch(&f,&rationalFactors,&padicFactors,P,11,e,1,mpz_sizeinbase(P,2));
	fprintf(stderr,"f:\n");
	printPoly_DUZP(f,xvar);
	if (padicFactors.size >= 2)
		NaiveFactorSearch(&f,&rationalFactors,&padicFactors,P,11,e,2,mpz_sizeinbase(P,2));
	if (padicFactors.size >= 3)
		NaiveFactorSearch(&f,&rationalFactors,&padicFactors,P,11,e,3,mpz_sizeinbase(P,2));
		
	// Current issues to fix:
	//   - degree six factor found
	
	fprintf(stderr,"\n%ld rationalFactors:\n",rationalFactors.size);
	for (int i=0; i<rationalFactors.size; ++i) {
		printPoly_DUZP(rationalFactors.polys[i],xvar);
	}
//	fprintf(stderr,"alloc = %ld:\n",rationalFactors.alloc);
//	fprintf(stderr,"size = %ld:\n",rationalFactors.size);
	fprintf(stderr,"\n%ld padicFactors:\n",padicFactors.size);
	for (int i=0; i<padicFactors.size; ++i) {
		printPoly_DUZP(padicFactors.polys[i],xvar);
	}
//	free(f);
//	free(liftedModularFactors);
	
}
