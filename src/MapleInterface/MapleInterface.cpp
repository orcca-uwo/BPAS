


#include "MapleInterface/MapleInterface.hpp"
#include "IntegerPolynomial/mzpolynomial.hpp"

// #include "Utils/Unix_Timer.h"
// #include <atomic>

// std::atomic<long> mapleGCDTimeSec;
// std::atomic<long> mapleGCDTimeUSec;

// std::atomic<long> mapleFactorTimeSec;
// std::atomic<long> mapleFactorTimeUSec;

SparseMultivariateIntegerPolynomial MapleInterface::gcd(const SparseMultivariateIntegerPolynomial& a, const SparseMultivariateIntegerPolynomial& b) {

	// timer_id id = start_timer();
	std::vector<std::string> request, response;
	request.emplace_back("gcd");
	request.emplace_back(a.toString());
	request.emplace_back(b.toString());

	// std::cerr << std::this_thread::get_id() << " is about to wait for CPP GCD" << std::endl;
	unsigned long long start;
	startTimer(&start);

	sendCommand(request, response);

	stopTimerAddElapsed(&start, &g_mapleWaitTime);
	// std::cerr << std::this_thread::get_id() << " finished waiting for CPP GCD" << std::endl;

	char* c_names[a.nvar];
    for (int i = 0; i < a.nvar; ++i) {
        std::string str = a.names[i+1].toString();
        c_names[i] = (char*) malloc(sizeof(char)*str.length()+1);
        strcpy(c_names[i], str.c_str());
    }
	AltArr_t* opQ = generate_altarr_var_defined(response[0].c_str(), (const char**)c_names, a.nvar);

    // timer_time elapsed1 = elapsed_time(&id);
    // mapleGCDTimeSec += elapsed1.tv_sec;
    // mapleGCDTimeUSec += elapsed1.tv_usec;
    // double time = (mapleGCDTimeSec + ((double)mapleGCDTimeUSec / 1000000));
    // fprintf(stderr, "\nmapleGCDTime: %f\n", time);
	
	return SparseMultivariateIntegerPolynomial (RationalNumber(1), opQ, a.nvar, a.names);

}

Factors<SparseMultivariateIntegerPolynomial> MapleInterface::factor(const SparseMultivariateIntegerPolynomial& p) {
	// timer_id id = start_timer();

	std::vector<std::string> request, response;
	request.emplace_back("factors");
	request.emplace_back(p.toString());

	// std::cerr << std::this_thread::get_id() << " is about to wait for factorization" << std::endl;
	unsigned long long start;
	startTimer(&start);
	sendCommand(request, response);
	stopTimerAddElapsed(&start, &g_mapleWaitTime);
	// std::cerr << std::this_thread::get_id() << " finished waiting for factorization" << std::endl;


	Factors<SparseMultivariateIntegerPolynomial> retFacts;


	//response comes as numericFact, fact1, exp1, ..., factn, expn,
	RationalNumber ringElem(response[0]);

	char* c_names[p.nvar];
    for (int i = 0; i < p.nvar; ++i) {
        std::string str = p.names[i+1].toString();
        c_names[i] = (char*) malloc(sizeof(char)*str.length()+1);
        strcpy(c_names[i], str.c_str());
    }

    mpq_t cont;
    mpq_init(cont);
	for (size_t i = 1; i < response.size(); i += 2) {	
		AltArr_t* opQ = generate_altarr_var_defined(response[i].c_str(), (const char**) c_names, p.nvar);
		AltArrZ_t* opZ = primitivePartAndContent_AAZFromAA(opQ, cont);

		// retFacts[i-1] = deepCopyPolynomial_AAZFromAA(opQ);
		
		int exp = atoi(response[i+1].c_str());
		for (int k = 0; k < exp; ++k) {
			mpq_mul(ringElem.get_mpq_t(), ringElem.get_mpq_t(), cont);
		}

		freePolynomial_AA(opQ);

		SparseMultivariateIntegerPolynomial zfact(opZ, p.nvar, p.names);
		retFacts.addFactor(std::move(zfact), exp);
	}
	mpq_clear(cont);

	retFacts.setRingElement(ringElem);

    // timer_time elapsed1 = elapsed_time(&id);
	// mapleFactorTimeSec += elapsed1.tv_sec;
 //    mapleFactorTimeUSec += elapsed1.tv_usec;
 //    double time = (mapleFactorTimeSec + ((double)mapleFactorTimeUSec / 1000000));
 //    fprintf(stderr, "\nmapleFactorTime: %f\n", time);
	return retFacts;
}


    
bool MapleInterface::TriangularizeValidation(bool isLazard, std::vector<std::string>& inputs) {

	std::vector<std::string> response;
	inputs.insert(inputs.begin(), std::to_string(isLazard));
	inputs.insert(inputs.begin(), "triangularizeValidate");
	sendCommand(inputs, response);

	if (response[0] == "1") {
		return true;
	} else {
		return false;
	}

}