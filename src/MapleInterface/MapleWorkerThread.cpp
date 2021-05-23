

#include "MapleInterface/MapleWorkerThread.hpp"
#include "Utils/Unix_Timer.h"


bool MapleWorkerThread::processTask(const std::vector<std::string>& request, std::vector<std::string>& response) {
	if (request.size() < 1) {
		return false;
	}

	bool ret;

	// unsigned long long start;
	// _startTimer(&start);
	if (request[0] == "factors") {
		ret = factors_mapleWorker(request, response);
	} else if (request[0] == "gcd") {
		ret = gcd_mapleWorker(request, response);
	} else if (request[0] == "triangularizeValidate") {
		ret = triangularizeValidate_mapleWorker(request, response);
	} else if (request[0] == "startKernel") {
		MKernelVector* kv_p = getMapleKVSingleton();
		response.push_back("started");
		ret = true;
	} else if (request[0] == "stopKernel") {
		restartMapleKernel(); //actually stopping the kernel makes Maple (2020) freak out
		response.push_back("stopped");
		ret = true;
	}
	// float time = 0;
	// _stopTimer(&start, &time);
	// std::cerr << "Actual kernal time: " << time << std::endl;
	return ret;
}

void MapleWorkerThread::threadCleanup() {
	std::vector<std::string> request;
	request.emplace_back("stopKernel");
	sendRequestAndWait(request);
}

bool MapleWorkerThread::factors_mapleWorker(const std::vector<std::string>& request, std::vector<std::string>& response) {

	if(request.size() < 2) {
		return false;
	}

	std::string poly = request[1] + ":";

	char** facts = NULL;
	char** exps = NULL;
	char* numericFact = NULL;
	int nfacts = 0;

	factorPolynomial_MplInt_string(poly.c_str(), &nfacts, &facts, &exps, &numericFact);

	std::string tmpStr(numericFact);
	tmpStr.erase(std::remove(tmpStr.begin(), tmpStr.end(), '\n'), tmpStr.end());
	tmpStr.erase(std::remove(tmpStr.begin(), tmpStr.end(), '\\'), tmpStr.end());

	response.emplace_back(std::move(tmpStr));
	free(numericFact);

	for(int i = 0; i < nfacts; ++i) {
		tmpStr = std::string(facts[i]);
		tmpStr.erase(std::remove(tmpStr.begin(), tmpStr.end(), '\n'), tmpStr.end());
		tmpStr.erase(std::remove(tmpStr.begin(), tmpStr.end(), '\\'), tmpStr.end());
		response.emplace_back(std::move(tmpStr));

		response.emplace_back(exps[i]);
		free(facts[i]);
		free(exps[i]);
	}

	return true;
}

bool MapleWorkerThread::gcd_mapleWorker(const std::vector<std::string>& request, std::vector<std::string>& response) {

	if(request.size() < 3) {
		return false;
	}

	std::string f = request[1] + ":";
	std::string g = request[2] + ":";

	char* gcd = gcd_MplInt_string(f.c_str(), g.c_str());

	std::string respStr(gcd);
	respStr.erase(std::remove(respStr.begin(), respStr.end(), '\n'), respStr.end());
	respStr.erase(std::remove(respStr.begin(), respStr.end(), '\\'), respStr.end());

	response.emplace_back(std::move(respStr));
	free(gcd);

	return true;
}

bool MapleWorkerThread::triangularizeValidate_mapleWorker(const std::vector<std::string>& request, std::vector<std::string>& response) {

	if (request.size() < 3) {
		 return false;
	}

	const char* inputs[request.size()-2];
	for (size_t i = 2; i < request.size(); ++i) {
		inputs[i-2] = request[i].c_str();
	}

	int isLazard = atoi(request[1].c_str());


	bool ret = triangularizeValidate_MplInt(inputs, request.size()-2, isLazard);

	response.push_back(std::to_string(ret));
	return true;

}
