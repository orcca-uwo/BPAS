


#include "MapleInterface/MapleInterfacePipe.hpp"
#include "MapleInterface/MapleInterfaceStream.hpp"

#include "IntegerPolynomial/mzpolynomial.hpp"
#include <unistd.h>
#include <sys/wait.h>

#define MAPLE_PIPE_END_TRANSACTION "__mapleInterfacePipe_End_Transaction"
#define MAPLE_PIPE_END_CONNECTION "__mapleInterfacePipe_End_Connection"

MapleInterfacePipe::MapleInterfacePipe() : child(0), childPID(-1)
{
	// free((int*) -1);

	int result1, result2;
	result1 = pipe(c2p);
	result2 = pipe(p2c);
	if (result1 == -1 || result2 == -1)
	{
		perror("MapleInterfacePipe could not create a pipe.\n");
		exit(1);
	}

	int pid = fork();
	if (pid == 0) {
		child = 1;
		//child
		outstr = new boost::fdostream(c2p[1]);
		instr = new boost::fdistream(p2c[0]);

		//close unused sides of pipe
		close(c2p[0]);
		close(p2c[1]);

		readLoop();

		this->~MapleInterfacePipe();
	}
	else {

		childPID = pid;

		//parent
		outstr = new boost::fdostream(p2c[1]);
		instr = new boost::fdistream(c2p[0]);

		//close unused sides of pipe
		close(p2c[0]);
		close(c2p[1]);

		std::vector<std::string> firstCmd, firstResponse;
		firstCmd.push_back("startKernel");
		sendCommand(firstCmd, firstResponse);
	}

}

MapleInterfacePipe::~MapleInterfacePipe()
{
	if (child) {
		delete outstr;
		delete instr;
		close(c2p[1]);
		close(p2c[0]);
		exit(0); //This fork is done
	} else {
		*outstr << MAPLE_PIPE_END_CONNECTION << std::endl;
		outstr->flush();
		delete outstr;
		delete instr;

		close(p2c[1]);
		close(c2p[0]);

		waitpid(childPID, NULL, 0);
	}
}


void MapleInterfacePipe::readLoop() {
	if (!child) {
		fprintf(stderr, "MapleInterfacePipe parent in read loop!\n");
		exit(1);
	}

	MapleWorkerThread* worker = NULL;
	worker = new MapleWorkerThread();

	if (worker == NULL) {
		fprintf(stderr, "MapleInterfacePipe could not create MapleWorkerThread\n");
		exit(1);
	}

	std::vector<std::string> command, response;

	std::string line;
	while (!(*instr).eof()) {
		getline(*instr, line);

		if (line == MAPLE_PIPE_END_TRANSACTION) {
			response = worker->sendRequestAndWait(command);

			for (std::string s : response) {
				*outstr << s << std::endl;
			}

			*outstr << MAPLE_PIPE_END_TRANSACTION << std::endl;

			command.clear();
			response.clear();

		} else if(line == MAPLE_PIPE_END_CONNECTION) {
			break;
		} else {
			command.push_back(line);	
		}
	}

	delete worker;
}

void MapleInterfacePipe::sendCommand(const std::vector<std::string>& request, std::vector<std::string>& outvec) {

	if (child) {
		fprintf(stderr, "MapleInterfacePipe child in calling sendTransactionAndWait.\n");
		exit(1);
	}

	synchronized(m_mutex) {
		for (std::string s : request) {
			*outstr << s << std::endl;
		}
		*outstr << MAPLE_PIPE_END_TRANSACTION << std::endl;

		std::string line;
		while(!(*instr).eof()) {
			getline(*instr, line);

			if (line == MAPLE_PIPE_END_TRANSACTION) {
				break;
			}

			outvec.push_back(line);
		}
	}

}
