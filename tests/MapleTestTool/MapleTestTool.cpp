

#include "MapleTestTool.hpp"

MapleTestTool* MapleTestTool::_mapleInstance = NULL;

void MapleTestTool::mapleKernelTextCB(void* data, int tag, const char* output) {
	std::cout << "Maple output: " << output << std::endl;
}

void MapleTestTool::mapleKernelErrorCB(void* data, long int tag, const char* output) {
	std::cerr << std::endl << std::endl << "MAPLE ERROR: " << output << std::endl;
	exit(1);
}

/**
 * Construct a maple test tool. This starts up the maple kernel.
 */
MapleTestTool::MapleTestTool() {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	char err[2048];
	MCallBackVectorDesc cb = { mapleKernelTextCB,
								mapleKernelErrorCB,
								0,
								0,
								0,
								0,
								0,
								0
							};
	if ( (kv=StartMaple(0, NULL, &cb, NULL, NULL, err)) == NULL ) {
		std::cerr << "MapleTestTool ERROR: Cannot start maple kernel " << err << std::endl;
		exit(1);
	}
#else
	kv = NULL;
#endif
}

/**
 * Construct a maple test tool using argv parameters.
 */
MapleTestTool::MapleTestTool(int argc, char* argv[]) {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	char err[2048];
	MCallBackVectorDesc cb = { mapleKernelTextCB,
								mapleKernelErrorCB,
								0,
								0,
								0,
								0,
								0,
								0
							};
	if ( (kv=StartMaple(argc, argv, &cb, NULL, NULL, err)) == NULL ) {
		std::cerr << "MapleTestTool ERROR: Cannot start maple kernel " << err << std::endl;
		exit(1);
	}	
#else 
	kv = NULL;
#endif
}

/**
 * Destructor. Closes the connection to the maple kernel.
 */
MapleTestTool::~MapleTestTool() {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	StopMaple(kv);
#endif
}

/**
 * Get a pointer to the MapleTestTool object.
 */
MapleTestTool* MapleTestTool::getMapleTestTool() {
	if (_mapleInstance == NULL) {
		_mapleInstance = new MapleTestTool();
	}
	return _mapleInstance;
}

MKernelVector MapleTestTool::getMKernelVector() {
	return kv;
}

/**
 * Restart the maple session.
 */
void MapleTestTool::restartMapleKernel() {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	char err[2048];
	bool restart = RestartMaple(kv, err);
	if (RestartMaple(kv, err) == 0) {
		std::cerr << "MapleTestTool ERROR: Cannot restart the marple kernel " << err << std::endl;
		exit(1);
	}
#endif
}

/**
 * Helper function for ALGEB to string.
 */
std::string MapleTestTool::algebToString(MKernelVector kv, ALGEB in) const {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	char printStr[] = "sprintf:";
	ALGEB f = EvalMapleStatement(kv, printStr);
	ALGEB resStrObj = EvalMapleProc(kv, f, 2, ToMapleString(kv, "%a"), in);
	char* resStr = MapleToString(kv, resStrObj);
	return std::string(resStr);
#else
	return "";
#endif
}

ALGEB MapleTestTool::expressionTreeToAlgeb(const ExpressionTree& tree) const {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	std::string str = tree.toMapleString();
	str += ":";
	ALGEB f = EvalMapleStatement(kv, const_cast<char*>(str.c_str()));
	return f;
#else
	return NULL;
#endif
}

bool MapleTestTool::testEquality(ALGEB test, const std::string& other, const std::string& verifyMethod) const {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	std::string otherColon = other + ":";
	ALGEB otherALGEB = EvalMapleStatement(kv, const_cast<char*>(otherColon.c_str()));
	return testEquality(test, otherALGEB, verifyMethod);
#else
	return 1;
#endif
}

bool MapleTestTool::testEquality(ALGEB test, ALGEB other, const std::string& verifyMethod) const {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	char expandFunc[] = "expand:";
	ALGEB expF = EvalMapleStatement(kv, expandFunc);
	ALGEB expTest = EvalMapleProc(kv, expF, 1, test);
	ALGEB expOther = EvalMapleProc(kv, expF, 1, other);

	char compareFunc[] = "verify:";
	ALGEB cmpF = EvalMapleStatement(kv, compareFunc);
	

	if (verifyMethod != "") {
		std::string verStr = verifyMethod + ":";
		ALGEB verMthd = EvalMapleStatement(kv, const_cast<char*>(verStr.c_str()));
		ALGEB comp = EvalMapleProc(kv, cmpF, 3, expTest, expOther, verMthd);
		M_BOOL compBool = MapleToM_BOOL(kv, comp);
		return (bool) compBool;
	
	} else {
		ALGEB comp = EvalMapleProc(kv, cmpF, 2, expTest, expOther);
		M_BOOL compBool = MapleToM_BOOL(kv, comp);
		return (bool) compBool;
		
	}
#else
	return 1;
#endif
}


/**
 * Determine if an expression tree evaluates to 0.
 */
bool MapleTestTool::testIfZero(const ExpressionTree& testTree, std::string* retIfNotZero) const {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	std::string exprStr = testTree.toMapleString();
	exprStr += ":"; //to end the statement;
	char* cstr = new char[exprStr.length()+1];
  	std::strcpy (cstr, exprStr.c_str());
	ALGEB result = EvalMapleStatement(kv, cstr);

  	char expandStr[] = "expand:";
  	ALGEB f = EvalMapleStatement(kv, expandStr);
	result = EvalMapleProc(kv, f, 1, result);

	char simplifyStr[] = "simplify:";
	f = EvalMapleStatement(kv, simplifyStr);
	result = EvalMapleProc(kv, f, 1, result);

	delete[] cstr;
	if (IsMapleUnnamedZero(kv, result)) {
		*retIfNotZero = "";
		return true;
	} else {
		*retIfNotZero = algebToString(kv, result);
		return false;
	}
#else
	return 1;
#endif
}

/**
 * Test the provided proc with the supplied inputs and compare them to the 
 * expectedResult. 
 * All strings must be in maple format. This proc automatically ends maple statements
 * with ":" as appropriate.
 * Maximum length of the inputs vector is 6.
 * If the result is not as expected, retIfNotExpected is set to the value of 
 * the proc evaluted with the inputs specified. Otherwise it is the empty string.
 *
 * returns true iff the maple proc evalutes the inputs to be equal to expectedResult.
 */
bool MapleTestTool::testProcReturn(const std::string& inProcStr, const std::vector<std::string>& inputs, const std::string& expectedResult, const std::string& verifyMethod, std::string* retIfNotExpected) const {
#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)
	std::string procStr = inProcStr + ":";
	char* cstr = new char[procStr.length()+1];
	std::strcpy (cstr, procStr.c_str());
	ALGEB testProc = EvalMapleStatement(kv, cstr);
	delete[] cstr;

	std::vector<ALGEB> algebList;
	for (int i = 0; i < inputs.size(); ++i) {
		cstr = new char[inputs[i].length()+2];
		std::string tempStr = inputs[i] + ":";
		std::strcpy(cstr, tempStr.c_str());
		ALGEB res = EvalMapleStatement(kv, cstr);
		algebList.push_back(res);
		delete[] cstr;
	}

	//TODO is there an elegant way to do this? 
	ALGEB result;
	switch(inputs.size()) {
		case 1: {
			result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0]);
			break;
		}
		case 2: {
			result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1]);
			break;
		}
		case 3: {
			result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[2]);
			break;
		}
		case 4: {
			result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[2], algebList[3]);
			break;
		}
		case 5: {
			result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[2], algebList[3], algebList[4]);
			break;
		}
		case 6: {
			result = EvalMapleProc(kv, testProc, inputs.size(), algebList[0], algebList[1], algebList[2], algebList[3], algebList[4], algebList[5]);
			break;
		}
	}

	bool compBool = testEquality(result, expectedResult, verifyMethod);
	if (compBool) {
		*retIfNotExpected = "";
		return true;
	}

	*retIfNotExpected = algebToString(kv, result);
	return false;
#else
	return 1;	
#endif
}



#if (!defined(MAPLE_VALIDATE) || MAPLE_VALIDATE==0) && \
    (!defined(WITH_MAPLE) || WITH_MAPLE==0)

EXT_DECL ALGEB M_DECL EvalMapleStatement( MKernelVector kv, const char *statement ) {
	return NULL;
}

EXT_DECL ALGEB M_DECL ToMapleName( MKernelVector kv, const char *n, M_BOOL is_global ) {
	return NULL;
}

EXT_DECL ALGEB M_CDECL EvalMapleProc( MKernelVector kv, ALGEB fn, int nargs, /* ALGEB arg1, ALGEB arg2, */ ... ) {
	return NULL;
}

EXT_DECL M_BOOL M_DECL MapleToM_BOOL( MKernelVector kv, ALGEB s ) {
	return 1;
}


EXT_DECL MKernelVector M_DECL StartMaple( int argc, char *argv[],
                MCallBackVector cb, void *user_data, void *info, char *errstr ) {
	return NULL;
}

EXT_DECL void M_DECL StopMaple( MKernelVector kv ) {

}

EXT_DECL M_BOOL M_DECL RestartMaple( MKernelVector kv, char *errstr ) {
	return 1;
}

EXT_DECL M_BOOL M_DECL IsMapleUnnamedZero( MKernelVector kv, ALGEB s ) {
	return 1;
}

EXT_DECL char* M_DECL MapleToString( MKernelVector kv, ALGEB s ) {
	return NULL;
}

EXT_DECL ALGEB M_DECL ToMapleInteger( MKernelVector kv, M_INT i ) {
	return NULL;
}


#endif
