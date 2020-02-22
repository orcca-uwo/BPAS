
#ifndef _MAPLE_TEST_TOOL_HPP_
#define _MAPLE_TEST_TOOL_HPP_

#include <iostream>
#include "../../include/ExpressionTree/ExpressionTree.hpp"
#include <vector>

#if (defined(MAPLE_VALIDATE) && MAPLE_VALIDATE) || (defined(WITH_MAPLE) && WITH_MAPLE)

// #if defined(SERIAL) && SERIAL
#include <maplec.h>
// #else 
// #error "Temporarily, cannot compile maple test tool unless building serially: BPAS_BUILD_SERIAL is ON."
// #include <FakeHeaderToStopCompilaton>
// #endif

#else
#include "fakemaplec.h"
#endif

/**
 * MapleTestTool is a class to help automate testing. 
 * By using expression trees, a generic representation of various mathematical
 * expressions, and access to the maple kernal this class facilitates checking
 * the accuracy of various methods.
 *
 * MapleTestTool is a singleton. This is due to the limitations of Maple itself.
 * Please be aware that all references to the MapleTestTool object are the same, 
 * thus, the underlying maple session is shared. 
 */
class MapleTestTool {
private:
	static MapleTestTool* _mapleInstance;

	MKernelVector kv;

	static void mapleKernelTextCB(void* data, int tag, const char* output);

	static void mapleKernelErrorCB(void* data, long int tag, const char* output);

	/**
	 * Construct a maple test tool. This starts up the maple kernel.
	 */
	MapleTestTool();

	/**
	 * Construct a maple test tool using argv parameters.
	 */
	MapleTestTool(int argc, char* argv[]);

	/**
	 * Destructor. Closes the connection to the maple kernel.
	 */
	~MapleTestTool();

public:

	/**
	 * Get a pointer to the MapleTestTool object.
	 */
	static MapleTestTool* getMapleTestTool();
	
	/**
	 * Get the pointer to the maple kernel vector.
	 */
	MKernelVector getMKernelVector();

	/**
	 * Helper function for ALGEB to string.
	 */
	std::string algebToString(MKernelVector kv, ALGEB in) const;

	/**
	 * Helper function to convert an ExpressionTree to an ALGEB object.
	 */
	ALGEB expressionTreeToAlgeb(const ExpressionTree& tree) const;

	/**
	 * Compare an ALGEB and a string for equality.
	 */
	bool testEquality(ALGEB test, const std::string& other, const std::string& verifyMethod) const;
	
	inline bool testEquality(ALGEB test, const std::string& other) const {
		return testEquality(test, other, "");
	}

	/**
	 * Compare two ALGEBS for equality.
	 */
	inline bool testEquality(ALGEB test, ALGEB other) const {
		return testEquality(test, other, "");
	}

	/**
	 * Compare two ALGEBS for equality using some verification methods.
	 */
	bool testEquality(ALGEB test, ALGEB other, const std::string& verifyMethod) const;

	/**
	 * Restart the maple session.
	 * NOTE: In Maple 2017 this is broken. Variable names will not be reset,
	 * BUT memory will indeed be freed in the kernel. 
	 */
	void restartMapleKernel();

	/**
	 * Determine if an expression tree evaluates to 0.
	 */
	bool testIfZero(const ExpressionTree& testTree, std::string* retIfNotZero) const;

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
	bool testProcReturn(const std::string& inProcStr, const std::vector<std::string>& inputs, const std::string& expectedResult, const std::string& verifyMethod, std::string* retIfNotExpected) const;
	
	inline bool testProcReturn(const std::string& inProcStr, const std::vector<std::string>& inputs, const std::string& expectedResult, std::string* retIfNotExpected) const {
		return testProcReturn(inProcStr, inputs, expectedResult, "", retIfNotExpected);
	}


};






#endif

