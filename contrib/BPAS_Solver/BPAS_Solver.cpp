
#include <bpas.h>

#include <regex>
#include <string>
#include <iostream>
#include <fstream>
#include <ios>

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}
 
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}
 
std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    return ltrim(rtrim(str, chars), chars);
}

void printVariableList(std::vector<Symbol> vars) {
	if (vars.size() == 0) {
		std::cout << "[ ]";
		return; 
	}

	std::cout << "[" << vars[0];
	for (size_t i = 1; i < vars.size(); ++i) {
		std::cout << ", " << vars[i];
	}
	std::cout << "]";
}
void printPolyList(std::vector<SparseMultivariateRationalPolynomial> polys) {
	if (polys.size() == 0) {
		std::cout << "[ ]";
		return;
	}

	std::cout << "[" << polys[0];
	for (size_t i = 1; i < polys.size(); ++i) {
		std::cout << ", " << polys[i];
	}
	std::cout << "]";
}

int main(int argc, char** argv) {

	if (argc < 2) {
		std::cerr << "USAGE: BPAS_Solver <system_file> <Lazard_or_Kalkbrener> <Level_or_Bubble> <MapleValidate>" << std::endl;
		std::cerr << "       Lazard_or_Kalkbrener: 1 for Lazard, 0 for Kalkbrener (default)" << std::endl;
		std::cerr << "       Level_or_Bubble: 1 level-wise solving, 0 for bubble (default)" << std::endl;
		std::cerr << "       MapleValidate: 1 for validation against maple, 0 (default) no validation" << std::endl;
		return 1;
	}

	bool Lazard = 0;
	if (argc > 2 && atoi(argv[2])) {
		Lazard = atoi(argv[2]);
	}
	bool levelwise = 0;
	if (argc > 3 && atoi(argv[3])) {
		levelwise = atoi(argv[3]);
	}
	bool mapleValidate = 0;
	if (argc > 4 && atoi(argv[4])) {
		mapleValidate = atoi(argv[4]);
	}

	std::string fname = std::string(argv[1]);
	std::ifstream sysFile;
	sysFile.open(fname); 

	if (!sysFile.is_open() || !sysFile.good()) {
		std::cerr << "Failed to open the file: " << fname << std::endl;
		exit(1);
	}

	std::string varsString;
	while (std::getline(sysFile, varsString)) {
		if (varsString.length() > 0 && varsString.find("#") == std::string::npos) {
			break; 
		}
	}

	std::vector<std::string> polyStrings;
	std::string polyString;
	while (std::getline(sysFile, polyString)) {
		if (polyString.length() == 0 || polyString.find("#") != std::string::npos) {
			continue;
		}

		trim(polyString);
		polyStrings.push_back(polyString);
	}
	
	sysFile.close();
	
	if (polyStrings.size() == 0) {
		std::cerr << "ERROR: No polynomials supplied." << std::endl;
		return 1;
	}

	std::cerr << "lazard: " << Lazard << ", levelwise: " << levelwise << std::endl;


	trim(varsString);
	std::vector<Symbol> inputVars;
	size_t pos = 0;
	std::string delimiter = ",";
	std::string token;
	std::regex varRegex("[0-9a-zA-Z_]*");
	while ((pos = varsString.find(delimiter)) != std::string::npos) {
	    token = varsString.substr(0, pos);
	    trim(token);

	    if(!std::regex_match(token, varRegex)) {
	    	std::cerr << "ERROR: Invalid variable. Only the characters a-z, A-Z, 0-9, and _ are supported. Received: " << token << std::endl;
	    	return 1;
	    }

	    inputVars.emplace_back(token);
	    varsString.erase(0, pos + delimiter.length());
	}
	trim(varsString);
	if (varsString.length() != 0) {
		if(!std::regex_match(varsString, varRegex)) {
	    	std::cerr << "ERROR: Invalid variable. Only the characters a-z, A-Z, 0-9, and _ are supported. Received: " << varsString << std::endl;
	    	return 1;
	    }
		inputVars.emplace_back(varsString);
	}

	std::vector<SparseMultivariateRationalPolynomial> inputPolys;
	for (auto s : polyStrings) {
		inputPolys.emplace_back(s);
		inputPolys.back().setRingVariables(inputVars);
	}

	if (inputPolys.size() == 0) {
		std::cout << "Variables: ";
		printVariableList(inputVars);
		std::cout << "Result: " << std::endl;
		printPolyList(inputPolys);
		return 0;
	}

	RegularChain<RationalNumber, SparseMultivariateRationalPolynomial> rc;
	std::vector<RegularChain<RationalNumber,SparseMultivariateRationalPolynomial>> results;

	unsigned long long startTime = 0;
	startTimer(&startTime);
	results = rc.triangularize(inputPolys, Lazard, levelwise);
	float time = 0;
	stopTimer(&startTime, &time);


	std::cout << "Variables: ";
	printVariableList(inputVars);
	std::cout << std::endl;
	std::cout << "Result:" << std::endl;
	std::cout << "["; 
	for (size_t k=0; k<results.size(); ++k) {
		std::cout << results[k];
		if (k!=results.size()-1) {
			std::cout << ", ";
		}
	}
	std::cout << "]" << std::endl;
	std::cerr << "Triangularize Time: " << time << std::endl;

	if (mapleValidate) {

            ExpressionTree FTree;
            FTree.fromVector<SMQP>(inputPolys);
            ExpressionTree rcTree = rc.convertToExpressionTree();
            ExpressionTree RTree;
            RTree.fromVector<Symbol>(inputVars);
            ExpressionTree resultTrees;
            resultTrees.fromVector<RegularChain<RN,SMQP>>(results);
            
            std::vector<std::string> inputs;
            inputs.push_back(FTree.toMapleString() + ":");
            inputs.push_back(rcTree.toMapleString() + ":");
            inputs.push_back(RTree.toMapleString() + ":");
            inputs.push_back(resultTrees.toMapleString() + ":");
            
            MapleInterfaceStream& mis = MapleInterfaceStream::instance();
            bool passed = mis.TriangularizeValidation(Lazard, inputs);
            
			std::cerr << "BPAS Triangularize Time: " << time << std::endl;

            if (passed) {
                std::cerr << "\t\t\t\t\t\t\t PASSED\n";
            } else {
                std::cerr << "\t\t\t\t\t\t\t FAILED\n";
                exit(1);
            }

        }	    
    return 0;
}
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


