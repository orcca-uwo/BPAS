#include "../../include/polynomial.h"
#include "../../include/ring.h"
#include "../../include/RationalNumberPolynomial/urpolynomial.h"
#include "../../include/RingPolynomial/upolynomial.h"
#include "RationalFunction/rationalfunction_integrationprinting.h"


template <class RealField>
std::string fieldToFPS(RealField &a){
	float b = (float) a.get_d();
	std::ostringstream buff;
    buff << b;
    if (buff.str() == "1")
    	return "";
    else
		return buff.str();
}

float absval(float a){
	if (a >= 0)
		return a;
	else
		return -a;
}

/*template <>
std::string fieldToFPS<ComplexRationalNumber>(ComplexRationalNumber &c){
	float a = (float) c.realPart().get_d();
	float b = (float) c.imaginaryPart().get_d();
	std::ostringstream buff;	
	if (b == 0) {
		buff << a;
	}
	else if (a == 0) {
		if (absval(b) == 1) {
			if (b == 1)
				buff << "I";
			else
				buff << "-I";
		}
		else {
			if (b < 0)
				buff << "-I*" << absval(b);
			else
				buff << "I*" << b;
		}
	}
	else {
		buff << "(" << a;
		if (absval(b) == 1) {
			if (b == 1)
				buff << "+I";
			else
				buff << "-I";
		}
		else {
			if (b < 0)
				buff << "-I*" << absval(b);
			else
				buff << "+I*" << b;
		}
		buff << ")";
	}
	return buff.str();
}*/

template <class UnivariatePolynomialOverField, class Field>
std::string UPoFtoFPS(UnivariatePolynomialOverField &A, bool prettyPrinting){
	if (A.isZero())
    	return "";
    else {
		int lowDeg = 0;
		Integer deg = A.degree();
		if (deg > 0 && A.coefficient(0) == 0){
			int j = 1;
			while (lowDeg == 0){
				if (!(A.coefficient(j) == 0))
					lowDeg = j;
				j++;
			}
		}
		Field a;
		//float b;
		std::string variable = A.variable().toString();
		std::ostringstream buff;
		a = A.coefficient(0);
		if (lowDeg == 0){
			//b = (float) a.get_d();
    		//buff << b;
    		std::string sTemp = fieldToFPS<Field>(a);
    		if (sTemp == "")
    			buff << "1";
    		else
    			buff << fieldToFPS<Field>(a);
    	}
    	if (deg > 0){
    		a = A.coefficient(1);
			if (!(a == 0)){
    			if (lowDeg == 0 && a > 0)
    				buff << "+";
				//b = (float) a.get_d();
				std::string sTemp = fieldToFPS<Field>(a);
				if (sTemp == "-1")
					buff << "-";
				else if (sTemp != ""){
					//buff << b;
					buff << sTemp;
					if (!prettyPrinting)
						buff << "*";
				}
    			buff << variable;
    		}
    	}
    	if (deg > 1){
    		for (int i=2; i<deg+1; i++){
    			a = A.coefficient(i);
				if (!(a == 0)){
    				if (lowDeg < i && a > 0 && i < deg+1)
    					buff << "+";
					//b = (float) a.get_d();
					std::string sTemp = fieldToFPS<Field>(a);
					if (sTemp == "-1")
						buff << "-";
					else if (sTemp != ""){
						//buff << b;
    					buff << sTemp;
						if (!prettyPrinting)
							buff << "*";
					}
    				buff << variable << "^" << i;
  	  			}
    		}
    	}
		return buff.str();
	}
}

template <class UnivariatePolynomial>
std::string UPtoString(UnivariatePolynomial &A, bool prettyPrinting, bool vectorDotPrinting){
	std::stringstream buffer;
	std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
	std::cout << A;
	std::string text = buffer.str();
	std::cout.rdbuf( old );
	if (vectorDotPrinting){
		// strings to find in output
		std::string s1 = "^";
		// string to insert into output
		std::string s2 = ".^";
		int index;
		while( (index = text.find(s1,index) ) != std::string::npos){
			text.replace(index,1,s2);
			index += 2;
		}
	}
	else if (prettyPrinting){
		// string to remove from output
		std::string s1 = "*";
		int index;
		while( (index = text.find(s1)) != std::string::npos)
			text.erase(index,s1.length());
	}
	return text;
}

template <class UnivariatePolynomial>
std::string SUPtoString(SparseUnivariatePolynomial<UnivariatePolynomial> &A, bool prettyPrinting){
	std::stringstream buffer;
	std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
	std::cout << A;
	std::string text = buffer.str();
	std::cout.rdbuf( old );
	if (prettyPrinting){
		// strings to remove from output
		std::string s1 = "*";
		
		int index;
		while( (index = text.find(s1)) != std::string::npos)
			text.erase(index,s1.length());
	}
	return text;
}

void _centreEqualSpacing(std::string *a, std::string *b){
	int total = b->length();
	std::string larger;
	std::string smaller;
	bool isSwitched;
	if (a->length() < b->length()){
		total = b->length();
		larger = *b;
		smaller = *a;
		isSwitched = false;
	}
	else {
		total = a->length();
		larger = *a;
		smaller = *b;
		isSwitched = true;
	}
	int spacer = larger.length() - smaller.length();		
	spacer = floor(spacer/2.0);
	for (int i=0; i<spacer; i++){
		smaller.insert(smaller.begin(),' ');
	}
	spacer = total - smaller.length();
	for (int i=0; i<spacer; i++){
		smaller.push_back(' ');
	}
	if (!isSwitched)
		*a = smaller;
	else
		*b = smaller;	
}

template <class UnivariatePolynomialOverField, class Field>
void _printFormalIntegral(UnivariatePolynomialOverField &A, UnivariatePolynomialOverField &D, UnivariatePolynomialOverField &P, std::vector<UnivariatePolynomialOverField> &G, std::vector<UnivariatePolynomialOverField> &U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverField> > &S, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting){
	int defaultWindowSize = 80;
	
	std::string mainVariable(A.variable().toString());
	
	// control variables:
	// for printing linear logarithmic term
	bool existsLinearLog = false;
	for (int i=0; i<S.size(); i++){
		if (S.at(i).degree()==1)
			existsLinearLog = true;
	}
	// for printing plus/minus signs
	bool arePreviousTerms = false;
		
	// code to determine how to centre integrand
	int i = 0, j = 0;
	int numRationalTerms = G.size()/2;
	std::string dText = UPtoString<UnivariatePolynomialOverField>(D,prettyPrinting,vectorDotPrinting);
	std::string aText = UPtoString<UnivariatePolynomialOverField>(A,prettyPrinting,vectorDotPrinting);
	std::string sTemp;
	std::vector<std::string> gNumText;
	std::vector<std::string> gDenText;
	for (i=0; i<numRationalTerms; i++){
		if (floatingPointPrinting)
			sTemp = UPoFtoFPS<UnivariatePolynomialOverField,Field>(G.at(2*i),prettyPrinting);
		else
			sTemp = UPtoString<UnivariatePolynomialOverField>(G.at(2*i),prettyPrinting,vectorDotPrinting);
		gNumText.push_back(sTemp);
		if (floatingPointPrinting)
			sTemp = UPoFtoFPS<UnivariatePolynomialOverField,Field>(G.at(2*i+1),prettyPrinting);
		else
			sTemp = UPtoString<UnivariatePolynomialOverField>(G.at(2*i+1),prettyPrinting,vectorDotPrinting);
		gDenText.push_back(sTemp);
	}
	//std::string c,d;
	
	_centreEqualSpacing(&aText,&dText);
	for (i=0; i<numRationalTerms; i++){
		_centreEqualSpacing(&gNumText[i],&gDenText[i]);
	}
	
	
	// code to display the result of the integration
	std::cout << "  -" << std::endl;
	std::cout << " / " << aText;
	// print numerators of rational terms if any //
	if (G.size() > 0){
		arePreviousTerms = true;
		std::cout << "      ";
		for (j=0; j<numRationalTerms; j++){
			std::cout << gNumText[j];
			if ( !(j == (numRationalTerms-1)) ){
				std::cout << "   ";
			}
		}
	}
	std::cout << std::endl;
	std::cout << " | ";
	for (i=0; i<aText.length(); i++)
		std::cout << "-";
	std::cout << " d" << mainVariable << " =";
	// print division lines of rational terms if any //
	if (G.size() > 0){
		std::cout << " ";
		for (j=0; j<numRationalTerms; j++){
			for (i=0; i<gNumText[j].length(); i++)
				std::cout << "-";
			if ( !(j == (numRationalTerms-1)) ){
				std::cout << " + ";
			}
		}
	}
	// print polynomial term is it exists //
	if (!(P.isZero())){
		if (arePreviousTerms)
			std::cout << " + ";
		if (floatingPointPrinting)
			std::cout << UPoFtoFPS<UnivariatePolynomialOverField,Field>(P,prettyPrinting);
		else
			std::cout << UPtoString<UnivariatePolynomialOverField>(P,prettyPrinting,vectorDotPrinting);
		arePreviousTerms = true;
	}
	// print simple linear log terms if any //
	/*if (existsLinearLog){  // NOT SURE THIS CONDITION IS REQUIRED
		if (U.at(0).degree() == 1){ 
			RationalNumber qTemp = -U.at(0).coefficient(0);
			if (qTemp < 0){
				if (arePreviousTerms)
					std::cout << " -";
			}
			else{
				if (arePreviousTerms)
				std::cout << " +";
			}
			std::cout << " ";
			if (qTemp < 0){
				if (arePreviousTerms){
					qTemp = -qTemp;
					if (qTemp != 1)
						std::cout << qTemp;
				}
				else{
					if (qTemp != -1)
						std::cout << qTemp;
					else
						std::cout << "-";
				}
			}
			else {
				if (qTemp != 1)
					std::cout << qTemp;
			}
			std::cout << "ln(" << S.at(0) << ")";
		}
		arePreviousTerms = true;
	}*/
	// print plus sign if there's something preceding and more to come //
	if (arePreviousTerms && U.size()>0) {
	// ((!existsLinearLog && U.size()>0) || (existsLinearLog && U.size()>1))){
		std::cout << " +";
	}
	std::cout << std::endl;
	std::cout << " / " << dText;
	// print denominators of rational terms if any //
	if (G.size() > 0){
		std::cout << "      ";
		for (j=0; j<numRationalTerms; j++){
			std::cout << gDenText[j];
			if ( !(j == (numRationalTerms-1)) ){
				std::cout << "   ";
			}
		}
	}
	std::cout << std::endl;
	std::cout << "-  " << std::endl;
	// print complex logarithmic terms if any //
	//if ((!existsLinearLog && U.size()>0) || (existsLinearLog && U.at(0).degree()>1)){
	if (U.size()>0){
		std::string indeterminate(U.at(0).variable().toString());
		for (i=0; i<U.size(); i++){
			if (U.at(i).degree() > 0){
				std::cout << "     ----" << std::endl;
				std::cout << "     \\    " << indeterminate << "*ln(" << SUPtoString<UnivariatePolynomialOverField>(S.at(i),prettyPrinting) << ")";
				if (!(i == U.size()-1))
					std::cout << " +";
				std::cout << std::endl;
				std::cout << "     /   " << std::endl;
				std::cout << "     ----" << std::endl;
				std::cout << indeterminate << "|" << UPtoString<UnivariatePolynomialOverField>(U.at(i),prettyPrinting,vectorDotPrinting) << "=0" << std::endl << std::endl;
			}
		}
	}
	else
		std::cout << std::endl;
}

template <class UnivariatePolynomialOverField, class Field>
void _printIntegral(UnivariatePolynomialOverField A, UnivariatePolynomialOverField D, UnivariatePolynomialOverField P, std::vector<UnivariatePolynomialOverField> G, std::vector<Field> h, std::vector<UnivariatePolynomialOverField> H, std::vector<Field> k, std::vector<UnivariatePolynomialOverField> K1, std::vector<UnivariatePolynomialOverField> K2, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting, std::string outputFormatting){
	int defaultWindowSize = 80;
		
	bool exteRationalNumberalOutput = true;
	//bool vectorDot = false;
	// formatting defined globally in global.cpp //
		
	std::string mainVariable(A.variable().toString());
	
	// control variables:
	// for printing linear logarithmic term
	bool existsLogTerm = false;
	bool existsArctanTerm = false;
	bool existsUnwindingNumber = false;
	if (H.size() > 0)
		existsLogTerm = true;
	if (K1.size() > 0)
		existsArctanTerm = true;
	// for printing plus/minus signs
	bool arePreviousTerms = false;
		
	// code to determine how to centre integrand
	int i = 0, j = 0;
	int numRationalTerms = G.size()/2;
	std::string dText(UPtoString<UnivariatePolynomialOverField>(D,prettyPrinting,vectorDotPrinting));
	std::string aText(UPtoString<UnivariatePolynomialOverField>(A,prettyPrinting,vectorDotPrinting));
	std::string sTemp("");
	std::vector<std::string> gNumText;
	std::vector<std::string> gDenText;
	for (i=0; i<numRationalTerms; i++){
		if (floatingPointPrinting)
			sTemp = UPoFtoFPS<UnivariatePolynomialOverField,Field>(G.at(2*i),prettyPrinting);
		else
			sTemp = UPtoString<UnivariatePolynomialOverField>(G.at(2*i),prettyPrinting,vectorDotPrinting);
		if (sTemp == "")
			sTemp = "1";
		gNumText.push_back(sTemp);
		if (floatingPointPrinting)
			sTemp = UPoFtoFPS<UnivariatePolynomialOverField,Field>(G.at(2*i+1),prettyPrinting);
		else
			sTemp = UPtoString<UnivariatePolynomialOverField>(G.at(2*i+1),prettyPrinting,vectorDotPrinting);
		if (sTemp == "")
			sTemp = "1";
		gDenText.push_back(sTemp);
	}
	
	_centreEqualSpacing(&aText,&dText);
	if (!exteRationalNumberalOutput){
		for (i=0; i<numRationalTerms; i++){
			_centreEqualSpacing(&gNumText[i],&gDenText[i]);
		}
	}
	
	
	// code to display the result of the integration
	std::cout << "  -" << std::endl;
	std::cout << " / " << aText;
	if (!exteRationalNumberalOutput) {
		// print numerators of rational terms if any //
		if (G.size() > 0){
			arePreviousTerms = true;
			std::cout << "      ";
			for (j=0; j<numRationalTerms; j++){
				std::cout << gNumText[j];
				if ( !(j == (numRationalTerms-1)) ){
					std::cout << "   ";
				}
			}
		}
	}
	std::cout << std::endl;
	std::cout << " | ";
	for (i=0; i<aText.length(); i++)
		std::cout << "-";
	std::cout << " d" << mainVariable << " =";
	if (!exteRationalNumberalOutput) {
		// print division lines of rational terms if any //
		if (G.size() > 0){
			std::cout << " ";
			for (j=0; j<numRationalTerms; j++){
				for (i=0; i<gNumText[j].length(); i++)
					std::cout << "-";
				if ( !(j == (numRationalTerms-1)) ){
					std::cout << " + ";
				}
			}
		}
		// print polynomial term is it exists //
		if (!(P.isZero())){
			if (arePreviousTerms)
				std::cout << " + ";
			if (floatingPointPrinting)
				std::cout << UPoFtoFPS<UnivariatePolynomialOverField,Field>(P,prettyPrinting);
			else
				std::cout << UPtoString<UnivariatePolynomialOverField>(P,prettyPrinting,vectorDotPrinting);
			arePreviousTerms = true;
		}
	}
	// print plus sign if there's something preceding and more to come //
	//if (arePreviousTerms){
	//	std::cout << " +";
	//}
	std::cout << std::endl;
	std::cout << " / " << dText;
	if (!exteRationalNumberalOutput) {
		// print denominators of rational terms if any //
		if (G.size() > 0){
			std::cout << "      ";
			for (j=0; j<numRationalTerms; j++){
				std::cout << gDenText[j];
				if ( !(j == (numRationalTerms-1)) ){
					std::cout << "   ";
				}
			}
		}
	}
	std::cout << std::endl;
	std::cout << "-  " << std::endl;
	if (exteRationalNumberalOutput) {
		// print polynomial term is it exists //
		if (!(P.isZero())){
			std::cout << " " << UPtoString<UnivariatePolynomialOverField>(P,prettyPrinting,vectorDotPrinting);
			arePreviousTerms = true;
		}
		// print rational terms if any //
		if (G.size() > 0){
			if (arePreviousTerms)
				std::cout << " + ";
			else
				std::cout << " ";
			std::cout << "(" << gNumText[0] << ")";
			if (outputFormatting == "MATLAB_OUTPUT")
				std::cout << ".";
			std::cout << "/(" << gDenText[0] << ")";
			arePreviousTerms = true;
			if (G.size() > 1){
				for (j=1; j<numRationalTerms; j++) {
					std::cout << " + (" << gNumText[j] << ")";
					if (outputFormatting == "MATLAB_OUTPUT")
						std::cout << ".";
					std::cout << "/(" << gDenText[j] << ")";
				}
			}
		}
	}
	// print log terms if any //
	if (existsLogTerm){
		for (i=0; i<h.size(); i++) {
			if (h.at(i) < 0){
				std::cout << " -";
				if (arePreviousTerms)
					std::cout << " ";
				if (h.at(i) != -1) {
					h.at(i) = -h.at(i);
					if (floatingPointPrinting)
						std::cout << fieldToFPS<Field>(h.at(i));
					else
						std::cout << h.at(i);
					if (!prettyPrinting && fieldToFPS<Field>(h.at(i)) != "")
						std::cout << "*";
				}
			}
			else {
				if (arePreviousTerms)
					std::cout << " + ";
				else
					std::cout << " ";
				if (h.at(i) != 1) {
					if (floatingPointPrinting)
						std::cout << fieldToFPS<Field>(h.at(i));
					else
						std::cout << h.at(i);
					if (!prettyPrinting && fieldToFPS<Field>(h.at(i)) != "")
						std::cout << "*";
				}
			}
			if (outputFormatting == "MATLAB_OUTPUT")
				std::cout << "log(abs(";
			else
				std::cout << "ln(";
			if (floatingPointPrinting)
				std::cout << UPoFtoFPS<UnivariatePolynomialOverField,Field>(H.at(i),prettyPrinting) << ")";
			else
				std::cout << UPtoString<UnivariatePolynomialOverField>(H.at(i),prettyPrinting,vectorDotPrinting) << ")";
			if (outputFormatting == "MATLAB_OUTPUT")
				std::cout << ")";
			arePreviousTerms = true;
		}
	}
	// print arctan terms if any //
	if (existsArctanTerm){
		for (i=0; i<k.size(); i++) {
			if (k.at(i) < 0){
				std::cout << " -";
				if (arePreviousTerms)
					std::cout << " ";
				if (k.at(i) != -1) {
					k.at(i) = -k.at(i);
					if (floatingPointPrinting)
						std::cout << fieldToFPS<Field>(k.at(i));
					else
						std::cout << k.at(i);
					if (!prettyPrinting && fieldToFPS<Field>(k.at(i)) != "")
						std::cout << "*";
				}
			}
			else {
				if (arePreviousTerms)
					std::cout << " + ";
				else
					std::cout << " ";
				if (k.at(i) != 1) {
					if (floatingPointPrinting)
						std::cout << fieldToFPS<Field>(k.at(i));
					else
						std::cout << k.at(i);
					if (!prettyPrinting && fieldToFPS<Field>(k.at(i)) != "")
						std::cout << "*";
				}
			}
			if (K2.size() == 0) {
				if (outputFormatting == "MATLAB_OUTPUT")
					std::cout << "atan(";
				else
					std::cout << "arctan(";
			if (floatingPointPrinting)
				std::cout << UPoFtoFPS<UnivariatePolynomialOverField,Field>(K1.at(i),prettyPrinting) << ")";
			else
				std::cout << UPtoString<UnivariatePolynomialOverField>(K1.at(i),prettyPrinting,vectorDotPrinting) << ")";
			}
			else {
			if (outputFormatting == "MATLAB_OUTPUT")
				std::cout << "atan2(";
			else
				std::cout << "arctan(";
			if (floatingPointPrinting)
				std::cout << UPoFtoFPS<UnivariatePolynomialOverField,Field>(K1.at(i),prettyPrinting) << "," << UPoFtoFPS<UnivariatePolynomialOverField,Field>(K2.at(i),prettyPrinting) << ")";
			else
				std::cout << UPtoString<UnivariatePolynomialOverField>(K1.at(i),prettyPrinting,vectorDotPrinting) << "," << UPtoString<UnivariatePolynomialOverField>(K2.at(i),prettyPrinting,vectorDotPrinting) << ")";
			}
			arePreviousTerms = true;
		}
	}
	std::cout << std::endl;
}


// to avoid linking errors
//template std::string fieldToFPS<RationalNumber>(RationalNumber &a);
//template std::string fieldToFPS<RationalNumber>(RationalNumber &a);

template void _printFormalIntegral<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, DenseUnivariateRationalPolynomial &P, std::vector<DenseUnivariateRationalPolynomial> &G, std::vector<DenseUnivariateRationalPolynomial> &U, std::vector< SparseUnivariatePolynomial<DenseUnivariateRationalPolynomial> > &S, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting);
template void _printFormalIntegral<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, SparseUnivariatePolynomial<RationalNumber> &P, std::vector< SparseUnivariatePolynomial<RationalNumber> > &G, std::vector< SparseUnivariatePolynomial<RationalNumber> > &U, std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<RationalNumber> > > &S, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting);
//template void _printFormalIntegral<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(SparseUnivariatePolynomial<ComplexRationalNumber> &A, SparseUnivariatePolynomial<ComplexRationalNumber> &D, SparseUnivariatePolynomial<ComplexRationalNumber> &P, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > &G, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > &U, std::vector< SparseUnivariatePolynomial< SparseUnivariatePolynomial<ComplexRationalNumber> > > &S, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting);


template void _printIntegral<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial A, DenseUnivariateRationalPolynomial D, DenseUnivariateRationalPolynomial P, std::vector<DenseUnivariateRationalPolynomial> G, std::vector<RationalNumber> h, std::vector<DenseUnivariateRationalPolynomial> H, std::vector<RationalNumber> k, std::vector<DenseUnivariateRationalPolynomial> K1, std::vector<DenseUnivariateRationalPolynomial> K2, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting, std::string outputFormatting);
template void _printIntegral<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> A, SparseUnivariatePolynomial<RationalNumber> D, SparseUnivariatePolynomial<RationalNumber> P, std::vector< SparseUnivariatePolynomial<RationalNumber> > G, std::vector<RationalNumber> h, std::vector< SparseUnivariatePolynomial<RationalNumber> > H, std::vector<RationalNumber> k, std::vector< SparseUnivariatePolynomial<RationalNumber> > K1, std::vector< SparseUnivariatePolynomial<RationalNumber> > K2, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting, std::string outputFormatting);
//template void _printIntegral<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(SparseUnivariatePolynomial<ComplexRationalNumber> A, SparseUnivariatePolynomial<ComplexRationalNumber> D, SparseUnivariatePolynomial<ComplexRationalNumber> P, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > G, std::vector<ComplexRationalNumber> h, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > H, std::vector<ComplexRationalNumber> k, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > K1, std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > K2, bool prettyPrinting, bool floatingPointPrinting, bool vectorDotPrinting, std::string outputFormatting);
