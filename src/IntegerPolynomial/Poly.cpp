/**
	Implementation of 'Poly.h'	
	
	@author Farnam Mansouri
*/

#include "../../include/IntegerPolynomial/Poly.h"

static const std::string zero = "0";
static const int DEFAULT_SIZE = 64;

void UnivariateIntegerPolynomial::taylorShift(){
	//TODO implement
}

mpz_class UnivariateIntegerPolynomial::evaluate(mpz_class value){
	if (exposedSize == 0) return 0; 

	mpz_class eval = coefficients[exposedSize - 1];
	for (int i = exposedSize-2; i > -1; --i)
		eval = eval * value + coefficients[i];
	return eval;
}

void UnivariateIntegerPolynomial::addInPlace(UnivariateIntegerPolynomial* p){
	int minSize = MAX(size, p->getSize());
	for (int i = 0; i < minSize; ++i)
		coefficients[i] += p->getCoefficient(i);
}

UnivariateIntegerPolynomial UnivariateIntegerPolynomial::add(UnivariateIntegerPolynomial* p){

	int resultSize = MAX(size, p->getSize());

	UnivariateIntegerPolynomial result(resultSize);

	for (int i = 0; i < resultSize; ++i) {
		mpz_class value = 0;

		if (i < size)
			value += coefficients[i];
		if (i < p->getSize())
			value += p->getCoefficient(i);

		result.setCoefficient(i, value);
	}

	return result;
}

void UnivariateIntegerPolynomial::substractInPlace(UnivariateIntegerPolynomial* p){
        int minSize = MAX(size, p->getSize());
        for (int i = 0; i < minSize; ++i)
                coefficients[i] -= p->getCoefficient(i);
}

UnivariateIntegerPolynomial UnivariateIntegerPolynomial::substract(UnivariateIntegerPolynomial* p){

        int resultSize = MAX(size, p->getSize());

        UnivariateIntegerPolynomial result(resultSize);

        for (int i = 0; i < resultSize; ++i) {
                mpz_class value = 0;

                if (i < size)
                        value += coefficients[i];
                if (i < p->getSize())
                        value -= p->getCoefficient(i);

                result.setCoefficient(i, value);
        }

        return result;
}




bool UnivariateIntegerPolynomial::isEqual(UnivariateIntegerPolynomial* p) {
	if (size != p->getSize())
		return false;

	for (int i = 0; i < size; i++) {
		if (coefficients[i] != p->getCoefficient(i))
				return false;
	}

	return true;
}

void UnivariateIntegerPolynomial::assign(UnivariateIntegerPolynomial* p){
	if(size > 0) 
		delete [] coefficients;
	variate = p->variate;
	size = p->size;
	exposedSize = p->exposedSize;
	std::cout << size << "\t" << exposedSize << std::endl;
	coefficients = new mpz_class[size];
	std::copy(p->coefficients, p->coefficients + exposedSize, coefficients);
}

void UnivariateIntegerPolynomial::copy(UnivariateIntegerPolynomial* p) {
	for (int i = 0; i < size; i++) {
		if (i < p->getSize())
			coefficients[i] = p->getCoefficient(i);
		else
			coefficients[i] = 0;
	}
}

void UnivariateIntegerPolynomial::print(){
	if(exposedSize > 0)
		std::cout << coefficients[0];
	for(int i=1; i < exposedSize; ++i)
		std::cout << " + " << coefficients[i] << variate << "^" << i;
	std::cout << std::endl;
}


void UnivariateIntegerPolynomial::write(std::string name) {
	std::ofstream ofs(name.c_str(), std::ofstream::out);
	
	for(int j = 0; j < size; ++j) 
		ofs << coefficients[j] << "\n";

	ofs.close();	
}


void UnivariateIntegerPolynomial::read(std::string name){
	std::ifstream ifs(name.c_str(), std::ifstream::in); 

	for(int i = 0; i < size; ++i)
		ifs >> coefficients[i];

	ifs.close();
}



void UnivariateIntegerPolynomial::generateRandom(){

	gmp_randclass rr(gmp_randinit_default);
        rr.seed(time(0));

	//debug
	//rr.seed(1);
	//mpz_class v = abs(rr.get_z_bits(coefficientDigits));
	//mpz_class v = 1;
	for(int i=0; i < size; ++i){
		coefficients[i] = pow(-1, i) * rr.get_z_bits(coefficientDigits);
		if(coefficients[i] != 0) exposedSize++;
	}
	//coefficients[i] = abs(pow(-1, i) * rr.get_z_bits(coefficientDigits));//nonneg
	//coefficients[i] = v;

	/* 
	rr.seed(time(NULL));
	srand (time(NULL));
	for(int i=0; i < size; ++i) {
	  int k = (int) rand() % 2;
	  if (k)
	  	coefficients[i] = rr.get_z_bits(coefficientDigits);
	  else
		coefficients[i] = -rr.get_z_bits(coefficientDigits);
	}
	*/
}


void UnivariateIntegerPolynomial::setToZero(){
	for(int i=0; i < exposedSize; ++i)
                coefficients[i] = 0;
	exposedSize = 0;
}

void UnivariateIntegerPolynomial::setToOne(){
	if(size == 0){
		size = DEFAULT_SIZE;
		coefficients = new mpz_class[size];
		coefficients[0] = 1;
		exposedSize = 1;
		return;
	}
	
	coefficients[0] = 1;	
        for(int i=1; i < exposedSize; ++i)
                coefficients[i] = 0;
        exposedSize = 1;
}


void UnivariateIntegerPolynomial::copyFromCoefficients(mpz_class *result){
	for(int i = 0; i < size; i++)
                result[i] = coefficients[i];
}


void UnivariateIntegerPolynomial::addCoefficients(UnivariateIntegerPolynomial *p, int startIndex){
	for(int i = 0; i < p->size; ++i)
                coefficients[startIndex + i] += p->getCoefficient(i);
}


mpz_class UnivariateIntegerPolynomial::getBigIntegerUnsigned(int s){
	/*size++;
	string resultStr;
        resultStr.reserve(size * size);
        for(int i = 0; i < size; ++i){
		mpz_class tmp;
		cout << "coefficientDigits: " << coefficientDigits << endl;
		if(coefficients[i] < 0)
			tmp = (1 << size) / 2 + (1 << coefficientDigits) + coefficients[i];
		else
	                tmp = coefficients[i];
		//cout << "representationBase: " << representationBase << endl;
		string str = tmp.get_str(representationBase);
		//cout << "str " << i << ": " << str << endl;
                int zeroCount = size - str.size();
                while(zeroCount > 0){
                        zeroCount--;
                        resultStr.append(zero);
                }
                resultStr.append(str);

        }
        return mpz_class(resultStr, representationBase);*/

	std::string resultStr;
	resultStr.reserve(s * size);
	for(int i = 0; i < size; ++i){
		std::string str = coefficients[i].get_str(representationBase);
		int zeroCount = s - str.size();
		resultStr.append(zeroCount, '0');
		resultStr.append(str);

	}
	return mpz_class(resultStr, representationBase);
}


mpz_class UnivariateIntegerPolynomial::getBigIntegerSigned(int s){
	if(size == 0) return mpz_class(0);

	int index = 1;
	while(coefficients[size - index] == 0) index++;

	mpz_class * tmp;
	if(coefficients[size - index] < 0){
                negativeLC = true;
		tmp = negate();
        } else{
		tmp = new mpz_class[size];
		copyFromCoefficients(tmp);
	}


	std::string fString;
	fString.append("1");
	fString.append(s, '0');
	mpz_class F =  mpz_class(fString, representationBase);

	for(int i = 0; i < size - 1; ++i){
                if(tmp[i] < 0){
			tmp[i + 1]--;
			tmp[i] += F;
                }
	}
	
	std::string resultStr;
        resultStr.reserve(s * size);
        for(int i = size; i > 0; i--){
                std::string str = tmp[i - 1].get_str(representationBase);
		//string str = coefficients[i - 1].get_str(representationBase);
                int zeroCount = s - str.size();
                resultStr.append(zeroCount, '0');
                resultStr.append(str);
        }
	delete[] tmp;
        return mpz_class(resultStr, representationBase);
}

void UnivariateIntegerPolynomial::getBigIntegerSigned(mpz_class *r, int s){
	*r = getBigIntegerSigned(s);
}


void UnivariateIntegerPolynomial::convertFromBigIntegerUnsigned(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount){
	/*digitCount++;
        string str = integerRepresentation.get_str(representationBase);
        int start = str.size() % digitCount;
        coefficients[startIndex] += mpz_class(str.substr(0, start), representationBase);

        int range = lastIndex - startIndex;
        for(int i = 1; i < range; ++i){
                string tmp = str.substr((i - 1) * digitCount + start, digitCount);
		mpz_class tmpInt = mpz_class(tmp, representationBase);
		cout << "tmp " << i << ": " << tmp << endl;
		if(tmpInt > (1 << (digitCount - 1)))
			coefficients[startIndex + i] += tmpInt;// - (1 << digitCount) / 2;
		else
                	coefficients[startIndex + i] += tmpInt;
        }*/
	std::string str = integerRepresentation.get_str(representationBase);
	int start = str.size() % digitCount;
	int count = ceil((float)str.size() / digitCount);
	int range = lastIndex - startIndex;
	int stIndex = (range <= count) ? 1 : range - count + 1;
	coefficients[startIndex + stIndex - 1] += mpz_class(str.substr(0, start), representationBase);
	for(int i = stIndex; i < range; ++i)
		coefficients[startIndex + i] += mpz_class(str.substr((i - stIndex) * digitCount + start, digitCount), representationBase);
}


void UnivariateIntegerPolynomial::convertFromBigIntegerSigned(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount){
	std::string fString;
        fString.append("1");
        fString.append(digitCount, '0');
        mpz_class F =  mpz_class(fString, representationBase);
        mpz_class HALF_F = F / 2;

	std::string str = integerRepresentation.get_str(representationBase);
        int start = str.size() % digitCount;
        int count = 	str.size() / digitCount + 1;
        int range = lastIndex - startIndex;

	if(range > count){
		lastIndex -= (range - count);
		range = count;
	}

	mpz_class tmp;
	if(lastIndex > 0 && start > 0){
		tmp = mpz_class(str.substr(0, start), representationBase);
		coefficients[lastIndex - 1] += tmp;
		if(tmp > HALF_F )	coefficients[lastIndex - 1] -= F;
	}

        for(int i = 1; i < range; ++i){
		tmp = mpz_class(str.substr((i - 1) * digitCount + start, digitCount), representationBase);
		coefficients[lastIndex - i - 1] += tmp;
		if(tmp > HALF_F){
			coefficients[lastIndex - i - 1] -= F;
			coefficients[lastIndex - i]++;
		}
	}
}


void UnivariateIntegerPolynomial::convertFromBigIntegerSignedNegative(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount){
        std::string fString;
        fString.append("1");
        fString.append(digitCount, '0');
        mpz_class F =  mpz_class(fString, representationBase);
        mpz_class HALF_F = F / 2;

        std::string str = integerRepresentation.get_str(representationBase);
        int start = str.size() % digitCount;
	int count = str.size() / digitCount + 1;
        int range = lastIndex - startIndex;

        if(range > count){
                lastIndex -= (range - count);
                range = count;
        }

	mpz_class tmp;
        if(lastIndex > 0 && start > 0){
                tmp = mpz_class(str.substr(0, start), representationBase);
                coefficients[lastIndex - 1] -= tmp;
                if(tmp > HALF_F )       coefficients[lastIndex - 1] += F;
        }

        for(int i = 1; i < range; ++i){
                tmp = mpz_class(str.substr((i - 1) * digitCount + start, digitCount), representationBase);
		coefficients[lastIndex - i - 1] -= tmp;
                if(tmp > HALF_F){
                        coefficients[lastIndex - i - 1] += F;
			coefficients[lastIndex - i]--;
                }
        }

}


void UnivariateIntegerPolynomial::convertFromBigInteger(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount){
        if(integerRepresentation > 0){
                convertFromBigIntegerSigned(integerRepresentation, startIndex, lastIndex, digitCount);
	} else{
		integerRepresentation = -integerRepresentation;
                convertFromBigIntegerSignedNegative(integerRepresentation, startIndex, lastIndex, digitCount);
	}
}


mpz_class UnivariateIntegerPolynomial::getReverseBigIntegerUnsigned(int s){
        std::string resultStr;
        resultStr.reserve(s * size);
        for(int i = size; i > 0; --i){
                std::string str = coefficients[i - 1].get_str(representationBase);
                int zeroCount = s - str.size();
		resultStr.append(zeroCount, '0');
                resultStr.append(str);
        }
        return mpz_class(resultStr, representationBase);
}


void UnivariateIntegerPolynomial::extractCoeffsicients(mpz_class c1Int, mpz_class c2Int, int startIndex, int lastIndex, int digitCount){
	int range = lastIndex - startIndex;
	
	mpz_class * h1 = new mpz_class[range + 1];
	mpz_class * h2 = new mpz_class[range + 1];

	std::string h1Str = c1Int.get_str(representationBase);
	std::string h2Str = c2Int.get_str(representationBase);

        int h1Start = h1Str.size() % digitCount;
        h1[0] += mpz_class(h1Str.substr(0, h1Start), representationBase);
	
	int h2Start = h2Str.size() % digitCount;
        h2[0] += mpz_class(h2Str.substr(0, h2Start), representationBase);

        for(int i = 1; i <= range; ++i){
                std::string tmp = h1Str.substr((i - 1) * digitCount + h1Start, digitCount);
                h1[i] += mpz_class(tmp, representationBase);
		
		tmp = h2Str.substr((i - 1) * digitCount + h2Start, digitCount);
                h2[i] += mpz_class(tmp, representationBase);
        }
	
	mpz_class * a = new mpz_class[range + 1];
	mpz_class * b = new mpz_class[range + 1];
	mpz_class * c = new mpz_class[range + 1];
	mpz_class * d = new mpz_class[range + 1];
	a[0] = 0; b[0] = 0; c[0] = 0; d[0] = 0;

	std::string shStr; shStr.append("1");
	for(int i = 0; i < digitCount; ++i)
		shStr.append(zero);
	mpz_class shift = mpz_class(shStr, representationBase);

	for(int i = 0; i < range; ++i){
		c[i + 1] = (h1[range - i ] < b[i] + c[i]);
		a[i + 1] = (h1[range - i ] - b[i] - c[i]) % (shift);
		d[i + 1] = (h2[i + 1] < a[i + 1]);
		b[i + 1] = (h2[i] - a[i] - d[i + 1]) % (shift);
		coefficients[startIndex + range - i - 1] += a[i + 1] + (b[i + 1] * shift);
	}

	delete[] h1; delete[] h2;
	delete[] a; delete[] b; delete[] c; delete[] d;
}

int UnivariateIntegerPolynomial::maxCoefficientSize()
{//-----------------------------------
  size_t max=0, tmp=0;

  for(int i = 0; i < size; ++i){
    tmp = mpz_sizeinbase((coefficients[i]).get_mpz_t(), 2);
    if (tmp>max)
      max = tmp;  
  }
  return max;
}


/*int UnivariateIntegerPolynomial::getBitCount0(){
  //while(size > 1 && coefficients[size - 1] == 0) size--;
//ToDo, find max abs coefficient, then get bit count
  mpz_class max = 0, tmp;
  for(int i = 0; i < size; ++i){
    tmp = abs(coefficients[i]);
    if (tmp>max)
      max = tmp;    
  }
  return mpz_sizeinbase(max.get_mpz_t(), 2);
}*/

void UnivariateIntegerPolynomial::fixDegree(){
	while(size > 1 && coefficients[size - 1] == 0) size--;
}


UnivariateIntegerPolynomial UnivariateIntegerPolynomial::multiply(mpz_class factor){
	UnivariateIntegerPolynomial result(size);
        result.variate = variate;
        for (int i = 0; i < size; ++i)
                result.coefficients[i] = coefficients[i] * factor;
        return result;
}

UnivariateIntegerPolynomial UnivariateIntegerPolynomial::negateOutOfPlace(){
	UnivariateIntegerPolynomial result(size);
        result.variate = variate;
        for (int i = 0; i < size; ++i)
        	result.coefficients[i] = 0 - coefficients[i];
        return result;
}

mpz_class * UnivariateIntegerPolynomial::negate(){
	mpz_class * result = new mpz_class[size];
	negate(result);
	return result;

}


void UnivariateIntegerPolynomial::negate(mpz_class * result){
        for(int i = 0; i < size; ++i)
                        result[i] = 0 - coefficients[i];
}


void UnivariateIntegerPolynomial::negateInPlace(){
        for(int i = 0; i < size; ++i)
                        coefficients[i] = 0 - coefficients[i];
}


void UnivariateIntegerPolynomial::negateInPlace(int startIndex, int lastIndex){
        for(int i = startIndex; i < lastIndex; ++i)
                        coefficients[i] = 0 - coefficients[i];
}


void UnivariateIntegerPolynomial::grow(int newSize){
	if(newSize <= size) return;
	mpz_class * newCoefficients = new mpz_class[newSize];
	for(int i = 0; i < size; ++i)
		newCoefficients[i] = coefficients[i];
	freeHeap();
	coefficients = newCoefficients;
	size = newSize;
}


void UnivariateIntegerPolynomial::shrink(int newSize){
	if(newSize >= size) return;
	mpz_class * newCoefficients = new mpz_class[newSize];
	for(int i = 0; i < newSize; ++i)
		newCoefficients[i] = coefficients[i];
	freeHeap();
	coefficients = newCoefficients;
	size = newSize;
}
