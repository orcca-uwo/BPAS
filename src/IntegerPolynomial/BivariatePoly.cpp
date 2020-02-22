/**
	 Implementation of 'BivariatePoly.h'      
          
	 @author Farnam Mansouri
 */

#include "../../include/IntegerPolynomial/BivariatePoly.h"

BivariatePolynomial::BivariatePolynomial(int newD, int newK, int newM){
	d = newD;
	K = newK;
	M = newM;
	size = K * d;

	//initialize to zero
	coefficients = new sfixn[size]();
}

BivariatePolynomial::BivariatePolynomial(int newD, int newK, int newM, int n){
	d = newD;
	K = newK;
	M = newM;
	size = K * d;

	//initialize to zero
	coefficients = new sfixn[size*n]();
}

void BivariatePolynomial::freeHeap(){
        delete[] coefficients;
}

void BivariatePolynomial::print(){
        std::cout << "\t Size: " << size << std::endl;
	std::cout << "\t d: " << d << std::endl;
	std::cout << "\t K: " << K << std::endl;
	std::cout << "\t M: " << M << std::endl;
	for(int i = 0; i < d; ++i){
		std::cout << "\t\t Partial-Y: " << i << std::endl;
		for(int j = 0; j < K; ++j)
			std::cout << "\t\t\t Partial-X: " << j << "\t Coefficient: " << coefficients[i * K +j] << std::endl;
		std::cout << std::endl;
	}
        std::cout << std::endl;
}

void BivariatePolynomial::setToZero(){
        for(int i=0; i < size; ++i)
                coefficients[i] = 0;
}

void BivariatePolynomial::convertFromBigIntegerSigned(mpz_class coeff, int index){

        mpz_class absValue = abs(coeff);
        bool isNegative = coeff < 0;

        std::string str = absValue.get_str(2);
        int start = str.size() % M;
        int count =  ceil((float)str.size() / M);
        int stIndex = index * K;
        stIndex += (K <= count) ? 1 : K - count + 1;

	if(start == 0) start = M;

        int tmp = mpz_class(str.substr(0, start), 2).get_ui();
	coefficients[stIndex - 1] = isNegative ? -tmp : tmp;

        for(int i = stIndex; i < K * (index + 1); ++i) {
                tmp = mpz_class(str.substr((i - stIndex) * M + start, M), 2).get_ui();
                coefficients[i] = isNegative ? -tmp : tmp;
        }

}

mpz_class BivariatePolynomial::coefficientReconstruction(int index){
	int startIndex = K * index;
	mpz_class result = coefficients[startIndex];

	for(int i = 0; i < K; i++){
		result <<= M;
		result += coefficients[startIndex + i];
	}

	return result;
}

int * BivariatePolynomial::adapt(int prime){

  int * adapted = new int[size]; //reuse coefficients?

	//TODO check wheather "2^M < prime" is always true!
	// if not, we have to change the modularization...
	for(int i = 0; i < size; i++)
		adapted[i] = coefficients[i] >= 0 ? coefficients[i] 
						: (prime + coefficients[i]);

	return adapted;
}

void BivariatePolynomial::writeToFile(const char* name) {
	std::ofstream ofs(name, std::ofstream::out);
	for(int j = 0; j < size; ++j)
		ofs << coefficients[j] << " ";
	ofs.close();
}

void BivariatePolynomial::writeToFile(const char* name, int n) {
	std::ofstream ofs(name, std::ofstream::out);
	for(int i = 0; i < n; ++i) {
		for (int j = 0; j < size; ++j)
			ofs << coefficients[i*size+j] << " ";
		ofs << "\n";
	}
	ofs.close();	
}
