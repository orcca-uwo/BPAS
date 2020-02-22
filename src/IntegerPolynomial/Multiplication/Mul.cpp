/**
	Implementation of 'Mul.h'

	@author Farnam Mansouri
*/

#include "../../../include/IntegerPolynomial/Multiplication/Mul.h"

void Mul::setBitCounts(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b, UnivariateIntegerPolynomial *c){
	int aBitCount = a->maxCoefficientSize();
	int bBitCount = b->maxCoefficientSize();
	if(aBitCount < 32 || bBitCount < 32){
                a->setBitPackage(1);
		b->setBitPackage(1);
		c->setBitPackage(1);
                a->setRepresentationBase(2);
		b->setRepresentationBase(2);
		c->setRepresentationBase(2);
        }
	a->setCoefficientDigits(aBitCount / a->getBitPackage());
	b->setCoefficientDigits(bBitCount / b->getBitPackage());
}

void Mul::setBitCounts(UnivariateIntegerPolynomial *a, UnivariateIntegerPolynomial *b){
        int aBitCount = a->maxCoefficientSize();
        int bBitCount = b->maxCoefficientSize();
        if(aBitCount < 32 || bBitCount < 32){
                a->setBitPackage(1);
                b->setBitPackage(1);
                a->setRepresentationBase(2);
                b->setRepresentationBase(2);
        }
        a->setCoefficientDigits(aBitCount / a->getBitPackage());
        b->setCoefficientDigits(bBitCount / b->getBitPackage());
}


void Mul::setNumberOfWorkers(int nworkers){
	std::stringstream strs;
	strs << nworkers;
        std::string temp_str = strs.str();
        char* char_type = (char*) temp_str.c_str();
        __cilkrts_set_param("nworkers", char_type);
}
