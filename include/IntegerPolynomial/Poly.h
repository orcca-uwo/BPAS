#pragma once
/**
	Polynomial class for representing Univariate Integer Polynomial in the densed-fashion which is an
	array of integers using GMP Integer Objects.
	
	This class provide basic functions for polynomials, plus functions for doing the Kronecker Substitution.

	@author Farnam Mansouri
*/

#include "../global.h"

class UnivariateIntegerPolynomial{

	#define MAX(a, b) (a > b) ? a : b
	#define MIN(a, b) (a < b) ? a : b

	#define RESULT_SIZE(aSize, bSize) (aSize == 0 || bSize == 0) ? 0 : (aSize + bSize - 1)

	private:
	
		int size;// Degree + 1, The size of the allocated array.
		int exposedSize; // The real size of the polynomial, not the size of the allocated array.
                mpz_class * coefficients;
		int coefficientDigits;
		int bitPackage;
		int representationBase;
		bool negativeLC;
		char variate;

		/**
		 * Initialize some default values for the polynomial.
		 */
		void init(){
			bitPackage = 5;
		        representationBase = 1 << bitPackage;
        		negativeLC = false;
			variate = 'x';
		}

		/**
		 * Convert from the GMP's big integer to the polynomial representation.
		 *
		 * @param integerRepresentation The GPM integer corresponding to the result polynomial
		 * @param startIndex starting index for the coefficients of the polynomial
		 * @param lastIndex Last index for the coefficients of the polynomial
		 * @param digitCount maximum number of the digits for the coefficients
		 */
		void convertFromBigIntegerSigned(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount);

		/**
		 * Convert from the GMP's big integer to the polynomial representation. This is for when the LC was negative. 
		 *
		 * @param integerRepresentation The GPM integer corresponding to the result polynomial
		 * @param startIndex starting index for the coefficients of the polynomial
		 * @param lastIndex Last index for the coefficients of the polynomial
		 * @param digitCount maximum number of the digits for the coefficients
		 */
        void convertFromBigIntegerSignedNegative(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount);

        public:


		/******************************* CONSTRUCTORS **************************************/

		/**
		 * Polynomial Constructor.
		 *
		 * @param s size
		 * @param cDigits Maximum # of bits for the coefficients
		 * @param c the pointer to the coefficients
		 */
		UnivariateIntegerPolynomial(int s, int cDigits, mpz_class *c){
			init();
			size = s;
		        coefficientDigits = cDigits;
			coefficients = c;
			computeExposedSize();
		}

		/**
		 * Polynomial Constructor.
		 * This will instantiate new array for the coefficients...
		 *
		 * @param s size
		 * @param cDigits Maximum # of bits for the coefficients
		 */
		UnivariateIntegerPolynomial(int s, int cDigits){
			init();
		        size = s;
        		coefficientDigits = cDigits;
		        coefficients = new mpz_class[size](); //initialize to zero
			exposedSize = 0;
		}

		/**
		 * Polynomial Constructor.
		 * This will instantiate new array for the coefficients...
		 *
		 * @param s size
		 */
		UnivariateIntegerPolynomial(int s){
			init();
        		size = s;
		        coefficients = new mpz_class[size](); //initialize to zero
			exposedSize = 0;
		}

		/**
		 * Polynomial Constructor.
		 *
		 * @param s size
		 * @param c the pointer to the coefficients
		 */
		UnivariateIntegerPolynomial(int s, mpz_class *c){
			init();
			size = s;
			coefficients = c;
			computeExposedSize();
		}

		/**
		 * Polynomial Constructor.
		 *
		 * @param c the pointer to the coefficients
		 */
		UnivariateIntegerPolynomial(mpz_class *c){
			init();
			coefficients = c;
		}

		/**
		 * Default Constructor
		 */
		UnivariateIntegerPolynomial(){
			init();
			exposedSize = 0;
		}

		/*inline UnivariateIntegerPolynomial& operator= (UnivariateIntegerPolynomial q) {
                        if (this != &q) {
				size = q.size;
				exposedSize = q.exposedSize;
				coefficients = q.coefficients;
			}
                        return *this;
                }*/

		/**
		 * Free the heap memory used by the coefficients.
		 * To be called manually when using the object is finished!
		 */
		inline void freeHeap(){
			delete[] coefficients;
		}

		/******************************* Access to Properties **************************************/

		/**
		 * Get the real size of the polynomial.
		 */
		inline int getExposedSize() {
			return exposedSize;
		}

		/**
		 * Get the address of a specific element of the polynomial.
		 *
		 * @param index starting index of the coefficient
		 */
		inline mpz_class * getPointer(int index){
			return &coefficients[index];
		}

		/**
		 * Get the degree of the polynomial.
		 */
		inline int getDegree(){
			return size - 1;
		}

		/**
		 * Set the degree of the polynomial. 
		 * It also reallocate the array of coefficients based on the passed degree.
		 *
		 * @param d degree
		 */
		inline void setDegree(int d){
			size = d + 1;
			coefficients = new mpz_class[size]();
		}

		/**
		 * Get the size of the polynomial.
		 */
		inline int getSize(){
			return size;
		}

		/**
		 * Set the degree of the polynomial.
		 * It also reallocate the array of coefficients based on the passed size.
		 *
		 * @param s size
		 */
		inline void setSize(int s){
			size = s;
			coefficients = new mpz_class[size]();
		}

		/**
		 * Get the leading coefficient of the polynomial.
		 */
		inline mpz_class getLeadingCoefficient(){
			return coefficients[exposedSize - 1];
		}

		/**
		 * Get a coefficient of the polynomial.
		 *
		 * @param index index of the coefficient.
		 */
		inline mpz_class getCoefficient(int index){
			return coefficients[index];
		}

		/**
		 * Set the value of a coefficient of the polynomial.
		 * It also recompute the exposed size.
		 *
		 * @param index index of the coefficient.
		 * @param value value of the coefficient.
		 */
		inline void setCoefficient(int index, mpz_class value){
			coefficients[index] = value;

			if(value == 0) computeExposedSize();
			else if(index > exposedSize - 1) exposedSize = index + 1;
		}

		/**
		 * Get the maximum number of digits per coefficients.
		 */
		inline int getCoefficientDigits(){
			return coefficientDigits;
		}

		/**
		 * Set the maximum number of digits per coefficients.
		 *
		 * @param cDigits Maximum # of bits for the coefficients
		 */
		inline void setCoefficientDigits(int cDigits){
			coefficientDigits = cDigits;
		}

		/**
		 * Get the bit package representing the coefficient in.
		 */
		inline int getBitPackage(){
			return bitPackage;
		}

		/**
		 * Set the bit package representing the coefficient in.
		 *
		 * @param bPack bit package
		 */
		inline void setBitPackage(int bPack){
			bitPackage = bPack;
		}

		/**
		 * Get the base where the coefficients are represented in.
		 */
		inline int getRepresentationBase(){
			return representationBase;
		}

		/**
		 * Set the base where the coefficients are represented in.
		 *
		 * @param repBase representation base
		 */
		inline void setRepresentationBase(int repBase){
			representationBase = repBase;
		}

		/**
		 * Is the Leading Coefficient negative.
		 */
		inline bool isNegativeLC(){
			return negativeLC;
		}

		/**
		 * Is the Leading Coefficient negative.
		 */
		inline mpz_class* getCoefficients(){
			return coefficients;
		}

		/**
		 * Is the Leading Coefficient negative.
		 *
		 * @param coeff pointer to the coefficients.
		 */
		inline void setCoefficients(mpz_class * coeff){
			coefficients = coeff;
		}

		/**
		 * Evaluate the polynomial at an arbitrary point.
		 *
		 * @param value the point
		 */
		mpz_class evaluate(mpz_class value);

		/**
		 * Set the variable's name for the polynomial.
		 *
		 * @param v variable's name
		 */
		inline void setVariableName(char v) {
			variate = v;
		}

		/**
		 * Get the variable's name of the polynomial.
		 */
		inline char getVariableName(){
			return variate;
		}

		/**
		 * Compute the Taylor-Shift.
		 * This means incrementing the variable (replace 'x' by 'x+1') and computing the coefficients.
		 */
		void taylorShift();

		/******************* Polynomial Arithmetic ****************************************/

		/**
		 * Add the polynomial to another polynomial and return the result polynomial.
		 *
		 * @param p the polynomial to be added to
		 */
		UnivariateIntegerPolynomial add(UnivariateIntegerPolynomial* p);

		/**
		 * Add another polynomial to the polynomial in-place.
		 *
		 * @param p the polynomial which will be added from
		 */
		void addInPlace(UnivariateIntegerPolynomial* p);

		/**
		 * Subtract the polynomial from another polynomial and return the result polynomial.
		 *
		 * @param p the polynomial which will be substracted.
		 */
		UnivariateIntegerPolynomial substract(UnivariateIntegerPolynomial* p);
	
		/**
		 * Substract the polynomial from another polynomial in-place.
		 *
		 * @param p the polynomial which will be substrated
		 */
		void substractInPlace(UnivariateIntegerPolynomial* p);

		/**
		 * Multiply the polynomial to another polynomial and return the result polynomial.
		 *
		 * @param p the polynomial to be multiplied with
		 */
		inline UnivariateIntegerPolynomial multiply(UnivariateIntegerPolynomial* p);
		
		/**
		 * Multiply the polynomial with an integer.
		 *
		 * @param factor the integer
		 */
		UnivariateIntegerPolynomial multiply(mpz_class factor);

		/**
		 * Check if the polynomial is equivalant to another polynomial.
		 *
		 * @param p the other polynomial
		 */
		bool isEqual(UnivariateIntegerPolynomial* p);

		/**
		 * Check if the polynomial has the root of 0. 
		 */
		inline bool isTrailingCoefficientZero() {
			return (coefficients[0] == 0);
		}
		
		/**
		 * Copy the coefficients from another polynomial.
		 *
		 * @param p the polynomial to be copied from
		 */
		void copy(UnivariateIntegerPolynomial* p);

		/**
		 * Assign the polynomial to the given polynomial.
		 *
		 * @param p the given polynomial
		 */
		void assign(UnivariateIntegerPolynomial* p);

		/******************* Kronecker Substitution Related *******************************/

		/**
		 * Convert the Polynomial representation to the big integer object in GMP library.
		 * This function assumes that all of the coefficients are positive.
		 *
		 * @param size Maximum number of digits (in binary) of the coefficients in the result polynomial
		 */
		mpz_class getBigIntegerUnsigned(int size);
		/**
		 * Convert the Polynomial representation to the big integer object in GMP library.
		 *
		 * @param size Maximum number of digits (in binary) of the coefficients in the result polynomial
		 */
		mpz_class getBigIntegerSigned(int size);
		void getBigIntegerSigned(mpz_class *r, int size);

		/**
		 * To be documented ""2^{N(L-1)}a(2^{-N})"".
		 *
		 * @param size Maximum number of digits (in binary) of the coefficients in the result polynomial
		 */
		mpz_class getReverseBigIntegerUnsigned(int size);
		
		/**
		 * Convert from the GMP's big integer to the polynomial representation.
		 * This function assumes that all of the coefficients are positive.
		 * 
		 * @param integerRepresentation The GPM integer corresponding to the result polynomial
		 * @param startIndex starting index for the coefficients of the polynomial
		 * @param lastIndex Last index for the coefficients of the polynomial
		 * @param digitCount maximum number of the digits for the coefficients
		 */
		void convertFromBigIntegerUnsigned(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount);
		/**
		 * Convert from the GMP's big integer to the polynomial representation.
		 *
		 * @param integerRepresentation The GPM integer corresponding to the result polynomial
		 * @param startIndex starting index for the coefficients of the polynomial
		 * @param lastIndex Last index for the coefficients of the polynomial
		 * @param digitCount maximum number of the digits for the coefficients
		 */
		void convertFromBigInteger(mpz_class& integerRepresentation, int startIndex, int lastIndex, int digitCount);

		/**
		 * To be documented ""extract from 2^{2N(L-1)}c(2^{-N}) and c(2^N)""
		 *
		 * @param c1Int The first intermediate result
		 * @param c2Int The second intermediate result
		 * @param startIndex starting index for the coefficients of the polynomial
		 * @param lastIndex Last index for the coefficients of the polynomial
		 * @param digitCount maximum number of the digits for the coefficients
		 */
		void extractCoeffsicients(mpz_class c1Int, mpz_class c2Int, int startIndex, int lastIndex, int digitCount);
		
		/******************************* UTIL ********************************************/

		/**
		 * Prints all of the coefficients of the polynomial.
		 */
		void print();
		/**
		 * Write all of the coefficients to a file
		 *
		 * @param fileName Name of the file
		 */
		void write(std::string fileName);
		/**
		 * Read all of the coefficients from a file
		 *
		 * @param fileName Name of the file
		 */
		void read(std::string fileName);

		/**
		 * Generate random coefficients for the polynomial.
		 */
		void generateRandom();

		/**
		 * Set all of the coeeficients of the polynomial to 0.
		 */
		void setToZero();
		/**
		 * Set all of the coefficients of the polynomial to 0 except the first one which will be set to 1.
		 */
		void setToOne();
		/**
		 * Copy the coefficients of the polynomial to the input array.
		 *
		 * @param result the target array.
		 */
		void copyFromCoefficients(mpz_class *result);

		/**
		 * Add coefficients of a polynomial from a specific index.
		 *
		 * @param p The polynomial from which the coefficients are going to be added.
		 * @param startIndex The index that the coefficients will be added from.
		 */
		void addCoefficients(UnivariateIntegerPolynomial *p, int startIndex);

		/**
		 * Get the maximum bits required for representing the polynomial's coefficients
		 */
		int maxCoefficientSize();

		//int getBitCount0();

		/**
		 * Negate all of the coefficients of the polynomial in the returning polynomial object.
		 */
		UnivariateIntegerPolynomial negateOutOfPlace();

		/**
		 * Negate all of the coefficients of the polynomial in the returning array.
		 */
		mpz_class * negate();
		/**
		 * Negate all of the coefficients of the polynomial in-place.
		 */
		void negateInPlace();
		/**
		 * Negate the coefficients of the polynomial (in a specific interval) in-place.
		 *
		 * @param startIndex starting index of the interval
		 * @param lastIndex last index of the interval
		 */
		void negateInPlace(int startIndex, int lastIndex);
		/**
		 * Negate all of the coefficients of the polynomial.
		 *
		 * @param result The array including the negated coefficients
		 */
		void negate(mpz_class * result);

		/**
		 * Reduce the degree of the polynomial based on the number of zeros in the most significant coefficients.
		 */
		void fixDegree();
		/**
		 * Grow the coefficients array.
		 *
		 * @param newSize The new size for the polynomial
		 */
		void grow(int newSize);
		/**
		 * Shrink the coefficients array.
		 *
		 * @param newSize The new size for the polynomial
		 */
		void shrink(int newSize);

		/**
		 * Compute the real size of the polynomial and set it to the exposedSize.
		 */
		void computeExposedSize(){
			exposedSize = size;
			while(exposedSize > 1 && coefficients[exposedSize - 1] == 0) exposedSize--;
		}
	

};
