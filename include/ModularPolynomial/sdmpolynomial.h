#ifndef _SMALLMMPOLYNOMIAL_H_
#define _SMALLMMPOLYNOMIAL_H_

#include "../Polynomial/BPASMultivarPolynomial.hpp"
#include "../FFT/src/modpn.h"
#include "../FFT/src/general_routine.h"
#include "../FFT/src/basic_routine.h"
#include "../FFT/src/DDMP.h"
//#include "../../include/FFT/src/modpn_hfiles/Types.h"
//#if FURER
//#include "src/fft_furer1.h"
//#else
//#include "src/fft_iter1.h"
//#endif


#define SDMPBASESIZE 1024

/**
 * A multivariate polynomial with coefficients in a small prime field using a dense representation.
 */
class SmallPrimeFieldDistributedDenseMultivariateModularPolynomial : public BPASGCDDomain<SmallPrimeFieldDistributedDenseMultivariateModularPolynomial> {
//TODO wrap sfixn in BPASRing class.
//public BPASMultivariatePolynomial<sfixn,SmallPrimeFieldDistributedDenseMultivariateModularPolynomial>,

	private:
		Symbol* names;     // Variable names
		int var;		// Number of variables
		int n;			// Number of terms
		sfixn* coefs;		// Coefficients in Prime field
		int* degs;		// Partial size
		sfixn p;		// Prime
		int* acc;

		void zeros();
		bool isSameRing(const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const;
		bool shrink();
		void plain_multiplication(SmallPrimeFieldDistributedDenseMultivariateModularPolynomial* c, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& a, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const;
		bool fftbased_multiplication(SmallPrimeFieldDistributedDenseMultivariateModularPolynomial* c, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& a, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const;

	public:
		static mpz_class characteristic;
		// static bool isPrimeField;
		// static bool isSmallPrimeField;
        // static bool isComplexField;

		/**
		 * Constructor using a default prime
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial () : var(0), n(0) {
			p = BPASPRIME;
			coefs = new sfixn[1];
			coefs[0] = 0;
			degs = new int[1];
			degs[0] = 0;
			acc = new int[1];
			acc[0] = 1;
			names = new Symbol[1];
			names[0] = "1";
		}

		/**
		 * Constructor with the field
		 *
		 * @param m: The prime
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial (sfixn m) : var(0), n(0), p(m) {
			coefs = new sfixn[1];
			coefs[0] = 0;
			degs = new int[1];
			degs[0] = 0;
			acc = new int[1];
			acc[0] = 1;
			names = new Symbol[1];
			names[0] = "1";
		}

		/**
		 * Constructor with number of variables and terms
		 *
		 * @param v: Number of variables
		 * @param ds: Partial degrees
		 * @param m: The prime
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial (int v, int* ds, sfixn m) : var(v), p(m) {
			degs = new int[var];
			acc = new int[var];
			acc[0] = 1;
			for (int i = 0; i < var; ++i) {
				degs[i] = ds[i];
				if (i < var - 1)
					acc[i+1] = acc[i] * (degs[i] + 1);
			}
			n = acc[var-1] * (degs[var-1] + 1);
			coefs = new sfixn[n];
			zeros();
			names = new Symbol[var+1];
			names[0] = "1";
			for (int i = 1; i <= var; ++i) {
				std::ostringstream convert;
				convert << var - i + 1;
				names[i] = "_";
				names[i] += convert.str();
			}
		}

		/**
		 * Construct with a variable name
		 * such that f(x) = x;
		 *
		 * @param x: The variable name
		 * @param m: The prime
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial (Symbol x, sfixn m) : var(1), n(2), p(m) {
			names = new Symbol[2];
			names[0] = "9";
			names[1] = x;
			degs = new int[1];
			degs[0] = 1;
			acc = new int[1];
			acc[0] = 1;
			coefs = new sfixn[2];
			coefs[0] = 0;
			coefs[1] = 1;
		}

		/**
		 * Copy constructor
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) : var(b.var), n(b.n), p(b.p) {
			degs = new int[var];
			std::copy(b.degs, b.degs+var, degs);
			acc = new int[var];
			std::copy(b.acc, b.acc+var, acc);
			coefs = new sfixn[n];
			std::copy(b.coefs, b.coefs+n, coefs);
			names = new Symbol[var+1];
			std::copy(b.names, b.names+var+1, names);
		}

		/**
		 * Deconstructor
		 *
		 * @param
		 **/
		~SmallPrimeFieldDistributedDenseMultivariateModularPolynomial() {
			delete [] coefs;
			delete [] degs;
			delete [] names;
			delete [] acc;
		}

		/**
		 * Overload operator =
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator= (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) {
			if (this != &b) {
				delete [] coefs;
				delete [] degs;
				delete [] names;
				delete [] acc;

				var = b.var;
				degs = new int[var];
				std::copy(b.degs, b.degs+var, degs);
				acc = new int[var];
				std::copy(b.acc, b.acc+var, acc);
				n = b.n;
				coefs = new sfixn[n];
				std::copy(b.coefs, b.coefs+n, coefs);
				names = new Symbol[var+1];
				std::copy(b.names, b.names+var+1, names);
				p = b.p;
			}
			return *this;
		}

		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator= (const sfixn& b) {
			return *this = SmallPrimeFieldDistributedDenseMultivariateModularPolynomial(b);
		}

		/**
	     * The characteristic of this ring class.
		 */
		mpz_class getCharacteristic() const override {
			return SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::characteristic;
		}

		/**
		 * Is a zero polynomial
		 *
		 **/
		inline bool isZero() const {
			for (int i = 0; i < n; ++i) {
				if (coefs[i] != 0)
					return 0;
			}
			return 1;
		}

		/**
		 * Zero polynomial
		 *
		 * @param
		 **/
		inline void zero() {
			zeros();
		}

		/**
		 * Is polynomial 1
		 *
		 * @param
		 **/
		inline bool isOne() const {
			for (int i = 1; i < n; ++i) {
				if (coefs[i] != 0)
					return 0;
			}
			return (coefs[0] == 1);
		}

		/**
		 * Set polynomial to 1
		 *
		 * @param
		 **/
		inline void one() {
			coefs[0] = 1;
			for (int i = 1; i < n; ++i)
				coefs[i] = 0;
		}

		/**
		 * Is polynomial -1
		 *
		 * @param
		 **/
		inline bool isNegativeOne() const {
			for (int i = 1; i < n; ++i) {
				if (coefs[i] != 0)
					return 0;
			}
			return (coefs[0] == p - 1);
		}

		/**
		 * Set polynomial to -1
		 *
		 * @param
		 **/
		inline void negativeOne() {
			coefs[0] = p - 1;
			for (int i = 1; i < n; ++i)
				coefs[i] = 0;
		}

		/**
		 * Is a constant
		 *
		 * @param
		 **/
		inline int isConstant() const {
			for (int i = 1; i < n; ++i) {
				if (coefs[i] != 0)
					return 0;
			}
			if (coefs[0] < (p >> 1)) { return 1; }
			else { return -1; }
		}

		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial unitCanonical(SmallPrimeFieldDistributedDenseMultivariateModularPolynomial* u = NULL,
																						  SmallPrimeFieldDistributedDenseMultivariateModularPolynomial* v = NULL) const {
			sfixn leadCoef = this->leadingCoefficient();
			mpz_class mpzLead(leadCoef);
			mpz_class mpzPrime(p);

			mpz_t temp;
			mpz_init(temp);
			if (!mpz_invert(temp, mpzLead.get_mpz_t(), mpzPrime.get_mpz_t()))
				mpz_set_si(temp, 0);
			sfixn leadInv = mpz_get_si(temp);
			mpz_clear(temp);

			SmallPrimeFieldDistributedDenseMultivariateModularPolynomial ret = *this * leadInv;
			if (u != NULL) {
				*u = leadCoef;
			}
			if (v != NULL) {
				*v = leadInv;
			}
			return ret;
		}

		/**
		 * Get the number of variables
		 *
		 * @param
		 **/
		inline int numberOfVariables() const {
			return variables().size();
		}

		/**
		 * Get the number of variables in this polynomial ring.
		 */
		inline int numberOfRingVariables() const {
			return ringVariables().size();
		}

		/**
		 * Get the number of non-zero terms
		 *
		 * @param
		 **/
		inline Integer numberOfTerms() const {
			int k = 0;
			for (int i = 0; i < n; ++i)
				if (coefs[i] != 0) { k++; }
			return k;
		}

		/**
		 * Get the size of the polynomial
		 *
		 * @param
		 **/
		inline int size() const {
			return n;
		}

		/**
		 * Total degree.
		 */
		inline Integer degree() const {
			std::cerr << "SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::degree() NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return 0;
		}

		/**
		 * Get a partial degree of variable x
		 *
		 * @param x: The variable name
		 **/
		inline Integer degree(const Symbol& x) const {
			int k = 0, d = 0;
			for (int i = 1; i <= var; ++i) {
				if (names[i] == x) {
					k = i;
					break;
				}
			}
			if (k) {
				k--;
				for (int i = 0; i < n; ++i) {
					int e = (i / acc[k]) % (degs[k] + 1);
					if (coefs[i] != 0 && e > d)
						d = e;
				}
			}
			return d;
		}

		/**
		 * Get the leading coefficient
		 *
		 * @param
		 **/
		inline sfixn leadingCoefficient() const {
			for (int i = n-1; i > -1; --i) {
				if (coefs[i] != 0)
					return coefs[i];
			}
			return 0;
		}

		inline sfixn trailingCoefficient() const {
			for (int i = 0; i < n; ++i) {
				if (coefs[i] != 0) {
					return coefs[i];
				}
			}
			return 0;
		}

		inline bool isConstantTermZero() const {
			return (coefs[0] == 0);
		}

		/**
		 * Get a coefficient
		 *
		 * @param v: Number of variables
		 * @param d: The exponent of each variable
		 **/
		inline sfixn coefficient(int v, const int* d) const {
			int k = 0;
			for (int i = var-1, j = 0; i > -1 && j < v; --i, ++j) {
				if (d[j] <= degs[i])
					k += d[j] * acc[i];
				else { return 0; }
			}
			for (int i = v-var-1; i > -1; --i)
				if (d[i]) { return 0; }
			return coefs[k];
		}

		inline sfixn coefficient(const std::vector<int>& v) const {
			return coefficient(v.size(), v.data());
		}

		/**
		 * Set a coefficient
		 *
		 * @param v: Number of variables
		 * @param d: The exponent of each variable
		 * @param val: Value of the coefficient
		 **/
		inline void setCoefficient(int v, const int* d, const sfixn& val) {
			if (v != var) {
				std::cout << "BPAS: error, SFDDMMP(" << var << "), but trying to setCoefficient with " << v << " variables." << std::endl;
				exit(1);
			}
			int k = 0;
			for (int i = var-1, j = 0; i > -1 && j < v; --i, ++j) {
				if (d[j] <= degs[i])
					k += d[j] * acc[i];
				else {
					std::cout << "BPAS: error, the degree of " << names[i+1] << " in SFDDMMP is " << degs[i] << "." << std::endl;
					exit(1);
				}
			}
			coefs[k] = val;
		}

		inline void setCoefficient(const std::vector<int>& v, const sfixn& val) {
			setCoefficient(v.size(), v.data(), val);
		}

		/**
		 * Set a coefficient
		 *
		 * @param k: The offset in the coefficient array
		 * @param val: Value of the coefficient
		 **/
		inline void setCoefficient(int k, const sfixn& val) {
			if (k >= n || k < 0) {
				std::cout << "BPAS: error, trying to access a non-exist coefficient in SFDDMMP<Field>." << std::endl;
				exit(1);
			}
			coefs[k] = val;
		}

		/**
		 * Get variable names
		 *
		 * @param
		 **/
		inline std::vector<Symbol> ringVariables() const {
			std::vector<Symbol> xs;
			for (int i = var; i > 0; --i)
				xs.push_back(names[i]);
			return xs;
		}

		/**
		 * Set variable names
		 *
		 * @param xs: Variable names
		 **/
		inline void setRingVariables(const std::vector<Symbol>& xs) {
			int ns = xs.size();
			if (ns != var) {
				std::cerr << "BPAS ERROR: SDMP shrinking and expanding polynomial ring NOT YET IMPLEMENTED" << std::endl;
				return;
			}
			names[0] = "9";
			for (int i = var, j = 0; i > 0 && j < ns; --i, ++j)
				names[i] = xs[j];
		}

		inline std::vector<Symbol> variables() const {
			std::cerr << "BPAS ERROR: SDMP::variables() NOT YET IMPLEMENTED" << std::endl;
			return ringVariables();
		}

		/**
		 * Convert current object to its k-th derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the derivative, k > 0
		 **/
    	void differentiate(const Symbol& s, int k) {}

		/**
		 * Convert current object to its derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/
    	inline void differentiate(const Symbol& s) {
    		this->differentiate(s,0);
    	}

		/**
		 * Return k-th derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the k-th derivative, k > 0
		 **/
		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial derivative(const Symbol& s, int k) const {
    		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial a(*this);
    		a.differentiate(s,k);
	 	return a;
        	 }
		/**
		 * Compute derivative
		 *
		 * @param s: Symbol to differentiate with respect to
		 **/
    	inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial derivative(const Symbol& s) const {
    	 	return this->derivative(s,0);
    	 }

		/**
		 * Evaluate f(x)
		 *
		 * @param syms: Array of Symbols to evaluate at corresponding xs
		 * @param xs: Evaluation points
		 **/
    	inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial evaluate(int, const Symbol* syms, const sfixn* xs) const {
    		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial a;
    		return a;
    	}

		/**
		 * Evaluate f(x)
		 *
		 * @param syms: Vector of Symbols to evaluate at corresponding xs
		 * @param xs: Corresponding evaluation points
		 **/
    	inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial evaluate(const std::vector<Symbol>& syms, const std::vector<sfixn>& xs) const {
    		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial a;
    		return a;
    	}

		/**
		 * Overload operator ==
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		bool operator== (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const;

		/**
		 * Overload operator !=
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline bool operator!= (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const {
			return !(*this == b);
		}

		/**
		 * Overload operator +
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator+ (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const;

		/**
		 * Overload operator +=
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator+= (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) {
			*this = *this + b;
			return *this;
		}

		/**
		 * Overload operator +
		 *
		 * @param e: A constant
		 **/
		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator+ (const sfixn& e) const {
			SmallPrimeFieldDistributedDenseMultivariateModularPolynomial r (*this);
			return (r += e);
		}

		inline friend SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator+ (const sfixn& e, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& f) {
			return (f + e);
		}

		/**
		 * Overload operator +=
		 *
		 * @param e: A constant
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator+= (const sfixn& e);

		/**
		 * Overload operator -
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator- (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const;

		/**
		 * Overload operator -=
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator-= (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) {
			*this = *this - b;
			return *this;
		}

		/**
		 * Overload operator -
		 *
		 * @param e: A constant
		 **/
		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator- (const sfixn& e) const {
			SmallPrimeFieldDistributedDenseMultivariateModularPolynomial r (*this);
			return (r -= e);
		}

		inline friend SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator- (const sfixn& e, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& r) {
			return (-r + e);
		}

		/**
		 * Overload operator -=
		 *
		 * @param e: A constant
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator-= (const sfixn& e);

		/**
		 * Overload operator -, negate
		 *
		 * @param
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator- () const;

		/**
		 * Negate, f(-x)
		 *
		 * @param
		 **/
		void negate();

		/**
		 * Overload operator *
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator* (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) const;

		/**
		 * Overload operator *=
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator*= (const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& b) {
			*this = *this * b;
			return *this;
		}

		/**
		 * Overload operator *
		 *
		 * @param e: A constant
		 **/
		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator* (const sfixn& e) const {
			SmallPrimeFieldDistributedDenseMultivariateModularPolynomial r(*this);
			return (r *= e);
		}

		inline friend SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator* (const sfixn& e, const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& f) {
			return (f * e);
		}

		/**
		 * Overload operator *=
		 *
		 * @param e: A constant
		 **/
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator*= (const sfixn& e);

		/**
		 * Overload operator ^ for exponentiation.
		 *
		 *
		 */
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator^ (long long int e) const {
			//TODO
			std::cerr << "BPAS ERROR: SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::opeartor^ NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		/**
		 * Overload operator ^ for exponentiation.
		 *
		 *
		 */
		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator^= (long long int e) {
			*this = *this ^ e;
			return *this;
		}

		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator/ (const sfixn& e) const {
			std::cerr << "BPAS ERROR: SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator/(sfixn) NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator/= (const sfixn& e) {
			*this = *this / e;
			return *this;
		}
		sfixn content() const  {
			std::cerr << "BPAS ERROR: SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::content NOT YET IMPLEMENTED" << std::endl;
			return sfixn(1);
		}

		SmallPrimeFieldDistributedDenseMultivariateModularPolynomial primitivePart() const {
			std::cerr << "BPAS ERROR: SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::primitivePart NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		void print(std::ostream& out) const;

		/**
		 * Convert *this to an expression tree
		 *
		 */
		ExpressionTree convertToExpressionTree() const;

		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial operator/(const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& p) const {
			std::cerr << "SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator/ NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& operator/=(const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& p) {
			std::cerr << "SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::operator/= NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		inline SmallPrimeFieldDistributedDenseMultivariateModularPolynomial gcd(const SmallPrimeFieldDistributedDenseMultivariateModularPolynomial& p) const {
			std::cerr << "SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::gcd NOT YET IMPLEMENTED" << std::endl;
			//TODO
			return *this;
		}

		/**
		 * Compute squarefree factorization of *this
		 */
		inline Factors<SmallPrimeFieldDistributedDenseMultivariateModularPolynomial> squareFree() const {
			std::cerr << "SmallPrimeFieldDistributedDenseMultivariateModularPolynomial::squareFree NOT YET IMPLEMENTED" << std::endl;
			//TODO
			std::vector<SmallPrimeFieldDistributedDenseMultivariateModularPolynomial> ret;
			ret.push_back(*this);
			return ret;
		}

};

#endif

