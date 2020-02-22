#ifndef _MMPolynomial_H_
#define _MMPolynomial_H_


#include "../polynomial.h"
#include "../Utils/TemplateHelpers.hpp"
#include "../Ring/BPASField.hpp"

template <class Field>
inline void coef_add_mod(Field* c, Field a, Field b, Field p) {
	*c = (a + b) % p;
}
template <class Field>
inline void coef_sub_mod(Field* c, Field a, Field b, Field p) {
	*c = (a - b) % p;
	if (*c < 0) { *c += p; }
}
template <class Field>
inline void coef_neg_mod(Field* c, Field a, Field p) {
	*c = p - a;
}
template <class Field>
inline void coef_mul_mod(Field* c, Field a, Field b, Field p) {
	*c = (a * b) % p;
}

/**
 * A multivariate polynomial with coefficients in an arbitrary finite field represented densely.
 * The class is templated by a Field which should be a BPASField.
 */
template <class Field>
class DistributedDenseMultivariateModularPolynomial : public BPASMultivariatePolynomial<Field,DistributedDenseMultivariateModularPolynomial<Field>>, 
												      private Derived_from<Field, BPASField<Field>> {
	private:
		Symbol* names;     // Variable names
		int var;		// Number of variables
		int n;			// Number of terms
		Field* coefs;		// Coefficients in Prime field
		int* degs;		// Partial size
		Field p;		// Prime


		int* acc;

		inline void zeros() {
			for (int i = 0; i < n; ++i)
				coefs[i] = 0;
		}

		inline bool isSameRing(const DistributedDenseMultivariateModularPolynomial& b) const {
			if (var != b.var || names[0] != b.names[0])
				return 0;
			if (names[0] == "9") {
				for (int i = 1; i <= var; ++i) {
					if (names[i] != b.names[i])
						return 0;
				}
			}
			return 1;
		}

	public:
		static mpz_class characteristic;
		// static bool isPrimeField;
		// static bool isSmallPrimeField;
		// static bool isComplexField;
		/**
		 * Constructor using a default field
		 **/
		DistributedDenseMultivariateModularPolynomial () : var (0),  n(0), p() {
			coefs = new Field[1];
			coefs[0] = 0;
			degs = new int[1];
            degs[0] = 0;
            acc = new int[1];
            acc[0] = 1;
            names = new Symbol[1];
            names[0] = "1";
		}

		/**
		 * Constructor with the Field
		 *
		 * @param m: The prime
		 **/
		DistributedDenseMultivariateModularPolynomial (const Field& m) : var(0), n(0), p(m) {
			coefs = new Field[1];
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
		DistributedDenseMultivariateModularPolynomial(int v, int* ds, Field m) : var(v), p(m) {
			degs = new int[var];
			acc = new int[var];
			acc[0] = 1;
			for (int i = 0; i < var; ++i) {
				degs[i] = ds[i];
				if (i < var - 1)
					acc[i+1] = acc[i] * (degs[i] + 1);
			}
			n = acc[var-1] * (degs[var-1] + 1);
			coefs = new Field[n];
			zeros();
			//coefs[0] = 0;
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
		DistributedDenseMultivariateModularPolynomial (const Symbol& x, const Field& m) : var(1), n(2), p(m) {
                        names = new Symbol[2];
                        names[0] = "9";
                        names[1] = x;
			degs = new int[1];
			degs[0] = 1;
			acc = new int[1];
			acc[0] = 1;
			coefs = new Field[2];
			coefs[0] = 0;
			coefs[1] = 1;
		}
		/**
		 * Copy constructor
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		DistributedDenseMultivariateModularPolynomial(const DistributedDenseMultivariateModularPolynomial<Field>& b) : var(b.var), n(b.n), p(b.p) {
			degs = new int[var];
			std::copy(b.degs, b.degs+var, degs);
			acc = new int[var];
			std::copy(b.acc, b.acc+var, acc);
			coefs = new Field[n];
			std::copy(b.coefs, b.coefs+n, coefs);
			names = new Symbol[var+1];
			std::copy(b.names, b.names+var+1, names);
		}
		/**
		 * Deconstructor
		 *
		 * @param
		 **/
		~DistributedDenseMultivariateModularPolynomial() {
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
		DistributedDenseMultivariateModularPolynomial<Field>& operator= (const DistributedDenseMultivariateModularPolynomial<Field>& b) {
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
				coefs = new Field[n];
				std::copy(b.coefs, b.coefs+n, coefs);
				names = new Symbol[var+1];
				std::copy(b.names, b.names+var+1, names);
				p = b.p;
			}
			return *this;
		}

		/**
		 * Overload operator =. Assign this to a be a base Field element. 
		 */		
		DistributedDenseMultivariateModularPolynomial<Field>& operator= (const Field& f) {
			*this = DistributedDenseMultivariateModularPolynomial(f);
			return *this;
		}

		/**
		 * Is a zero polynomial
		 *
		 * @param
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
			if (coefs[0] < (p / 2)) { return 1; }
			else { return -1; }
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
		 * Get the total degree.
		 */
		inline Integer degree() const {
			std::cerr << "DistributedDenseMultivariateModularPolynomial::degree() NOT YET IMPLEMENTED" << std::endl;
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
				if (names[i] == x)
					k = i;
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
		inline Field leadingCoefficient() const {
			for (int i = n-1; i > -1; --i) {
				if (coefs[i] != 0)
					return coefs[i];
			}
			return 0;
		}

		inline Field trailingCoefficient() const {
			for (int i = 0; i < n; ++i) {
				if (coefs[i] != 0) {
					return coefs[i];
				}
			}
			return 0;
		}

		bool isConstantTermZero() const {
			int d[var];
			Field f = coefficient(var, d);
			return (f == 0);
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> unitCanonical(DistributedDenseMultivariateModularPolynomial* u = NULL, 
																			      DistributedDenseMultivariateModularPolynomial* v = NULL) const {
			Field lead = leadingCoefficient();
			Field leadInv = lead.inverse();
			DistributedDenseMultivariateModularPolynomial ret = *this * leadInv;
			if (u != NULL) {
				*u = lead;
			} 
			if (v != NULL) {
				*v = leadInv;
			}
			return ret;
		}

		/**
		 * Get a coefficient
		 *
		 * @param v: Number of variables
		 * @param d: The exponent of each variable
		 **/
		inline Field coefficient(int v, const int* d) const {
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

		inline Field coefficient(const std::vector<int>& v) const {
			return coefficient(v.size(), v.data());
		}

		/**
		 * Set a coefficient
		 *
		 * @param v: Number of variables
		 * @param d: The exponent of each variable
		 * @param val: Value of the coefficient
		 **/
		inline void setCoefficient(int v, const int* d, const Field& val) {
			if (v != var) {
				std::cout << "BPAS: error, DDMMP(" << var << "), but trying to setCoefficient with " << v << " variables." << std::endl;
				exit(1);
			}
			int k = 0;
			for (int i = var-1, j = 0; i > -1 && j < v; --i, ++j) {
				if (d[j] <= degs[i])
					k += d[j] * acc[i];
				else {
					std::cout << "BPAS: error, the degree of " << names[i+1] << " in DDMMP is " << degs[i] << "." << std::endl;
					exit(1);
				}
			}
			coefs[k] = val;
		}

		inline void setCoefficient(const std::vector<int>& v, const Field& val) {
			setCoefficient(v.size(), v.data(), val);
		}

		/**
		 * Set a coefficient
		 *
		 * @param k: The offset in the coefficient array
		 * @param val: Value of the coefficient
		 **/
		inline void setCoefficient(int k, const Field& val) {
			if (k >= n || k < 0) {
				std::cout << "BPAS: error, trying to access a non-exist coefficient in DDMMP<Field>." << std::endl;
				exit(1);
			}
			coefs[k] = val;
		}

		inline Field content() const {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::content() NOT YET IMPLEMENTED" << std::endl;
			return Field(1); 
		}

		inline DistributedDenseMultivariateModularPolynomial primitivePart() const {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::primitivePart() NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		/**
		 * Overload operator ==
		 *
		 * @param b: A multivariate modular polynomial
		 **/ 
		inline bool operator== (const DistributedDenseMultivariateModularPolynomial<Field>& b) const {
			if (var != b.var || p != b.p)
				return 0;

			int prev = 0;
			for (int i = 0; i < n; ++i) {
				int k = 0;
				for (int j = 0; j < var; ++j) {
					int e = (i / acc[j]) % (degs[j] + 1);
					if (e <= b.degs[j])
						k += e * b.acc[j];
					else if (coefs[i] != 0)
						return 0;
				}
				for (int j = prev+1; j < k; ++j) {
					if (b.coefs[j] != 0)
						return 0;
				}
				if (coefs[i] != b.coefs[k])
					return 0;
				prev = k;
			}
			for (int i = n; i < b.n; ++i) {
				if (b.coefs[i] != 0)
					return 0;
			}
			return 1;
		}
		/**
		 * Overload operator !=
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline bool operator!= (const DistributedDenseMultivariateModularPolynomial<Field>& b) const  {
			return !(*this == b);
		}

		/**
		 * Overload operator +
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field> operator+ (const DistributedDenseMultivariateModularPolynomial<Field>& b) const {
			if (p != b.p) {
				std::cout << "BPAS: error, trying to add between Z/" << p << "Z and Z/" << b.p << "Z from DDMMP<Field>." << std::endl;
				exit(1);
			}
			if (isConstant()) { return b + coefs[0]; }
			if (b.isConstant()) { return *this + b.coefs[0]; }

			bool isSame = isSameRing(b);
			if (!isSame) {
				std::cout << "BPAS: error, trying to add between Z/" << p << "Z[";
				for (int i = 1; i <= var; ++i) {
					std::cout << names[i];
					if (i < var)
						std::cout << ", ";
				}
				std::cout << "] and Z/" << b.p << "Z[";
				for (int i = 1; i <= b.var; ++i) {
					std::cout << b.names[i];
					if (i < b.var)
						std::cout << ", ";
				}
				std::cout << "] from DDMMP<Field>." << std::endl;
				exit(1);
			}

			int* ds = new int[var];
			for (int i = 0; i < var; ++i)
				ds[i] = (degs[i] >= b.degs[i])? degs[i] : b.degs[i];
			DistributedDenseMultivariateModularPolynomial<Field> res (var, ds, p);
			std::copy(names, names+var+1, res.names);

			//#pragma cilk_grainsize = 1024;
			for (int i = 0; i < res.n; ++i) {
				Field elem = 0;
				int offseta = 0, offsetb = 0;
				for (int j = 0; j < var; ++j) {
					int k = (i / res.acc[j]) % (res.degs[j] + 1);
					if (offseta >= 0 && k <= degs[j])
						offseta += k * acc[j];
					else
						offseta = -1;
					if (offsetb >= 0 && k <= b.degs[j])
						offsetb += k * b.acc[j];
					else
						offsetb = -1;
				}
				if (offseta >= 0 && offsetb >= 0)
					coef_add_mod(&elem, coefs[offseta], b.coefs[offsetb], p);
				else if (offseta >= 0)
					elem = coefs[offseta];
				else if (offsetb >= 0)
					elem = b.coefs[offsetb];
				res.coefs[i] = elem;
			}

			delete [] ds;
			return res;
		}

		/**
		 * Overload operator +=
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field>& operator+= (const DistributedDenseMultivariateModularPolynomial<Field>& b) {
			*this = *this + b;
			return *this;
		}
		
		/**
		 * Overload operator +
		 *
		 * @param e: A constant
		 **/ 
		inline DistributedDenseMultivariateModularPolynomial<Field> operator+ (const Field& e) const {
			DistributedDenseMultivariateModularPolynomial<Field> r (*this);
			return (r += e);
		}
		
		inline friend DistributedDenseMultivariateModularPolynomial<Field> operator+ (const Field& e, const DistributedDenseMultivariateModularPolynomial<Field>& f) {
			return (f + e);
		}

		/**
		 * Overload operator +=
		 *
		 * @param e: A constant
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field>& operator+= (const Field& e) {
			coef_add_mod(&coefs[0], coefs[0], e, p);
			return *this;
		}
		
		/**
		 * Overload operator -
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field> operator- (const DistributedDenseMultivariateModularPolynomial<Field>& b) const {
			if (p != b.p) {
				std::cout << "BPAS: error, trying to subtract between Z/" << p << "Z and Z/" << b.p << "Z from DDMMP<Field>." << std::endl;
				exit(1);
			}
			if (isConstant()) { return -b + coefs[0]; }
			if (b.isConstant()) { return *this - b.coefs[0]; }

			bool isSame = isSameRing(b);
			if (!isSame) {
				std::cout << "BPAS: error, trying to subtract between Z/" << p << "Z[";
				for (int i = 1; i <= var; ++i) {
					std::cout << names[i];
					if (i < var)
						std::cout << ", ";
				}
				std::cout << "] and Z/" << b.p << "Z[";
				for (int i = 1; i <= b.var; ++i) {
					std::cout << b.names[i];
					if (i < b.var)
						std::cout << ", ";
				}
				std::cout << "] from DDMMP<Field>." << std::endl;
				exit(1);
			}

			int* ds = new int[var];
			for (int i = 0; i < var; ++i)
				ds[i] = (degs[i] >= b.degs[i])? degs[i] : b.degs[i];
			DistributedDenseMultivariateModularPolynomial<Field> res (var, ds, p);
			std::copy(names, names+var+1, res.names);

			for (int i = 0; i < res.n; ++i) {
				Field elem = 0;
				int offseta = 0, offsetb = 0;
				for (int j = 0; j < var; ++j) {
					int k = (i / res.acc[j]) % (res.degs[j] + 1);
					if (offseta >= 0 && k <= degs[j])
						offseta += k * acc[j];
					else { offseta = -1; }
					if (offsetb >= 0 && k <= b.degs[j])
						offsetb += k * b.acc[j];
					else { offsetb = -1; }
				}
				if (offseta >= 0 && offsetb >= 0)
					coef_sub_mod(&elem, coefs[offseta], b.coefs[offsetb], p);
				else if (offseta >= 0)
					elem = coefs[offseta];
				else if (offsetb >= 0)
					coef_neg_mod(&elem, b.coefs[offsetb], p);
				res.coefs[i] = elem;
			}
			delete [] ds;
			return res;
		}

		/**
		 * Overload operator -=
		 * 
		 * @param b: A multivariate modular polynomial
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field>& operator-= (const DistributedDenseMultivariateModularPolynomial<Field>& b) {
			*this = *this - b;
			return *this;
		}

		/**
		 * Overload operator -, negate
		 *
		 * @param
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field> operator- () const {
			DistributedDenseMultivariateModularPolynomial<Field> res (var, degs, p);
			std::copy(names, names+var+1, res.names);
			for (int i = 0; i < res.n; ++i)
				coef_neg_mod(&res.coefs[i], coefs[i], p);
			return res;
		}

		/**
		 * Overload operator -
		 *
		 * @param e: A constant
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field> operator- (const Field& e) const {
			DistributedDenseMultivariateModularPolynomial<Field> r (*this);
			return (r -= e);
		}

		inline friend DistributedDenseMultivariateModularPolynomial<Field> operator- (const Field& e, const DistributedDenseMultivariateModularPolynomial<Field>& f) {
			return (-f + e);
		}

		/**
		 * Overload operator -=
		 *
		 * @param e: A constant
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field>& operator-= (const Field& e) {
			coef_sub_mod(&coefs[0], coefs[0], e, p);
			return *this;
		}

		/**
		 * Negate, f(-x)
		 *
		 * @param
		 **/
		inline void negate() {
			for (int i = 0; i < n; ++i)
				coef_neg_mod(&coefs[i], coefs[i], p);
		}

		/**
		 * Overload operator *
		 *
		 * @param b: A multivariate modular polynomial
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field> operator* (const DistributedDenseMultivariateModularPolynomial<Field>& b) const {
			if (p != b.p) {
				std::cout << "BPAS: error, trying to multiply between Z/" << p << "Z and Z/" << b.p << "Z from DDMMP<Field>." << std::endl;
				exit(1);
			}
			if (isConstant()) { return b * coefs[0]; }
			if (b.isConstant()) { return *this * b.coefs[0]; }

			bool isSame = isSameRing(b);
			if (!isSame) {
				std::cout << "BPAS: error, trying to multiply between Z/" << p << "Z[";
				for (int i = 1; i <= var; ++i) {
					std::cout << names[i];
					if (i < var)
						std::cout << ", ";
				}
				std::cout << "] and Z/" << b.p << "Z[";
				for (int i = 1; i <= b.var; ++i) {
					std::cout << b.names[i];
					if (i < b.var)
						std::cout << ", ";
				}
				std::cout << "] from DDMMP<Field>." << std::endl;
				exit(1);
			}

			int* ds = new int[var];
			for (int i = 0; i < var; ++i)
				ds[i] = degs[i] + b.degs[i];
			DistributedDenseMultivariateModularPolynomial<Field> res (var, ds, p);
			std::copy(names, names+var+1, res.names);

			for (int i = 0; i < n; ++i) {
				for (int v = 0; v < var; ++v)
					ds[v] = (i / acc[v]) % (degs[v] + 1);
				for (int j = 0; j < b.n; ++j) {
					int k = 0;
					for (int v = 0; v < b.var; ++v) {
						int e = (j / b.acc[v]) % (b.degs[v] + 1);
						if (v < var)
							k += (ds[v] + e) * res.acc[v];
						else 
							k += e * res.acc[v];
					}
					for (int v = b.var; v < var; ++v)
						k += ds[v] * res.acc[v];
					// res.coefs[k] += coefs[i] * b.coefs[j];
					Field t;
					coef_mul_mod(&t, coefs[i], b.coefs[j], p);
					coef_add_mod(&res.coefs[k], res.coefs[k], t, p);
				}
			}

			delete [] ds;
			return res;
		}

		/**
		 * Overload operator *=
		 * 
		 * @param b: A multivariate modular polynomial
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field>& operator*= (const DistributedDenseMultivariateModularPolynomial<Field>& b) {
			*this = *this * b;
			return *this;
		}

		/**
		 * Overload operator *
		 *
		 * @param e: A constant
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field> operator* (const Field& e) const {
			DistributedDenseMultivariateModularPolynomial<Field> r (*this);
			return (r *= e);
		}

		inline friend DistributedDenseMultivariateModularPolynomial<Field> operator* (const Field& e, const DistributedDenseMultivariateModularPolynomial<Field>& f) {
			return (f * e);
		}

		/**
		 * Overload operator *=
		 *
		 * @param e: A constant
		 **/
		inline DistributedDenseMultivariateModularPolynomial<Field>& operator*= (const Field& f) {
			if (f != 0 && f != 1) {
				Field e(f);
				if (e < 0) { e = e % p + p; }
				for (int i = 0; i < n; ++i)
					coef_mul_mod(&coefs[i], coefs[i], e, p);
			}
			else if (f == 0) { zero(); }
			return *this;
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> operator^ (long long int e) const {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::operator^ NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		inline DistributedDenseMultivariateModularPolynomial<Field>& operator^= (long long int e) {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::operator^= NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> operator/ (const DistributedDenseMultivariateModularPolynomial<Field>& p) const {
			DistributedDenseMultivariateModularPolynomial ret(*this);
			ret /= p;
			return ret;
		}

		inline DistributedDenseMultivariateModularPolynomial<Field>& operator/= (const DistributedDenseMultivariateModularPolynomial<Field>& p) {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::operator/=(DistributedDenseMultivariateModularPolynomial) NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> operator/ (const Field& e) const {
			DistributedDenseMultivariateModularPolynomial ret(*this);
			ret /= e;
			return ret;
		}		

		inline DistributedDenseMultivariateModularPolynomial<Field>& operator/= (const Field& e) {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::operator/= NOT YET IMPLEMENTED" << std::endl;
			return *this;
		}

		/**
		 * Set variable names
		 *
		 * @param xs: Variable names
		 **/ 
		inline void setRingVariables (const std::vector<Symbol>& xs) {
			int ns = xs.size();
			if (ns != var) {
				std::cerr << "BPAS ERROR: DDMMP shrinking and expanding polynomial ring NOT YET IMPLEMENTED" << std::endl;
				return;
			}
			names[0] = "9";
			for (int i = var, j = 0; i > 0 && j < ns; --i, ++j)
				names[i] = xs[j];
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

		inline std::vector<Symbol> variables() const {
			std::cerr << "BPAS ERROR: DDMMP::variables() NOT YET IMPLEMENTED" << std::endl;
			return ringVariables();
		}


		/**
		 * Overload stream operator <<
		 *
		 * @param out: Stream object
		 * @param b: The multivariate modular polynomial
		 **/ 
		inline void print(std::ostream &out) const {
			bool isFirst = 0;
			for (int i = 0; i < this->n; ++i) {
				if (this->coefs[i] != 0) {
					if (isFirst) {
						if (this->coefs[i] >= 0)
							out << "+";
						else if (this->coefs[i] == -1)
							out << "-";
						if (this->coefs[i] != 1 && this->coefs[i] != -1)
							out << this->coefs[i];
						bool isIt = 1;
						for (int j = 0; j < this->var; ++j) {
							int exp = (i / this->acc[j]) % (this->degs[j] + 1);
							if (exp) {
								if ((this->coefs[i] != 1 && this->coefs[i] != -1 && isIt) || !isIt)
									out << "*";
								out << this->names[j+1];
								if (exp > 1)
									out << "^" << exp;
								isIt = 0;
							}
						}
					}
					else { out << this->coefs[i]; }
					isFirst = 1;
				}
			}
			if (!isFirst) { out << "0"; }
		}

		/**
		 * Convert *this to an expression tree
		 */
		inline ExpressionTree convertToExpressionTree() const {
			//TODO
			std::cerr << "DistributedDenseMultivariateModularPolynomial::convertToExpressionTree() NOT YET IMPLEMENTED" << std::endl;
			return ExpressionTree();
		}

		inline void differentiate(const Symbol& s, int k) {
			std::cerr << "DistributedDenseMultivariateModularPolynomial::differentiate NOT YET IMPLEMENTED" << std::endl;
			//TODO
		}

		inline void differentiate(const Symbol& s) {
			differentiate(s, 1);
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> derivative(const Symbol& s, int k) const {
			DistributedDenseMultivariateModularPolynomial ret(*this);
			ret.differentiate(s, k);
			return ret;
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> derivative(const Symbol& s) const {
			return derivative(s, 1);
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> evaluate(const std::vector<Symbol>& syms, const std::vector<Field>& vals) const {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::evaluate NOT YET IMPLEMENTED" << std::endl;
			return *this; 
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> evaluate(int n, const Symbol* syms, const Field* vals) const {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::evaluate NOT YET IMPLEMENTED" << std::endl;
			return *this; 
		}

		inline DistributedDenseMultivariateModularPolynomial<Field> gcd(const DistributedDenseMultivariateModularPolynomial<Field>& p) const {
			//TODO
			std::cerr << "BPAS ERROR: DistributedDenseMultivariateModularPolynomial::gcd NOT YET IMPLEMENTED" << std::endl;
			return *this;
		} 
		
		/**
		 * Compute squarefree factorization of *this
		 */
		inline Factors<DistributedDenseMultivariateModularPolynomial> squareFree() const {
			std::cerr << "DistributedDenseMultivariateModularPolynomial::squareFree NOT YET IMPLEMENTED" << std::endl;
			//TODO
			std::vector<DistributedDenseMultivariateModularPolynomial> ret;
			ret.push_back(*this);
			return ret;
		}







};

#endif
