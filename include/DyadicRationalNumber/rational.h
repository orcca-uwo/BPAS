#ifndef _RATIONAL_H_
#define _RATIONAL_H_

#include <gmpxx.h>

/**
 * A Dyadic rational number.
 */
class DyadicRationalNumber {
	private:
		mpz_class num;			// numinator 
                unsigned long long den; 	// the power of denominator based on 2

	public:
		DyadicRationalNumber() {
			num = 0;
			den = 0;
		}
		DyadicRationalNumber(mpz_class a) {
			num = a;
			den = 0;
		}
		DyadicRationalNumber(mpz_class a, unsigned long long b) {
			num = a;
			den = b;
		}
		inline mpz_class getNum() {
			return num;
		}
		inline unsigned long long getDen() {
			return den;
		}
		inline void setNum(mpz_class a) {
			num = a;
		}
		inline void setDen(unsigned long long a) {
			den = a;
		}
		inline mpq_class getGMP() {
			mpq_class elem;
			elem.get_num() = num;
			elem.get_den() = 1;
			elem.get_den() <<= den;
			return elem;
		}

		/* Overwrite assign operatiors */
		inline DyadicRationalNumber& operator= (DyadicRationalNumber rhs) {
			if (this != &rhs) {
				num = rhs.getNum();
				den = rhs.getDen();
			}
			return *this;
		}
		inline DyadicRationalNumber& operator= (mpz_class a) {
			num = a;
			den = 0;
			return *this;
		}
		inline DyadicRationalNumber& operator+= (DyadicRationalNumber rhs) {
			if (rhs.getNum() != 0) {
				if (num == 0) {
					num = rhs.getNum();
					den = rhs.getDen();
				}
				else if (!den && !rhs.getDen())
					num += rhs.getNum();
				else if (den > rhs.getDen())
					num += rhs.getNum() << (den - rhs.getDen());
				else if (den < rhs.getDen()) {
					num = (num << (rhs.getDen() - den)) + rhs.getNum();
					den = rhs.getDen();
				}
				else {
					mpz_class sum = num + rhs.getNum();
					if (sum == 0) {
						num = 0;
						den = 0;
					}
					else {
						mpz_class com = sum & -sum;
						unsigned long long cbits = mpz_sizeinbase(com.get_mpz_t(), 2) - 1;
						if (cbits >= den) {
							num = sum >> den;
							den = 0;
						}
						else {
							num = sum >> cbits;
							den -= cbits;
						}
					}
				}
			}
			return *this;
		}
		inline DyadicRationalNumber& operator+= (mpz_class a) {
			DyadicRationalNumber rhs(a);
			*this += rhs;
			return *this;
		}
		inline DyadicRationalNumber& operator-= (DyadicRationalNumber rhs) {
			rhs = -rhs;
			*this += rhs;
			return *this;
		}
		inline DyadicRationalNumber& operator-= (mpz_class a) {
			DyadicRationalNumber rhs(a);
			*this -= rhs;
			return *this;
		}
		inline DyadicRationalNumber& operator*= (DyadicRationalNumber rhs) {
			if (num == 0 || rhs.getNum() == 0) {
				num = 0;
				den = 0;
			}
			else if (!den && !rhs.getDen())
				num *= rhs.getNum();
			else if (!den) {
				mpz_class com = num & -num;
				unsigned long long cbits = mpz_sizeinbase(com.get_mpz_t(), 2) - 1;
				if (cbits < rhs.getDen()) {
					num = (num >> cbits) * rhs.getNum();
					den = rhs.getDen() - cbits;
				}
				else {
					num = (num >> rhs.getDen()) * rhs.getNum();
					den = 0;
				}
			}
			else if (!rhs.getDen()) {
				mpz_class com = rhs.getNum() & -rhs.getNum();
				unsigned long long cbits = mpz_sizeinbase(com.get_mpz_t(), 2) - 1;
				if (cbits < den) {
					num *= rhs.getNum() >> cbits;
					den -= cbits;
				}
				else {
					num *= rhs.getNum() >> den;
					den = 0;
				}
			}
			else {
				num *= rhs.getNum();
				den += rhs.getDen();
			}
			return *this;
		}
		inline DyadicRationalNumber& operator*= (mpz_class a) {
			DyadicRationalNumber rhs(a);
			*this *= rhs;
			return *this;
		}
		/* Overwrite arithmetic operators */
		inline DyadicRationalNumber operator+ (DyadicRationalNumber rhs) {
			DyadicRationalNumber res = *this;
			res += rhs;
			return res;
		}
		inline DyadicRationalNumber operator+ (mpz_class a) {
			DyadicRationalNumber rhs(a);
			DyadicRationalNumber res = *this;
			res += rhs;
			return res;
		}
		inline DyadicRationalNumber operator++ (int a) {
			if (!den)
				num++;
			else {
				mpz_class one = 1;
				one <<= den;
				num += one;
			}
			return *this;
		}
		inline DyadicRationalNumber operator- () const {
			DyadicRationalNumber res;
			res.setNum(-num);
			res.setDen(den);
			return res;
		}
		inline DyadicRationalNumber operator- (DyadicRationalNumber rhs) {
			DyadicRationalNumber res = *this;
			res -= rhs;
			return res;
		}
		inline DyadicRationalNumber operator- (mpz_class a) {
			DyadicRationalNumber rhs(a);
			DyadicRationalNumber res = *this;
			res -= rhs;
			return res;
		}
		inline DyadicRationalNumber operator* (DyadicRationalNumber rhs) {
			DyadicRationalNumber res = *this;
			res *= rhs;
			return res;
		}
		inline DyadicRationalNumber operator* (mpz_class a) {
			DyadicRationalNumber rhs(a);
			DyadicRationalNumber res = *this;
			res *= rhs;
			return res;
		}
		inline mpq_class operator/ (DyadicRationalNumber rhs) {
			mpq_class res = num;
			res /= rhs.getNum();
			if (den < rhs.getDen())
				res <<= (rhs.getDen() - den);
			else if (den > rhs.getDen())
				res >>= (den - rhs.getDen());
			return res;
		}
		/* Overwrite bit shift operators */
		inline DyadicRationalNumber& operator<<= (unsigned long long k) {
			if (num != 0) {
				if (!den)
					num <<= k;
				else if (den < k) {
					num <<= k - den;
					den = 0;
				}
				else
					den -= k;
			}
			return *this;
		}
		inline DyadicRationalNumber& operator>>= (unsigned long long k) {
			if (num != 0) {
				if (!den) {
					mpz_class com = num & -num;
					unsigned long long cbits = mpz_sizeinbase(com.get_mpz_t(), 2) - 1;
					if (cbits < k) {
						num >>= cbits;
						den = k - cbits;
					}
					else
						num >>= k;
				}
				else
					den += k;
			}
			return *this;
		}
		inline DyadicRationalNumber operator<< (unsigned long long k) {
			DyadicRationalNumber res = *this;
			res <<= k;
			return res;
		}
		inline DyadicRationalNumber operator>> (unsigned long long k) {
			DyadicRationalNumber res = *this;
			res >>= k;
			return res;
		}
		/* Overwrite comparison operators */
		inline bool operator== (DyadicRationalNumber rhs) {
			if (rhs.getNum() == num && rhs.getDen() == den)
				return 1;
			else
				return 0;
		}
		inline bool operator== (mpz_class a) {
			if (num == a && !den)
				return 1;
			else
				return 0;
		}
		inline bool operator!= (DyadicRationalNumber rhs) {
			if (rhs.getNum() != num || rhs.getDen() != den)
				return 1;
			else
				return 0;
		}
		inline bool operator!= (mpz_class a) {
			if (num != a || den)
				return 1;
			else
				return 0;
		}
		inline bool operator< (DyadicRationalNumber rhs) {
			if (!den && !rhs.getDen())
				return (num < rhs.getNum())? 1 : 0;
			else {
				unsigned long long min = (den <= rhs.getDen())? den : rhs.getDen();
				mpz_class a = num << (rhs.getDen() - min);
				mpz_class b = rhs.getNum() << (den - min);
				return (a < b)? 1 : 0;
			}
		}
		inline bool operator< (mpz_class a) {
			if (a == 0) {
				if (num < 0)
					return 1;
				else
					return 0;
			}
			if (!den)
				return (num < a)? 1 : 0;
			else {
				a <<= den;
				return (num < a)? 1 : 0;
			}
		}
		inline bool operator<= (DyadicRationalNumber rhs) {
			if (*this == rhs || *this < rhs)
				return 1;
			else
				return 0;
		}
		inline bool operator<= (mpz_class a) {
			if (a == 0) {
				if (num <= 0)
					return 1;
				else
					return 0;
			}
			if (*this == a || *this < a)
				return 1;
			else
				return 0;
		}
		inline bool operator> (DyadicRationalNumber rhs) {
			if (rhs < *this)
				return 1;
			else
				return 0;
		}
		inline bool operator> (mpz_class a) {
			if (a == 0) {
				if (num > 0)
					return 1;
				else
					return 0;
			}
			DyadicRationalNumber rhs(a);
			if (rhs < *this)
				return 1;
			else
				return 0;
		}
		inline bool operator>= (DyadicRationalNumber rhs) {
			if (rhs == *this || rhs < *this)
				return 1;
			else
				return 0;
		}
		inline bool operator>= (mpz_class a) {
			if (a == 0) {
				if (num >= 0)
					return 1;
				else
					return 0;
			}
			DyadicRationalNumber rhs(a);
			if (rhs == *this || rhs < *this)
				return 1;
			else
				return 0;
		}
};

/**
 * Compute the sign of mpq_class a
 *
 * If a > 0, return 1
 * Elif a < 0, return -1
 * Else, return 0
 **/
inline int sign(mpq_class a) {
	if (a > 0)
		return 1;
	else if (a < 0)
		return -1;
	else
		return 0;
}

/**
 * Compute the sign of DyadicRationalNumber a
 *
 * If a > 0, return 1
 * Elif a < 0, return -1
 * Else, return 0
 **/
inline int sign(DyadicRationalNumber a) {
	if (a.getNum() > 0)
		return 1;
	else if (a.getNum() < 0)
		return -1;
	else
		return 0;
}

/**
 * Overwrite abs function
 * Return the absulote value of a rational number
 **/
inline DyadicRationalNumber abs (DyadicRationalNumber a) {
	DyadicRationalNumber res;
	mpz_class num = a.getNum();
	if (num >= 0)
		res.setNum(num);
	else
		res.setNum(-num);
	res.setDen(a.getDen());
	return res;
}

/**
 * Overwrite std::ostream& operator<<
 * Use GMP mpq_class operator<<
 **/
std::ostream& operator<<(std::ostream& os, DyadicRationalNumber obj);

#endif
