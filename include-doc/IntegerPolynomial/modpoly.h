#ifndef _MODPOLY_H_
#define _MODPOLY_H_

#include "Multiplication/Multiplication.h"

sfixn nextprime(int* s, mpz_class m);
void monicGCD(sfixn* gp, int* gd, sfixn* fp, int fd, sfixn p, sfixn lc);
//sfixn nextprime(int*, mpz_class);
//void monicGCD(sfixn*, int*, sfixn*, int, sfixn, sfixn);

class UniModularSubresultantChain {
	public:
		int n;		// Number of polynomials, at least 2
		int* deg;	// Degree of each polynomial
		int* size;	// Size of each polynomial array
		sfixn** coef;	// Coefficients

		UniModularSubresultantChain () : n (0) {}		
		UniModularSubresultantChain (int _n) : n (_n) {
			if (n > 0) {
				deg = new int[n];
				size = new int[n];
				for (int i = 0; i < n; ++i)
					size[i] = deg[i] = 0;
				coef = new sfixn*[n];
			}
			else { n = 0; }
		}
		UniModularSubresultantChain (const UniModularSubresultantChain& s) : n (s.n) {
			deg = new int[n];
			size = new int[n];
			std::copy(s.deg, s.deg+n, deg);
			std::copy(s.size, s.size+n, size);
			coef = new sfixn*[n];
			for (int i = 0; i < n; ++i) {
				if (size[i] > 0) {
					coef[i] = new sfixn[size[i]];
					std::copy(s.coef[i], s.coef[i]+size[i], coef[i]);
				}
			}
		}
		~UniModularSubresultantChain() {
			if (n) {
				for (int i = 0; i < n; ++i)
					if (size[i]) { delete [] coef[i]; }
				delete [] coef;
				delete [] deg;
				delete [] size;
			}
		}

		inline UniModularSubresultantChain& operator= (UniModularSubresultantChain s) {
			if (this != &s) {
				if (n) {
					for (int i = 0; i < n; ++i)
						if (size[i]) { delete [] coef[i]; }
					delete [] coef;
					delete [] deg;
					delete [] size;
				}
				n = s.n;
				deg = new int[n];
				size = new int[n];
				std::copy(s.deg, s.deg+n, deg);
				std::copy(s.size, s.size+n, size);
				if (n > 1) { coef = new sfixn*[n]; }
				for (int i = 0; i < n; ++i) {
					if (size[i] > 0) {
						coef[i] = new sfixn[size[i]];
						std::copy(s.coef[i], s.coef[i]+size[i], coef[i]);
					}
				}
			}
			return *this;
		}
		inline friend std::ostream& operator<< (std::ostream &out, UniModularSubresultantChain s) {
			if (!s.n) { out << "0"; }
			else {
				bool isFirst = 0, isZero = 1;
				for (int i = 0; i < s.n; ++i) {
					if (s.size[i]) {
						isZero = 1;
						if (isFirst) { out << "\n"; }
						for (int j = 0; j <= s.deg[i]; ++j) {
							if (s.coef[i][j]) {
								out << s.coef[i][j];
								if (j == 1) { out << "*x"; }
								else if (j) { out << "*x^" << j; }
								if (j < s.deg[i]) { out << "+"; }
								isZero = 0;
							}
						}
						if (isZero) { out << "0"; }
						isFirst = 1;
					}
				}
			}
			return out;
		}
};

template <class Ring>
class BiModularSubresultantChain {
	public:
		int n;		// Number of polynomials, at least 2
		int** deg;	// Partial degree of each polynomial
		int* size;	// Size of each polynomial array
		Ring** coef;	// Coefficients of a polynomial stored in a distributed way

		BiModularSubresultantChain<Ring> () : n (0) {}
		BiModularSubresultantChain<Ring> (int _n) : n (_n) {
			if (n > 0) {
				if (n > 1) {
					deg = new int*[n];
					coef = new Ring*[n];
				}
				size = new int[n];
				for (int i = 0; i < n; ++i) {
					deg[i] = new int[2];
					deg[i][1] = deg[i][0] = 0;
					size[i] = 0; 
				}
			}
			else { n = 0; }
		}
		BiModularSubresultantChain<Ring> (const BiModularSubresultantChain<Ring>& s) : n (s.n) {
			size = new int[n];
			std::copy(s.size, s.size+n, size);
			deg = new int*[n];
			coef = new Ring*[n];
			for (int i = 0; i < n; ++i) {
				if (size[i] > 0) {
					deg[i] = new int[2];
					std::copy(s.deg[i], s.deg[i]+2, deg[i]);
					coef[i] = new Ring[size[i]];
					std::copy(s.coef[i], s.coef[i]+size[i], coef[i]);
				}
			}
		}
		~BiModularSubresultantChain<Ring> () {
			if (n) {
				for (int i = 0; i < n; ++i) {
					if (size[i]) {
						delete [] deg[i];
						delete [] coef[i];
					}
				}
				delete [] coef;
				delete [] deg;
				delete [] size;
			}
		}

		inline BiModularSubresultantChain<Ring>& operator= (BiModularSubresultantChain<Ring> s) {
			if (this != &s) {
				if (n) {
					for (int i = 0; i < n; ++i) {
						if (size[i]) {
							delete [] deg[i];
							delete [] coef[i];
						}
					}
					if (n > 1) {
						delete [] coef;
						delete [] deg;
					}
					delete [] size;
				}
				size = new int[n];
				std::copy(s.size, s.size+n, size);
				if (n > 1) {
					deg = new int*[n];
					coef = new Ring*[n];
				}
				for (int i = 0; i < n; ++i) {
					if (size[i] > 0) {
						deg[i] = new int[2];
						std::copy(s.deg[i], s.deg[i]+2, deg[i]);
						coef[i] = new Ring[size[i]];
						std::copy(s.coef[i], s.coef[i]+size[i], coef[i]);
					}
				}
			}
			return *this;
		}
		inline bool operator== (BiModularSubresultantChain<Ring>& s ) {
			for (int i = 0; i < n && i < s.n; ++i) {
				if (deg[i][1] != s.deg[i][1] || deg[i][0] != s.deg[i][0])
					return 0;
				for (int j = 0; j < size[i]; ++j) {
					if (coef[i][j] != s.coef[i][j])
						return 0;
				}
			}
			for (int i = n; i < s.n; ++i) {
				for (int j = 0; j < s.size[i]; ++j)
					if (s.coef[i][j] != 0) { return 0; }
			}
			for (int i = s.n; i < n; ++i) {
				for (int j = 0; j < size[i]; ++j)
					if (coef[i][j] != 0) { return 0; }
			}
			return 1;
		}
		inline friend std::ostream& operator<< (std::ostream &out, BiModularSubresultantChain<Ring> s) {
			bool isZero = 1;
			for (int i = 0; i < s.n; ++i) {
				if (s.size[i] && (s.deg[i][1] || s.deg[i][0] || s.coef[i][0] != 0)) {
					if (!isZero) { out << "\n"; }
					out << i << ": ";
					int k = 1 + s.deg[i][0];
					for (int j = 0; j <= s.deg[i][1]; ++j) {
						out << "(";
						for (int l = 0; l <= s.deg[i][0]; ++l) {
							out << s.coef[i][j*k+l];
							if (l == 1) { out << "*t"; }
							else if (l) { out << "*t^"<< l; }
							if (l < s.deg[i][0] && s.coef[i][j*k+l] >= 0) { out << "+"; }
						}
						out << ")";
						if (j == 1) { out << "*x"; }
						else if (j) { out << "*x^" << j; }
						if (j < s.deg[i][1]) { out << "+"; }
					}
					isZero = 0;
				}
			}
			if (isZero) { out << "0"; }
			return out;
		}
};

//void prem (sfixn*, int*, sfixn*, int, sfixn);
//void brownSubresultant(UniModularSubresultantChain*, sfixn*, int, sfixn*, int, sfixn);
//void lagrangeBasis (sfixn*, int, sfixn);
//sfixn longDivision(sfixn*, sfixn, sfixn*, int, sfixn);
//void interpolateSubresultant (BiModularSubresultantChain<sfixn>*, UniModularSubresultantChain*, int, sfixn, int);
BiModularSubresultantChain<sfixn> modularSubresultantChain(mpz_class*, int*, mpz_class*, int*, sfixn);
BiModularSubresultantChain<mpz_class> modularSetSubresultantChain(mpz_class*, int*, mpz_class*, int*, sfixn*, int, mpz_class*);
BiModularSubresultantChain<mpz_class> modularStableSubresultantChain (mpz_class*, int*, mpz_class*, int*, int v=2);

#endif
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


