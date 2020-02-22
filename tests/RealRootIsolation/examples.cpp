#include "examples.h"

void assignRCs(std::vector<RationalRegularChain>& rcs, int choice, std::string filepath, Symbol* xs) {
	std::ostringstream convert;
	convert << choice;
	std::string filename = filepath + "multivariate_" + convert.str() + ".dat";
        std::ifstream inputfile(filename.c_str(), std::ifstream::in);

	int m;
	inputfile >> m;
        int vars;
        inputfile >> vars;

	for (int s = 0; s < m; ++s) {
		for (int i = 0; i < vars; ++i) {
			int n;
			inputfile >> n;
			int v;
			inputfile >> v;
			if (v == 1) {
				mpz_class elem;
				inputfile >> elem;
				int* curd = new int[vars];
				for (int j = 0; j < vars; ++j)
					inputfile >> curd[j];
				int d = (curd[vars-1] > n-1)? curd[vars-1] : n-1;
				DUQP p(d+1);
				p.setVariableName(xs[vars-1]);
				p.setCoefficient(d, elem);
				while (n > 1) {
					--d, --n;
					inputfile >> elem;
					for (int j = 0; j < vars; ++j) {
						inputfile >> curd[j];
					}
					p.setCoefficient(curd[vars-1], elem);
					if (curd[vars-1] < d) {
						for (int j = d; j > curd[vars-1]; --j ) {
							p.setCoefficient(j, 0);
						}
					}
					d = curd[vars-1];
				}
				delete [] curd;
				rcs[s] += p;
			}
			else if (v > 1) {
				std::vector<Symbol> ns;
				for (int i = 0; i < v; ++i)
					ns.push_back(xs[vars-v+i]);
				SMQP mp(v);
				mp.setRingVariables(ns);
				for (int j = 0; j < n; ++j) {
					RationalNumberTerm t;;
					inputfile >> t.coef;
					t.v = v;
					t.degs = new int[v];
					int* curd = new int[vars];
					for (int k = 0; k < vars; ++k)
						inputfile >> curd[k];
					for (int k = 0; k < v; ++k)
						t.degs[v-k-1] = curd[vars-1-k];
					delete [] curd;
					mp.setCoefficient(v, t.degs, t.coef);
				}
				rcs[s] += mp;
			}
		}
	}
        inputfile.close();
}

void assignPoly(DUQP* p, int n, int choice, std::string filepath) {
        if (choice) {
                if (choice == 1) {         // Bn,d(x) = 2^degree*x^n + ... + 2^degree
			mpz_class elem = 1;
			elem <<= n-1;
                        for (int i = 0; i < n; i++)
                                p->setCoefficient(i, elem);
                }
                else if (choice == 2) {         // Cn,d(x) = x^n + degree
                        for (int i = 0; i < n; i++) {
                                if (!i)
                                        p->setCoefficient(i, n-1);
                                else if (i == n-1)
                                        p->setCoefficient(i, 1);
                                else
                                        p->setCoefficient(i, 0);
                        }
                }
		else if (choice == 3) {         // Mignotte polynomial
                        for (int i = 0; i < n; i++) {
                                if (!i)
                                        p->setCoefficient(i, -2);
                                else if (i == 1)
                                        p->setCoefficient(i, 20);
                                else if (i == 2)
                                        p->setCoefficient(i, -50);
                                else if (i == n-1)
                                        p->setCoefficient(i, 1);
                                else
                                        p->setCoefficient(i, 0);
                        }
                }
                else if (choice == 4) {         // Chebycheff polynomial
			std::ostringstream convert;
			if (n == 128 || n == 256 || n == 512 || n == 1024)
				convert << n-1;
			else {
				std::cout << "ERROR: no input file!" << std::endl;
				return;
			}
			std::string filename = filepath + "chebycheff_" + convert.str() + ".dat";
			
                        std::ifstream inputfile(filename.c_str(), std::ifstream::in);
                        mpz_class elem;
                        inputfile >> elem;
                        for (int i = 0; i < n; i++) {
                                inputfile >> elem;
                                p->setCoefficient(i, elem);
                        }
                        inputfile.close();
                }
                else if (choice == 5) {         // Laguerre polynomial
			std::ostringstream convert;
			if (n == 128 || n == 256 || n == 512 || n == 1024)
				convert << n-1;
			else {
				std::cout << "ERROR: no input file!" << std::endl;
				return;
			}
			std::string filename = filepath + "laguerre_" + convert.str() + ".dat";

                        std::ifstream inputfile(filename.c_str(), std::ifstream::in);
                        mpz_class elem;
                        inputfile >> elem;
                        for (int i = 0; i < n; i++) {
                                inputfile >> elem;
                                p->setCoefficient(i, elem);
                        }
                        inputfile.close();
                }
                else if (choice == 6) {         // Wilkinson polynomial
			std::ostringstream convert;
			if (n == 128 || n == 256 || n == 512 || n == 1024)
				convert << n-1;
			else {
				std::cout << "ERROR: no input file!" << std::endl;
				return;
			}
			std::string filename = filepath + "wilkinson_" + convert.str() + ".dat";

                        std::ifstream inputfile(filename.c_str(), std::ifstream::in);
                        mpz_class elem;
                        inputfile >> elem;
                        for (int i = 0; i < n; i++) {
                                inputfile >> elem;
                                p->setCoefficient(i, elem);
                        }
                        inputfile.close();
                }
		else if (choice == 7) {		// Random polynomial with expected number of roots
			std::ostringstream convert;
			if (n == 128 || n == 256 || n == 512 || n == 1024)
				convert << n-1;
			else {
				std::cout << "ERROR: no input file!" << std::endl;
				return;
			}
			std::string filename = filepath + "random2_" + convert.str() + ".dat";

			std::ifstream inputfile(filename.c_str(), std::ifstream::in);
			mpz_class elem;
			inputfile >> elem;
			for (int i = 0; i < n; i++) {
				inputfile >> elem;
				p->setCoefficient(i, elem);
			}
			inputfile.close();
		}
		else if (choice == 8) {         // Random GMP polynomial
			std::ostringstream convert;
			if (n == 128 || n == 256 || n == 512 || n == 1024 || n == 2048 || n == 4096 || n == 8192 || n == 16384 || n == 32768 || n == 65536)
				convert << n-1;
			else {
				std::cout << "ERROR: no input file!" << std::endl;
				return;
			}
			std::string filename = filepath + "randomGMP_" + convert.str() + ".dat";

                        std::ifstream inputfile(filename.c_str(), std::ifstream::in);
                        mpz_class elem;
                        inputfile >> elem;
                        for (int i = 0; i < n; i++) {
                                inputfile >> elem;
                                p->setCoefficient(i, elem);
                        }
                        inputfile.close();
                }
                else {			// Generate random GMP
                        gmp_randclass rr (gmp_randinit_default);
                        rr.seed(time(NULL));
                        srand (time(NULL));
                        for(int i = 0; i < n; i++) {
                                mpz_class elem = rr.get_z_bits(n-1);
                                int k = (int) rand() % 2;
                                if (k)
                                        p->setCoefficient(i, elem);
                                else {
					mpz_class tmp = -elem;
                                        p->setCoefficient(i, tmp);
				}
                        }
                }
        }
        else {
		std::string filename = filepath + "poly_input.dat";
                std::ifstream inputfile(filename.c_str(), std::ofstream::in);
                Integer n = p->degree()+1;
                mpz_class elem;
                inputfile >> elem;
                for (int i = 0; i < n; i++) {
                        inputfile >> elem;
                        p->setCoefficient(i, elem);
                }
                inputfile.close();
        }
}
