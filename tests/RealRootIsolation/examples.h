#include <bpas.h>

/**
 * Assign regular chains
 *
 * Ouput:
 * @rcs: Assigned regular chains
 *
 * All regular chains written in a file
 * File name begins with "multivariate_" and ends with a @choice number
 *
 * Input file format:
 * Line 0: <Number of RC> <Number of variables, v>
 * For each polynomial (univariate or multivariate):
 * Line 1: <Number of terms> <Number of variables, x <= v>
 * Line i: <Coefficient in a descending order> <Each degree in v, in a descending order>
 *
 * @choice:
 * 1-5: 4-body-homog, Arnborg-Lazard, Arnborg-Lazard-rev, Barry, Caprasse
 * 6-10: Caprasse-Li, chemical-reaction, circles, cyclic-5, Czapor-Geddes-Wang
 * 11-15: d2v10, d2v15, d4v5, d4v10, easy1
 * 16-20: easy2 - easy6
 * 21-25: easy7 - easy9, fabfaux, geometric-constraints
 * 26-30: GonzalezGonzalez, Katsura-4, l-3, lllp1, lllp2
 * 31-35: lllp3 - lllp6, neural-network
 * 36-40: PlateForme2d-easy, r-5, r-6, Rose, Takeuchi-Lu
 * 41-45: Trinks-difficult, trivial-5, wilkinsonxy, 5-body-homog, simple-nql-20-30
 * 46-50: nld-9-3, nld-10-3, nql-10-4, nql-15-2, 13_sings_9
 * 51-53: challenge_12, SA_2_4_eps, ten_circles
 * @filepath: Directory of the file stored
 **/
void assignRCs(std::vector<RationalRegularChain>& rcs, int choice, std::string filepath, Symbol*);


/**
 * Assign univariate polynomial
 *
 * Output:
 * @p: Assigned univariate rational polynomial given degree
 * 
 * Examples 1, 2, 3 are genereated, 0, 4, 5, 6, 7, 8, 9 are read from a file
 * any other number are to generate polynomial with randon GMP coefficient
 * File name for custom example: poly_input.dat
 *
 * Input file format:
 * Line 0: <Number of terms>
 * Line i: Coefficient in a ascending order
 *
 * @choice:
 * 0: Read from a file: poly_input.dat
 * 1: Bn,d(x) = 2^degree*x^n + ... + 2^degree
 * 2: Cn,d(x) = x^n + degree
 * 3: Mignotte x^n-50x^2+20x-2
 * 4: Chebycheff T_{n} = 2*x*T_{n-1}-T_{n-2}, T_1 = x, T_0 = 1
 * 5: Laguerre L_n = n! * sum_{i=0} binomial(n, i) * (-1)^i / i! * x^i
 * 6: Wilkinson W_n = mul_{i=0} (x - i - 1)
 * 7: Random polynomials with expected number of roots
 * 8: Polynomial coefficients: Random GMP integer, where bits based on degree
 * @filepath: Directory of the file stored
 **/
void assignPoly(DUQP* p, int n, int choice, std::string filepath);
