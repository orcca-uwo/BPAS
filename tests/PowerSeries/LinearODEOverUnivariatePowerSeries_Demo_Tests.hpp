#ifndef _LODEoUPS_DEMO_TESTS_HPP_
#define _LODEoUPS_DEMO_TESTS_HPP_

#include <sstream>
#include <string>

// LODEoUPS
#include "PowerSeries/LinearODEOverUnivariatePowerSeries.h"
#include "PowerSeries/LinearODEOverUnivariatePowerSeries_parallel.hpp"

/**
 * Runs timings on the polynomial family of test ODE's of order n
 *  The polynomial ode family is defined as this equation 
 * (written below in latex-style pseudocode as well as Maple ASCII  where n corresponds to the ode_order parameter). 
 * y^{(n)}(x) + \sum_{i=0}^{n-1} x^{n-i} y^{(i)}(x) = 0
 * 
 *                       /n - 1                    \
            / n      \   |----- / i      \         |
            |d       |   | \    |d       |  (n - i)|         n
            |--- y(x)| + |  )   |--- y(x)| x       | + y(x) x  = 0
            |  n     |   | /    |  i     |         |
            \dx      /   |----- \dx      /         |
                         \i = 1                    /

 * 
 * @param ode_order the order of the modified Euler ODE
 * @param precision the precision of the ODE solution power series
 */
void testFamily_poly_demo(int ode_order,  int precision);

/**
 * Runs timings on the exponential family of test ODE's of order n
 *  The exponential ode family is defined as this equation 
 * (written below in latex-style pseudocode as well as Maple ASCII  where n corresponds to the ode_order parameter). 
 * \sum_{i=0}^{n} (i+1) e^{(i+1)x} y^{(i)}(x)  = e^{nx}
 * 
 *    /  n                                    \
      |----- / i      \                       |
      | \    |d       |                       |
      |  )   |--- y(x)| (i + 1) exp((i + 1) x)| + y(x) exp(x) = exp(n x)
      | /    |  i     |                       |
      |----- \dx      /                       |
      \i = 1                                  /

 * 
 * @param ode_order the order of the modified Euler ODE
 * @param precision the precision of the ODE solution power series
 */
void testFamily_exp_demo(int ode_order,  int precision);

/**
 * Runs timings on the rational family of test ODE's of order n
 * The rational ode family is defined as this equation 
 * (written below in latex-style pseudocode as well as Maple ASCII where n corresponds to the ode_order parameter). 
 * \sum_{i=0}^{n} \frac{1}{1+\sum_{j=1}^{i+2} jx^j} y^{(i)}(x)  = \frac{1}{1+x}
 * 
*   /                          i                          \
    |                         d                           |
    |  n                      --- y(x)                    |
    |-----                      i                         |
    | \                       dx                          |       y(x)         1
    |  )   -----------------------------------------------| + ------------ = -----
    | /         (i + 3)                                   |      2           1 + x
    |-----     x        ((i + 3) x - i - 3 - x)      x    |   2 x  + x + 1
    |i = 1 1 + -------------------------------- + --------|
    |                             2                      2|
    \                      (x - 1)                (x - 1) /

 * 
 * @param ode_order the order of the modified Euler ODE
 *  @param precision the precision of the ODE solution power series
 */
void testFamily_ratio_demo(int ode_order,  int precision);




/**
 * Test header for printing
 */
void test_header (string test_string){
    std::cout << "\n===========================================================================================" << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
    std::cout << test_string << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "===========================================================================================\n" << std::endl;


}

#endif
