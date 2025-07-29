#include "LinearODEOverUnivariatePowerSeries_Demo_Tests.hpp"

int main(){

    /**
     * To Use this Demo file easily, simply change the testprecision parameter 
     * to select the number of power series terms to be calculated.
     * 
     * Comment terminators have partially been inserted below, so it is easy to comment 
     * out blocks which run different test families for different ODE orders
     */
    

    int testprecision = 256;

    /*
    testFamily_poly_demo(2, testprecision);
    testFamily_poly_demo(8, testprecision);
    testFamily_poly_demo(32, testprecision);
    //*/


    //*
    testFamily_exp_demo(2,  testprecision);
    testFamily_exp_demo(8,  testprecision);
    testFamily_exp_demo(32,  testprecision);
    //*/


    /*
    testFamily_ratio_demo(2, testprecision);
    testFamily_ratio_demo(8, testprecision);
    testFamily_ratio_demo(32, testprecision);
    //*/


    return 0;
}



/**
 * Runs timings on the polynomial family of test ODE's of order n
 * @param ode_order the order of the modified Euler ODE
 * @param precision the precision of the ODE solution power series
 */
void testFamily_poly_demo(int ode_order,  int precision){

    // For printing header
    std::ostringstream oss;
    oss << "Test Family, Polynomial, order " << ode_order << ", up to precision " << precision;
    std::string header_string = oss.str();
    test_header(header_string);

    // Helper mpq_t's
    mpq_t onempq; mpq_init(onempq); mpq_set_si(onempq,1,1);
    mpq_t zerompq; mpq_init(zerompq); mpq_set_si(zerompq,1,1);


    // The coefficients as polynomials, which are monomials
    Poly_ptr monomials[ode_order];
    for(int i=0; i<ode_order; i++){
        monomials[i] = make_univariate_monomial(i+1, onempq);
    }

    // The coefficient of the highest order term, which is 1
    PowerSeries_t * one_PSS = onePowerSeries_PS(1);

    // The ode coefficients as power series
    PowerSeries_t * ode_coefs[3][ode_order+1];
    for(int j=0; j<3; j++){
        for(int i=0; i<ode_order; i++){
            ode_coefs[j][i] = convertPolyToPowerSeries_PS(monomials[ode_order-i-1]);
        }
        ode_coefs[j][ode_order] = onePowerSeries_PS(1);
    }

    // The rhs of the equation is 0
    PowerSeries_t * zero_PSS = zeroPowerSeries_PS();

    PowerSeries_t * ode_rhs[3];
    for(int i=0; i<3; i++){
        ode_rhs[i] = zeroPowerSeries_PS();
    }

    // The initial conditions will all be (D^n)(y)(0) = 1
    mpq_t inits[ode_order];
    for (int i=0; i<ode_order; i++){
        mpq_init(inits[i]);
        mpq_set(inits[i], onempq);
    }

    /* ---------------------------------------------------------------------------------------------
    ** Using the power series solution generator, serially **
    --------------------------------------------------------------------------------------------- */

    // Construct a LODE structure which will use the serial generator
    LODEoUPS_t * lode_serial_generator = construct_LODEoUPS(ode_order, ode_coefs[0], ode_rhs[0], inits);

    // Calculate the power series solution up to a degree using the generator
    auto start = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_using_sol_generator(precision, lode_serial_generator);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Using solution generator, serial, time for calculating " << precision << " coefficients: " << duration_ms.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the solution generator serially
    mpq_t print_mpq_helper; mpq_init(print_mpq_helper); 
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_serial_generator->sol, precision, &zerompq));
    std::cout << "For using solution generator, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n\n" << std::endl;


    /* ---------------------------------------------------------------------------------------------
    ** Using the truncated recursive method, serially **
    --------------------------------------------------------------------------------------------- */

    
    // Construct a LODE structure which will use the truncated method
    LODEoUPS_t * lode_truncated_serial = construct_LODEoUPS(ode_order, ode_coefs[1], ode_rhs[1], inits);

    // Calculate the power series solution up to a degree using truncated serial construction
    start = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_truncated_construction(precision, lode_truncated_serial);
    end = std::chrono::high_resolution_clock::now();
    duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Using truncated construction, serial, time for calculating " << precision << " coefficients: " << duration_ms.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the truncated serial construction
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_truncated_serial->sol, precision, &zerompq));
    std::cout << "For using serial truncated construction, serial, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n\n" << std::endl;
    

    /* ---------------------------------------------------------------------------------------------
    ** Using the parallel method, truncated **
    --------------------------------------------------------------------------------------------- */
    
    // Construct a LODE structure which will use the parallel method
    LODEoUPS_t * lode_parallel = construct_LODEoUPS(ode_order, ode_coefs[2], ode_rhs[2], inits);

    // Calculate the power series solution up to a degree using truncated serial construction
    start = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_truncated_parallel(precision, lode_parallel);
    end = std::chrono::high_resolution_clock::now();
    auto duration_ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Using parallel construction, time for calculating " << precision << " coefficients: " << duration_ms2.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the truncated serial construction
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_parallel->sol, precision, &zerompq));
    std::cout << "For using truncated parallel, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n" << std::endl;

    // Parallel speedup
    std::cout << "Parallel speedup: " << (double)duration_ms.count()/duration_ms2.count() << "\n" << std::endl;
    


    // Free memory
    mpq_clear(print_mpq_helper);
    mpq_clear(onempq);
    mpq_clear(zerompq);
    destroyPowerSeries_PS(one_PSS);
    destroyPowerSeries_PS(zero_PSS);
    for(int i=0; i<ode_order; i++){
        freePolynomial_AA(monomials[i]);
    }
    for(int j=0; j<3; j++){
        for(int i=0; i<ode_order; i++){
            destroyPowerSeries_PS(ode_coefs[j][i]);
        }
    }
    for(int i=0; i<3; i++){
        destroyPowerSeries_PS(ode_rhs[i]);
    }
    for(int i=0; i<ode_order; i++){
        mpq_clear(inits[i]);
    }
    destroy_LODEoUPS(lode_serial_generator);
    destroy_LODEoUPS(lode_truncated_serial);
    destroy_LODEoUPS(lode_parallel);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/**
 * Runs timings on the rational family of test ODE's of order n
 * @param ode_order the order of the modified Euler ODE
 *  @param precision the precision of the ODE solution power series
 */
void testFamily_ratio_demo(int ode_order,  int precision){

    // For printing header
    std::ostringstream oss;
    oss << "Modified Euler, Rational, order " << ode_order << ", up to precision " << precision;
    std::string header_string = oss.str();
    test_header(header_string);

    // Helper mpq_t's
    mpq_t onempq; mpq_init(onempq); mpq_set_si(onempq,1,1);
    mpq_t zerompq; mpq_init(zerompq); mpq_set_si(zerompq,1,1);

    // A Power Series for the number 1;
    PowerSeries_t * one_PSS = onePowerSeries_PS(1);


    // The monomials which will be added to 1 in the denominator of the rational functions
    Poly_ptr monomials[ode_order];
    for(int i=0; i<ode_order; i++){
        monomials[i] = make_univariate_monomial(i+1, onempq);
    }

    // A univariate polynomial which is just the number 1
    Poly_ptr oneAA = make_univariate_monomial(0,onempq);

    // Construct the ode rhs
    Poly_ptr xx = make_univariate_monomial(1,onempq);
    Poly_ptr ode_rhs_denom_poly = addPolynomials_AA(oneAA, xx,1);
    PowerSeries_t * ode_rhs_denominator = convertPolyToPowerSeries_PS(ode_rhs_denom_poly);
    PowerSeries_t * ode_rhs[3];
    for(int i=0; i<3; i++){
        ode_rhs[i] = dividePowerSeries_PS(deepCopy_PS(one_PSS), deepCopy_PS(ode_rhs_denominator)); 
    }
    //PowerSeries_t * ode_rhs = dividePowerSeries_PS(deepCopy_PS(one_PSS), deepCopy_PS(ode_rhs_denominator)); 

    // Auxiliaries to create the ode coefficients
    mpq_t current_coef; mpq_init(current_coef); mpq_set(current_coef,onempq);
    Poly_ptr running_sum = make_univariate_monomial(0,onempq); // initialize with 1
    Poly_ptr running_sum_register = running_sum;
    Poly_ptr current_monomial = make_univariate_monomial(1,onempq);
    running_sum = addPolynomials_AA(running_sum, current_monomial,1); // Starts at 1+x
    freePolynomial_AA(running_sum_register);
    freePolynomial_AA(current_monomial);
    PowerSeries_t * running_sum_PS;

    // The ode coefficients as power series
    PowerSeries_t * ode_coefs[3][ode_order+1];
    for(int i=0; i<=ode_order; i++){

        mpq_set_si(current_coef, i+2,1);
        current_monomial = make_univariate_monomial(i+2,current_coef);

        running_sum_register = running_sum;
        running_sum = addPolynomials_AA(running_sum_register, current_monomial,1);

        freePolynomial_AA(running_sum_register);
        freePolynomial_AA(current_monomial);
        
        running_sum_PS = convertPolyToPowerSeries_PS(running_sum);

        for(int j = 0; j<3; j++){
            ode_coefs[j][i] = dividePowerSeries_PS(deepCopy_PS(one_PSS), deepCopy_PS(running_sum_PS));
        }
        
        destroyPowerSeries_PS(running_sum_PS);
    }

    freePolynomial_AA(running_sum);
    mpq_clear(current_coef);

    // The initial conditions will all be (D^n)(y)(0) = 1
    mpq_t inits[ode_order];
    for (int i=0; i<ode_order; i++){
        mpq_init(inits[i]);
        mpq_set(inits[i], onempq);
    }

    // Helpers for printing
    mpq_t print_mpq_helper; mpq_init(print_mpq_helper);

    /* ---------------------------------------------------------------------------------------------
    ** Using the power series solution generator, serially **
    --------------------------------------------------------------------------------------------- */
    //*
    // Construct a LODE structure which will use the serial generator
    LODEoUPS_t * lode_serial_generator = construct_LODEoUPS(ode_order, ode_coefs[0], ode_rhs[0], inits);

    //Calculate the power series solution up to a degree using the generator
    auto start = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_using_sol_generator(precision, lode_serial_generator);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Using solution generator, serial, time for calculating " << precision << " coefficients: " << duration_ms.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the solution generator serially 
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_serial_generator->sol, precision, &zerompq));
    std::cout << "For using solution generator, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n\n" << std::endl;
    //*/

    /* ---------------------------------------------------------------------------------------------
    ** Using the truncated recursive method, serially **
    --------------------------------------------------------------------------------------------- */
    //*
    // Construct a LODE structure which will use the truncated method
    LODEoUPS_t * lode_truncated_serial = construct_LODEoUPS(ode_order, ode_coefs[1], ode_rhs[1], inits);

    // Calculate the power series solution up to a degree using truncated serial construction
    auto start1 = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_truncated_construction(precision, lode_truncated_serial);
    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration_ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
    std::cout << "Using truncated construction, serial, time for calculating " << precision << " coefficients: " << duration_ms1.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the truncated serial construction
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_truncated_serial->sol, precision, &zerompq));
    std::cout << "For using serial truncated construction, serial, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n\n" << std::endl;
    //*/

    /* ---------------------------------------------------------------------------------------------
    ** Using the parallel method, truncated **
    --------------------------------------------------------------------------------------------- */
    //*
    // Construct a LODE structure which will use the parallel method
    LODEoUPS_t * lode_parallel = construct_LODEoUPS(ode_order, ode_coefs[2], ode_rhs[2], inits);

    // Calculate the power series solution up to a degree using truncated serial construction
    auto start2 = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_truncated_parallel(precision, lode_parallel);
    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration_ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
    std::cout << "Using parallel construction, time for calculating " << precision << " coefficients: " << duration_ms2.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the truncated serial construction
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_parallel->sol, precision, &zerompq));
    std::cout << "For using truncated parallel, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n" << std::endl;
    //*/

    // Parallel speedup
    std::cout << "Parallel speedup: " << (double)duration_ms1.count()/duration_ms2.count() << "\n" << std::endl;



    // Free memory
    mpq_clear(print_mpq_helper);
    mpq_clear(onempq);
    mpq_clear(zerompq);
    for(int i=0; i<ode_order; i++){
        freePolynomial_AA(monomials[i]);
    }
    for(int j=0; j<3; j++){
        for(int i=0; i<ode_order; i++){
            destroyPowerSeries_PS(ode_coefs[j][i]);
        }
    }
    for(int i=0; i<ode_order; i++){
        mpq_clear(inits[i]);
    }
    destroy_LODEoUPS(lode_serial_generator);
    destroy_LODEoUPS(lode_truncated_serial);
    destroy_LODEoUPS(lode_parallel);
    freePolynomial_AA(oneAA);

    freePolynomial_AA(xx);
    freePolynomial_AA(ode_rhs_denom_poly);
    destroyPowerSeries_PS(ode_rhs_denominator);
    for(int i=0; i<3; i++){
        destroyPowerSeries_PS(ode_rhs[i]);
    }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/**
 * Runs timings on the exponential family of test ODE's of order n
 * @param ode_order the order of the modified Euler ODE
 * @param precision the precision of the ODE solution power series
 */
void testFamily_exp_demo(int ode_order,  int precision){

    // For printing header
    std::ostringstream oss;
    oss << "Test Family, Exponential, order " << ode_order << ", up to precision " << precision;
    std::string header_string = oss.str();
    test_header(header_string);

    // Helper mpq_t's
    mpq_t onempq; mpq_init(onempq); mpq_set_si(onempq,1,1);
    mpq_t zerompq; mpq_init(zerompq); mpq_set_si(zerompq,1,1);

    // mpq's which are factors of the exponential functions
    mpq_t exp_factors[ode_order + 1];

    // Give a value to the exponential function factors
    for(int i=0; i<=ode_order; i++){
        mpq_init(exp_factors[i]);
        mpq_set_si(exp_factors[i], i+1, 1);
    }

    // Create the ODE rhs which is e^(ode_order*x)
    mpq_t ode_order_mpq; mpq_init(ode_order_mpq); mpq_set_ui(ode_order_mpq, ode_order, 1);
    PowerSeries_t * ode_rhs[3];
    for(int i=0; i<3; i++){
        ode_rhs[i] = exponential_univariate(&onempq, &ode_order_mpq);
    }

    // The ode coefficients as power series
    PowerSeries_t * ode_coefs[3][ode_order+1];
    for(int j=0; j<3; j++){
        for(int i=0; i<=ode_order; i++){
            ode_coefs[j][i] = exponential_univariate(&exp_factors[i], &exp_factors[i]);
        }
    }


    // The initial conditions will all be (D^n)(y)(0) = 1
    mpq_t inits[ode_order];
    for (int i=0; i<ode_order; i++){
        mpq_init(inits[i]);
        mpq_set(inits[i], onempq);
    }

    /* ---------------------------------------------------------------------------------------------
    ** Using the power series solution generator, serially **
    --------------------------------------------------------------------------------------------- */

    // Construct a LODE structure which will use the serial generator
    LODEoUPS_t * lode_serial_generator = construct_LODEoUPS(ode_order, ode_coefs[0], ode_rhs[0], inits);

    // Calculate the power series solution up to a degree using the generator
    auto start = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_using_sol_generator(precision, lode_serial_generator);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Using solution generator, serial, time for calculating " << precision << " coefficients: " << duration_ms.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the solution generator serially
    mpq_t print_mpq_helper; mpq_init(print_mpq_helper); 
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_serial_generator->sol, precision, &zerompq));
    std::cout << "For using solution generator, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n\n" << std::endl;


    /* ---------------------------------------------------------------------------------------------
    ** Using the truncated recursive method, serially **
    --------------------------------------------------------------------------------------------- */

    
    // Construct a LODE structure which will use the truncated method
    LODEoUPS_t * lode_truncated_serial = construct_LODEoUPS(ode_order, ode_coefs[1], ode_rhs[1], inits);

    // Calculate the power series solution up to a degree using truncated serial construction
    start = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_truncated_construction(precision, lode_truncated_serial);
    end = std::chrono::high_resolution_clock::now();
    duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Using truncated construction, serial, time for calculating " << precision << " coefficients: " << duration_ms.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the truncated serial construction
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_truncated_serial->sol, precision, &zerompq));
    std::cout << "For using serial truncated construction, serial, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n\n" << std::endl;
    

    /* ---------------------------------------------------------------------------------------------
    ** Using the parallel method, truncated **
    --------------------------------------------------------------------------------------------- */
    
    // Construct a LODE structure which will use the parallel method
    LODEoUPS_t * lode_parallel = construct_LODEoUPS(ode_order, ode_coefs[2], ode_rhs[2], inits);

    // Calculate the power series solution up to a degree using truncated serial construction
    start = std::chrono::high_resolution_clock::now();
    solve_LODE_upto_degree_truncated_parallel(precision, lode_parallel);
    end = std::chrono::high_resolution_clock::now();
    auto duration_ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Using parallel construction, time for calculating " << precision << " coefficients: " << duration_ms2.count() << " ms\n" << std::endl;

    // Print the last coefficient of using the truncated serial construction
    mpq_set(print_mpq_helper, *get_univariate_power_series_coefficient_mpq(lode_parallel->sol, precision, &zerompq));
    std::cout << "For using truncated parallel, coefficient # " << precision << " of the solution is: " << mpq_numref(print_mpq_helper) << "\n ";
    std::cout << "/\n" << mpq_denref(print_mpq_helper) <<" \n\n" << std::endl;

    // Parallel speedup
    std::cout << "Parallel speedup: " << (double)duration_ms.count()/duration_ms2.count() << "\n" << std::endl;
    


    // Free memory
    mpq_clear(print_mpq_helper);
    mpq_clear(onempq);
    mpq_clear(zerompq);
    for(int j=0; j<3; j++){
        destroyPowerSeries_PS(ode_rhs[j]);
        for(int i=0; i<ode_order; i++){
            destroyPowerSeries_PS(ode_coefs[j][i]);
        }
    }
    for(int i=0; i<ode_order; i++){
        mpq_clear(inits[i]);
    }
    for(int i=0; i<=ode_order; i++){
        mpq_clear(exp_factors[i]);
    }
    mpq_clear(ode_order_mpq);
    destroy_LODEoUPS(lode_serial_generator);
    destroy_LODEoUPS(lode_truncated_serial);
    destroy_LODEoUPS(lode_parallel);

}


