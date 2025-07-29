#include "PowerSeries/LinearODEOverUnivariatePowerSeries_parallel.hpp"


/**
 * Solves a linear ODE up to the given degree, using a method that treats the degree as a fixed precision
 * The method assumes that no prior coefficients for the power series solution have been computed
 *  This is the parallel version using ExecutorThread Pool
 * @param deg the precision for the power series solution to be calculated up to
 * @param lode the linear ODE structure to solve
 */
void solve_LODE_upto_degree_truncated_parallel(int deg, LODEoUPS_t * lode){

    // Return if the desired coefficient has already been calculated
    if(lode->update_deg >= deg || lode->sol->deg >= deg){
        return;
    }

    // These are used in the calculations below
    int d = lode->ord;
    int k = deg - lode->ord; 

    // Number of parallel cores available
    int number_of_cores = std::thread::hardware_concurrency();
    int iterations_per_core = k/number_of_cores;

    // Parallel thread pool
    ExecutorThreadPool& pool = ExecutorThreadPool::getThreadPool();

    // We will update the requirements for solving the ODE up to the given degree, in parallel
    //*
    int d_iterations_per_core = d/number_of_cores;
    if( (deg >= lode->ord) && (deg > lode->update_deg)){

        // Destroy the old rf table
        if (lode->sol->deg >= lode->ord){
            destroy_rising_factorial_table_mpq_t(lode->rf, lode->ord, lode->update_deg - lode->ord);
        }

        // Update the ODE coefficients, in parallel
        for (int i=d; i>=d_iterations_per_core; i--){

            
            // Lambda function for the inner loop
            std::function<void()>  inner_loop = [=]() { 
                updateToDeg_PS(k+1, lode->coefs[i]);
            };

            // Update coefficient in parallel, if it has a generator
            if(lode->coefs[i]->genOrder != -1){
                pool.addTask(inner_loop);
            }
            

        }

        // Update a few ODE coefficients on the parent thread
        for (int i=d_iterations_per_core-1; i>= 0; i--){

            updateToDeg_PS(k+1, lode->coefs[i]);

        }
        
        // Update the ODE rhs, done serially to use parent thread
        updateToDeg_PS(k+1, lode->forcing);

        // Update the rising factorial table, done serieally to use parent thread
        lode->rf = rising_factorial_table_recursive_mpq_t(lode->ord, k);

        lode->update_deg = deg;


    }
    //*/


    // Update the requirements for solving up to the given degree (done serially)
    //update_LODE_auxiliaries(deg, lode);

    // Update the power series solution to include the initial conditions (constant time operation)
    updateToDeg_PS(lode->ord-1, lode->sol);

    // This is used in the calculations below
    mpq_t ** rf = lode->rf;

    // Create a vector which will contain the power series solution coefficients
    mpq_t * c = (mpq_t*)malloc((d+k+1) * sizeof(mpq_t));

    // Copy in the initial conditions
    for(int i=0; i<d; i++){
        mpq_init(c[i]); mpq_set(c[i], GC_Y(i));
    }

    // Array with just the necessary terms of the forcing function
    mpq_t * b = (mpq_t *) malloc((k+1) * sizeof(mpq_t));
    for(int n=0; n<=k; n++){
        mpq_init(b[n]);
        mpq_set(b[n], GETCOEF_MPQ(lode->forcing, n));
    }

    // Parallel sync from updating auxiliaries
    pool.waitForAllThreads();

    // Array with just the necessary terms of the ODE coefficients
    mpq_t ** a = (mpq_t **) malloc((d+1) * sizeof(mpq_t *) );
    for(int i=0; i<=d; i++){
        a[i] = (mpq_t *) malloc((k+1) * sizeof(mpq_t));
    }
    for(int j=0; j<=d; j++){
        for(int n=0; n<=k; n++){
            mpq_init(a[j][n]);
            mpq_set(a[j][n], GC_A(j,n));
        }
    }

    // Set up and calculate the w vector
    mpq_t * w = (mpq_t *)malloc( (k+1) * sizeof(mpq_t));
    for(int i=0; i<=k; i++){
        mpq_init(w[i]); mpq_set(w[i], rf[d][i]);
        mpq_mul(w[i], w[i], a[d][0]); // w[i] = a[d][0] * rf[d][i]
        // mpq_inv(w[i], w[i]); // w[i] = 1/a[d][0] * 1/rf[d][i]
    }

    /*-----------------------------------------------------------------------------------------------
    Calculate table of g_{n,i's}
        g_{n,i} -> n ranges from 0 to k
                -> i ranges from 0 to d-1+n
    -----------------------------------------------------------------------------------------------*/

    // Set up the table g with nice indexing
    int size_of_g_table = d*(k+1) + (k)*(k+1)/2;
    mpq_t * g_data = (mpq_t*)malloc( (size_of_g_table )*sizeof(mpq_t) );
    mpq_t ** g = (mpq_t**)malloc((k+1) * sizeof(mpq_t*));
    g[0] = g_data;
    for(int i=1; i<=k; i++){
        g[i] = g[i-1] + d - 1 + i;
    }

    /*-----------------------------------------------------------------------------------------------*/
    // Calculation of g table coefficients is where most of the work is, so we will parallelize it
    /*-----------------------------------------------------------------------------------------------*/
    
    // We will try to estimate a big enough pre-allocation size for the necessary mpq_t's
    //*
    size_t limbs_in_biggest_rf = mpz_size(mpq_numref(rf[d][k-1]));
    size_t limbs_in_biggest_coefficient_numerator = mpz_size(mpq_numref(a[0][k]));
    size_t limbs_in_biggest_coefficient_denominator = mpz_size(mpq_denref(a[0][k]));
    for(int i=1; i<=d; i++){
        if(limbs_in_biggest_coefficient_numerator < mpz_size(mpq_numref(a[i][k])) ){
            limbs_in_biggest_coefficient_numerator = mpz_size(mpq_numref(a[i][k]));
        }
        if(limbs_in_biggest_coefficient_denominator < mpz_size(mpq_denref(a[i][k])) ){
            limbs_in_biggest_coefficient_denominator = mpz_size(mpq_denref(a[i][k]));
        }
    }
    limbs_in_biggest_coefficient_denominator*=2;
    limbs_in_biggest_coefficient_numerator*=2; limbs_in_biggest_coefficient_numerator+=limbs_in_biggest_rf;

    size_t bits_in_biggest_coefficient_denominator = limbs_in_biggest_coefficient_denominator * mp_bits_per_limb;
    size_t bits_in_biggest_coefficient_numerator = limbs_in_biggest_coefficient_numerator * mp_bits_per_limb;
    //*/

    
    // This array of temporary mpq_t's is made to attempt to minimize memory operations by threads
    mpq_t * iter_save_array = (mpq_t * )malloc((k+1) * sizeof(mpq_t));
    //mpq_t iter_save_array[k+1];
    for(int i=0; i<=k; i++){
        mpq_init(iter_save_array[i]);
        mpz_realloc2(mpq_numref(iter_save_array[i]), bits_in_biggest_coefficient_numerator);
        mpz_realloc2(mpq_denref(iter_save_array[i]), bits_in_biggest_coefficient_denominator);
    }

    // Pre-initialize the g table, to ensure memory safety
    for(int i=0; i< size_of_g_table; i++){
        mpq_init(g_data[i]);
        //mpz_realloc2(mpq_numref(g_data[i]), bits_in_biggest_coefficient_numerator);
        //mpz_realloc2(mpq_denref(g_data[i]), bits_in_biggest_coefficient_denominator);
    }

    /*
    // Array of indices for each entry in the g table, used to granularize the parallelization
    int g_n_indices[size_of_g_table], g_i_indices[size_of_g_table];
    int g_index_accumulated[k+1];
    for(int n=0; n<=k; n++){

        g_index_accumulated[n] = n*d - 1 + n*(n+1)/2;
        for(int i=0; i<n+d; i++){
            g_n_indices[g_index_accumulated[n]] = n;
            g_i_indices[g_index_accumulated[n]] = i;
            g_index_accumulated[n] ++;
        }

    }
    */


    // Calculate the g table in parallel
    for(int n=k; n>=iterations_per_core; n--){

        // Lambda function for the inner loop
        std::function<void()>  inner_loop = [=, &iter_save_array]() { 

            for(int i=0; i<n+d; i++){
                //mpq_init(g[n][i]); 
                mpq_set_ui(g[n][i],0,1);
                for(int m=d; m>=0; m--){
                    if( (i>=m) && (n>=(i-m)) ){
                        mpq_mul(iter_save_array[n], rf[m][i-m], a[m][n-i+m]);
                        mpq_sub(g[n][i], g[n][i], iter_save_array[n]);
                    }
                }
            }

        };

        // Add function call to the Executor pool
        pool.addTask(inner_loop);

    }

    // Calculate the first iterations_per_core iterations on the parent thread
    for(int n=iterations_per_core-1; n>=0; n--){

        for(int i=0; i<n+d; i++){
            //mpq_init(g[n][i]); 
            mpq_set_ui(g[n][i],0,1);
            for(int m=d; m>=0; m--){
                if( (i>=m) && (n>=(i-m)) ){
                    mpq_mul(iter_save_array[n], rf[m][i-m], a[m][n-i+m]);
                    mpq_sub(g[n][i], g[n][i], iter_save_array[n]);
                }
            }
        }

    }


    // Set up the T table
    int size_of_T_table = (k+1)*(k+2)/2;
    mpq_t * T_data = (mpq_t *)malloc( (size_of_T_table )*sizeof(mpq_t) );
    mpq_t ** T = (mpq_t**)malloc((k+1) * sizeof(mpq_t*));
    T[0] = T_data;
    for(int i=1; i<=k; i++){
        T[i] = T[i-1] + i;
    }

    // Pre-initialize the T table, to ensure memory safety
    for(int i=0; i< size_of_T_table; i++){
        mpq_init(T_data[i]);
    }

    // Synchronize the parallelization of the for loops from g table calculation
    pool.waitForAllThreads();

    /*-----------------------------------------------------------------------------------------------*/
    // Calculate the T table and the c[i]'s
    /*-----------------------------------------------------------------------------------------------*/

    // Calculate the T[i][0]'s - the entries in the leftmost column in the T table
    // There is not as much work here in other parts, but we will paralellize it anyway
    for(int i=k; i>=iterations_per_core; i--){

        // Lambda function for the inner loop
        std::function<void()>  inner_loop = [=, &iter_save_array]() { 

            mpq_init(T[i][0]); mpq_set(T[i][0], b[i]);
            for(int j=d-1; j>=0; j--){
                mpq_mul(iter_save_array[i], g[i][j], c[j]);
                mpq_add(T[i][0], T[i][0], iter_save_array[i]);
            }

        };

        pool.addTask(inner_loop);

    }

    // Calculate the first iterations_per_core iterations on the parent thread
    for(int i=iterations_per_core-1; i>=0; i--){

        mpq_init(T[i][0]); mpq_set(T[i][0], b[i]);
        for(int j=d-1; j>=0; j--){
            mpq_mul(iter_save_array[i], g[i][j], c[j]);
            mpq_add(T[i][0], T[i][0], iter_save_array[i]);
        }

    }

    // Synchronize the for loop from calculating first column of T table
    pool.waitForAllThreads();

    // Clear the temporary mpq_t's
    for(int i=0; i<=k; i++){
        mpq_clear(iter_save_array[i]);
    }
    free(iter_save_array);

    // Calculate T[0][0]
    mpq_div(T[0][0], T[0][0], w[0]);
    mpq_init(c[d]); mpq_set(c[d], T[0][0]);


    // Call the parallel recursive T table calculation
    calculate_T_table_triangular_parallel(1,k,1,k,d,w,g,c,T);

    /*-----------------------------------------------------------------------------------------------*/
    // Update the LODEoUPS object
    /*-----------------------------------------------------------------------------------------------*/

    lode->update_deg = deg;

    resizePowerSeries_PS(lode->sol, deg+1);

    // Copy in the calculated coefficients
    for(int i=d; i<=deg; i++){
        lode->sol->polys[i] = make_univariate_monomial(i, c[i]);
    }

    lode->sol->deg = deg;



    /*-----------------------------------------------------------------------------------------------*/
    // Free used Memory
    /*-----------------------------------------------------------------------------------------------*/

    
    for(int j=0; j<=d; j++){
        for(int n=0; n<=k; n++){
            mpq_clear(a[j][n]);
        }
        free(a[j]);
    }  
    free(a);
    for(int n=0; n<=k; n++){
        mpq_clear(b[n]);
    }
    free(b);
        for(int n=0; n<=k; n++){
        mpq_clear(w[n]);
    }
    free(w);
    for(int n=0; n<=d+k; n++){
        mpq_clear(c[n]);
    }
    free(c);
    for(int i=0; i< size_of_g_table; i++){
        mpq_clear(g_data[i]);
    }
    free(g);
    free(g_data);
    for(int i=0; i< size_of_T_table; i++){
        mpq_clear(T_data[i]);
    }
    free(T);
    free(T_data);


    // Wait for all threads doing deallocations
    pool.waitForAllThreads();


}


/**
 * An internal function - Parallel version using ExecutorThread Pool
 * Recursively calculates the "T" table which at position T[i][i] holds the coefficient  c[d+i] of the power series
 * Solution to an ODE
 * This is the triangular recursive call, for when the section of the T table being calculated is triangular
 * @param left the leftmost index in the T table being calculated in this recursive call
 * @param right the rightmost index in the T table being calculated in this recursive call
 * @param low the lowest vertical index in the T table being calculated in this recursive call
 * @param high the highest vertical index in the T table being calculated in this recursive call
 * @param d the order of the differential equation
 * @param w the w vector, which contains the constant terms 1/rf[d][i] * 1/a[d][0]
 * @param g the g table
 * @param c an mpq_t vector which contains the coefficients of the power series solution to the ODE
 * @param T the T table being calculated

 */
void calculate_T_table_triangular_parallel(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T){

    // Very basic work estimate slightly taking into account the number of digits
    double input_size_estimate = sqrt(high*right);
    double input_size_bits = log2(input_size_estimate);
    double bits_estimate = input_size_bits * log2(input_size_bits);
    double work_estimate = (high-low)*(right-left)* (bits_estimate) * log2(bits_estimate) * log2(bits_estimate);

    if( work_estimate < PARALLEL_WORK_BASE ){ // The base case

        // Call the Serial version
        calculate_T_table_triangular(left, right, low, high, d, w, g, c, T);

    }
    else{ // Recursive Calls

        int mid_v = (low+high)/2;
        int mid_h = (left+right)/2;

        // None of these recursive calls can be done in parallel, they must be done in this order

        calculate_T_table_triangular_parallel(left, mid_h, low, mid_v, d, w, g, c, T);

        calculate_T_table_rectangular_parallel(left, mid_h, mid_v+1, high, d, w, g, c, T);

        calculate_T_table_triangular_parallel(mid_h+1, right, mid_v+1, high, d, w, g, c, T);

    }

}

/**
 * An internal function - Parallel version using ExecutorThread Pool
 * Recursively calculates the "T" table which at position T[i][i] holds the coefficient  c[d+i] of the power series
 * Solution to an ODE
 * This is the rectangular recursive call, for when the section of the T table being calculated is a rectangle
 * @param left the leftmost index in the T table being calculated in this recursive call
 * @param right the rightmost index in the T table being calculated in this recursive call
 * @param low the lowest vertical index in the T table being calculated in this recursive call
 * @param high the highest vertical index in the T table being calculated in this recursive call
 * @param d the order of the differential equation
 * @param w the w vector, which contains the constant terms 1/rf[d][i] * 1/a[d][0]
 * @param g the g table
 * @param c an mpq_t vector which contains the coefficients of the power series solution to the ODE
 * @param T the T table being calculated

 */
void calculate_T_table_rectangular_parallel(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T){

    // Very basic work estimate slightly taking into account the number of digits
    double input_size_estimate = sqrt(high*right);
    double input_size_bits = log2(input_size_estimate);
    double bits_estimate = input_size_bits * log2(input_size_bits);
    double work_estimate = (high-low)*(right-left)* (bits_estimate) * log2(bits_estimate) * log2(bits_estimate) * 2;

    if( work_estimate < PARALLEL_WORK_BASE ){ // The base case

        // Call the Serial version
        calculate_T_table_rectangular(left, right, low, high, d, w, g, c, T);

    }
    else{ // Recursive Calls

        int mid_v = (low+high)/2;
        int mid_h = (left+right)/2;

        // Parallel thread pool setup
        ExecutorThreadPool::threadID id; 
        ExecutorThreadPool & pool = ExecutorThreadPool::getThreadPool();

        /*----------------------------------------------------------------------------------*/
        // The first two recursive calls can be done in parallel
        /*----------------------------------------------------------------------------------*/

        // Lambda function for the first recursive function call
        std::function<void()>  p_r1 = [=]() { 
            calculate_T_table_rectangular_parallel(left, mid_h, low, mid_v, d, w, g, c, T);
        };

        // Obtain a thread from the pool
        //id = pool.obtainThread(); 
        pool.obtainThread(id); 

        // Execute first function call in parallel on other thread
        pool.executeTask(id, p_r1); 

        // Do the second recursive call in parallel on the calling thread 
        calculate_T_table_rectangular_parallel(left, mid_h, mid_v+1, high, d, w, g, c, T);

        pool.returnThread(id); // Parallel sync
        //pool.waitForThread(id); // Parallel sync



        /*----------------------------------------------------------------------------------*/
        /* The next two recursive calls can also be done in parallel */
        /*----------------------------------------------------------------------------------*/

        //id = pool.obtainThread(); // Obtain a new thread
        pool.obtainThread(id);

        // Lambda function for the third recursive function call
        std::function<void()>  p_r2 = [=]() { 
            calculate_T_table_rectangular_parallel(mid_h+1, right, low, mid_v, d, w, g, c, T);
        };

        // Execute first function call in parallel on other thread
        pool.executeTask(id, p_r2);

        //calculate_T_table_rectangular_parallel(mid_h+1, right, low, mid_v, d, w, g, c, T);


        // Do the second recursive call on the calling thread in parallel
        calculate_T_table_rectangular_parallel(mid_h+1, right, mid_v+1, high, d, w, g, c, T);

        pool.returnThread(id); // Parallel sync


    }

}





