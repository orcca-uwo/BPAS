#include "PowerSeries/LinearODEOverUnivariatePowerSeries.h"


/**
 * Allocates necessary memory for a LODEoUPS
 *  @param ord : the order of the ODE
 */
LODEoUPS_t* allocate_LODEoUPS(int ord){

    LODEoUPS_t * lode = (LODEoUPS_t*) malloc(sizeof(LODEoUPS_t));

    if(ord < 1){
        fprintf(stderr, "ERROR: DE order must be at least 1\n");
        exit(1);
    }

    lode->coefs = (PowerSeries_t **) malloc( (ord+1) * sizeof(PowerSeries_t *));
    for(int i=0; i<=ord; i++){
        lode->coefs[i] = allocatePowerSeries_PS(1);
    }

    lode->sol = allocatePowerSeries_PS(1);

    lode->forcing = allocatePowerSeries_PS(1);

    lode->inits = (ratNum_t * ) malloc((ord+1) * sizeof(ratNum_t));

    mpq_init(lode->zero_mpq); mpq_set_si(lode->zero_mpq, 0, 1);

    return lode;

}

/**
 * Constructor for LODEoUPS structure
 * @param ord the order of the ODE
 * @param coefs an array of univariate power series which represent the coefficients of the ode. These will be deep copied
 * @param forcing the forcing function of the ODE (the LHS of the equation)
 * @param inits the initial conditions for the ODE, rational mpq_t numbers
*/
LODEoUPS_t * construct_LODEoUPS(int ord, PowerSeries_t ** coefs, PowerSeries_t * forcing, ratNum_t * inits){

    // Check that the first coefficient is a unit
    if (!isUnit_PS(coefs[ord])) {
        fprintf(stderr, "ERROR: Leading coefficient of the ODE must be a unit \n" );
        exit(1);
    }

    // Check that all ODE coefficients are univariate or zero
    for (int i=0; i<= ord; i++){
        if((coefs[i]->nvar != 1) && (coefs[i]->nvar != -1) ){
            fprintf(stderr, "ERROR: ODE coefficients must all be univariate power series\n");
            exit(1);
        }
    }

    // Check that all ODE rhs coefficient (forcing) is univariate or zero
    if((forcing->nvar != 1) && (forcing->nvar != -1) ){
        fprintf(stderr, "ERROR: ODE rhs must be univariate power series\n");
        exit(1);
    }

    LODEoUPS_t * lode = allocate_LODEoUPS(ord);

    lode->ord = ord;

    lode->forcing = deepCopy_PS(forcing);

    // Copy the initial conditions array (deep copy)
    for(int i=0; i<ord; i++){
        mpq_init(lode->inits[i]);
        mpq_set(lode->inits[i], inits[i]);
    }
    turn_LODE_initial_conditions_to_coefficients(lode);

    // Copy the power series coefficients (deep copy)
    for(int i=0; i<= ord; i++){
        lode->coefs[i] = deepCopy_PS(coefs[i]);
    }

    lode->update_deg = 0;


    /*************************************************************
     * Set up the power series solution sol
    *************************************************************/

    // This just helps to make the code below a bit cleaner
    PowerSeries_t * lodesol = lode->sol;

    lodesol->nvar = 1;

    // We will just pass the whole LODEoUPS structure as a single parameter to the generator
    lodesol->genParam1 = (void*) lode; 

    // Use this type so destroyPowerSeries_PS doesn't destroy the LODEoUPS
    lodesol->paramType1 = PLAIN_DATA; 

    lodesol->genOrder = 1;

    lodesol->gen.unaryGen = &(homogPartVoid_LODEoUPS_sol_PS);

    lodesol->polys[0] = make_univariate_monomial(0,lode->inits[0]);


    /* // CODE to Copy initial conditions during creation (obsolete)

    // Resize to include the initial conditions
    resizePowerSeries_PS(lodesol, lode->ord);

    // Copy the initial conditions into the power series solution
    for(int i=0; i<ord; i++){
        lodesol->polys[i] = make_univariate_monomial(i, inits[i]);
    }

    lodesol->deg = lode->ord - 1; // Is it ord or ord-1?

    */

    return lode;


}

/**
 * Function to deallocate and clear memory for LODEoUPS structure
 * @param lode the LODEoUPS structure
 */
void destroy_LODEoUPS(LODEoUPS_t * lode){

    if (lode->sol->deg >= lode->ord){
        destroy_rising_factorial_table_mpq_t(lode->rf, lode->ord, lode->sol->deg - lode->ord);
    }

    for(int i=0; i<=lode->ord; i++){
        destroyPowerSeries_PS(lode->coefs[i]);
    }
    free(lode->coefs);

    destroyPowerSeries_PS(lode->sol);

    destroyPowerSeries_PS(lode->forcing);

    for(int i=0; i<lode->ord; i++){
        mpq_clear(lode->inits[i]);
    }
    free(lode->inits);

    mpq_clear(lode->zero_mpq);

    free(lode);

}

/**
 * Function to return the coefficient of a power series
 * Had to make this instead of a macro for the case where power series do not have a generator that returns 0
 * @param ps the power series to return the coefficient from
 * @param n the degree of the power series coefficient to return
 * @param zero_mpq an mpq_t which has the value 0 (0/1). The user has to manually manage (clear) the memory used by this mpq
 */
mpq_t * get_univariate_power_series_coefficient_mpq(PowerSeries_t * ps, int n, mpq_t * zero_mpq){


    if(ps->genOrder == -1){ // This means there is essentially no generator

        if (ps->deg == -1){ // This means that the power series is the zero power series
            return zero_mpq;
        }
        else{ // We have a fixed number of terms

            if(n <= ps->deg){ // The nth degree term in the power series exists
                if (ps->polys[n] == NULL){ // This happens when a polynomial term's coefficient is 0
                    return zero_mpq;
                }
                else{
                    return &(ps->polys[n]->elems[0].coef);
                }              
            }
            else{ // The nth term doesn't exist, so return 0
                return zero_mpq;
            }

        }

    }
    else{ // There is a generator
        if (ps->polys[n] == NULL){ // Haven't checked closely if this makes sense, copied from above
            return zero_mpq;
        }
        else{
            return &(ps->polys[n]->elems[0].coef);
        }
    }

}

/**
 * Internal function
 * Turns the initial conditions given to the ODE which are y(0), y'(0), y''(0), etc
 * Into coefficients in the power series solution, that is, divide the initial condition by a factorial
 * @param lode the LODE structure for which we are correcting the intial conditions
 */
void turn_LODE_initial_conditions_to_coefficients(LODEoUPS_t * lode){

    mpq_t fact_ratio; mpq_init(fact_ratio); mpq_set_ui(fact_ratio,1,1);
    mpq_t iter_save; mpq_init(iter_save); mpq_set_ui(iter_save,1,1);

    for(unsigned long i=2; i< lode->ord; i++){

        mpq_set_ui(iter_save, 1, i); // 1/i
        mpq_mul(fact_ratio, iter_save, fact_ratio); // Saves 1/i!
        mpq_mul(iter_save, lode->inits[i], fact_ratio); // Multiply initial condition by 1/i!

        mpq_set(lode->inits[i], iter_save);

    }

    mpq_clear(fact_ratio);
    mpq_clear(iter_save);

}


/**
 * Helper function to make a univariate monomial of degree d with a specific rational coefficient
 * @param deg the degree of the univariate monomial
 * @param coef the rational coefficient of the univariate monomial
 */
Poly_ptr make_univariate_monomial(int deg, mpq_t coef){

    Poly_ptr unimon = makeConstPolynomial_AA(1, 1, coef);

    
    unimon->elems->degs = deg;


    return unimon;

}



/**
 * An internal function.
 * A void generator wrapper for the generator of a
 * power series solution for a Linear ODE with Univariate Power Series coefficients
 * @param d: the requested degree to generate
 * @param param1: the LODEoUPS structure
 * @return the homogeneous part of degree d of the power series solution
 */
Poly_ptr homogPartVoid_LODEoUPS_sol_PS(int d, void* param1) {
    return homogPart_LODEoUPS_sol_PS(d, (LODEoUPS_t*) param1);
}


/**
 * An internal function.
 * Computes the homogeneous part of the
 * power series solution for a Linear ODE with Univariate Power Series coefficients
 * @param deg: the requested degree to generate
 * @param lode: the LODEoUPS structure
 * @return the homogeneous part of degree d of the power series solution
 */
Poly_ptr homogPart_LODEoUPS_sol_PS(int deg, LODEoUPS_t * lode) {

    // If the requested coefficient is within the initial conditions, we will return immediately
    if(deg<lode->ord){
        return make_univariate_monomial(deg, lode->inits[deg]);
    }

    if (lode->update_deg < deg){
        update_LODE_auxiliaries(deg, lode);
    }

    int k = deg - lode->ord; // The "k" value of the of the new coefficients to be calculated
    int d = lode->ord;
    mpq_t ** rf = lode->rf;


    // Auxiliaries
    mpq_t sum_store; mpq_init(sum_store);
    mpq_t inner_store; mpq_init(inner_store);
    mpq_t product_store; mpq_init(product_store);

    mpq_set(product_store, GETCOEF_MPQ(lode->forcing, k)); // product_store = b_k

    // The coefficient to be returned, as a rational number
    mpq_t cdplusk; mpq_init(cdplusk);
    mpq_t cdplusk_temp; mpq_init(cdplusk_temp);

    // cdplusk = 1/rf[d][k] * 1/ a[d][0]
    mpq_mul(cdplusk_temp, rf[d][k], GC_A(d,0));
    mpq_inv(cdplusk_temp, cdplusk_temp);

    
    for(int i=d+k-1; i>=0; i--){
        mpq_set_ui(sum_store,0,1); // reset sum_store = 0

        for(int m=d; m>=0; m--){

            if( (i>=m) && (k>=(i-m)) ){
                mpq_mul(inner_store, rf[m][i-m], GC_A(m,k-i+m)); 
                mpq_add(sum_store, sum_store, inner_store); // sum_store += rf[m][i-m] * A[m][k-i+m]
            }

        }

        mpq_mul(sum_store, sum_store, GC_Y(i)); // sum_store *= C[i]
        mpq_sub(product_store, product_store, sum_store); // product_store -= sum_store

    }
    
    
    mpq_mul(cdplusk, cdplusk_temp, product_store); // c[d+k] *= product_store

    // Convert the mpq_k coefficient to a Polynomial to be a coefficient of the power series
    Poly_ptr cdplusk_poly = make_univariate_monomial(deg, cdplusk);

    mpq_clear(cdplusk); mpq_clear(product_store); mpq_clear(sum_store); mpq_clear(inner_store); mpq_clear(cdplusk_temp);

    return cdplusk_poly; 

}


/**
 * An internal function.
 * Computes the homogeneous part of the
 * power series solution for a Linear ODE with Univariate Power Series coefficients
 * Uses the first non-optimized formula
 * @param deg: the requested degree to generate
 * @param lode: the LODEoUPS structure
 * @return the homogeneous part of degree d of the power series solution
 */
Poly_ptr homogPart_LODEoUPS_sol_PS_first(int deg, LODEoUPS_t * lode) {

    // If the requested coefficient is within the initial conditions, we will return immediately
    if(deg<lode->ord){
        return make_univariate_monomial(deg, lode->inits[deg]);
    }

    if (lode->update_deg < deg){
        update_LODE_auxiliaries(deg, lode);
    }

    int k = deg - lode->ord; // The "k" value of the of the new coefficients to be calculated
    int d = lode->ord;
    mpq_t ** rf = lode->rf;

    // Auxiliary
    mpq_t sum_store; mpq_init(sum_store);

    // The coefficient to be returned, as a rational number
    mpq_t cdplusk; mpq_init(cdplusk);

    // cdplusk = 1/rf[d][k] * 1/ a[d][0]
    mpq_set_ui(sum_store,1,1);
    mpq_mul(cdplusk, rf[d][k], GC_A(d,0));
    mpq_inv(cdplusk, cdplusk);

    mpq_t product_store; mpq_init(product_store);
    mpq_set(product_store, GETCOEF_MPQ(lode->forcing, k)); // product_store = b_k


    for(int j=0; j<d; j++){
        for(int n=0; n<=k; n++){
            //product_store -= rf[j][n] * A[j][k-n] * C[j+n]
            mpq_mul(sum_store, rf[j][n], GC_A(j,k-n)); // sum_store = rf[j][n] * A[j][k-n]
            mpq_mul(sum_store, sum_store, GC_Y(j+n)); // sum_store *= C[j+n]
            mpq_sub(product_store, product_store, sum_store); // product_store -= sum_store
        }
    }
    for(int n=0; n<k; n++){
        //product_store -= rf[d][n] * A[d][k-n] * C[d+n]
        mpq_mul(sum_store, rf[d][n], GC_A(d,k-n)); // sum_store = rf[d][n] * A[d][k-n]
        mpq_mul(sum_store, sum_store, GC_Y(d+n)); // sum_store *= C[d+n]
        mpq_sub(product_store, product_store, sum_store); // product_store -= sum_store
    }

    mpq_mul(cdplusk, cdplusk, product_store); // c[d+k] *= product_store

    // Convert the mpq_k coefficient to a Polynomial to be a coefficient of the power series
    Poly_ptr cdplusk_poly = make_univariate_monomial(deg, cdplusk);

    mpq_clear(cdplusk); mpq_clear(product_store); mpq_clear(sum_store);

    return cdplusk_poly; 

}



/**
 * Solves a linear ODE up to the given degree using the generator for the solution power series
 * @param deg the precision for the power series solution to be calculated up to
 * @param lode the linear ODE structure to solve
 */
void solve_LODE_upto_degree_using_sol_generator(int deg, LODEoUPS_t * lode){

    // Update all LODE auxiliaries to the necessary degree
    update_LODE_auxiliaries(deg, lode);

    //Now call the generator
    updateToDeg_PS(deg, lode->sol);

}

/**
 * Internal function
 * This function will update the LODE coefficients and rising factorial table to 
 * The depth required to calculate a desired degree of the power series solution
 * This helps so these operations are only done once, and not many times in the generator
 * @param deg the desired degree of the power series solution
 * @param lode the LODE structure
 */
void update_LODE_auxiliaries(int deg, LODEoUPS_t * lode){

    //printf("updating_aux");

    if (deg < lode->ord ){
        return;
    }
    else if (deg <= lode->update_deg) {
        return;
    }

    // Update all ODE coefficients and forcing function b to the necessary precision
    int k = deg -lode->ord;
    for (int i=0; i<= lode->ord; i++){
        updateToDeg_PS(k+1, lode->coefs[i]);
    }
    updateToDeg_PS(k+1, lode->forcing);

    // Update the rf table
    if (lode->sol->deg >= lode->ord){
        destroy_rising_factorial_table_mpq_t(lode->rf, lode->ord, lode->update_deg - lode->ord);
    }
    lode->rf = rising_factorial_table_recursive_mpq_t(lode->ord, k);

    lode->update_deg = deg;

}


/**
 * Solves a linear ODE up to the given degree, using a method that treats the degree as a fixed precision
 * The method assumes that no prior coefficients for the power series solution have been computed
 * @param deg the precision for the power series solution to be calculated up to
 * @param lode the linear ODE structure to solve
 */
void solve_LODE_upto_degree_truncated_construction(int deg, LODEoUPS_t * lode){

    // Return if the desired coefficient has already been calculated
    if(lode->update_deg >= deg || lode->sol->deg >= deg){
        return;
    }

    // Update the requirements for solving up to the given degree (done serially)
    update_LODE_auxiliaries(deg, lode);

    // Update the power series solution to include the initial conditions (constant time operation)
    updateToDeg_PS(lode->ord-1, lode->sol);

    mpq_t ** rf = lode->rf;
    int d = lode->ord;
    int k = deg - lode->ord; 

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

    /* 
    Calculate table of g_{n,i's}
        g_{n,i} -> n ranges from 0 to k
                -> i ranges from 0 to d-1+n
    */

    // Set up the table g with nice indexing
    int size_of_g_table = d*(k+1) + (k)*(k+1)/2;
    mpq_t * g_data = malloc( (size_of_g_table )*sizeof(mpq_t) );
    mpq_t ** g = malloc((k+1) * sizeof(mpq_t*));
    g[0] = g_data;
    for(int i=1; i<=k; i++){
        g[i] = g[i-1] + d - 1 + i;
    }

    // Pre-initialize the g table, to ensure memory safety
    for(int i=0; i< size_of_g_table; i++){
        mpq_init(g_data[i]);
    }

    //mpq_t iter_save; mpq_init(iter_save);

    // Calculate the g table
    for(int n=0; n<=k; n++){

        for(int i=0; i<n+d; i++){
            mpq_t iter_save; mpq_init(iter_save);
            //mpq_init(g[n][i]); 
            mpq_set_ui(g[n][i],0,1);
            for(int m=d; m>=0; m--){
                if( (i>=m) && (n>=(i-m)) ){
                    mpq_mul(iter_save, rf[m][i-m], a[m][n-i+m]);
                    mpq_sub(g[n][i], g[n][i], iter_save);
                }
            }
            mpq_clear(iter_save);           
        }

    }

    // Set up the T table
    int size_of_T_table = (k+1)*(k+2)/2;
    mpq_t * T_data = (mpq_t *)malloc( (size_of_T_table )*sizeof(mpq_t) );
    mpq_t ** T = (mpq_t * )malloc((k+1) * sizeof(mpq_t*));
    T[0] = T_data;
    for(int i=1; i<=k; i++){
        T[i] = T[i-1] + i;
    }

    // Pre-initialize the T table, to ensure memory safety
    for(int i=0; i< size_of_T_table; i++){
        mpq_init(T_data[i]);
    }

    // Calculate the T[i][0]'s - the entries in the leftmost column in the T table
    for(int i=0; i<=k; i++){

        mpq_t iter_save; mpq_init(iter_save);
        //mpq_init(T[i][0]); 
        mpq_set(T[i][0], b[i]);
        for(int j=d-1; j>=0; j--){
            mpq_mul(iter_save, g[i][j], c[j]);
            mpq_add(T[i][0], T[i][0], iter_save);
        }
        mpq_clear(iter_save);

    }
    // Calculate T[0][0]
    mpq_div(T[0][0], T[0][0], w[0]);
    mpq_init(c[d]); mpq_set(c[d], T[0][0]);

    // Call the recursive T table calculation
    calculate_T_table_triangular(1,k,1,k,d,w,g,c,T);


    // Update the LODEoUPS object
    //------------------------------------------------------------------------------------------------

    lode->update_deg = deg;

    resizePowerSeries_PS(lode->sol, deg+1);

    for(int i=d; i<=deg; i++){
        lode->sol->polys[i] = make_univariate_monomial(i, c[i]);
    }

    lode->sol->deg = deg;



    // Free used memory
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
    for(int n=0; n<=k; n++){
        mpq_clear(w[n]);
    }
    free(w);
    for(int n=0; n<=d+k; n++){
        mpq_clear(c[n]);
    }
    free(c);

}

/**
 * An internal function
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
void calculate_T_table_triangular(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T){

    // Very basic memory use estimate slightly taking into account the number of digits
    double input_size_estimate = sqrt(high*right);
    double input_size_bits = log2(input_size_estimate);
    double bytes_estimate = input_size_bits * log2(input_size_bits) / 8;
    double memory_use_estimate = (high-low)*(right-left)*bytes_estimate /2;

    if( memory_use_estimate < SERIAL_SIZE_BASE ){ // The base case

        mpq_t calc_temp; mpq_init(calc_temp);
        for(int i=low; i<=high; i++){
            for(int j=left; j<=right; j++){

                if(j<=i){
                    //mpq_init(T[i][j]);
                    mpq_mul(calc_temp, g[i][d+j-1], c[d+j-1]);
                    mpq_add(T[i][j], calc_temp, T[i][j-1]);
                }

                if(i==j){ // The entry at T[i][i] will be the d+i th coefficient of the solution
                    mpq_set(calc_temp, T[i][j]);                                        
                    mpq_div(T[i][j], calc_temp, w[i]); // divide by a[d][0] * rf[d][i]
                    mpq_init(c[d+i]); mpq_set(c[d+i], T[i][j]); // Copy the table value into the c array for easier access
                }

            }
        }
        mpq_clear(calc_temp);



    }
    else{ // Recursive Calls

        int mid_v = (low+high)/2;
        int mid_h = (left+right)/2;

        calculate_T_table_triangular(left, mid_h, low, mid_v, d, w, g, c, T);

        calculate_T_table_rectangular(left, mid_h, mid_v+1, high, d, w, g, c, T);

        calculate_T_table_triangular(mid_h+1, right, mid_v+1, high, d, w, g, c, T);

    }

}

/**
 * An internal function
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
void calculate_T_table_rectangular(int left, int right, int low, int high, int d, mpq_t * w, mpq_t ** g, mpq_t * c,  mpq_t ** T){

    // Very basic work estimate slightly taking into account the number of digits
    double input_size_estimate = sqrt(high*right);
    double input_size_bits = log2(input_size_estimate);
    double bytes_estimate = input_size_bits * log2(input_size_bits) / 8;
    double memory_use_estimate = (high-low)*(right-left)*bytes_estimate;

    if( memory_use_estimate < SERIAL_SIZE_BASE ){ // The base case

        mpq_t calc_temp; mpq_init(calc_temp);
        for(int i=high; i>=low; i--){
            for(int j=left; j<=right; j++){
      
                //mpq_init(T[i][j]);
                mpq_mul(calc_temp, g[i][d+j-1], c[d+j-1]);
                mpq_add(T[i][j], calc_temp, T[i][j-1]);
                
                if(i==j){ // The entry at T[i][i] will be the d+i th coefficient of the solution
                    mpq_set(calc_temp, T[i][j]);                    
                    mpq_div(T[i][j], calc_temp, w[i]); // divide by a[d][0] * rf[d][i]
                    mpq_init(c[d+i]); mpq_set(c[d+i], T[i][j]); // Copy the table value into the c array for easier access
                }

            }
        }
        mpq_clear(calc_temp);


    }
    else{ // Recursive Calls

        int mid_v = (low+high)/2;
        int mid_h = (left+right)/2;

        calculate_T_table_rectangular(left, mid_h, low, mid_v, d, w, g, c, T);
        calculate_T_table_rectangular(left, mid_h, mid_v+1, high, d, w, g, c, T);

        calculate_T_table_rectangular(mid_h+1, right, low, mid_v, d, w, g, c, T);
        calculate_T_table_rectangular(mid_h+1, right, mid_v+1, high, d, w, g, c, T);

    }

}



















/**
* Deallocate the memory allocated for the rising factorial table
* @param d the order of the DE for which the rising factorial table is being created
* @param k the precision (minus order) of the coefficient of the power series solution being calculated
*/
void destroy_rising_factorial_table_mpq_t (mpq_t ** rf, int d, int k){

    for (int i = 0; i <= d; i++) {

        // Clear the mpz_t structures
        for (int j = 0; j <= k; j++) {
            mpq_clear(rf[i][j]);
        }

        // Free each row array
        free (rf[i]);
    }
    
    free(rf);
    
}

/**
Create the table (2-D array) of rising factorials using simple iterative construction
*/
mpq_t ** rising_factorial_table_mpq_t (int dd, int kk){

    // Cast conversion for using GMP functions
    unsigned long d = (unsigned long) dd;
    unsigned long k = (unsigned long) kk;

    // Dynamically allocate space to 2-D table and initialize the structures
    mpq_t ** rf = (mpq_t **) malloc(sizeof (mpq_t *) * (d+1) );
    for (unsigned long i=0; i<=d; i++){
        rf[i] = (mpq_t * ) malloc (sizeof(mpq_t) * (k+1) );

        // Initialize the mpz_t structures
        for (unsigned long j = 0; j <= k; j++) {
            mpq_init(rf[i][j]);
        }
    }   

    // The 0th row of the table is all 1's
    for(unsigned long n=0; n<=k; n++){
        mpq_set_ui(rf[0][n],1,1); 
    }

    /**
     * We will do the calculations using integers, then convert to rational as a last step
    */

   // Helper for first column
   mpz_t first_col_save; mpz_init_set_ui(first_col_save,1);

   // Registers to save values, used for iteration and also to improve performance
   mpz_t prod_tmp; mpz_init(prod_tmp);
   mpz_t rf_j_n_1; mpz_init(rf_j_n_1); // Saves value rf[j][n-1]

   // Build table size of size k*d, use formulas from notes
    for (unsigned long j=1; j<=d; j++){

        mpz_mul_ui(first_col_save, first_col_save, j); // first_col_save *= j
        mpq_set_z(rf[j][0], first_col_save); // rf[j][0] = first_col_save
        mpz_set(rf_j_n_1, first_col_save); // save in a register as mpz

        for (unsigned long n=1; n<=k; n++){

            // rf[j][n] = rf[j][n-1] * (j+n) / 
            mpz_mul_ui(prod_tmp, rf_j_n_1, (j+n)); 
            mpz_tdiv_q_ui(rf_j_n_1, prod_tmp, n);
            mpq_set_z(rf[j][n], rf_j_n_1);

       }
   }

   mpz_clear(first_col_save);
   mpz_clear(prod_tmp);
   mpz_clear(rf_j_n_1);

   return rf;

}

/*
Create the table (2-D array) of rising factorials using recursive construction
*/
mpq_t ** rising_factorial_table_recursive_mpq_t (int dd, int kk){


    // Cast conversion for using GMP functions
    unsigned long d = (unsigned long) dd;
    unsigned long k = (unsigned long) kk;

    // Dynamically allocate space to 2-D table and initialize the structures
    mpq_t ** rf = (mpq_t **) malloc(sizeof (mpq_t *) * (d+1) );
    for (unsigned long i=0; i<=d; i++){
        rf[i] = (mpq_t * ) malloc (sizeof(mpq_t) * (k+1) );

        // Initialize the mpz_t structures
        for (unsigned long j = 0; j <= k; j++) {
            mpq_init(rf[i][j]);
        }
    }

    // Call recursive function to build the table
    rising_factorial_mpq_t(0,kk, rf, dd);

    return rf;

}

/*
Recursive function to call to create rising factorial table recursively dividing only along witdh (treat height as constant)
*/
void rising_factorial_mpq_t (int left_i, int right_i, mpq_t ** rf, int height_i){


    unsigned long BASE = 2048;
    unsigned long height = (unsigned long) height_i;

    // If height >= BASE it would be possible to never terminate
    if (height >= BASE){
        printf(stderr, "Table height is too high\n");
        exit(1);
    }
    // BASE *= height_i;

    // Base Case
    if ( (right_i-left_i)*height_i < (int)BASE){

        unsigned long left = (unsigned long) left_i;
        unsigned long right = (unsigned long) right_i;

        // The 0th row of the table will all be 1's
        for (unsigned long n=left; n<=right; n++){
            mpq_set_ui(rf[0][n],1,1);
        }

        // Register to save the first term from the previous row to avoid new cache misses
        mpz_t first_col_save; mpz_init_set_ui(first_col_save,1);

        // Registers to save values, used for iteration and also to improve performance
        mpz_t prod_tmp; mpz_init(prod_tmp);
        mpz_t rf_j_n_1; mpz_init(rf_j_n_1); // Saves value rf[j][n-1]

        // Calculate the table entries
        for (unsigned long j=1; j<=height; j++){

            mpz_mul_ui(first_col_save, first_col_save, j+left); // first_col_save *= j
            mpq_set_z(rf[j][left], first_col_save); // rf[j][left] = first_col_save
            mpz_set(rf_j_n_1, first_col_save); // save in a register as mpz

            for (int n=left+1; n<=right; n++){

                // rf[j][n] = rf[j][n-1] * (j+n) / n
                mpz_mul_ui(prod_tmp, rf_j_n_1, (j+n)); 
                mpz_tdiv_q_ui(rf_j_n_1, prod_tmp, n);
                mpq_set_z(rf[j][n], rf_j_n_1);

            }

        }

        mpz_clear(first_col_save);
        mpz_clear(prod_tmp);
        mpz_clear(rf_j_n_1);

    }

    // Recursive calls
    else{
        int mid = (right_i+left_i)/2;
        rising_factorial_mpq_t(left_i, mid, rf, height_i);
        rising_factorial_mpq_t(mid+1, right_i, rf, height_i);
    }

}











/**
 Function to generate the mclaurin series for a univariate exponential function
 Of the form A*e^(B*x) --> front_factor * e ^ (inner_factor * x)
 @param front_factor the front factor of the exponential function so front_factor * e^x
 @param inner_factor the inner factor of the exponential function so e^(inner_factor * x)
*/
PowerSeries_t * exponential_univariate(mpq_t* front_factor, mpq_t* inner_factor){

        PowerSeries_t * expuni = allocatePowerSeries_PS(1);

        expuni->deg = 0;

        expuni->nvar = 1;
    
        expuni->genParam1 = (void*) front_factor; 
        expuni->genParam2 = (void*) inner_factor; 

        expuni->paramType1 = PLAIN_DATA;
        expuni->paramType2 = PLAIN_DATA; 
    
        expuni->genOrder = 2;
    
        expuni->gen.binaryGen = &(homogPartVoid_expUnivariate_PS);
    
        expuni->polys[0] = makeConstPolynomial_AA(1,1,*front_factor);

        return expuni;

}


/**
 * An internal function.
 * A void generator wrapper for the generator of a power series of the exponential function A * e ^ Bx
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPartVoid_expUnivariate_PS(int d, void* param1, void* param2) {
    return homogPart_expUnivariate_PS(d, (mpq_t*) param1, (mpq_t*) param2);
}


/**
 * An internal function.
 * Computes the homogeneous part of a power series of the exponential funciton A * e ^ Bx
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPart_expUnivariate_PS(int d, mpq_t* front_factor , mpq_t* inner_factor) {

    unsigned long n = (unsigned long) d;

    mpq_t front, inner; mpq_init(front); mpq_init(inner);
    mpq_set(front, *front_factor); mpq_set(inner, *inner_factor);

    mpz_t inner_num_exp, inner_den_exp;
    mpz_init(inner_num_exp); mpz_init(inner_den_exp);

    mpz_pow_ui(inner_num_exp, mpq_numref(inner), n);
    mpz_pow_ui(inner_den_exp, mpq_denref(inner), n);

    mpq_t inner_exp; mpq_init(inner_exp);
    mpq_set_num(inner_exp, inner_num_exp);
    mpq_set_den(inner_exp, inner_den_exp);
    mpq_canonicalize(inner_exp); // inner_exp = inner_factor ^ d

    mpz_t fact; mpz_init(fact);
    mpz_fac_ui(fact, n); // d!

    mpq_t fact_mpq; mpq_init(fact_mpq);
    mpq_set_z(fact_mpq, fact);
    mpq_inv(fact_mpq, fact_mpq); // 1/d!

    mpq_t coeffic; mpq_init(coeffic);
    mpq_mul(coeffic,inner_exp,fact_mpq); 
    mpq_mul(coeffic, coeffic, front); // 1/d! * (inner^d) * front 

    Poly_ptr return_coeff = make_univariate_monomial(d,coeffic);

    mpz_clear(inner_num_exp);
    mpz_clear(inner_den_exp);
    mpq_clear(inner_exp);
    mpz_clear(fact);
    mpq_clear(fact_mpq);
    mpq_clear(coeffic);
    mpq_clear(front);
    mpq_clear(inner);

    return return_coeff;
    

}




/**
 Function to generate the mclaurin series for a univariate exponential function with x to the exponent of some natural number
 Of the form A*e^(B*(x^n)) --> front_factor * e ^ (inner_factor * x^n)
 @param front_factor the front factor of the exponential function so front_factor * e^(x^n)
 @param inner_factor the inner factor of the exponential function so e^(inner_factor * x^n)
 @param nat the exponent of x within the exponential function, so e^(inner_factor * x^n)
*/
PowerSeries_t * exponential_univariate_natural_power(mpq_t* front_factor, mpq_t* inner_factor, int nat){

        if(nat<=0){
            fprintf(stderr, "ERROR: (In exponential_univariate_natural_power) Parameter nat must be a natural number (integer >=0)\n");
            exit(1);
        }

        PowerSeries_t * expuni = allocatePowerSeries_PS(1);

        expuni->deg = 0;

        expuni->nvar = 1;
    
        expuni->genParam1 = (void*) front_factor; 
        expuni->genParam2 = (void*) inner_factor;
        expuni->genParam3 = (void*) nat; 
 

        expuni->paramType1 = PLAIN_DATA;
        expuni->paramType2 = PLAIN_DATA;
        expuni->paramType3 = PLAIN_DATA; 
 
    
        expuni->genOrder = 3;
    
        expuni->gen.tertiaryGen = &(homogPartVoid_expUnivariate_natural_power_PS);
    
        expuni->polys[0] = makeConstPolynomial_AA(1,1,*front_factor);

        return expuni;

}


/**
 * An internal function.
 * A void generator wrapper for the generator of a power series of the exponential function
 * With x to the exponent of some natural number
 *  A * e ^ (B*(x^n))
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 * @param nat the exponent of x within the exponential function
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPartVoid_expUnivariate_natural_power_PS(int d, void* param1, void* param2, void* param3) {
    return homogPart_expUnivariate_natural_power_PS(d, (mpq_t*) param1, (mpq_t*) param2, (int) param3);
}


/**
 * An internal function.
 * Computes the homogeneous part of a power series of the exponential funciton
 * With x to the exponent of some natural number
 *  A * e ^ (B*(x^n))
 * @param d: the requested degree to generate
 * @param param1 the front factor A
 * @param param2 the inner factor B
 *  @param nat the exponent of x within the exponential function
 * @return the homogeneous part of degree d
 */
Poly_ptr homogPart_expUnivariate_natural_power_PS(int d, mpq_t* front_factor , mpq_t* inner_factor, int nat) {


    if(d%nat != 0){ // if the desired degree is not divisible by n, we return nothing
        return NULL;        
    }

    unsigned long small_index = (unsigned long) d/nat;

    mpq_t inner_exp; mpq_init(inner_exp);
    mpz_pow_ui(mpq_numref(inner_exp), mpq_numref(*inner_factor), small_index);
    mpz_pow_ui(mpq_denref(inner_exp), mpq_denref(*inner_factor), small_index);
    mpq_canonicalize(inner_exp); // inner_exp = inner_factor ^ (d/nat)

    mpq_t fact_mpq; mpq_init(fact_mpq);
    mpz_set_ui(mpq_numref(fact_mpq),1);
    mpz_fac_ui(mpq_denref(fact_mpq), small_index); // 1/(d/nat)!

    mpq_t coeffic; mpq_init(coeffic);
    mpq_mul(coeffic,inner_exp,fact_mpq); 
    mpq_mul(coeffic, coeffic, *front_factor); // 1/d! * (inner^d) * front 

    Poly_ptr return_coeff = make_univariate_monomial(d,coeffic);

    mpq_clear(inner_exp);
    mpq_clear(fact_mpq);
    mpq_clear(coeffic);

    return return_coeff;
    

}











