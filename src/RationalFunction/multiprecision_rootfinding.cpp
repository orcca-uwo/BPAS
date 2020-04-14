#include "../../include/ring.h"
#include "../../include/RationalNumberPolynomial/urpolynomial.h"
#include "../../include/RingPolynomial/upolynomial.h"
#include <mps/mps.h>
#include "RationalFunction/multiprecision_rootfinding.h"


double _checkInt_d(const double val_in, const double radius){
	/* _checkInt(val_in,val_out,radius):                      */
	/* Check for an integer                                  */
	/*                                                       */
	/* If the input value is within the radius of an integer */
	/* return that integer in the output value.              */

	double val_out;
	double diff;

	diff = val_in - floor(val_in);
  	if (diff < radius) {
  		val_out = floor(val_in);
  	}
  	else {
  		diff = ceil(val_in) - val_in;
  		if (diff < radius){
  			val_out = ceil(val_in);
  		}
  		else {
  			val_out = val_in;
  		}
  	}

  	return val_out;

  	/*
    mpf_get_2dl (&val_d, &exp_d, val_in);
  	rdpe_set_2dl (val_dpe, val_d, exp_d);
  	rdpe_abs (val_mrdpe,val_mrdpe);
  	rdpe_lt(val_mrdpe,rad[j])
	*/
}

void _checkInt_m(mpf_t val_out, const mpf_t val_in, const rdpe_t radius, int prec){
	/* _checkInt(val_in,val_out,radius):                      */
	/* Check for an integer                                  */
	/*                                                       */
	/* If the input value is within the radius of an integer */
	/* return that integer in the output value.              */

	double diff_d;
	long int exp_li;
	rdpe_t diff_dpe;
	mpf_t diff;
	mpf_init2(diff,prec);

	mpf_floor(diff,val_in);
	mpf_sub(diff,val_in,diff);
    mpf_get_2dl (&diff_d, &exp_li, diff);
  	rdpe_set_2dl (diff_dpe, diff_d, exp_li);
  	if (rdpe_lt(diff_dpe,radius))
  		mpf_floor(val_out,val_in);
  	else {
  		mpf_ceil(diff,val_in);
		mpf_sub(diff,diff,val_in);
    	mpf_get_2dl (&diff_d, &exp_li, diff);
  		rdpe_set_2dl (diff_dpe, diff_d, exp_li);
  		if (rdpe_lt(diff_dpe,radius))
  			mpf_ceil(val_out,val_in);
  		else
  			mpf_set(val_out,val_in);
  	}
}

template <class UnivariatePolynomialOverField, class Field>
std::vector< std::vector<ComplexRationalNumber> > _rootsDoublePrecision(std::vector<UnivariatePolynomialOverField> &U){
	// find roots of elements of U in double precision //
	int i,j;
	int n;
	Field coeff;
	cplx_t *roots;
	double *radii;
	ComplexRationalNumber ComplexRationalNumberTemp;
	std::vector<ComplexRationalNumber> r;
	std::vector< std::vector<ComplexRationalNumber> > E;
	mpq_t zero;
	double val;

	mpq_init (zero);
	mpq_set_si (zero, 0, 1);

	mps_monomial_poly *p;
	mps_context *s;

	s = mps_context_new ();

	mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);

	mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);

  	mps_context_set_input_prec (s, 0);
	mps_context_set_output_prec (s, 53);

	for (i=0; i<U.size(); i++){

		n = U.at(i).degree().get_si();

		roots = new cplx_t[n];
		radii = new double[n];

		p = mps_monomial_poly_new (s, n);

		for (j=0; j<n+1; j++){
			coeff = U.at(i).coefficient(j);
			if (!(coeff == 0))
				mps_monomial_poly_set_coefficient_q (s, p, j, coeff.get_mpq_t(), zero);
		}

		mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));

		mps_mpsolve (s);

		mps_context_get_roots_d (s, &roots, &radii);

		for (j=0; j<n; j++)
    		{
    			val = cplx_Re(roots[j]);
    			val = _checkInt_d(val,radii[j]);
    			ComplexRationalNumberTemp.setRealPart(Field(mpq_class(val)));
    			val = cplx_Im(roots[j]);
    			val = _checkInt_d(val,radii[j]);
    			ComplexRationalNumberTemp.setImaginaryPart(Field(mpq_class(val)));
    			r.push_back(ComplexRationalNumberTemp);
    		}

		E.push_back(r);
		r.clear();

		delete [] roots;
		delete [] radii;
	}

	/*vector<DenseUnivariateRationalPolynomial> Z;
	double real,imag;
	cout << endl << "Roots in double precision:" << endl;
	for (i=0; i<E.size(); i++){
		Z = E.at(i);
		for (j=0; j<Z.size(); j++){
			real = Z.at(j).coefficient(0).get_d();
			imag = Z.at(j).coefficient(1).get_d();
			cout << Z.at(j) << " = " << real << " + " << imag << "i" << endl;
		}
	}*/

	return(E);
}

template <class UnivariatePolynomialOverField, class Field>
std::vector< std::vector<ComplexRationalNumber> > _rootsMultiprecision(std::vector<UnivariatePolynomialOverField> &U, int prec){
	/* roots(U):                                                 */
	/* Multiprecision complex roots of a polynomial              */
	/*                                                           */
	/* return the roots t=a+ib of U(t) as ComplexRationalNumbers */

	int precCorr;
	int i,j;
	int n;
	Field coeff;
	ComplexRationalNumber ComplexRationalNumberTemp;
	std::vector<ComplexRationalNumber> r;
	std::vector< std::vector<ComplexRationalNumber> > F;
	mpq_t zero;
	//Field mpqTemp;
	double val;
	//double radius;
	mpf_t val_mf;
	mpf_init2(val_mf,prec);
	rdpe_t radius_m;
	rdpe_t val_mrdpe;
	long int expo;
	mpq_t val_mq;
	mpq_init(val_mq);

	mpq_init (zero);
	mpq_set_si (zero, 0, 1);

	mps_monomial_poly *p;
	mps_context *s;

	s = mps_context_new ();

	mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA); // alternative MPS_ALGORITHM_STANDARD_MPSOLVE

	mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE); // alternatives MPS_OUTPUT_GOAL_ISOLATE, MPS_OUTPUT_GOAL_COUNT

  	mps_context_set_input_prec (s, 0);
	mps_context_set_output_prec (s, prec);

	for (i=0; i<U.size(); i++){

		n = U.at(i).degree().get_si();
		// Adjust precision to ensure that the relative error in the root is < 2^{-prec}
		//precCorr = (int) ceil(log2(n));
		//precCorr += prec;

		//cout << "Computing roots for U[" << i << "] = " << U->at(i) << endl;

		/* Allocate space to hold the results. We are requesting multiprecision
		 * roots here, which requires the following initialization scheme */
		mpc_t *rts;
		rdpe_t *rad;
		rts = mpc_valloc(n);
		mpc_vinit2(rts,n,prec);
		rad = rdpe_valloc(n);
		rdpe_vinit(rad,n);

		mpc_t mpc_zero;
		mpc_set_d (mpc_zero, 0.0, 0);

		for (int l = 0; l < n; ++l) {
			rdpe_set (rad[l], rdpe_zero);
			mpc_set (rts[l], mpc_zero);
		}

		p = mps_monomial_poly_new (s, n);

		for (j=0; j<n+1; j++){
			coeff = U.at(i).coefficient(j);
			if (!(coeff == 0))
				mps_monomial_poly_set_coefficient_q (s, p, j, coeff.get_mpq_t(), zero);
		}

		/* Set the input polynomial */
		mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));

		/* Actually solve the polynomial */
		mps_mpsolve (s);

		//cout << "data precision max: " << mps_context_get_data_prec_max(s) << endl;\
		//cout << "over max = " << mps_context_get_over_max(s) << endl;

		/* Save roots computed in the vector results */
		mps_context_get_roots_m (s, &rts, &rad);

		for (j=0; j<n; j++)
    		{
    			/* Get the radius of the circle containing the root for zero check */
    			//rdpe_set(radius_m,rad[j]);
    			//radius = rdpe_get_d(radius_m);

    			/* Get real part as mpf_t and as double for zero check */
    			mpf_set(val_mf,mpc_Re(rts[j]));
    			//gmp_printf ("real part mpf %.*Ff with %d digits\n", prec, val_mf, prec);
    			_checkInt_m(val_mf,val_mf,rad[j],prec);

    			/* convert from mpf_t -> mpq_t -> mpq_class (using constructor) */
    			mpq_set_f(val_mq,val_mf);
    			mpq_class val_mxxq_r(val_mq);
    			ComplexRationalNumberTemp.setRealPart(val_mxxq_r);

    			/* Get imaginary part as mpf_t and as double for zero check */
    			mpf_set(val_mf,mpc_Im(rts[j]));
    			//gmp_printf ("imaginary part mpf %.*Ff with %d digits\n", prec, val_mf, prec);
    			_checkInt_m(val_mf,val_mf,rad[j],prec);

    			/* convert from mpf_t -> mpq_t -> mpq_class (using constructor) */
    			mpq_set_f(val_mq,val_mf);
    			mpq_class val_mxxq_i(val_mq);
    			ComplexRationalNumberTemp.setImaginaryPart(val_mxxq_i);

    			r.push_back(ComplexRationalNumberTemp);
    		}

    		F.push_back(r);
    		r.clear();

		free(rts);
	}

	mpf_clear(val_mf);
	mpq_clear(val_mq);

	//cout << endl << "Roots in multiprecision:" << endl;
	/*std::vector<DenseUnivariateRationalPolynomial> Z;
	double real,imag;
	for (i=0; i<F.size(); i++){
		Z = F.at(i);
		for (j=0; j<Z.size(); j++){
			real = Z.at(j).coefficient(0).get_d();
			imag = Z.at(j).coefficient(1).get_d();
			//cout << Z.at(j) << " = " << real << " + " << imag << "i" << endl;
		}
	}*/

	return(F);
}

// to avoid linking errors
template std::vector< std::vector<ComplexRationalNumber> > _rootsMultiprecision<DenseUnivariateRationalPolynomial,RationalNumber>(std::vector<DenseUnivariateRationalPolynomial> &U, int prec);
template std::vector< std::vector<ComplexRationalNumber> > _rootsMultiprecision<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(std::vector< SparseUnivariatePolynomial<RationalNumber> > &U, int prec);
// WILL NEED A SWITCH OR OTHER CONTROL FOR COMPLEX COEFFICIENTS
//template std::vector< std::vector<ComplexRationalNumber> > _rootsMultiprecision<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *U, int prec);

template std::vector< std::vector<ComplexRationalNumber> > _rootsDoublePrecision<DenseUnivariateRationalPolynomial,RationalNumber>(std::vector<DenseUnivariateRationalPolynomial> &U);
template std::vector< std::vector<ComplexRationalNumber> > _rootsDoublePrecision<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(std::vector< SparseUnivariatePolynomial<RationalNumber> > &U);
// WILL NEED A SWITCH OR OTHER CONTROL FOR COMPLEX COEFFICIENTS
//template std::vector< std::vector<ComplexRationalNumber> > _rootsDoublePrecision<SparseUnivariatePolynomial<ComplexRationalNumber>,ComplexRationalNumber>(std::vector< SparseUnivariatePolynomial<ComplexRationalNumber> > *U);

