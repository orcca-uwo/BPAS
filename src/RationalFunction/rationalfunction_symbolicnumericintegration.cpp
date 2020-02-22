//#include "bpas.h"
#include <mps/mps.h>
#include <mps/mpc.h>
#include <complex>
#include <algorithm>
#include "../../include/polynomial.h"
#include "../../include/ring.h"
#include "../../include/RationalNumberPolynomial/urpolynomial.h"
#include "../../include/RingPolynomial/upolynomial.h"
//#include "../../include/RationalFunction/rationalfunction_euclideanmethods.h"
#include "../../include/RationalFunction/urationalfunction.h"
#include "../../include/RationalFunction/rationalfunction_symbolicintegration.h"
#include "../../include/RationalFunction/multiprecision_rootfinding.h"
#include "../../include/RationalFunction/rationalfunction_integrationpostprocessing.h"
#include "../../include/RationalFunction/rationalfunction_symbolicnumericintegration.h"
#include "../../include/RationalFunction/rationalfunction_complexrationalnumberordering.h"
#include "../../include/RationalFunction/complexdouble.h"
#include "../../include/RationalFunction/unumericalpolynomial.h"

template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _computePartialDerivatives(std::vector<SparseUnivariatePolynomial<UnivariatePolynomialOverRealField>>& S_in,std::vector<UnivariatePolynomialOverRealField>& U_in,std::vector<SparseUnivariatePolynomial<UnivariatePolynomialOverRealField>> *dSdx_out,std::vector<SparseUnivariatePolynomial<UnivariatePolynomialOverRealField>> *dSdt_out,std::vector<SparseUnivariatePolynomial<UnivariatePolynomialOverRealField>> *d2Sdxdt_out,std::vector<UnivariatePolynomialOverRealField> *dU_out, float *ea_total, std::ostringstream *outStr,bool PROFILING){

	UnivariatePolynomialOverRealField uporfTemp;
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> supTemp;
	std::vector<SparseUnivariatePolynomial<UnivariatePolynomialOverRealField>> S(S_in),dSdx,dSdt,d2Sdxdt;
	std::vector<UnivariatePolynomialOverRealField> U(U_in),dU;
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING)
		startTimer(&start);
	
	for (int i=0; i < U.size(); i++){
		uporfTemp = U.at(i);
		uporfTemp.differentiate(1);
		dU.push_back(uporfTemp);
		supTemp = S.at(i);
//		std::cout << "S = " << supTemp << std::endl;
		supTemp.differentiate(1);
		dSdx.push_back(supTemp);
//		std::cout << "dSdx = " << supTemp << std::endl;
		supTemp.zero();
		for (int j=0; j<S.at(i).degree()+1; j++){
			uporfTemp = S.at(i).coefficient(j);
			if (!uporfTemp.isZero()){
				uporfTemp.differentiate(1);
				supTemp.setCoefficient(j,uporfTemp);
			}
		}
		dSdt.push_back(supTemp);
//		std::cout << "dSdt = " << supTemp << std::endl;
		supTemp.differentiate(1);
		d2Sdxdt.push_back(supTemp);
//		std::cout << "d2Sdxdt = " << supTemp << std::endl;
	}
	
	*dSdx_out = dSdx;
	*dSdt_out = dSdt;
	*d2Sdxdt_out = d2Sdxdt;
	*dU_out = dU;
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*ea_total += elapsed;
		*outStr << "\t\tPartial derivatives\t" << elapsed << "\n";
		startTimer(&start);
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_rootfinding(UnivariatePolynomialOverRealField &H, int prec, std::vector<ComplexRationalNumber> *denomRoots, float *ea_total, std::ostringstream *outStr,bool PROFILING,bool PFD){
	std::vector<UnivariatePolynomialOverRealField> U2;
	std::vector< std::vector<ComplexRationalNumber> > roots;
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		*outStr << "------------------------------\n";
		*outStr << "\tError Analysis\n";
		*outStr << "\t------------------------------\n";
		startTimer(&start);
	}
			
	U2.push_back(H);
	roots = _rootsMultiprecision<UnivariatePolynomialOverRealField,RealField>(U2,prec);
	*denomRoots = roots.at(0);

	if (PROFILING){
		stopTimer(&start,&elapsed);
		if (PFD)
			std::cout << "\t\tdenom Rootfinding\t" << elapsed << std::endl;
		else {
			*ea_total += elapsed;
			*outStr << "\t\tdenom Rootfinding\t" << elapsed << std::endl;
		}
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_computeReside(UnivariatePolynomialOverRealField &G, UnivariatePolynomialOverRealField &dH, ComplexRationalNumber root, ComplexRationalNumber *residue){
	*residue = G.template evaluate<ComplexRationalNumber>(root);
	*residue /= dH.template evaluate<ComplexRationalNumber>(root);
}

bool epsilonEqual(ComplexRationalNumber a, ComplexRationalNumber b, int prec) {
	if (ComplexRationalNumberOrdering::epsilonEqual(a.realPart(),b.realPart(),prec) && ComplexRationalNumberOrdering::epsilonEqual(a.imaginaryPart(),b.imaginaryPart(),prec))
		return true;
	else
		return false;
}

bool epsilonEqual(RationalNumber a, RationalNumber b, int prec) {
	if (ComplexRationalNumberOrdering::epsilonEqual(a,b,prec))
		return true;
	else
		return false;
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_rootsAndResidues(UnivariatePolynomialOverRealField &G, UnivariatePolynomialOverRealField &H, int prec, UnivariatePolynomialOverRealField *dH, std::vector<ComplexRationalNumber> *denomRoots, std::vector<ComplexRationalNumber> *residues, residueRootIndexMap *rrim, float *ea_total, std::ostringstream *outStr,bool PROFILING,bool PFD){
	std::map< int, std::vector<int> > bundle;
	UnivariatePolynomialOverRealField ddH;
	UnivariatePolynomialOverRealField dG;
	ComplexRationalNumber residue;
	int localPrec(prec); // it increases the cost significantly to increase from the default prec
	int i(0), j(0);
	bool duplicate;
	//std::map<int, ComplexRationalNumber> rootIndexResidueMap;
	residueRootIndexMap resRoot;
	//std::vector<ComplexRationalNumber> res;
		
	// Profiling variables
	unsigned long long start;
	float elapsed(0);

	// This routine will need to check for close roots, and recompute at higher precision if necessary.
	// This will only happen with large coefficients or high degree, so there is a way to control for likelihood.
	// From experiments, it appears that the error in the redidue is of order of the difference
	//   in precision of the rootfinding and the distance between nearby roots, so this provides a
	//   means of computing the updated precision to ensure the residue error is less than eps.
	// Given that the error in the roots is typically MANY orders of magnitude smaller than eps, this does
	//   not appear to be a practical issue at this point.
	
	_errorAnalysisLRT_rootfinding<UnivariatePolynomialOverRealField,RealField>(H,localPrec,denomRoots,ea_total,outStr,PROFILING,PFD);
	
	if (PROFILING){
		startTimer(&start);
	}
	
	// sort the roots of H(x)
	sort(denomRoots->begin(),denomRoots->end(),ComplexRationalNumberOrdering(localPrec));
	
	// remove conjugate with positive imaginary part, gives us residues with positive imaginary part.
	for (i=denomRoots->size()-1; i>0; --i) {
		if (denomRoots->at(i).imaginaryPart() != 0) {
			denomRoots->erase(denomRoots->begin()+i-1);
			--i;
		}
	}

	if (PROFILING){
		stopTimer(&start,&elapsed);
		if (PFD)
			std::cout << "\t\tsort and edit roots\t" << elapsed << std::endl;
		else {
			*ea_total += elapsed;
			*outStr << "\t\tsort and edit roots\t" << elapsed << std::endl;
		}
		startTimer(&start);
	}
	
	*dH = H;
	dH->differentiate(1);
	
	// fill array of roots
	int size = denomRoots->size();
	//ComplexRationalNumber* rootArray = new ComplexRationalNumber[size];
	//for (i=0; i<size; ++i) {
	//	rootArray[i] = denomRoots->at(i);
	//}
	
	ComplexRationalNumber* residueArray = new ComplexRationalNumber[size];

	// compute residues accurately in parallel
	cilk_for (int k=0; k<size; ++k) {
	//for (int k=0; k<size; ++k) {
		ComplexRationalNumber c,root;
		root = denomRoots->at(k);
		c = G.template evaluate<ComplexRationalNumber>(root);
		c /= dH->template evaluate<ComplexRationalNumber>(root);
		residueArray[k] = c;
	}
	
	ComplexRationalNumber a;
	a.zero();
	for (i=0; i<size; i++) {
		if (epsilonEqual(residueArray[i],a,prec))
			denomRoots->erase(denomRoots->begin()+i);
		else
			residues->push_back(residueArray[i]);
	}
	//for (i=0; i<residues->size(); ++i)
	//	std::cout << "residues[" << i << "] = (" << residues->at(i).realPart().get_d() << "," << residues->at(i).imaginaryPart().get_d() << ")" << std::endl;
	
	delete [] residueArray;
	//delete [] rootArray;

	if (PROFILING){
		stopTimer(&start,&elapsed);
		*ea_total += elapsed;
		if (PFD)
			std::cout << "\t\tcompute residues\t" << elapsed << std::endl;
		else {
			*ea_total += elapsed;
			*outStr << "\t\tcompute residues\t" << elapsed << std::endl;
		}
		startTimer(&start);
	}
	
	for (i=0; i<residues->size(); ++i) {
		rrim->insert(std::make_pair(residues->at(i),i));
	}

	if (PROFILING){
		stopTimer(&start,&elapsed);
		*ea_total += elapsed;
		if (PFD)
			std::cout << "\t\tresidue-root map\t" << elapsed << std::endl;
		else {
			*ea_total += elapsed;
			*outStr << "\t\tresidue-root map\t" << elapsed << std::endl;
		}
		startTimer(&start);
	}
	
	
	
	/*if (PROFILING){
		startTimer(&start);
	}
	
	complexMPF* resArray = new complexMPF[size];
	mpq_t zero;
	mpq_init(zero);
	
	SparseUnivariateMPComplexPolynomial nG(G,53);
	SparseUnivariateMPComplexPolynomial nH(H,53);
	SparseUnivariateMPComplexPolynomial ndH(*dH,53);
	// compute residues accurately in parallel
	cilk_for (int m=0; m<size; ++m) {
		complexMPF c,root;
		mpc_init2(c.c,4*53);
		mpc_init2(root.c,4*53);
		mpc_set_q(root.c,denomRoots->at(m).realPart().get_mpq_t(),denomRoots->at(m).imaginaryPart().get_mpq_t());
		//c = nG.evaluate(root,53);
		c = nH.evaluate(root,53);
		//mpc_div(c.c,c.c,ndH.evaluate(root,53).c);
		resArray[m] = c;
	}
	
	
		
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*ea_total += elapsed;
		if (PFD)
			std::cout << "\t\tnumerical residues\t" << elapsed << std::endl;
		else {
			*ea_total += elapsed;
			*outStr << "\t\tnumerical residues\t" << elapsed << std::endl;
		}
	}
	
	for (int k=0; k<size; ++k) {
		*outStr << "resArray[" << k << "] = (" << mpc_Re(resArray[k].c) << "," << mpc_Im(resArray[k].c) << ")" << std::endl;
	}
	delete [] resArray;*/
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_computeNumericalPolynomial(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &A, std::vector< std::vector< std::complex<double> > > &roots, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > >  *A_a) {

	std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > sunpVec;
	
	for (int i=0; i<roots.size(); i++) {
		sunpVec.clear();
		for (int j=0; j<roots.at(i).size(); j++) {
			SparseUnivariateDoublePolynomial< std::complex<double> > sunp(A.at(i),roots.at(i).at(j));
			sunpVec.push_back(sunp);
		}
		A_a->push_back(sunpVec);
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_computeNumericalPolynomial(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &A, std::vector< std::vector<ComplexRationalNumber> > &roots, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > >  *A_a) {

	std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > sunpVec;
	
	for (int i=0; i<roots.size(); i++) {
		sunpVec.clear();
		for (int j=0; j<roots.at(i).size(); j++) {
			std::complex<double> cd(roots.at(i).at(j).realPart().get_d(),roots.at(i).at(j).imaginaryPart().get_d());
			SparseUnivariateDoublePolynomial< std::complex<double> > sunp(A.at(i),cd);
			sunpVec.push_back(sunp);
		}
		A_a->push_back(sunpVec);
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_computeNumericalPolynomial(std::vector<UnivariatePolynomialOverRealField> &A, std::vector< std::vector<ComplexRationalNumber> > &roots, std::vector< std::vector< std::complex<double> > >  *A_a) {

	std::vector< std::complex<double> > cdVec;
	std::complex<double> cdTemp;
	
	for (int i=0; i<roots.size(); i++) {
		SparseUnivariateDoublePolynomial< std::complex<double> > sunp(A.at(i));
		cdVec.clear();
		for (int j=0; j<roots.at(i).size(); j++) {
			std::complex<double> cd(roots.at(i).at(j).realPart().get_d(),roots.at(i).at(j).imaginaryPart().get_d());
			cdTemp = sunp.evaluate(cd);
			cdVec.push_back(cdTemp);
		}
		A_a->push_back(cdVec);
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_computeNumericalPolynomials(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &Stx, std::vector<UnivariatePolynomialOverRealField> &U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdx, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdt, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &d2Sdxdt, std::vector<UnivariatePolynomialOverRealField> &dU, std::vector< std::vector<ComplexRationalNumber> > &resultantRoots, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > >  *S_t, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > > *dSdx_t, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > > *dSdt_t, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > > *d2Sdxdt_t, std::vector< std::vector< std::complex<double> > > *t, std::vector< std::vector< std::complex<double> > > *dt, std::ostringstream *outStr, bool PROFILING){
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}

	// this computation can be split up and done in parallel, and maybe with exact computation
	// to a UPoF before conversion to SUNP
	std::vector< std::complex<double> > tVec,dtVec;
	
	int size(0),startIndex(0);
	for (int i=0; i<resultantRoots.size(); ++i) {
		size += resultantRoots.at(i).size();
	}
	
	std::complex<double>* tArray = new std::complex<double>[size];
	std::complex<double>* dtArray = new std::complex<double>[size];
	
	cilk_for (int l=0; l<size; ++l) {
	//for (int l=0; l<size; ++l) {
		int i=0;
		int j=0;
		int k=0;
		while (k != l) {
			if (l-k < resultantRoots.at(i).size()){
				j = l-k;
				k = l;
			}
			else {
				k += resultantRoots.at(i).size();
				i++;
			}
		}
		ComplexRationalNumber crn(resultantRoots.at(i).at(j));
		std::complex<double> cd(crn.realPart().get_d(),crn.imaginaryPart().get_d());
		tArray[l] = cd;
		ComplexRationalNumber crn2(-U.at(i).template evaluate<ComplexRationalNumber>(crn));
		crn2 /= dU.at(i).template evaluate<ComplexRationalNumber>(crn);
		cd = std::complex<double>(crn2.realPart().get_d(),crn2.imaginaryPart().get_d());
		dtArray[l] = cd;
	}
	startIndex = 0;
	for (int i=0; i<resultantRoots.size(); ++i) {
		for (int j=0; j<resultantRoots.at(i).size(); ++j) {
			tVec.push_back(tArray[startIndex+j]);
			dtVec.push_back(dtArray[startIndex+j]);
		}
		t->push_back(tVec);
		dt->push_back(dtVec);
		startIndex += resultantRoots.at(i).size();
	}
	_errorAnalysisLRT_computeNumericalPolynomial<UnivariatePolynomialOverRealField,RealField>(Stx,*t,S_t);
	_errorAnalysisLRT_computeNumericalPolynomial<UnivariatePolynomialOverRealField,RealField>(dSdx,*t,dSdx_t);
	_errorAnalysisLRT_computeNumericalPolynomial<UnivariatePolynomialOverRealField,RealField>(dSdt,*t,dSdt_t);
	_errorAnalysisLRT_computeNumericalPolynomial<UnivariatePolynomialOverRealField,RealField>(d2Sdxdt,*t,d2Sdxdt_t);
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*outStr << "\t\t\tNumerical polys\t" << elapsed << "\n";
	}	
	
	/*for (int i=0; i<t->size(); i++) {
		for (int j=0; j<t->at(i).size(); j++) {
			*outStr << "t[" << i << "][" << j << "] = " << t->at(i).at(j) << " +/ " << dt->at(i).at(j) << std::endl;
		}
	}*/
}

void _errorAnalysisLRT_computeIntegralTermRootMap(std::vector<ComplexRationalNumber>& denomRoots, std::vector< std::vector<ComplexRationalNumber> > &resultantRoots, int prec, residueRootIndexMap *rrim, std::vector<ComplexRationalNumber> *residues, std::vector<int> *realRootIndices, std::vector<int> *complexRootIndices, IntegralTermRootMap* itrm_real, IntegralTermRootMap* itrm_complex, std::ostringstream *outStr, bool PROFILING){
	// This function establishes a correspondence between the terms of the integral, expressed as an IntegralTerm
	// object, and the roots of H(x).
	
	int i,j,k,l;
	RootIntegralTermMap residueITM,rootITM;
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}
	
	realRootIndices->clear();
	complexRootIndices->clear();
	
	// find indices of real and complex roots of H(x)	
	for (i=0; i<denomRoots.size(); i++) {
		if (denomRoots.at(i).imaginaryPart() == 0)
			realRootIndices->push_back(i);
		else
			complexRootIndices->push_back(i);
	}
	
	IntegralTermRootMap itrm2;
	
	// temporary recomputing of residues vector
	for (i=0; i<residues->size(); i++) {
		if (residues->at(i).imaginaryPart() < 0)
			residues->at(i) = residues->at(i).conjugate();
	}
	
	// create arrays for parallel computation
	int size(0),startIndex(0),r2size(residues->size());
	for (i=0; i<resultantRoots.size(); ++i) {
		size += resultantRoots.at(i).size();
	}
	std::vector< std::pair<IntegralTerm,int> >* pii = new std::vector< std::pair<IntegralTerm,int> >[size];
	ComplexRationalNumber* resroot = new ComplexRationalNumber[r2size];
	std::pair<ComplexRationalNumber,IntegralTerm>* ci = new std::pair<ComplexRationalNumber,IntegralTerm>[size];
	
	// compute correspondence between computed residues and term of integral
	for (i=0; i<resultantRoots.size(); ++i) {
		for (j=0; j<resultantRoots.at(i).size(); j++) {
			IntegralTerm a;
			a.Sindex = i;
			a.tindex = j;
			std::pair<ComplexRationalNumber,IntegralTerm> b(std::make_pair(resultantRoots.at(i).at(j),a));
			ci[startIndex+j] = b;
		}
		startIndex += resultantRoots.at(i).size();
	}
	
	// test for equality of computed residues and roots of resultant and compute integral term-root correspondence
	cilk_for (int l=0; l<size; ++l) {
	//for (int l=0; l<size; ++l) {
		ComplexRationalNumber crn(ci[l].first);
		IntegralTerm iterm(ci[l].second);
		std::vector< std::pair<IntegralTerm,int> > p;
		for (int k=0; k<r2size; ++k) {
			ComplexRationalNumber crn2(residues->at(k));
			if (epsilonEqual(crn,crn2,prec)) {
				p.push_back(std::make_pair(iterm,k));
			}
		}
		pii[l] = p;
	}
	
	// fill IntegralTermRootMap object
	for (j=0; j<size; ++j) {
		for (k=0; k<pii[j].size(); ++k)
			itrm2.insert(pii[j].at(k));
	}
	
	delete [] pii;
	delete [] ci;
	
	IntegralTermRootMap realITRM,complexITRM;
	
	ITRMIter it2;
	for (it2 = itrm2.begin(); it2 != itrm2.end(); ++it2) {
		if (std::binary_search(realRootIndices->begin(), realRootIndices->end(), (*it2).second))
			realITRM.insert(*it2);
		else
			complexITRM.insert(*it2);
	}
	/*for (it2 = rootITRM.begin(); it2 != rootITRM.end(); ++it2) {
		*outStr << "itrm" << (*it2).first << " = " << (*it2).second << std::endl;
	}
	for (it2 = realITRM.begin(); it2 != realITRM.end(); ++it2) {
		*outStr << "realITRM" << (*it2).first << " = " << (*it2).second << std::endl;
	}
	for (it2 = complexITRM.begin(); it2 != complexITRM.end(); ++it2) {
		*outStr << "complexITRM" << (*it2).first << " = " << (*it2).second << std::endl;
	}*/
	
	*itrm_real = realITRM;
	*itrm_complex = complexITRM;
	
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*outStr << "\t\t\tterm-root map\t" << elapsed << "\n";
	}	
	
}

template <class UnivariatePolynomialOverRealField, class RealField>
SparseUnivariatePolynomial<ComplexRationalNumber> _evaluateCoefficients(SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> &A, ComplexRationalNumber c) {
	SparseUnivariatePolynomial<ComplexRationalNumber> out;
	ComplexRationalNumber crnTemp;
	
	for (int i=0; i<=A.degree(); ++i) {
		if (!A.coefficient(i).isZero()) {
			crnTemp = A.coefficient(i).template evaluate(c);
			out.setCoefficient(i,crnTemp);
		}
	}	
	return out;
}

template <class UnivariatePolynomialOverRealField, class RealField>
SparseUnivariatePolynomial<RationalNumber> _evaluateCoefficients(SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> &A, RationalNumber c) {
	SparseUnivariatePolynomial<RationalNumber> out;
	RationalNumber rnTemp;
	
	for (int i=0; i<=A.degree(); ++i) {
		if (!A.coefficient(i).isZero()) {
			rnTemp = A.coefficient(i).evaluate(c);
			out.setCoefficient(i,rnTemp);
		}
	}	
	return out;
}

void _computeBackwardErrorTerm(double S_tx, double dSdx_tx, double dSdt_tx, double d2Sdxdt_tx, double t, double dt, double *bwdErrTerm) {
	// Backward error term for LRT at real roots
	*bwdErrTerm = (dSdx_tx + (t*(d2Sdxdt_tx - ((dSdx_tx*dSdt_tx)/S_tx))))/S_tx;
	*bwdErrTerm *= dt;
}

void _computeBackwardErrorTerm(std::complex<double> S_tx, std::complex<double> dSdx_tx, std::complex<double> dSdt_tx, std::complex<double> d2Sdxdt_tx, std::complex<double> t, std::complex<double> dt, std::complex<double> *bwdErrTerm) {
	// Backward error term for LRT at complex roots
	*bwdErrTerm = (dSdx_tx + (t*(d2Sdxdt_tx - ((dSdx_tx*dSdt_tx)/S_tx))))/S_tx;
	*bwdErrTerm *= dt;
}

void _computeForwardErrorTerm(double S_tx, double dSdt_tx, double t, double dt, double *fwdErrTerm) {
	// Forward error term for LRT at real roots
	*fwdErrTerm = (t*(dSdt_tx/S_tx) + log(S_tx));
	*fwdErrTerm *= dt;
}

void _computeForwardErrorTerm(std::complex<double> S_tx, std::complex<double> dSdt_tx, std::complex<double> t, std::complex<double> dt, double *fwdErrTerm) {
	// Forward error term for LRT at complex roots
	double a,b,da,db,dErrTerm(0);
	a = t.real();
	da = dt.real();
	b = t.imag();
	db = dt.imag();
	*fwdErrTerm = 0.0;
	if (b == 0.0) {
		dErrTerm = (a*(dSdt_tx.real()/S_tx.real()) + log(S_tx.real()))*da;
		*fwdErrTerm += dErrTerm;
	}
	else {
		double A(S_tx.real());
		double B(S_tx.imag());
		double dA(dSdt_tx.real());
		double dB(dSdt_tx.imag());
		double C(A*A + B*B);
		if (a != 0.0) {
			dErrTerm = (a*(2.0*(A*dA + B*dB))/C + log(C))*da;
			dErrTerm += (a*(2.0*(B*dA - A*dB))/C)*db;
			dErrTerm *= 2.0; // since real and imaginary t give same value
			*fwdErrTerm += dErrTerm;
		}
		dErrTerm = -4.0*b*(A*(dB*da + dA*db))/C;
		*fwdErrTerm += dErrTerm;			
	}
}

void _computeBackwardErrorTerm(std::complex<double> x, std::complex<double> dx, std::complex<double> r, std::complex<double> dr, double xVal, std::complex<double> *bwdErrTerm) {
	// Backward error term for PFD at complex roots
	std::complex<double> xMinusGamma(xVal-x);
	*bwdErrTerm = ((dr + (r/xMinusGamma))/xMinusGamma)*dx;
}

void _computeForwardErrorTerm(std::complex<double> x, std::complex<double> dx, std::complex<double> r, std::complex<double> dr, double xVal, double *fwdErrTerm) {
	// Forward error term for PFD at complex roots
	double a,b,da,db,alpha,dalpha,beta,dbeta,logArg,alphaMinusX,dErrTerm(0);
	a = r.real();
	da = dr.real(); // da/dalpha = db/dbeta
	b = r.imag();
	db = dr.imag(); // db/dalpha = -da/dbeta
	alpha = x.real();
	dalpha = dx.real();
	beta = x.imag();
	dbeta = dx.imag();
	logArg = xVal*xVal - 2*alpha*xVal + alpha*alpha + beta*beta;
	alphaMinusX = alpha - xVal;
	*fwdErrTerm = 0.0;
	if (beta == 0.0) {
		dErrTerm = (da*log(abs(alphaMinusX)) + a/(alphaMinusX))*dalpha;
		dErrTerm -= (db*log(abs(alphaMinusX)))*dbeta;
		*fwdErrTerm += dErrTerm;
	}
	else {
		if (alpha != 0.0) {
			double denom = alphaMinusX*alphaMinusX + beta*beta;
			dErrTerm = 2*(db*atan(alphaMinusX/beta) + b*beta/denom)*dalpha;
			dErrTerm += 2*(da*atan(alphaMinusX/beta) - b*alphaMinusX/denom)*dbeta;
			*fwdErrTerm += dErrTerm;
		}
		dErrTerm = (da*log(logArg) + 2*a*alphaMinusX/logArg)*dalpha;
		dErrTerm += (-db*log(logArg) + 2*a*beta/logArg)*dbeta;
		*fwdErrTerm += dErrTerm;			
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_realRoots(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &Stx, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdx, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdt, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &d2Sdxdt, std::vector< std::vector< std::complex<double> > > &dt, std::vector<ComplexRationalNumber> &denomRoots, std::vector< std::vector<ComplexRationalNumber> > &resultantRoots, IntegralTermRootMap &itrm, int prec, std::vector<int> *realSingIndices, std::vector<RationalNumber> *BEIntervalWidths, std::vector<RationalNumber> *FEIntervalWidths, std::vector<double> *BEBoundaryError, std::vector<double> *FEBoundaryError, std::ostringstream *outStr, bool PROFILING) {
	RationalNumber x,dx,BEdx,FEdx,y,rnErrTerm;
	double S_tx,dSdx_tx,dSdt_tx,d2Sdxdt_tx,t,deltat,dBwdErrTerm,dFwdErrTerm;
	ITRMIter it;
	int iTemp,iTemp2,corr;
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}
	
	SparseUnivariatePolynomial<RationalNumber> uS,udSdx,udSdt,ud2Sdxdt;
	it = itrm.begin();
	int i(0),j(0),icurr(0),l(0);
	bool termChange(true);
	while(it !=itrm.end()) {
		i = (*it).first.Sindex;
		j = (*it).first.tindex;
		l = (*it).second;
		if (i > icurr)
			termChange = true;
		icurr = i;
		if (termChange) {
			RationalNumber c(resultantRoots.at(i).at(j).realPart());
			uS = _evaluateCoefficients<UnivariatePolynomialOverRealField,RealField>(Stx.at(i),c);
			udSdx = _evaluateCoefficients<UnivariatePolynomialOverRealField,RealField>(dSdx.at(i),c);
			udSdt = _evaluateCoefficients<UnivariatePolynomialOverRealField,RealField>(dSdt.at(i),c);
			ud2Sdxdt = _evaluateCoefficients<UnivariatePolynomialOverRealField,RealField>(d2Sdxdt.at(i),c);
		}
		t = resultantRoots.at(i).at(j).realPart().get_d();
		deltat = dt.at(i).at(j).real();
		for (int k=0; k<3; k++) {
			if (k == 0) {
				x = denomRoots.at(l).realPart();
				dx = RationalNumber(1);
				dx /= pow(2,prec+10);
				BEdx = dx;
				FEdx = dx;
			}
			else if (k == 1) { // correct backward error interval witdth (quadratic error variation)
				if (dBwdErrTerm < 0)
					dBwdErrTerm = -dBwdErrTerm;
				if (dBwdErrTerm > pow(2.0,-prec)) {
					dx *= pow(2,log2(dBwdErrTerm*pow(2,prec))/2);
					BEdx = dx;
				}
				else {
					dx /= pow(2,-log2(dBwdErrTerm*pow(2,prec))/2);
					BEdx = dx;
				}
			}
			else if (k == 2) { // correct forward error interval width (linear error variation)
				if (dFwdErrTerm < 0)
					dFwdErrTerm = -dFwdErrTerm;
				if (dFwdErrTerm > pow(2.0,-prec)) {
					dx = FEdx;
					dx *= pow(2,log2(dFwdErrTerm*pow(2,prec)));
				}
				else {
					dx = FEdx;
					dx /= pow(2,-log2(dFwdErrTerm*pow(2,prec)));
				}
			}
			//else
			//	dx = -dx;
			//*outStr << "dx = " << dx.get_d() << std::endl;
			y = x + dx;
			RationalNumber cuS(uS.evaluate(y));
			RationalNumber cudSdx(udSdx.evaluate(y));
			RationalNumber cudSdt(udSdt.evaluate(y));
			RationalNumber cud2Sdxdt(ud2Sdxdt.evaluate(y));
			S_tx = cuS.get_d();
			dSdx_tx = cudSdx.get_d();
			dSdt_tx = cudSdt.get_d();
			d2Sdxdt_tx = cud2Sdxdt.get_d();
			if (k != 2) {
				_computeBackwardErrorTerm(S_tx,dSdx_tx,dSdt_tx,d2Sdxdt_tx,t,deltat,&dBwdErrTerm);
			}
			if (k != 1) {
				if (S_tx < 0)
					S_tx = -S_tx;
				_computeForwardErrorTerm(S_tx,dSdt_tx,t,deltat,&dFwdErrTerm);
			}
			if (dBwdErrTerm == 0.0)
				k = 3;
		}
		realSingIndices->push_back(l);
		BEIntervalWidths->push_back(BEdx);
		FEIntervalWidths->push_back(dx);
		if (dBwdErrTerm < 0)
			dBwdErrTerm = -dBwdErrTerm;
		if (dFwdErrTerm < 0)
			dFwdErrTerm = -dFwdErrTerm;
		BEBoundaryError->push_back(dBwdErrTerm);
		FEBoundaryError->push_back(dFwdErrTerm);
		++it;
	}
	
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*outStr << "\t\t\tReal roots\t" << elapsed << "\n";
	}	
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT_complexRoots(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &Stx, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdx, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdt, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &d2Sdxdt, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > >  &S_t, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > > &dSdx_t, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > > &dSdt_t, std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > > &d2Sdxdt_t, std::vector< std::vector< std::complex<double> > > &t, std::vector< std::vector< std::complex<double> > > &dt, std::vector<ComplexRationalNumber> &denomRoots, std::vector< std::vector<ComplexRationalNumber> > &resultantRoots, IntegralTermRootMap &itrm_real, IntegralTermRootMap &itrm_complex, std::vector<int> &complexRootIndices, double *maxBwdErrVal, double *maxFwdErrVal, std::ostringstream *outStr, bool PROFILING) {
	ComplexRationalNumber crnTemp;
	std::complex<double> S_t_x,dSdx_t_x,dSdt_t_x,d2Sdxdt_t_x,tt,dtt,cdErrTerm;
	double a,b,x,da,db,A,B,C,dA,dB,bwdErr,fwdErr,dErrTerm;
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}
	
	*maxBwdErrVal = 0.0;
	*maxFwdErrVal = 0.0;
	
	for (int i=0; i<complexRootIndices.size(); i++) {
		crnTemp = denomRoots.at(complexRootIndices.at(i));
		x = crnTemp.realPart().get_d();
		for (int j=0; j<t.size(); j++) {
			for (int k=0; k<t.at(j).size(); k++) {
				S_t_x = S_t.at(j).at(k).evaluate(x);
				dSdx_t_x = dSdx_t.at(j).at(k).evaluate(x);
				dSdt_t_x = dSdt_t.at(j).at(k).evaluate(x);
				d2Sdxdt_t_x = d2Sdxdt_t.at(j).at(k).evaluate(x);
				tt = t.at(j).at(k);
				dtt = dt.at(j).at(k);
				// backward error computation
				_computeBackwardErrorTerm(S_t_x,dSdx_t_x,dSdt_t_x,d2Sdxdt_t_x,tt,dtt,&cdErrTerm);
				if (tt.imag() != 0.0)
					cdErrTerm *= 2.0;
				bwdErr += cdErrTerm.real();
				// forward error computation
				_computeForwardErrorTerm(S_t_x,dSdt_t_x,tt,dtt,&fwdErr);
			}
		}
		if (bwdErr < 0.0)
			bwdErr = -bwdErr;
		if (bwdErr > *maxBwdErrVal)
			*maxBwdErrVal = bwdErr;
		if (fwdErr < 0.0)
			fwdErr = -fwdErr;
		if (fwdErr > *maxFwdErrVal)
			*maxFwdErrVal = fwdErr;
	}
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*outStr << "\t\t\tComplex roots\t" << elapsed << "\n";
	}
}

void _errorAnalysisPFD_realRoots(std::vector< std::complex<double> > &r, std::vector< std::complex<double> > &dr, std::vector< std::complex<double> > &x, std::vector< std::complex<double> > &dx, std::vector<int> &realRootIndices, int prec, std::vector<double> *BEIntervalWidths, std::vector<double> *FEIntervalWidths, std::vector<double> *BEBoundaryError, std::vector<double> *FEBoundaryError, std::ostringstream *outStr, bool PROFILING) {
	//RationalNumber x,dx,BEdx,FEdx,y,rnErrTerm;
	//double S_tx,dSdx_tx,dSdt_tx,d2Sdxdt_tx,t,deltat,dBwdErrTerm,dFwdErrTerm;
	//ITRMIter it;
	//int iTemp,iTemp2,corr;
	double eps,dist,alpha,dalpha;
	int index;
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}
		
	for (int i=0; i<realRootIndices.size(); ++i) {
		index = realRootIndices.at(i);
		alpha = x.at(index).real();
		dalpha = dx.at(index).real();
		eps = pow(2.0,-prec);
		// forward error distance (x - alpha)
		dist = dalpha/eps;
		dist *= r.at(index).real();
		if (dist < 0.0)
			dist = -dist;
		FEIntervalWidths->push_back(dist);
		FEBoundaryError->push_back(eps);
		// backward error distance (x - alpha)
		dist = sqrt(dist);
		BEIntervalWidths->push_back(dist);
		BEBoundaryError->push_back(eps);
	}	
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*outStr << "\t\t\tReal roots\t" << elapsed << "\n";
	}	
}

void _errorAnalysisPFD_complexRoots(std::vector< std::complex<double> > &r, std::vector< std::complex<double> > &dr, std::vector< std::complex<double> > &x, std::vector< std::complex<double> > &dx, std::vector<int> &complexRootIndices, double *maxBwdErrVal, double *maxFwdErrVal, std::ostringstream *outStr, bool PROFILING) {
	//ComplexRationalNumber crnTemp;
	std::complex<double> root,cdErrTerm;
	//double a,b,x,da,db,A,B,C,dA,dB,bwdErr,fwdErr,dErrTerm;
	double bwdErr,fwdErr,xVal;
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}
	
	*maxBwdErrVal = 0.0;
	*maxFwdErrVal = 0.0;
	
	for (int i=0; i<complexRootIndices.size(); i++) {
		root = x.at(complexRootIndices.at(i));
		xVal = root.real();
		for (int j=0; j<x.size(); j++) {
			// backward error computation
			_computeBackwardErrorTerm(x.at(j),dx.at(j),r.at(j),dr.at(j),xVal,&cdErrTerm);
			if (root.imag() != 0.0)
				cdErrTerm *= 2.0;
			bwdErr += cdErrTerm.real();
			// forward error computation
			_computeForwardErrorTerm(x.at(j),dx.at(j),r.at(j),dr.at(j),xVal,&fwdErr);
		}
		if (bwdErr < 0.0)
			bwdErr = -bwdErr;
		if (bwdErr > *maxBwdErrVal)
			*maxBwdErrVal = bwdErr;
		if (fwdErr < 0.0)
			fwdErr = -fwdErr;
		if (fwdErr > *maxFwdErrVal)
			*maxFwdErrVal = fwdErr;
	}
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*outStr << "\t\t\tComplex roots\t" << elapsed << "\n";
	}
}



template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisLRT(std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &Stx, std::vector<UnivariatePolynomialOverRealField> &U, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdx, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &dSdt, std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > &d2Sdxdt, std::vector<UnivariatePolynomialOverRealField> &dU, std::vector< std::vector<ComplexRationalNumber> > &resultantRoots, std::vector<ComplexRationalNumber> &denomRoots, int prec, residueRootIndexMap *rrim, std::vector<ComplexRationalNumber> *residues, std::vector<int> *realSingIndices, std::vector<RationalNumber> *BEIntervalWidths, std::vector<RationalNumber> *FEIntervalWidths, std::vector<double> *BEBoundaryError, std::vector<double> *FEBoundaryError, double *maxBwdErr, double *maxFwdErr, float *ea_total, std::ostringstream *outStr, bool PROFILING){
	std::vector< std::vector< SparseUnivariateDoublePolynomial< std::complex<double> > > > S_t,dSdx_t,dSdt_t,d2Sdxdt_t;
	std::vector< std::vector< std::complex<double> > > t,dt;
	std::vector<double> realSings,complexSings,intervalRadii; // x values of integrand singularities
	std::vector<int> transitionIndices;
	ComplexRationalNumber root,cTemp;
	double* maxErr;
	double dTemp(0);
	int size_of_t(0),nworkers(1),dataLength(0);
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		*outStr << "\t\tLRT EA\n";
		*outStr << "\t\t------------------------------\n";
		startTimer(&start);
	}
	
	// Two branches here:
	//   1) compute the numerical polynomials for all the residues and resultant root error
	//   2) compute the integralTerm and denomRoot index correspondence
	// The resutant root error is currently the major bottleneck
	
	cilk_spawn
		_errorAnalysisLRT_computeNumericalPolynomials<UnivariatePolynomialOverRealField,RealField>(Stx,U,dSdx,dSdt,d2Sdxdt,dU,resultantRoots,&S_t,&dSdx_t,&dSdt_t,&d2Sdxdt_t,&t,&dt,outStr,PROFILING);
	
	std::vector<int> realRootIndices,complexRootIndices;
	IntegralTermRootMap realITRM,complexITRM;
	_errorAnalysisLRT_computeIntegralTermRootMap(denomRoots,resultantRoots,prec,rrim,residues,&realRootIndices,&complexRootIndices,
	&realITRM,&complexITRM,outStr,PROFILING);
		
	cilk_sync;
	
	// Two branches here, but their cost is low:
	//   1) compute the error bound over the complex roots
	//   2) compute the interval bounds where the error exceeds the tolerance (real singularities)
	 
	if (complexRootIndices.size() > 0) {
		//cilk_spawn
			_errorAnalysisLRT_complexRoots<UnivariatePolynomialOverRealField,RealField>(Stx,dSdx,dSdt,d2Sdxdt,S_t,dSdx_t,dSdt_t,d2Sdxdt_t,t,dt,denomRoots,resultantRoots,realITRM,complexITRM, complexRootIndices,maxBwdErr,maxFwdErr,outStr,PROFILING);
	}
	
	if (realRootIndices.size() > 0)
		_errorAnalysisLRT_realRoots<UnivariatePolynomialOverRealField,RealField>(Stx,dSdx,dSdt,d2Sdxdt,dt,denomRoots,resultantRoots,realITRM,prec,realSingIndices,BEIntervalWidths,FEIntervalWidths, BEBoundaryError,FEBoundaryError,outStr,PROFILING);
	
	//if (complexRootIndices.size() > 0)
	//	cilk_sync;
	
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*ea_total += elapsed;
		*outStr << "\t\t------------------------------\n";
		*outStr << "\t\t\t\t\t" << elapsed << "\n";
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisPFD_computeNumericalValues(UnivariatePolynomialOverRealField &G_in, UnivariatePolynomialOverRealField &H_in, UnivariatePolynomialOverRealField &dH_in, std::vector<ComplexRationalNumber> &denomRoots, std::vector<ComplexRationalNumber> &residues, residueRootIndexMap& rrim_in, std::vector< std::complex<double> > *r_out, std::vector< std::complex<double> > *dr_out, std::vector< std::complex<double> > *x_out, std::vector< std::complex<double> > *dx_out, std::ostringstream *outStr, bool PROFILING){
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}
	
	UnivariatePolynomialOverRealField G(G_in),H(H_in),dH(dH_in);
	std::complex<double> dTemp,xVal,dHVal;
	std::vector< std::complex<double> > r,dr,x,dx;
	residueRootIndexMap rrim(rrim_in);
	
	UnivariatePolynomialOverRealField dG(G),ddH(dH);
	dG.differentiate(1);
	ddH.differentiate(1);
	
	// this all can be made more efficient when there is residue structure by only computing
	// the residues and residue error for each one that is repeated. The rest of the error
	// analysis would need to change accordingly.
	int size(denomRoots.size());
	
	std::complex<double>* rArray = new std::complex<double>[size];
	std::complex<double>* xArray = new std::complex<double>[size];
	std::complex<double>* dxArray = new std::complex<double>[size];
	
	cilk_for (int i=0; i<size; ++i) {
		ComplexRationalNumber crn(denomRoots[i]);
		std::complex<double> cd(crn.realPart().get_d(),crn.imaginaryPart().get_d());
		xArray[i] = cd;
		ComplexRationalNumber crn2(-H.template evaluate<ComplexRationalNumber>(crn));
		crn2 /= dH.template evaluate<ComplexRationalNumber>(crn);
		cd = std::complex<double>(crn2.realPart().get_d(),crn2.imaginaryPart().get_d());
		dxArray[i] = cd;
	}
	
	RRIMIter it;
	for (it=rrim.begin(); it!=rrim.end(); ++it) {
		ComplexRationalNumber crn((*it).first);
		std::complex<double> cd(crn.realPart().get_d(),crn.imaginaryPart().get_d());
		rArray[(*it).second] = cd;
	}
	
	for (int i=0; i<denomRoots.size(); ++i) {
		x.push_back(xArray[i]);
		dx.push_back(dxArray[i]);
		r.push_back(rArray[i]);
	}
	SparseUnivariateDoublePolynomial<double> uG(G),udG(dG),udH(dH),uddH(ddH);
	
	for (int i=0; i<denomRoots.size(); ++i) {
		xVal = x[i];
		dHVal = udH.evaluate(xVal);
		dTemp = -uG.evaluate(xVal)*uddH.evaluate(xVal);
		dTemp /= dHVal;
		dTemp += udG.evaluate(xVal);
		dTemp /= dHVal;
		dr.push_back(dTemp);
	}
	
	delete [] rArray;
	delete [] xArray;
	delete [] dxArray;
	
	*r_out = r;
	*dr_out = dr;
	*x_out = x;
	*dx_out = dx;
	
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*outStr << "\t\t\tnumerical vals\t" << elapsed << "\n";
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _errorAnalysisPFD(UnivariatePolynomialOverRealField &G, UnivariatePolynomialOverRealField &H, UnivariatePolynomialOverRealField &dH, std::vector<ComplexRationalNumber> &denomRoots, std::vector<ComplexRationalNumber> &residues, int prec, std::vector<int> &realRootIndices, std::vector<int> &complexRootIndices, residueRootIndexMap& rrim, double *maxBwdErrVal, double *maxFwdErrVal, std::vector<double> *BEIntervalWidths, std::vector<double> *FEIntervalWidths, std::vector<double> *BEBoundaryError, std::vector<double> *FEBoundaryError, float *ea_total, std::ostringstream *outStr, bool PROFILING){
	std::vector< std::complex<double> > r,dr,x,dx;
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		*outStr << "\t\tPFD EA\n";
		*outStr << "\t\t------------------------------\n";
		startTimer(&start);
	}
	
	_errorAnalysisPFD_computeNumericalValues<UnivariatePolynomialOverRealField,RealField>(G,H,dH,denomRoots,residues,rrim,&r,&dr,&x,&dx,outStr,PROFILING);
	/*for (int i=0; i<denomRoots.size(); ++i) {
		*outStr << "x = " << x.at(i) << " +/- " << dx.at(i) << std::endl;
		*outStr << "r = " << r.at(i) << ", dr = " << dr.at(i) << std::endl;
	}*/
	
	_errorAnalysisPFD_complexRoots(r,dr,x,dx,complexRootIndices,maxBwdErrVal,maxFwdErrVal,outStr,PROFILING);
	_errorAnalysisPFD_realRoots(r,dr,x,dx,realRootIndices,prec,BEIntervalWidths,FEIntervalWidths,BEBoundaryError,  FEBoundaryError,outStr,PROFILING);
	
		
	if (PROFILING){
		stopTimer(&start,&elapsed);
		*ea_total += elapsed;
		*outStr << "\t\t------------------------------\n";
		*outStr << "\t\t\t\t\t" << elapsed << "\n";
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void reverseVariableOrder(SparseUnivariatePolynomial<UnivariatePolynomialOverRealField>* A){
	Symbol variable;
	RealField rfTemp;
	UnivariatePolynomialOverRealField coef;
	UnivariatePolynomialOverRealField upTemp;
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> out,supTemp;
	variable = A->variable();
	upTemp.setVariableName(variable);
	variable = A->leadingCoefficient().variable();
	supTemp.setVariableName(variable);
	out.setVariableName(variable);
	for (int i=A->degree(); i>-1; i--) {
		coef = A->coefficient(i);
		if (!coef.isZero()){
			for (int j=0; j<coef.degree()+1; j++){
				rfTemp = coef.coefficient(j);
				if (!(rfTemp == 0)) {
					upTemp.setCoefficient(i,rfTemp);
					supTemp.setCoefficient(j,upTemp);
				}
			}
			out += supTemp;
			upTemp.setCoefficient(i,RealField(RationalNumber(0)));
		}
	}
	*A = out;
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _realSNIntegrate(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS){
	
	UnivariatePolynomialOverRealField Num(A);
	UnivariatePolynomialOverRealField Den(D);
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> supTemp;
	
	std::vector<UnivariatePolynomialOverRealField> U,dU;
	std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > Stx,Sxt,dSdx,dSdt,d2Sdxdt;
	std::vector< std::vector<ComplexRationalNumber> > resultantRoots;
	std::vector<ComplexRationalNumber> denomRoots;
	std::vector<ComplexRationalNumber> residues;
	std::vector<RationalNumber> BEIntervalWidths,FEIntervalWidths;
	std::vector<double> BEBoundaryError,FEBoundaryError;
	std::vector<int> realSingIndices;
	double maxBwdErr,maxFwdErr;
	residueRootIndexMap rrim;
	
	// variable checking and canonicalization
	_initializeRationalIntegration<UnivariatePolynomialOverRealField,RealField>(A,D,PROFILING);
	
	// set main variable
	Symbol mainVariable;
	mainVariable = A.variable();
	
	// Profiling variables
	unsigned long long start, start2, start3;
	float elapsed(0), elapsed2(0), elapsed3(0), ea_total(0);
	std::ostringstream outStr;
	
	// compute the rational part of the integral of R/D by Hermite reduction, which is 
	// returned as G.at(0)/G.at(1), and the polynomial part, returned in P, along with the 
	// derivative of the rest of the integral, viz., the integrand of the logarithmic 
	// part, which is returned as a single rational function as a pair of polynomials 
	// R and H
	
	UnivariatePolynomialOverRealField R;  // stores numerator of log part of integral of A/D
	UnivariatePolynomialOverRealField H;  // stores squarefee denominator of the log part
	UnivariatePolynomialOverRealField dH; // stores derivative of H
	R.setVariableName(mainVariable);
	H.setVariableName(mainVariable);
	
	if (PROFILING) {
		std::cout << "\tRational function rational part" << std::endl;
		std::cout << "\t------------------------------" << std::endl;
		startTimer(&start);
	}
		
	_integrateRationalFunctionRationalPart<UnivariatePolynomialOverRealField,RealField>(A,D,P,G,&R,&H,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}
	
	// compute integral of R/H
	
	if (!R.isZero()) {
		
		if (ERROR_ANALYSIS) {
			
			cilk_spawn
				_errorAnalysisLRT_rootsAndResidues<UnivariatePolynomialOverRealField,RealField>(R,H,prec,&dH,&denomRoots,&residues,&rrim,&ea_total,&outStr,PROFILING,false);
				
		}
		
		if (PROFILING){
			std::cout << "\tRational function log part" << std::endl;
			std::cout << "\t------------------------------" << std::endl;
			startTimer(&start);
		}
		
		_integrateRationalFunctionLogPart<UnivariatePolynomialOverRealField,RealField>(&Stx,&U,R,H,PROFILING);
		
		if (PROFILING){
			stopTimer(&start,&elapsed);
			std::cout << "\t------------------------------" << std::endl;
			std::cout << "\t\t\t" << elapsed << std::endl;
		}
	}
	
	if (ERROR_ANALYSIS){
		cilk_spawn 
			_computePartialDerivatives<UnivariatePolynomialOverRealField,RealField>(Stx,U,&dSdx,&dSdt,&d2Sdxdt,&dU,&ea_total,&outStr,PROFILING);
	}

	
	if (PROFILING)
		startTimer(&start);
	
	resultantRoots = _rootsMultiprecision<UnivariatePolynomialOverRealField,RealField>(U,prec);
	
	for (int i=0; i<resultantRoots.size(); i++)
		sort(resultantRoots.at(i).begin(),resultantRoots.at(i).end(),ComplexRationalNumberOrdering(prec));
	
	std::vector< std::vector<ComplexRationalNumber> > resultantRootsTemp(resultantRoots);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\tRoots\t\t" << elapsed << std::endl;
	}
		
	// This needs to be moved to the optimal position
	if (ERROR_ANALYSIS){
		cilk_sync;
		cilk_spawn
			_errorAnalysisLRT<UnivariatePolynomialOverRealField,RealField>(Stx,U,dSdx,dSdt,d2Sdxdt,dU,resultantRootsTemp,denomRoots,prec,&rrim,&residues,&realSingIndices, &BEIntervalWidths,&FEIntervalWidths,&BEBoundaryError,&FEBoundaryError,&maxBwdErr,&maxFwdErr, &ea_total,&outStr,PROFILING);
	}
	
	if (PROFILING){
		std::cout << "\tLog to real" << std::endl;
		std::cout << "\t--------------------------------------" << std::endl;
		startTimer(&start);
	}
	
	_logToReal<UnivariatePolynomialOverRealField,RealField>(resultantRoots,Stx,lg,Lg,atn,Atn,prec,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t--------------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}
	
	if (ERROR_ANALYSIS) {
		cilk_sync;
		if (PROFILING){
			outStr << "\t------------------------------\n";
			outStr << "\t\t\t\t\t" << ea_total << "\n";
			outStr << "\t\tBackward error away from sings = " << maxBwdErr << std::endl;
			if (BEIntervalWidths.size() > 0) {
				outStr << "\t\tBackward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << BEIntervalWidths.at(i).get_d() << "  (" << BEBoundaryError.at(i) << ")" << std::endl;
				}
			}
			outStr << "\t\tForward Error away from sings = " << maxFwdErr << std::endl;
			if (FEIntervalWidths.size() > 0) {
				outStr << "\t\tForward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << FEIntervalWidths.at(i).get_d() << "  (" << FEBoundaryError.at(i) << ")" << std::endl;
				}
			}
			std::cout << outStr.str();
		}
	}

	// divide by the Davenport divisor
	//for (int i=0; i<lg.size(); i++)
	//	lg.at(i) /= D.degree();
	//for (int i=0; i<atn.size(); i++)
	//	atn.at(i) /= D.degree();
	
}

//template <class UnivariatePolynomialOverRealField, class RealField>
//void _realSNIntegrate(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn1, std::vector<UnivariatePolynomialOverRealField> *Atn2, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS){
//	
//	UnivariatePolynomialOverRealField Num(A);
//	UnivariatePolynomialOverRealField Den(D);
//	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> supTemp;
//	
//	std::vector<UnivariatePolynomialOverRealField> U,dU;
//	std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > Stx,Sxt,dSdx,dSdt,d2Sdxdt;
//	std::vector< std::vector<ComplexRationalNumber> > resultantRoots;
//	std::vector<ComplexRationalNumber> denomRoots;
//	std::vector<ComplexRationalNumber> residues;
//	std::vector<RationalNumber> BEIntervalWidths,FEIntervalWidths;
//	std::vector<double> BEBoundaryError,FEBoundaryError;
//	std::vector<int> realSingIndices;
//	double maxBwdErr,maxFwdErr;
//	residueRootIndexMap rrim;
//	
//	// variable checking and canonicalization
//	_initializeRationalIntegration<UnivariatePolynomialOverRealField,RealField>(A,D,PROFILING);
//	
//	// set main variable
//	Symbol mainVariable;
//	mainVariable = A.variable();
//	
//	// Profiling variables
//	unsigned long long start, start2, start3;
//	float elapsed(0), elapsed2(0), elapsed3(0), ea_total(0);
//	std::ostringstream outStr;
//	
//	// compute the rational part of the integral of R/D by Hermite reduction, which is 
//	// returned as G.at(0)/G.at(1), and the polynomial part, returned in P, along with the 
//	// derivative of the rest of the integral, viz., the integrand of the logarithmic 
//	// part, which is returned as a single rational function as a pair of polynomials 
//	// R and H
//	
//	UnivariatePolynomialOverRealField R;  // stores numerator of log part of integral of A/D
//	UnivariatePolynomialOverRealField H;  // stores squarefee denominator of the log part
//	UnivariatePolynomialOverRealField dH; // stores derivative of H
//	R.setVariableName(mainVariable);
//	H.setVariableName(mainVariable);
//	
//	if (PROFILING) {
//		std::cout << "\tRational function rational part" << std::endl;
//		std::cout << "\t------------------------------" << std::endl;
//		startTimer(&start);
//	}
//		
//	_integrateRationalFunctionRationalPart<UnivariatePolynomialOverRealField,RealField>(A,D,P,G,&R,&H,PROFILING);
//	
//	if (PROFILING){
//		stopTimer(&start,&elapsed);
//		std::cout << "\t------------------------------" << std::endl;
//		std::cout << "\t\t\t" << elapsed << std::endl;
//	}
//	
//	// compute integral of R/H
//	
//	if (!R.isZero()) {
//			
//		cilk_spawn
//			_errorAnalysisLRT_rootsAndResidues<UnivariatePolynomialOverRealField,RealField>(R,H,prec,&dH,&denomRoots,&residues,&rrim,&ea_total,&outStr,PROFILING,false);
//			
//		
//		if (PROFILING){
//			std::cout << "\tRational function log part" << std::endl;
//			std::cout << "\t------------------------------" << std::endl;
//			startTimer(&start);
//		}
//		
//		_integrateRationalFunctionLogPart<UnivariatePolynomialOverRealField,RealField>(&Stx,&U,R,H,PROFILING);
//		
//		if (PROFILING){
//			stopTimer(&start,&elapsed);
//			std::cout << "\t------------------------------" << std::endl;
//			std::cout << "\t\t\t" << elapsed << std::endl;
//		}
//	}
//	
//	cilk_spawn 
//		_computePartialDerivatives<UnivariatePolynomialOverRealField,RealField>(Stx,U,&dSdx,&dSdt,&d2Sdxdt,&dU,&ea_total,&outStr,PROFILING);

//	
//	if (PROFILING)
//		startTimer(&start);
//	
//	resultantRoots = _rootsMultiprecision<UnivariatePolynomialOverRealField,RealField>(U,prec);
//	
//	for (int i=0; i<resultantRoots.size(); i++)
//		sort(resultantRoots.at(i).begin(),resultantRoots.at(i).end(),ComplexRationalNumberOrdering(prec));
//	
//	std::vector< std::vector<ComplexRationalNumber> > resultantRootsTemp(resultantRoots);
//	
//	if (PROFILING){
//		stopTimer(&start,&elapsed);
//		std::cout << "\tRoots\t\t" << elapsed << std::endl;
//	}
//		
//	// This may need to be moved to an optimal position
//	cilk_sync;
//	cilk_spawn
//		_errorAnalysisLRT<UnivariatePolynomialOverRealField,RealField>(Stx,U,dSdx,dSdt,d2Sdxdt,dU,resultantRootsTemp,denomRoots,prec,&rrim,&residues,&realSingIndices, &BEIntervalWidths,&FEIntervalWidths,&BEBoundaryError,&FEBoundaryError,&maxBwdErr,&maxFwdErr, &ea_total,&outStr,PROFILING);
//	
//	if (PROFILING){
//		std::cout << "\tLog to real" << std::endl;
//		std::cout << "\t--------------------------------------" << std::endl;
//		startTimer(&start);
//	}
//	
//	_logToReal<UnivariatePolynomialOverRealField,RealField>(resultantRoots,Stx,lg,Lg,atn,Atn1,prec,PROFILING);
//	
//	if (PROFILING){
//		stopTimer(&start,&elapsed);
//		std::cout << "\t--------------------------------------" << std::endl;
//		std::cout << "\t\t\t" << elapsed << std::endl;
//	}
//	
//	cilk_sync;
//	if (PROFILING){
//		outStr << "\t------------------------------\n";
//		outStr << "\t\t\t\t\t" << ea_total << "\n";
//		outStr << "\t\tBackward error away from sings = " << maxBwdErr << std::endl;
//		if (BEIntervalWidths.size() > 0) {
//			outStr << "\t\tBackward error exceeds tolerance for x =" << std::endl;
//		
//			for (int i=0; i<BEIntervalWidths.size(); i++) {
//				outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << BEIntervalWidths.at(i).get_d() << "  (" << BEBoundaryError.at(i) << ")" << std::endl;
//			}
//		}
//		outStr << "\t\tForward Error away from sings = " << maxFwdErr << std::endl;
//		if (FEIntervalWidths.size() > 0) {
//			outStr << "\t\tForward error exceeds tolerance for x =" << std::endl;
//		
//			for (int i=0; i<BEIntervalWidths.size(); i++) {
//				outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << FEIntervalWidths.at(i).get_d() << "  (" << FEBoundaryError.at(i) << ")" << std::endl;
//			}
//		}
//		std::cout << outStr.str();
//	}

//	// divide by the Davenport divisor
//	//for (int i=0; i<lg.size(); i++)
//	//	lg.at(i) /= D.degree();
//	//for (int i=0; i<atn.size(); i++)
//	//	atn.at(i) /= D.degree();
//	
//}

template <class UnivariatePolynomialOverRealField, class RealField>
void _realSNIntegrate(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn1, std::vector<UnivariatePolynomialOverRealField> *Atn2, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS){
	
	UnivariatePolynomialOverRealField Num(A);
	UnivariatePolynomialOverRealField Den(D);
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> supTemp;
	
	std::vector<UnivariatePolynomialOverRealField> U,dU;
	std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > Stx,Sxt,dSdx,dSdt,d2Sdxdt;
	std::vector< std::vector<ComplexRationalNumber> > resultantRoots;
	std::vector<ComplexRationalNumber> denomRoots;
	std::vector<ComplexRationalNumber> residues;
	std::vector<RationalNumber> BEIntervalWidths,FEIntervalWidths;
	std::vector<double> BEBoundaryError,FEBoundaryError;
	std::vector<int> realSingIndices;
	double maxBwdErr,maxFwdErr;
	residueRootIndexMap rrim;
	
	// variable checking and canonicalization
	_initializeRationalIntegration<UnivariatePolynomialOverRealField,RealField>(A,D,PROFILING);
	
	// set main variable
	Symbol mainVariable;
	mainVariable = A.variable();
	
	// Profiling variables
	unsigned long long start, start2, start3;
	float elapsed(0), elapsed2(0), elapsed3(0), ea_total(0);
	std::ostringstream outStr;
	
	// compute the rational part of the integral of R/D by Hermite reduction, which is 
	// returned as G.at(0)/G.at(1), and the polynomial part, returned in P, along with the 
	// derivative of the rest of the integral, viz., the integrand of the logarithmic 
	// part, which is returned as a single rational function as a pair of polynomials 
	// R and H
	
	UnivariatePolynomialOverRealField R;  // stores numerator of log part of integral of A/D
	UnivariatePolynomialOverRealField H;  // stores squarefee denominator of the log part
	UnivariatePolynomialOverRealField dH; // stores derivative of H
	R.setVariableName(mainVariable);
	H.setVariableName(mainVariable);
	
	if (PROFILING) {
		std::cout << "\tRational function rational part" << std::endl;
		std::cout << "\t------------------------------" << std::endl;
		startTimer(&start);
	}
		
	_integrateRationalFunctionRationalPart<UnivariatePolynomialOverRealField,RealField>(A,D,P,G,&R,&H,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}
	
	// compute integral of R/H
	
	if (!R.isZero()) {
		
		if (ERROR_ANALYSIS) {
			
			cilk_spawn
				_errorAnalysisLRT_rootsAndResidues<UnivariatePolynomialOverRealField,RealField>(R,H,prec,&dH,&denomRoots,&residues,&rrim,&ea_total,&outStr,PROFILING,false);
				
		}
		
		if (PROFILING){
			std::cout << "\tRational function log part" << std::endl;
			std::cout << "\t------------------------------" << std::endl;
			startTimer(&start);
		}
		
		_integrateRationalFunctionLogPart<UnivariatePolynomialOverRealField,RealField>(&Stx,&U,R,H,PROFILING);
		
		if (PROFILING){
			stopTimer(&start,&elapsed);
			std::cout << "\t------------------------------" << std::endl;
			std::cout << "\t\t\t" << elapsed << std::endl;
		}
	}
	
	if (ERROR_ANALYSIS){
		cilk_spawn 
			_computePartialDerivatives<UnivariatePolynomialOverRealField,RealField>(Stx,U,&dSdx,&dSdt,&d2Sdxdt,&dU,&ea_total,&outStr,PROFILING);
	}

	
	if (PROFILING)
		startTimer(&start);
	
	resultantRoots = _rootsMultiprecision<UnivariatePolynomialOverRealField,RealField>(U,prec);
	
	for (int i=0; i<resultantRoots.size(); i++)
		sort(resultantRoots.at(i).begin(),resultantRoots.at(i).end(),ComplexRationalNumberOrdering(prec));
	
	std::vector< std::vector<ComplexRationalNumber> > resultantRootsTemp(resultantRoots);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\tRoots\t\t" << elapsed << std::endl;
	}
		
	// This needs to be moved to the optimal position
	if (ERROR_ANALYSIS){
		cilk_sync;
		cilk_spawn
			_errorAnalysisLRT<UnivariatePolynomialOverRealField,RealField>(Stx,U,dSdx,dSdt,d2Sdxdt,dU,resultantRootsTemp,denomRoots,prec,&rrim,&residues,&realSingIndices, &BEIntervalWidths,&FEIntervalWidths,&BEBoundaryError,&FEBoundaryError,&maxBwdErr,&maxFwdErr, &ea_total,&outStr,PROFILING);
	}
	
	if (PROFILING){
		std::cout << "\tLog to real" << std::endl;
		std::cout << "\t--------------------------------------" << std::endl;
		startTimer(&start);
	}
	
	_logToReal<UnivariatePolynomialOverRealField,RealField>(resultantRoots,Stx,lg,Lg,atn,Atn1,prec,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t--------------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}
	
	if (ERROR_ANALYSIS) {
		cilk_sync;
		if (PROFILING){
			outStr << "\t------------------------------\n";
			outStr << "\t\t\t\t\t" << ea_total << "\n";
			outStr << "\t\tBackward error away from sings = " << maxBwdErr << std::endl;
			if (BEIntervalWidths.size() > 0) {
				outStr << "\t\tBackward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << BEIntervalWidths.at(i).get_d() << "  (" << BEBoundaryError.at(i) << ")" << std::endl;
				}
			}
			outStr << "\t\tForward Error away from sings = " << maxFwdErr << std::endl;
			if (FEIntervalWidths.size() > 0) {
				outStr << "\t\tForward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << FEIntervalWidths.at(i).get_d() << "  (" << FEBoundaryError.at(i) << ")" << std::endl;
				}
			}
			std::cout << outStr.str();
		}
	}

	// divide by the Davenport divisor
	//for (int i=0; i<lg.size(); i++)
	//	lg.at(i) /= D.degree();
	//for (int i=0; i<atn.size(); i++)
	//	atn.at(i) /= D.degree();
	
}

SparseUnivariatePolynomial<RationalNumber> _realPart(SparseUnivariatePolynomial<ComplexRationalNumber> &A){
	SparseUnivariatePolynomial<RationalNumber> output;
	output.setVariableName(A.variable());
	RationalNumber rTemp;
	for (int i=0; i<A.degree()+1; i++){
		rTemp = A.coefficient(i).realPart();
		if (!rTemp.isZero())
			output.setCoefficient(i,rTemp);
	}
	return output;
}

SparseUnivariatePolynomial<RationalNumber> _imaginaryPart(SparseUnivariatePolynomial<ComplexRationalNumber> &A){
	SparseUnivariatePolynomial<RationalNumber> output;
	output.setVariableName(A.variable());
	RationalNumber rTemp;
	for (int i=0; i<A.degree()+1; i++){
		rTemp = A.coefficient(i).imaginaryPart();
		if (!rTemp.isZero())
			output.setCoefficient(i,rTemp);
	}
	return output;
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _sumArctanTerms(std::vector<UnivariatePolynomialOverRealField>& arctanTerms, UnivariatePolynomialOverRealField* A, UnivariatePolynomialOverRealField* B) {
	UnivariatePolynomialOverRealField termTwo;
	UnivariatePolynomialOverRealField one;
	one.one();
	if (arctanTerms.size() > 0) {
		if (arctanTerms.size() == 1) {
			*A = arctanTerms.at(0);
			*B = termTwo;
		}
		else {
			UnivariatePolynomialOverRealField termOne;
			termTwo = arctanTerms.at(0);
			for (int i=1; i<arctanTerms.size(); ++i) {
				termOne += arctanTerms.at(i);
				termTwo *= -arctanTerms.at(i);
				termTwo += one;
			}
			*A = termOne;
			*B = termTwo;
		}
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _computeIntegralTerm(std::vector<ComplexRationalNumber> &denomRoots, ComplexRationalNumber residue, std::vector<int> rootIndices, int prec, Symbol mainVariable, RealField *lg, UnivariatePolynomialOverRealField *Lg, RealField *atn, UnivariatePolynomialOverRealField *Atn1, UnivariatePolynomialOverRealField *Atn2) {

	UnivariatePolynomialOverRealField logTerm,arctanTerm1,arctanTerm2;
	ComplexRationalNumber root;
	RealField one(RationalNumber(1));
	RealField zero(RationalNumber(0));
	UnivariatePolynomialOverRealField LinPoly; // poly of the form a*x + b
	LinPoly.setVariableName(mainVariable);
	LinPoly = (LinPoly + one << 1);
	logTerm.one();
	arctanTerm1.one();
	arctanTerm2.one();
	
	for (int i=0; i<rootIndices.size(); ++i) {
		root = denomRoots.at(rootIndices.at(i));
		if (root.imaginaryPart() == 0) {
			if (!epsilonEqual(residue.realPart(),zero,prec)) {
				RealField rf(-root.realPart());
				LinPoly.setCoefficient(0,rf);
				logTerm *= LinPoly;
			}
		}
		else {
			RealField a(root.realPart());
			RealField b(root.imaginaryPart());
			if (!epsilonEqual(residue.imaginaryPart(),zero,prec)){
				RealField rTemp(-a);
				LinPoly.setCoefficient(0,rTemp); // LinPoly = (x - a)
				if (arctanTerm1.isOne()) {
					RealField rf(-b);
					arctanTerm1 = LinPoly;
					arctanTerm2 *= rf;
				}
				else {
					UnivariatePolynomialOverRealField ATemp(arctanTerm1),BTemp(arctanTerm2);
					ATemp *= b;
					BTemp *= b;
					arctanTerm1 *= LinPoly;
					arctanTerm2 *= LinPoly;
					arctanTerm1 += BTemp;
					arctanTerm2 -= ATemp;
				}
			}
			// (x - (a+ib))(x - (a-ib)) = (x^2 - 2ax + a^2+b^2)
			if (!epsilonEqual(residue.realPart(),zero,prec)){
				UnivariatePolynomialOverRealField QuadPoly;
				QuadPoly.setVariableName(mainVariable);
				QuadPoly = (QuadPoly + one << 2); // x^2
				RealField rTemp(-2*a);
				QuadPoly.setCoefficient(1,rTemp);
				rTemp = RealField(a*a+b*b);
				QuadPoly.setCoefficient(0,rTemp);
				logTerm *= QuadPoly;
			}
		}
	}
	
	// The rounding is turned off, since it can give incorrect results on some problems.
	// To incorporate it there must be checking that the rounding will not significantly affect the error.
	if (logTerm.degree() > 0) {
		RealField coef;
		for (int j=0; j<=logTerm.degree(); ++j) {
			coef = logTerm.coefficient(j);
			//if (epsilonEqual(coef,zero,prec))
			//	logTerm.setCoefficient(j,zero);
		}
		*lg = residue.realPart();
		*Lg = logTerm;
	}
	if (arctanTerm1.degree() > 0) {
		RealField coef;
		for (int j=0; j<=arctanTerm1.degree(); ++j) {
			coef = arctanTerm1.coefficient(j);
			//if (epsilonEqual(coef,zero,prec))
			//	arctanTerm1.setCoefficient(j,zero);
		}
		for (int j=0; j<=arctanTerm2.degree(); ++j) {
			coef = arctanTerm2.coefficient(j);
			//if (epsilonEqual(coef,zero,prec))
			//	arctanTerm2.setCoefficient(j,zero);
		}
		RealField rTemp(residue.imaginaryPart()*2);
		*atn = rTemp;
		*Atn1 = arctanTerm1;
		*Atn2 = arctanTerm2;
	}
}

template <class UnivariatePolynomialOverRealField, class RealField>
void _PFDStructuredIntegrate(std::vector<ComplexRationalNumber> &denomRoots, residueRootIndexMap &rrim, Symbol mainVariable, int prec, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn1, std::vector<UnivariatePolynomialOverRealField> *Atn2, bool PROFILING){
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}

	RRIMIter it;
	
	std::vector<ComplexRationalNumber> residues;
	std::vector< std::vector<int> > residueRootIndices;
	std::vector<int> ind;
	//TermRootMap resRootCor; // correspondence between residue and roots
	
	it=rrim.begin();
	residues.push_back((*it).first);
	//resRootCor.insert(std::make_pair(0,(*it).second));
	ind.push_back((*it).second);
	++it;
	for (; it!=rrim.end(); ++it) {
		if (!epsilonEqual(residues.back(),(*it).first,prec)) {
			residueRootIndices.push_back(ind);
			ind.clear();
			residues.push_back((*it).first);
		}
		//resRootCor.insert(std::make_pair(residues.size()-1,(*it).second));
		ind.push_back((*it).second);
	}
	residueRootIndices.push_back(ind);
	
	int size(residues.size());
	
	RealField* lgArray = new RealField[size];
	RealField* atnArray = new RealField[size];
	UnivariatePolynomialOverRealField* LgArray = new UnivariatePolynomialOverRealField[size];
	UnivariatePolynomialOverRealField* Atn1Array = new UnivariatePolynomialOverRealField[size];
	UnivariatePolynomialOverRealField* Atn2Array = new UnivariatePolynomialOverRealField[size];
	
	cilk_for (int i=0; i<size; ++i) {
		RealField log,arctan;
		UnivariatePolynomialOverRealField Log,Arctan1,Arctan2;
		_computeIntegralTerm<UnivariatePolynomialOverRealField,RealField>(denomRoots,residues.at(i),residueRootIndices.at(i),prec,mainVariable,&log,&Log,&arctan,&Arctan1,&Arctan2);
		lgArray[i] = log;
		LgArray[i] = Log;
		atnArray[i] = arctan;
		Atn1Array[i] = Arctan1;
		Atn2Array[i] = Arctan2;
	}
	
	for (int i=0; i<size; ++i) {
		if (LgArray[i].degree() > 0) {
			lg->push_back(lgArray[i]);
			Lg->push_back(LgArray[i]);
		}
		if (Atn1Array[i].degree() > 0) {
			atn->push_back(atnArray[i]);
			Atn1->push_back(Atn1Array[i]);
			Atn2->push_back(Atn2Array[i]);
		}
	}
	
	delete [] lgArray;
	delete [] LgArray;
	delete [] atnArray;
	delete [] Atn1Array;
	delete [] Atn2Array;
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t\tstruc integrate\t\t" << elapsed << "\n";
	}
}

/*template <class UnivariatePolynomialOverRealField, class RealField>
void _PFDStructuredIntegrate(std::vector<ComplexRationalNumber> &denomRoots, residueRootIndexMap &rrim, std::string mainVariable, int prec, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atnTemp, std::vector<UnivariatePolynomialOverRealField> *Atn1, std::vector<UnivariatePolynomialOverRealField> *Atn2, TermRootMap *logMap, TermRootMap *arctanMap, bool PROFILING){
	
	// Profiling variables
	unsigned long long start;
	float elapsed(0);
	
	if (PROFILING){
		startTimer(&start);
	}

	UnivariatePolynomialOverRealField logTerm,arctanTerm1,arctanTerm2;
	ComplexRationalNumber c,root,residue,previousResidue;
	RealField arctanCoef;
	RealField one(RationalNumber(1));
	RealField zero(RationalNumber(0));
	UnivariatePolynomialOverRealField LinPoly; // poly of the form a*x + b
	LinPoly.setVariableName(mainVariable);
	LinPoly = (LinPoly + one << 1);
	logTerm.one();
	arctanTerm1.one();
	arctanTerm2.one();
	RRIMIter it,itTemp,itEnd;
	itEnd = rrim.end();
	--itEnd;
	
	
	// this code should probably be rewritten to precompute residue multiplicity
	for (it=rrim.begin(); it!=rrim.end(); ++it){
		root = denomRoots.at((*it).second);
		residue = (*it).first;
		//std::cout << "iteration " << i << ": residue = (" << residue.realPart().get_d() << "," << residue.imaginaryPart().get_d() << ")" << std::endl;
		//std::cout << "denom root = (" << root.realPart().get_d() << "," << root.imaginaryPart().get_d() << ")" << std::endl;
		itTemp = it;
		--itTemp;
		if (it!=rrim.begin() && !epsilonEqual(residue,(*itTemp).first,prec)) {
			logTerm.one();
			arctanTerm1.one();
			arctanTerm2.one();
		}
		if (root.imaginaryPart() == 0) {
			if (!epsilonEqual(residue.realPart(),zero,prec)) {
				RealField rf(-root.realPart());
				LinPoly.setCoefficient(0,rf);
				logTerm *= LinPoly;
				logMap->insert(std::make_pair(Lg->size(),(*it).second));
			}
		}
		else {
			RealField a(root.realPart());
			RealField b(root.imaginaryPart());
			if (!epsilonEqual(residue.imaginaryPart(),zero,prec)){
				RealField rTemp(-a);
				LinPoly.setCoefficient(0,rTemp); // LinPoly = (x - a)
				if (arctanTerm1.isOne()) {
					RealField rf(-b);
					arctanTerm1 = LinPoly;
					arctanTerm2 *= rf;
				}
				else {
					UnivariatePolynomialOverRealField ATemp(arctanTerm1),BTemp(arctanTerm2);
					ATemp *= b;
					BTemp *= b;
					arctanTerm1 *= LinPoly;
					arctanTerm2 *= LinPoly;
					arctanTerm1 += BTemp;
					arctanTerm2 -= ATemp;
				}
				arctanMap->insert(std::make_pair(Atn1->size(),(*it).second));
			}
			// (x - (a+ib))(x - (a-ib)) = (x^2 - 2ax + a^2+b^2)
			if (!epsilonEqual(residue.realPart(),zero,prec)){
				UnivariatePolynomialOverRealField QuadPoly;
				QuadPoly.setVariableName(mainVariable);
				QuadPoly = (QuadPoly + one << 2); // x^2
				RealField rTemp(-2*a);
				QuadPoly.setCoefficient(1,rTemp);
				rTemp = RealField(a*a+b*b);
				QuadPoly.setCoefficient(0,rTemp);
				logTerm *= QuadPoly;
				logMap->insert(std::make_pair(Lg->size(),(*it).second));
			}
		}
		itTemp = it;
		++itTemp;
		if (it==itEnd || !epsilonEqual(residue,(*itTemp).first,prec)) {
			if (logTerm.degree() > 0) {
				RealField coef;
				for (int j=0; j<=logTerm.degree(); ++j) {
					coef = logTerm.coefficient(j);
					if (epsilonEqual(coef,zero,prec))
						logTerm.setCoefficient(j,zero);
				}
				lg->push_back(residue.realPart());
				Lg->push_back(logTerm);
			}
			if (arctanTerm1.degree() > 0) {
				RealField coef;
				for (int j=0; j<=arctanTerm1.degree(); ++j) {
					coef = arctanTerm1.coefficient(j);
					if (epsilonEqual(coef,zero,prec))
						arctanTerm1.setCoefficient(j,zero);
				}
				for (int j=0; j<=arctanTerm2.degree(); ++j) {
					coef = arctanTerm2.coefficient(j);
					if (epsilonEqual(coef,zero,prec))
						arctanTerm2.setCoefficient(j,zero);
				}
				RealField rTemp(residue.imaginaryPart()*2);
				atnTemp->push_back(rTemp);
				Atn1->push_back(arctanTerm1);
				Atn2->push_back(arctanTerm2);
			}
		}
	}
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t\tstruc integrate\t\t" << elapsed << "\n";
	}
}*/

template <class UnivariatePolynomialOverRealField, class RealField>
void _realSNIntegratePFD(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS){
	
	UnivariatePolynomialOverRealField Num(A);
	UnivariatePolynomialOverRealField Den(D);
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> supTemp;
	
	std::vector<RealField> atnTemp;
	std::vector<UnivariatePolynomialOverRealField> U,dU,Atn1,Atn2;
	std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > Stx,Sxt,dSdx,dSdt,d2Sdxdt;
	std::vector< std::vector<ComplexRationalNumber> > resultantRoots;
	std::vector<ComplexRationalNumber> denomRoots;
	std::vector<ComplexRationalNumber> residues;
	std::vector<double> BEIntervalWidths,FEIntervalWidths,BEBoundaryError,FEBoundaryError;
	std::vector<int> realSingIndices;
	std::vector<int> realRootIndices,complexRootIndices;
	double maxBwdErr,maxFwdErr;
	residueRootIndexMap rrim;
	
	// variable checking and canonicalization
	_initializeRationalIntegration<UnivariatePolynomialOverRealField,RealField>(A,D,PROFILING);
	
	// set main variable
	Symbol mainVariable;
	mainVariable = A.variable();
	
	// Profiling variables
	unsigned long long start, start2, start3;
	float elapsed(0), elapsed2(0), elapsed3(0), ea_total(0);
	std::ostringstream outStr;
	
	// compute the rational part of the integral of R/D by Hermite reduction, which is 
	// returned as G.at(0)/G.at(1), and the polynomial part, returned in P, along with the 
	// derivative of the rest of the integral, viz., the integrand of the logarithmic 
	// part, which is returned as a single rational function as a pair of polynomials 
	// R and H
	
	UnivariatePolynomialOverRealField R;  // stores numerator of log part of integral of A/D
	UnivariatePolynomialOverRealField H;  // stores squarefee denominator of the log part
	UnivariatePolynomialOverRealField dH; // stores derivative of H
	R.setVariableName(mainVariable);
	H.setVariableName(mainVariable);
	
	if (PROFILING) {
		std::cout << "\tRational function rational part" << std::endl;
		std::cout << "\t------------------------------" << std::endl;
		startTimer(&start);
	}
		
	_integrateRationalFunctionRationalPart<UnivariatePolynomialOverRealField,RealField>(A,D,P,G,&R,&H,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}
	
	// compute integral of R/H
	
	if (!R.isZero()) {
		
		if (PROFILING){
			std::cout << "\tRational function trancendental part" << std::endl;
			std::cout << "\t------------------------------" << std::endl;
			startTimer(&start);
		}
		
		// if I keep this routine here, I need to change the name of it to _rootsAndResidues
		_errorAnalysisLRT_rootsAndResidues<UnivariatePolynomialOverRealField,RealField>(R,H,prec,&dH,&denomRoots,&residues,&rrim,&ea_total,&outStr,PROFILING,true);
		
		if (ERROR_ANALYSIS) {	
			
			// find indices of real and complex roots of H(x)	
			for (int i=0; i<denomRoots.size(); i++) {
				if (denomRoots.at(i).imaginaryPart() == 0)
					realRootIndices.push_back(i);
				else
					complexRootIndices.push_back(i);
			}
			cilk_spawn
				_errorAnalysisPFD<UnivariatePolynomialOverRealField,RealField>(R,H,dH,denomRoots,residues,prec,realRootIndices,complexRootIndices,rrim,&maxBwdErr,&maxFwdErr, &BEIntervalWidths,&FEIntervalWidths,&BEBoundaryError,&FEBoundaryError,&ea_total,&outStr,PROFILING);
		}
		
		// this needs to also yield the needed permutation of the roots
		/*typedef std::pair<ComplexRationalNumber,ComplexRationalNumber> ResidueRootPair;
		std::vector<ResidueRootPair> rrp;
		for (int i=0; i<denomRoots.size(); ++i) {
			ResidueRootPair rp(std::make_pair(residues.at(i),denomRoots.at(i)));
			rrp.push_back(rp);
		}
		sort(rrp.begin(),rrp.end(),CompareByRealThenReverseImaginary());*/
		//sort(residues.begin(),residues.end(),CompareByRealThenReverseImaginary());
		
		_PFDStructuredIntegrate<UnivariatePolynomialOverRealField,RealField>(denomRoots,rrim,mainVariable,prec,lg,Lg,&atnTemp,&Atn1,&Atn2,PROFILING);
		
		if (PROFILING){
			stopTimer(&start,&elapsed);
			std::cout << "\t------------------------------" << std::endl;
			std::cout << "\t\t\t" << elapsed << std::endl;
			startTimer(&start);
		}
		
		// postprocessing, using Rioboo's method to remove spurious discontinuities
		for (int i=0; i<Atn1.size(); ++i) {
			_arctan2ToArctan<UnivariatePolynomialOverRealField,RealField>(Atn1.at(i),Atn2.at(i),atnTemp.at(i),atn,Atn);
		}
		
		if (PROFILING){
			stopTimer(&start,&elapsed);
			std::cout << "\t------------------------------" << std::endl;
			std::cout << "\tPost processing\t" << elapsed << std::endl;
		}
		
		if (ERROR_ANALYSIS)
			cilk_sync;
	}
		
	if (ERROR_ANALYSIS) {
		//cilk_sync;
		if (PROFILING){
			outStr << "\t------------------------------\n";
			outStr << "\t\t\t\t\t" << ea_total << "\n";
			outStr << "\t\tBackward error away from sings = " << maxBwdErr << std::endl;
			if (BEIntervalWidths.size() > 0) {
				outStr << "\t\tBackward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realRootIndices.at(i)).realPart().get_d() << " +/- " << BEIntervalWidths.at(i) << "  (" << BEBoundaryError.at(i) << ")" << std::endl;
				}
			}
			outStr << "\t\tForward Error away from sings = " << maxFwdErr << std::endl;
			if (FEIntervalWidths.size() > 0) {
				outStr << "\t\tForward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realRootIndices.at(i)).realPart().get_d() << " +/- " << FEIntervalWidths.at(i) << "  (" << FEBoundaryError.at(i) << ")" << std::endl;
				}
			}
			std::cout << outStr.str();
		}
	}
	
}

// keeping this in case it is worth having a non-structure-preserving version, though I think not
// once the mpc computations are working
template <class UnivariatePolynomialOverRealField, class RealField>
void _realSNIntegrateSimplePFD(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS){
	
	UnivariatePolynomialOverRealField Num(A);
	UnivariatePolynomialOverRealField Den(D);
	SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> supTemp;
	
	std::vector<UnivariatePolynomialOverRealField> U,dU;
	std::vector< SparseUnivariatePolynomial<UnivariatePolynomialOverRealField> > Stx,Sxt,dSdx,dSdt,d2Sdxdt;
	std::vector< std::vector<ComplexRationalNumber> > resultantRoots;
	std::vector<ComplexRationalNumber> denomRoots;
	std::vector<ComplexRationalNumber> residues;
	std::vector<RationalNumber> BEIntervalWidths,FEIntervalWidths;
	std::vector<double> BEBoundaryError,FEBoundaryError;
	std::vector<int> realSingIndices;
	double maxBwdErr,maxFwdErr;
	residueRootIndexMap rrim;
	
	// variable checking and canonicalization
	_initializeRationalIntegration<UnivariatePolynomialOverRealField,RealField>(A,D,PROFILING);
	
	// set main variable
	Symbol mainVariable;
	mainVariable = A.variable();
	
	// Profiling variables
	unsigned long long start, start2, start3;
	float elapsed(0), elapsed2(0), elapsed3(0), ea_total(0);
	std::ostringstream outStr;
	
	// compute the rational part of the integral of R/D by Hermite reduction, which is 
	// returned as G.at(0)/G.at(1), and the polynomial part, returned in P, along with the 
	// derivative of the rest of the integral, viz., the integrand of the logarithmic 
	// part, which is returned as a single rational function as a pair of polynomials 
	// R and H
	
	UnivariatePolynomialOverRealField R;  // stores numerator of log part of integral of A/D
	UnivariatePolynomialOverRealField H;  // stores squarefee denominator of the log part
	UnivariatePolynomialOverRealField dH; // stores derivative of H
	R.setVariableName(mainVariable);
	H.setVariableName(mainVariable);
	
	if (PROFILING) {
		std::cout << "\tRational function rational part" << std::endl;
		std::cout << "\t------------------------------" << std::endl;
		startTimer(&start);
	}
		
	_integrateRationalFunctionRationalPart<UnivariatePolynomialOverRealField,RealField>(A,D,P,G,&R,&H,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}
	
	// compute integral of R/H
	
	if (!R.isZero()) {
		
		if (PROFILING){
			std::cout << "\tRational function trancendental part" << std::endl;
			std::cout << "\t------------------------------" << std::endl;
			startTimer(&start);
		}
		
		// if I keep this routine here, I need to change the name of it to _rootsAndResidues
		_errorAnalysisLRT_rootsAndResidues<UnivariatePolynomialOverRealField,RealField>(R,H,prec,&dH,&denomRoots,&residues,&rrim,&ea_total,&outStr,PROFILING,true);
		
		// There are some 0 residues that are included in the list. This is only supposed to happen for non-zero
		// rational part, but they are also popping up in example 2. Genuine ones need to be removed from the list 
		// because they do not contribute to the integral, and the spurious ones should not be there at all.
		ComplexRationalNumber c,root,residue;
		RealField one(RationalNumber(1));
		RealField zero(RationalNumber(0));
		UnivariatePolynomialOverRealField LinPoly; // poly of the form a*x + b
		LinPoly.setVariableName(mainVariable);
		LinPoly = (LinPoly + one << 1);
		//LinPoly.setCoefficient(1,r);
		
		for (int i=0; i<denomRoots.size(); i++){
			root = denomRoots.at(i);
			residue = residues.at(i);
			if (denomRoots.at(i).imaginaryPart() == 0) {
				if (!epsilonEqual(residue.realPart(),zero,prec)) {
					RealField rf(-root.realPart());
					LinPoly.setCoefficient(0,rf);
					lg->push_back(residue.realPart());
					Lg->push_back(LinPoly);
				}
			}
			else {
				RealField a(root.realPart());
				RealField b(root.imaginaryPart());
				if (!epsilonEqual(residue.imaginaryPart(),zero,prec)){
					RealField rTemp(-a);
					LinPoly.setCoefficient(0,rTemp);
					rTemp = RealField(-b);
					LinPoly /= rTemp; // LinPoly = (x - a)/(-b)
					rTemp = RealField(residue.imaginaryPart()*2);
					atn->push_back(rTemp);
					Atn->push_back(LinPoly);
					LinPoly.setCoefficient(1,one);
				}
				// (x - (a+ib))(x - (a-ib)) = (x^2 - 2ax + a^2-b^2)
				if (!epsilonEqual(residue.realPart(),zero,prec)){
					UnivariatePolynomialOverRealField QuadPoly;
					QuadPoly.setVariableName(mainVariable);
					QuadPoly = (QuadPoly + one << 2); // x^2
					RealField rTemp(-2*a);
					QuadPoly.setCoefficient(1,rTemp);
					rTemp = RealField(a*a-b*b);
					QuadPoly.setCoefficient(0,rTemp);
					lg->push_back(residue.realPart());
					Lg->push_back(QuadPoly);
				}
			}
		}
		
		if (PROFILING){
			stopTimer(&start,&elapsed);
			std::cout << "\t------------------------------" << std::endl;
			std::cout << "\t\t\t" << elapsed << std::endl;
		}
	}
		
	// This needs to be moved to the optimal position
	/*if (ERROR_ANALYSIS){
		// to be replaced with _errorAnalysisPFD
		cilk_spawn
			_errorAnalysisLRT<UnivariatePolynomialOverRealField,RealField>(Stx,U,dSdx,dSdt,d2Sdxdt,dU,resultantRootsTemp,denomRoots,prec,&rrim,&residues,&realSingIndices, &BEIntervalWidths,&FEIntervalWidths,&BEBoundaryError,&FEBoundaryError,&maxBwdErr,&maxFwdErr, &ea_total,&outStr,PROFILING);
	}*/
	
	/*if (PROFILING){
		std::cout << "\tLog to real" << std::endl;
		std::cout << "\t--------------------------------------" << std::endl;
		startTimer(&start);
	}
	
	// this needs to be replaced with the PFD integration method
	//_logToReal<UnivariatePolynomialOverRealField,RealField>(resultantRoots,Stx,lg,Lg,atn,Atn,prec,PROFILING);
	
	if (PROFILING){
		stopTimer(&start,&elapsed);
		std::cout << "\t--------------------------------------" << std::endl;
		std::cout << "\t\t\t" << elapsed << std::endl;
	}*/
	
	/*if (ERROR_ANALYSIS) {
		//cilk_sync;
		if (PROFILING){
		/*	outStr << "\t------------------------------\n";
			outStr << "\t\t\t\t\t" << ea_total << "\n";
			outStr << "\t\tBackward error away from sings = " << maxBwdErr << std::endl;
			if (BEIntervalWidths.size() > 0) {
				outStr << "\t\tBackward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << BEIntervalWidths.at(i).get_d() << "  (" << BEBoundaryError.at(i) << ")" << std::endl;
				}
			}
			outStr << "\t\tForward Error away from sings = " << maxFwdErr << std::endl;
			if (FEIntervalWidths.size() > 0) {
				outStr << "\t\tForward error exceeds tolerance for x =" << std::endl;
			
				for (int i=0; i<BEIntervalWidths.size(); i++) {
					outStr << "\t\t\t" << denomRoots.at(realSingIndices.at(i)).realPart().get_d() << " +/- " << FEIntervalWidths.at(i).get_d() << "  (" << FEBoundaryError.at(i) << ")" << std::endl;
				}
			}*/
	//		std::cout << outStr.str();
	//	}
	//}
	
}

// to avoid linking errors
template void _realSNIntegrate<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, DenseUnivariateRationalPolynomial *P, std::vector<DenseUnivariateRationalPolynomial> *G, std::vector<RationalNumber> *lg, std::vector<DenseUnivariateRationalPolynomial> *Lg, std::vector<RationalNumber> *atn, std::vector<DenseUnivariateRationalPolynomial> *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template void _realSNIntegrate<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, SparseUnivariatePolynomial<RationalNumber> *P, std::vector< SparseUnivariatePolynomial<RationalNumber> > *G, std::vector<RationalNumber> *lg, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Lg, std::vector<RationalNumber> *atn, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template void _realSNIntegrate<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, DenseUnivariateRationalPolynomial *P, std::vector<DenseUnivariateRationalPolynomial> *G, std::vector<RationalNumber> *lg, std::vector<DenseUnivariateRationalPolynomial> *Lg, std::vector<RationalNumber> *atn, std::vector<DenseUnivariateRationalPolynomial> *Atn1, std::vector<DenseUnivariateRationalPolynomial> *Atn2, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template void _realSNIntegrate<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, SparseUnivariatePolynomial<RationalNumber> *P, std::vector< SparseUnivariatePolynomial<RationalNumber> > *G, std::vector<RationalNumber> *lg, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Lg, std::vector<RationalNumber> *atn, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn1, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn2, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template void _realSNIntegratePFD<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, DenseUnivariateRationalPolynomial *P, std::vector<DenseUnivariateRationalPolynomial> *G, std::vector<RationalNumber> *lg, std::vector<DenseUnivariateRationalPolynomial> *Lg, std::vector<RationalNumber> *atn, std::vector<DenseUnivariateRationalPolynomial> *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template void _realSNIntegratePFD<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, SparseUnivariatePolynomial<RationalNumber> *P, std::vector< SparseUnivariatePolynomial<RationalNumber> > *G, std::vector<RationalNumber> *lg, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Lg, std::vector<RationalNumber> *atn, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template void _realSNIntegrateSimplePFD<DenseUnivariateRationalPolynomial,RationalNumber>(DenseUnivariateRationalPolynomial &A, DenseUnivariateRationalPolynomial &D, DenseUnivariateRationalPolynomial *P, std::vector<DenseUnivariateRationalPolynomial> *G, std::vector<RationalNumber> *lg, std::vector<DenseUnivariateRationalPolynomial> *Lg, std::vector<RationalNumber> *atn, std::vector<DenseUnivariateRationalPolynomial> *Atn1, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template void _realSNIntegrateSimplePFD<SparseUnivariatePolynomial<RationalNumber>,RationalNumber>(SparseUnivariatePolynomial<RationalNumber> &A, SparseUnivariatePolynomial<RationalNumber> &D, SparseUnivariatePolynomial<RationalNumber> *P, std::vector< SparseUnivariatePolynomial<RationalNumber> > *G, std::vector<RationalNumber> *lg, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Lg, std::vector<RationalNumber> *atn, std::vector< SparseUnivariatePolynomial<RationalNumber> > *Atn1, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
