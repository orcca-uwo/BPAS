#ifndef _RATIONALFUNCTION_SYMBOLICNUMERICINTEGRATION_H
#define _RATIONALFUNCTION_SYMBOLICNUMERICINTEGRATION_H

template <class UnivariatePolynomialOverRealField, class RealField>
extern void _realSNIntegrate(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template <class UnivariatePolynomialOverRealField, class RealField>
extern void _realSNIntegrate(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn1, std::vector<UnivariatePolynomialOverRealField> *Atn2, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template <class UnivariatePolynomialOverRealField, class RealField>
extern void _realSNIntegratePFD(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);
template <class UnivariatePolynomialOverRealField, class RealField>
extern void _realSNIntegrateSimplePFD(UnivariatePolynomialOverRealField &A, UnivariatePolynomialOverRealField &D, UnivariatePolynomialOverRealField *P, std::vector<UnivariatePolynomialOverRealField> *G, std::vector<RealField> *lg, std::vector<UnivariatePolynomialOverRealField> *Lg, std::vector<RealField> *atn, std::vector<UnivariatePolynomialOverRealField> *Atn1, int prec, bool PROFILING, bool PFD_LOG_PART, bool ERROR_ANALYSIS);

#endif
