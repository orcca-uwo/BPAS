
#include "SubResultantChain/SRC_SMQP.hpp"
#include "Utils/SymbolHelpers.hpp"

float SMQP_SRC_TIME = .0f;

SMQPSubResultantChain::SubResultantChain (const SMQP& a, const SMQP& b, const Symbol& v) : P(a), Q(b), var(v) {
	if (a.leadingVariable() != v && b.leadingVariable() != v) {
		throw std::invalid_argument("In SubResultantChain::SubResultantChain at least one of the input polynomials must have " + v.toString() + " as leading variable.");
	}

	if (P.degree(v) < Q.degree(v)) {
		SMQP tmp = std::move(Q);
		Q = std::move(P);
		P = std::move(tmp);
	}

	if (Q.isZero()) {
		chain = std::vector<SMQP>(2);
		chain[0] = Q;
		chain[1] = P;
		valid = std::vector<bool>(2, true);
	} else {
		chain = std::vector<SMQP>(Q.degree(v).get_ui() + 2u);
		chain[chain.size()-1] = P;
		chain[chain.size()-2] = Q;

		valid = std::vector<bool>(chain.size(), false);

		valid[chain.size()-1] = true;
		valid[chain.size()-2] = true;
	}

	chainCoefs = std::vector<SMQP>(chain.size());
	validCoefs = std::vector<bool>(chain.size(), false);

	// std::cerr << "Computing Chain for " << var << " Between: " << std::endl << "P: " << P << std::endl << "Q: " << Q << std::endl;

	std::vector<Symbol> pVars = P.variables();
	std::vector<Symbol> qVars = Q.variables();
	std::vector<Symbol> both = orderPreservingSetUnion(pVars, qVars);
	if (both.size() > 1) {
		fillChain();
	}
}

SMQPSubResultantChain::SubResultantChain (const SMQP& a, const SMQP& b, const Symbol& v, std::vector<Symbol> ringVars) : SMQPSubResultantChain(a,b,v) {}

bool SMQPSubResultantChain::operator==(SMQPSubResultantChain& a) {
	return (P == a.P && Q == a.Q && var == a.var);
}

std::vector<SMQP> SMQPSubResultantChain::polynomials() const {

	fillChain();

	return chain;
}

void SMQPSubResultantChain::fillChain() const {

	bool allValid = true;
	for (int i = 0; i < valid.size() && allValid; ++i) {
			allValid = valid[i];
	}
	if (allValid) {
		return;
	}

	chain = P.subresultantChain(Q, var);
	valid = std::vector<bool>(chain.size(), true);
}

SMQP SMQPSubResultantChain::subResultantOfIndex(size_t i, bool lazy) const {

	if (!chain.empty() && i < valid.size() && valid[i]) {
		// std::cerr << "retruning subResultantOfIndex[" << i << "]: " << chain[i];
		return chain[i];
	}

	Integer idx(i);
	if (idx > P.degree(var) || i < 0) {
		std::cerr << "P : " << P << std::endl;
		std::cerr << "Q : " << Q << std::endl;
		std::cerr << "var : " << var << std::endl;
		std::cerr << "i : " << i << std::endl;
		std::cerr << *this << std::endl;
		std::cerr << "BPAS: SMQPSubResultantChain ERROR: requested subresultant does not exist!" << std::endl;
		exit(1);
	}

	if (idx == P.degree(var) || idx == Q.degree(var) + 1) {
		return P;
	}

	if (idx == Q.degree(var)) {
		return Q;
	}

	if (i >= chain.size()) {
		SMQP ret;
		ret.zero();
		return ret;
	}


	size_t mdegQ = Q.degree(var).get_ui();
	std::vector<SMQP> tmp;
	if (!lazy && !validCoefs[0]) {
		std::vector<SMQP> tmpCoefs;
		tmp = P.subresultantChainAtIdx(Q, var, i, &(tmpCoefs));
		for (size_t j = 0; j < tmpCoefs.size() && i+j < mdegQ; ++j) {
			//we don't store values between Q and P.
			// fprintf(stderr, "chainCoefs size: %d\n", chainCoefs.size());
			chainCoefs[i+j] = std::move(tmpCoefs[j]);
			validCoefs[i+j] = true;
			// std::cerr << "princoef[" << j << "]: " << chainCoefs[j] << std::endl;
		}
	} else {
		tmp = P.subresultantChainAtIdx(Q, var, i);
	}

	for (size_t j = 0; j < tmp.size(); ++j) {
		chain[i + j] = std::move(tmp[j]);
		valid[i+j] = true;
	}

	return chain[i];
}

SMQP SMQPSubResultantChain::principleSubResultantCoefficientOfIndex(int i) const {

	Integer idx(i);
	if (idx > P.degree(var) || i < 0) {
		std::cerr << "BPAS: SMQPSubResultantChain ERROR: requested subresultant does not exist!" << std::endl;
		exit(1);
	}

	if (idx == P.degree(var) || idx == Q.degree(var) + 1) {
		if (!validCoefs[i]) {
			chainCoefs[i] = P.leadingCoefficientInVariable(var);
			validCoefs[i] = true;
		}
		return chainCoefs[i];
	}

	if (idx == Q.degree(var)) {
		if (!validCoefs[i]) {
			chainCoefs[i] = Q.leadingCoefficientInVariable(var);
			validCoefs[i] = true;
		}
		return chainCoefs[i];
	}

	if (i >= chain.size()) {
		SMQP ret;
		ret.zero();
		return ret;
	}

	if (!chainCoefs.empty() && i < validCoefs.size() && validCoefs[i]) {
		return chainCoefs[i];
	}

	if (!chain.empty() && i < valid.size() && valid[i]) {
		chainCoefs[i] = chain[i].leadingCoefficientInVariable(var);
		validCoefs[i] = true;
		return chainCoefs[i];
	}

	this->subResultantOfIndex(i);

	if (!validCoefs[i]) {
		chainCoefs[i] = chain[i].leadingCoefficientInVariable(var);
	}
	return chainCoefs[i];
}



SMQP SMQPSubResultantChain::resultant(bool lazy) const {
	std::vector<bool> localValid = valid;
	std::vector<SMQP> localChain = chain;

	if (!chain.empty() && valid[0]) {
		return chain[0];
	}


	//TODO: remove check on nvars once we have better resultant methods
	std::vector<Symbol> pVars = P.variables();
	std::vector<Symbol> qVars = Q.variables();
	std::vector<Symbol> both = orderPreservingSetUnion(pVars, qVars);
	
	if (lazy && both.size() == 1) {
		chain[0] = P.resultant(Q, var);
		valid[0] = true;
	} else {
		//this actually computes chain[0] and chain[1].
		this->subResultantOfIndex(0);
	}


	return chain[0];
}
		

void SMQPSubResultantChain::print(std::ostream& out) const {
	// fillChain();

	bool isNotFirst = 0;
	out << "var: " << var << std::endl;
	out << "[" << P << ",\n" << Q;
	for (size_t i = chain.size()-2; i-- > 0; ) {
		out << ",\n ";
		out << chain[i];
	}
	out << "]";
}
	

/**
 * Convert subresultant chain to an expression tree.
 *
 * @param
 **/
ExpressionTree SMQPSubResultantChain::convertToExpressionTree() const {
	if (!chain.size()) {
		ExprTreeNode* node = new ExprTreeNode(EXPR_ARRAY, NULL, NULL, NULL);
		return ExpressionTree(node);
	}
	else {
		ExpressionTree t;
		t.fromVector<SMQP>(polynomials());
		return t;
	}
}