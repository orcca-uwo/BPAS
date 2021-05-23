
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
	chainDegs = std::vector<int>(chain.size(), 0);

	for (size_t i = 0; i < chain.size(); ++i) {
		if (valid[i]) {
			chainCoefs[i] = chain[i].leadingCoefficientInVariable(v);
			chainDegs[i] = chain[i].degree(v).get_si();
			validCoefs[i] = true;
		}
	}

	// std::cerr << "Computing Chain for " << var << " Between: " << std::endl << "P: " << P << std::endl << "Q: " << Q << std::endl;
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

	// bool allValid = true;
	// for (int i = 0; i < valid.size() && allValid; ++i) {
	// 		allValid = valid[i];
	// }
	// if (allValid) {
	// 	return;
	// }
	// chain = P.subresultantChain(Q, var);
	// valid = std::vector<bool>(chain.size(), true);

	for (size_t i = 0; i < chain.size(); ++i) {
		if (!valid[i]) {
			this->subResultantOfIndex(i); //compute and cache S_i and S_{i+1}.
		}
	}
}

SMQP SMQPSubResultantChain::subResultantOfIndex(size_t i) const {

	if (!chain.empty() && i < valid.size() && valid[i]) {
		// std::cerr << "retruning subResultantOfIndex[" << i << "]: " << chain[i] << std::endl;
		return chain[i];
	}

	Integer idx(i);
	Integer degP = P.degree(var);
	Integer degQ = Q.degree(var);
	if (idx > P.degree(var) || i < 0) {
		std::cerr << "P : " << P << std::endl;
		std::cerr << "Q : " << Q << std::endl;
		std::cerr << "var : " << var << std::endl;
		std::cerr << "i : " << i << std::endl;
		std::cerr << *this << std::endl;
		std::cerr << "BPAS: SMQPSubResultantChain ERROR: requested subresultant does not exist!" << std::endl;
		exit(1);
	}

	if (idx == degP || idx == degQ + 1) {
		return P;
	}

	if (idx == degQ) {
		return Q;
	}

	if (i >= chain.size()) {
		SMQP ret;
		ret.zero();
		return ret;
	}

	std::vector<SMQP> tmp = P.subresultantAtIdx(Q, var, i, &(this->lazyInfo));
	// std::cerr << "idx: " << i << " j: " << 0 << " chain[i+j] = " << tmp[0] << std::endl << std::endl;
	chain[i] = std::move(tmp[0]);
    if (!validCoefs[i]) {
		chainCoefs[i] = chain[i].leadingCoefficientInVariable(var);
		chainDegs[i] = chain[i].degree(var).get_si();
		validCoefs[i] = true;
	}

	size_t mdeg = tmp[1].degree(var).get_si();
    if (idx < degQ) {
        for (size_t j = 1; j < mdeg-i; ++j) {
            chain[i+j].zero();
            valid[i+j] = true;
            if (!validCoefs[i+j]) {
				chainCoefs[i+j].zero();
				chainDegs[i+j] = 0;
				validCoefs[i+j] = true;
			}
        }
    }
	// std::cerr << "idx: " << i << " j: " << mdeg-i << " chain[i+j] = " << tmp[1] << std::endl << std::endl;

    chain[mdeg] = std::move(tmp[1]);
    valid[i] = true;
    valid[mdeg] = true;
    chainCoefs[mdeg] = chain[mdeg].leadingCoefficientInVariable(var);
	chainDegs[mdeg] = chain[mdeg].degree(var).get_si();
	validCoefs[mdeg] = true;


	// for (size_t j = 0; j < tmp.size(); ++j) {
	// 	std::cerr << "idx: " << i << " j: " << j << " chain[i+j] = " << tmp[j] << std::endl << std::endl;
	// 	chain[i + j] = std::move(tmp[j]);
	// 	valid[i+j] = true;

		// if (!validCoefs[i+j]) {
		// 	fprintf(stderr, "i: %d, j: %d, chainDegs size: %d, chainCoefs size: %d, validCoefs size: %d", i, j, chainDegs.size(), chainCoefs.size(), validCoefs.size());
		// 	chainCoefs[i+j] = chain[i+j].leadingCoefficientInVariable(var);
		// 	chainDegs[i+j] = chain[i+j].degree(var).get_si();
		// 	validCoefs[i+j] = true;
		// }
	// }

	return chain[i];
}

SMQP SMQPSubResultantChain::principalSubResultantCoefficientOfIndex(size_t i) const {
	Integer idx(i);
	if (idx > P.degree(var) || i < 0) {
		std::cerr << "BPAS: SMQPSubResultantChain ERROR: requested subresultant does not exist!" << std::endl;
		exit(1);
	}

	if (i >= chain.size()) {
		SMQP ret;
		ret.zero();
		return ret;
	}

	//compute if not yet found
	if (!validCoefs[i]) {
		this->subResultantInitialOfIndex(i);
	}

	if ((int) i == chainDegs[i]) {
		return chainCoefs[i];
	} else {
		SMQP ret;
		ret.zero();
		return ret;
	}
}


SMQP SMQPSubResultantChain::subResultantInitialOfIndex(size_t i) const {
	Integer idx(i);
	if (idx > P.degree(var) || i < 0) {
		std::cerr << "BPAS: SMQPSubResultantChain ERROR: requested subresultant does not exist!" << std::endl;
		exit(1);
	}

	if (idx == P.degree(var) || idx >= Q.degree(var) + 1) {
		if (!validCoefs[i]) {
			chainCoefs[i] = P.leadingCoefficientInVariable(var);
			chainDegs[i] = P.degree(var).get_si();
			validCoefs[i] = true;
		}
		return chainCoefs[i];
	}

	if (idx == Q.degree(var)) {
		if (!validCoefs[i]) {
			chainCoefs[i] = Q.leadingCoefficientInVariable(var);
			chainDegs[i] = Q.degree(var).get_si();
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
		chainDegs[i] = chain[i].degree(var).get_si();
		validCoefs[i] = true;
		return chainCoefs[i];
	}

	Integer degi, degi1;
	std::vector<SMQP> tmp = P.subresultantInitialAtIdx(Q, var, i, degi, degi1, &(this->lazyInfo));

	chainCoefs[i] = std::move(tmp[0]);
	chainDegs[i] = degi.get_si();
	validCoefs[i] = true;

	size_t mdeg = degi1.get_si();
    if (idx < Q.degree(var)) {
        for (size_t j = 1; j < mdeg-i; ++j) {
            chainCoefs[i+j].zero();
			chainDegs[i+j] = 0;
            validCoefs[i+j] = true;
        }
    }
	// std::cerr << "idx: " << i << " j: " << mdeg-i << " chain[i+j] = " << tmp[1] << std::endl << std::endl;

    chainCoefs[mdeg] = std::move(tmp[1]);
    chainDegs[mdeg] = mdeg;
	validCoefs[mdeg] = true;

	return chainCoefs[i];
}



SMQP SMQPSubResultantChain::resultant(bool lazy) const {
	if (!chain.empty() && valid[0]) {
		return chain[0];
	}

	// fprintf(stderr, "inited3 the vectors: coefs: %d, valid: %d, degs: %d\n", chainCoefs.size(), validCoefs.size(), chainDegs.size());

	if (lazy) {
		chain[0] = P.resultant(Q, var);
		valid[0] = true;

		chainCoefs[0] = chain[0].leadingCoefficientInVariable(var);
		// chainDegs[0] = chain[0].degree(var).get_si();
		validCoefs[0] = true;
	} else {
		//this actually computes chain[0] and chain[1].
		this->subResultantOfIndex(0);
	}


	return chain[0];
}



void SMQPSubResultantChain::print(std::ostream& out) const {
	fillChain();

	bool isNotFirst = 0;
	out << "var: " << var << std::endl;
	out << "[" << P << ",\n" << Q;
	if (!Q.isZero()) {
		for (size_t i = chain.size()-2; i-- > 0; ) {
			out << ",\n ";
			out << chain[i];
		}
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
