
#include "ExpressionTree/ExpressionTree.hpp"
#include <sstream>

ExpressionTree::ExpressionTree() : root(NULL) {

}

ExpressionTree::ExpressionTree(ExprTreeNode* theRoot) : root(theRoot) {

}

ExpressionTree::ExpressionTree(const ExpressionTree& otherTree) {
	root = otherTree.root->deepCopy();
}

ExpressionTree::~ExpressionTree() {
	delete root;
}

ExpressionTree::ExpressionTree(const std::vector<ExpressionTree>& trees) : root(NULL) {
	for (size_t i = 0; i < trees.size(); ++i) {
		ExpressionTree temp = trees[i];
		root = ExprTreeNode::combineExprTreeNodes(root, temp.root, EXPR_ARRAY, NULL);
		temp.root = NULL;
	}
	root = ExprTreeNode::combineExprTreeNodes(root, NULL, EXPR_ARRAY, NULL);
}

ExpressionTree& ExpressionTree::operator= (const ExpressionTree& otherTree) {
	if (this != &otherTree) {
		delete root;
		this->root = otherTree.root->deepCopy();
	}
	return *this;
}

#if __GNUC__ >= 5

ExpressionTree::ExpressionTree(ExpressionTree&& otherTree) {
	root = otherTree.root;
	otherTree.root = NULL;
}

ExpressionTree& ExpressionTree::operator= (ExpressionTree&& otherTree) {
	if (this != &otherTree) {
		delete root;
		this->root = otherTree.root;
		otherTree.root = NULL;
	}
	return *this;
}

#endif

ExpressionTree ExpressionTree::combineExpressionTrees(const ExpressionTree& left, const ExpressionTree& right, ExprTreeType type, ExprTreeVal* val) {
	ExprTreeNode* leftTree = left.root->deepCopy();
	ExprTreeNode* rightTree = right.root->deepCopy();
	ExprTreeNode* newRoot = ExprTreeNode::combineExprTreeNodes(leftTree, rightTree, type, val);
	return ExpressionTree(newRoot);
}

ExpressionTree ExpressionTree::operator+ (const ExpressionTree& rightTree) {
	return ExpressionTree::combineExpressionTrees(*this, rightTree, EXPR_ADD);
}

ExpressionTree& ExpressionTree::operator+= (const ExpressionTree& rightTree) {
	ExprTreeNode* rightNode = rightTree.root->deepCopy();
	root = ExprTreeNode::combineExprTreeNodes(root, rightNode, EXPR_ADD);
	return *this;
}

ExpressionTree ExpressionTree::operator- (const ExpressionTree& rightTree) {
	return ExpressionTree::combineExpressionTrees(*this, rightTree, EXPR_SUB);
}

ExpressionTree& ExpressionTree::operator-= (const ExpressionTree& rightTree) {
	ExprTreeNode* rightNode = rightTree.root->deepCopy();
	root = ExprTreeNode::combineExprTreeNodes(root, rightNode, EXPR_SUB);
	return *this;
}

ExpressionTree ExpressionTree::operator* (const ExpressionTree& rightTree) {
	return ExpressionTree::combineExpressionTrees(*this, rightTree, EXPR_MULT);
}

ExpressionTree& ExpressionTree::operator*= (const ExpressionTree& rightTree) {
	ExprTreeNode* rightNode = rightTree.root->deepCopy();
	root = ExprTreeNode::combineExprTreeNodes(root, rightNode, EXPR_MULT);
	return *this;
}

ExpressionTree ExpressionTree::operator/ (const ExpressionTree& rightTree) {
	return ExpressionTree::combineExpressionTrees(*this, rightTree, EXPR_DIV);
}

ExpressionTree& ExpressionTree::operator/= (const ExpressionTree& rightTree) {
	ExprTreeNode* rightNode = rightTree.root->deepCopy();
	root = ExprTreeNode::combineExprTreeNodes(root, rightNode, EXPR_DIV);
	return *this;
}

ExpressionTree ExpressionTree::operator^(const ExpressionTree& rightTree) {
	return ExpressionTree::combineExpressionTrees(*this, rightTree, EXPR_EXP);
}

ExpressionTree& ExpressionTree::operator^=(const ExpressionTree& rightTree) {
	ExprTreeNode* rightNode = rightTree.root->deepCopy();
	root = ExprTreeNode::combineExprTreeNodes(root, rightNode, EXPR_EXP);
	return *this;
}

std::string ExpressionTree::toString() const {
	if (root == NULL) {
		return "";
	}
	return toMapleString();
}

//Forward declaration;
void nodeToMapleString(ExprTreeNode* node, std::stringstream* ss);

void arrayNodeToMapleString(ExprTreeNode* node, std::stringstream* ss) {
	if (node->type != EXPR_ARRAY) {
		return;
	}

	if (node->left == NULL) {
		(*ss) << "[";
		if (node->right == NULL) {
			(*ss) << "]";
		} else {
			nodeToMapleString(node->right, ss);
		}
		return;
	}

	arrayNodeToMapleString(node->left, ss);
	if (node->right != NULL) {
		(*ss) << ", ";
		nodeToMapleString(node->right, ss);
	} else {
		(*ss) << "]";
	}
}

void nodeToMapleString(ExprTreeNode* node, std::stringstream* ss) {
	if (node->type == EXPR_ARRAY) {
		arrayNodeToMapleString(node, ss);
		return;
	}

	if (node->left != NULL) {
		(*ss) << "(";
		nodeToMapleString(node->left, ss);
		(*ss) << ")";
	}
	(*ss) << node->toString();
	if (node->right != NULL) {
		(*ss) << "(";
		nodeToMapleString(node->right, ss);
		(*ss) << ")";
	}
}

std::string ExpressionTree::toMapleString() const {
	if (root == NULL) {
		return "";
	}
	std::stringstream ss;
	nodeToMapleString(root, &ss);
	std::string ret = ss.str();
	return ret;
}

//Forward declaration;
void nodeToLaTeXString(ExprTreeNode* node, std::stringstream* ss);

void arrayNodeToLaTeXString(ExprTreeNode* node, std::stringstream* ss) {
	if (node->type != EXPR_ARRAY) {
		return;
	}

	if (node->left == NULL) {
		(*ss) << "\\left[";
		if (node->right == NULL) {
			(*ss) << "\\right]";
		} else {
			nodeToLaTeXString(node->right, ss);
		}
		return;
	}

	arrayNodeToLaTeXString(node->left, ss);
	if (node->right != NULL) {
		(*ss) << ", ";
		nodeToLaTeXString(node->right, ss);
	} else {
		(*ss) << "\\right]";
	}
}


void nodeToLaTeXString(ExprTreeNode* node, std::stringstream* ss) {
	if (node->type == EXPR_ARRAY) {
		arrayNodeToLaTeXString(node, ss);
		return;
	}

	if(node->left != NULL) {
		nodeToLaTeXString(node->left, ss);
	}
	switch (node->type) {
		case EXPR_ADD: {
			(*ss) << "\\ " << node->toString() << "\\ ";
			if (node->right != NULL) {
				nodeToLaTeXString(node->right, ss);
			}
			break;
		}
		case EXPR_SUB:
		case EXPR_DIV: {
			(*ss) << "\\ " << node->toString() << "\\ \\left(";
			if (node->right != NULL) {
				nodeToLaTeXString(node->right, ss);
			}
			(*ss) << "\\right)";
			break;
		}
		case EXPR_MULT: {
			if (node->isPolynomialTerm()) {
				//if poly term, exclude multiplication symbol;
				nodeToLaTeXString(node->right, ss);
			} else if (node->right != NULL) {
				(*ss) << "\\ " << node->toString() << "\\ ";
				nodeToLaTeXString(node->right, ss);
			}
			break;
		}
		case EXPR_EXP: {
			(*ss) << node->toString();
			(*ss) << "{";
			if (node->right != NULL) {
				nodeToLaTeXString(node->right, ss);
				(*ss) << "}";
			}
			break;
		}
		default : {
			(*ss) << node->toString();
			if (node->right != NULL) {
				nodeToLaTeXString(node->right, ss);
			}
			break;
		}
	}
}

std::string ExpressionTree::toLaTeXString() const {
	if (root == NULL) {
		return "";
	}

	std::stringstream ss;
	nodeToLaTeXString(root, &ss);
	std::string ret = ss.str();
	return ret;
}
