
#ifndef _EXPRESSION_TREE_HPP_
#define _EXPRESSION_TREE_HPP_

#include <string>
#include <gmpxx.h>
#include <vector>
#include "ExprTreeNode.hpp"

/**
 * An ExpressionTree encompasses various forms of data that can be expressed generically
 * as a binary tree with data elements connected by operators. 
 *
 * Provides functions for converting this generic datatype into various external formats. 
 * For example, maple expression strings or latex mathmode strings.
 */
class ExpressionTree {
private:
	ExprTreeNode* root;

public: 

	/**
	 * Construct the default, empty tree.
	 */
	ExpressionTree();

	/** 
	 * Construct a tree given the root node.
	 */
	ExpressionTree(ExprTreeNode* theRoot);

	/**
	 * Copy constructor.
	 */
	ExpressionTree(const ExpressionTree& otherTree);

	/**
	 * Construct a stree from a vector of trees as an array of trees.
	 */
	ExpressionTree(const std::vector<ExpressionTree>& trees);

	/**
	 * Destructor.
	 */
	~ExpressionTree();

	/**
	 * Copy assignment.
	 */
	ExpressionTree& operator= (const ExpressionTree& otherTree);

#if __GNUC__ >= 5
	
	/**
	 * Move constructor.
	 */
	ExpressionTree(ExpressionTree&& otherTree); 

	/**
	 * Move assignment.
	 */
	ExpressionTree& operator= (ExpressionTree&& otherTree);

#endif

	/** 
	 * Combine two expression trees using the supplied type and value as their common 
	 * parent. This method will create copies of the input trees to use. 
	 * 
	 * returns the adjoined ExpressionTree.
	 */
	static ExpressionTree combineExpressionTrees(const ExpressionTree& left, const ExpressionTree& right, ExprTreeType type, ExprTreeVal* val = NULL);

	/** 
	 * Combine *this and the input rightTree into a new expression tree
	 * adjoining the two trees with addition.
	 *
	 * returns the new tree.
	 */
	ExpressionTree operator+(const ExpressionTree& rightTree);

	/** 
	 * Add the input rightTree to the right side of *this, adjoining the two
	 * trees with addition.
	 *
	 * returns a reference to *this, the updated tree.
	 */
	ExpressionTree& operator+=(const ExpressionTree& rightTree);

	/** 
	 * Combine *this and the input rightTree into a new expression tree
	 * adjoining the two trees with subtraction.
	 *
	 * returns the new tree.
	 */
	ExpressionTree operator-(const ExpressionTree& rightTree);

	/** 
	 * Add the input rightTree to the right side of *this, adjoining the two
	 * trees with subtraction.
	 *
	 * returns a reference to *this, the updated tree.
	 */
	ExpressionTree& operator-=(const ExpressionTree& rightTree);

	/** 
	 * Combine *this and the input rightTree into a new expression tree
	 * adjoining the two trees with multiplication.
	 *
	 * returns the new tree.
	 */
	ExpressionTree operator*(const ExpressionTree& rightTree);

	/** 
	 * Add the input rightTree to the right side of *this, adjoining the two
	 * trees with multiplication.
	 *
	 * returns a reference to *this, the updated tree.
	 */
	ExpressionTree& operator*=(const ExpressionTree& rightTree);

	/** 
	 * Combine *this and the input rightTree into a new expression tree
	 * adjoining the two trees with divison.
	 *
	 * returns the new tree.
	 */
	ExpressionTree operator/(const ExpressionTree& rightTree);

	/** 
	 * Add the input rightTree to the right side of *this, adjoining the two
	 * trees with division.
	 *
	 * returns a reference to *this, the updated tree.
	 */
	ExpressionTree& operator/=(const ExpressionTree& rightTree);

	/** 
	 * Combine *this and the input rightTree into a new expression tree
	 * adjoining the two trees with divison.
	 *
	 * returns the new tree.
	 */
	ExpressionTree operator^(const ExpressionTree& rightTree);

	/** 
	 * Add the input rightTree to the right side of *this, adjoining the two
	 * trees with division.
	 *
	 * returns a reference to *this, the updated tree.
	 */
	ExpressionTree& operator^=(const ExpressionTree& rightTree);

	/** 
	 * Convert *this to a generic string representation.
	 */
	std::string toString() const;

	/**
	 * Convert *this to a string in the format expected of a maple expression.
	 *
	 * returns a string in maple format
	 */
	std::string toMapleString() const;

	/**
	 * Convert *this to a string in the format expected of LaTeX.
	 *
	 * returns a string in latex format.
	 */
	std::string toLaTeXString() const;
	
	/**
	 * Fill this expression tree with the vector of BPASRing elements specified.
	 * This clears this tree's current contents.
	 * ringVec: a vector of BPASRing elements. 
	 */
	template <class ExpTreeConvert> 
	void fromVector(const std::vector<ExpTreeConvert>& ringVec) {
		delete root;
		root = NULL;
		for(size_t i = 0; i < ringVec.size(); ++i) {
			ExpressionTree elem = ringVec[i].convertToExpressionTree();
			root = ExprTreeNode::combineExprTreeNodes(root, elem.root, EXPR_ARRAY, NULL);
			elem.root = NULL;
		}
		root = ExprTreeNode::combineExprTreeNodes(root, NULL, EXPR_ARRAY, NULL);
	}
};

/**
 * An interface defining conversion of a class to an ExpressionTree. 
 */
class ExpressionTreeConvert {

public:
 	/**
     * Convert this to an expression tree.
     *
     * returns an expression tree describing *this.
     */
	virtual ExpressionTree convertToExpressionTree() const = 0;
};


#endif
