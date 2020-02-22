
#ifndef _EXPR_TREE_NODE_
#define _EXPR_TREE_NODE_

#include <string>
#include <gmpxx.h>

//foward declaration.
class Symbol;

/**
 * The ExprTreeType enum indentifies the type of an ExprTreeNode.
 * It is a way to differentiate between data-types and operators and so on. 
 *
 * The special value EXPR_NONE is a default value to indicate that the ExprTreeNode
 * is of no type and holds no value or meaning. Generally, this default value will 
 * occur during the building of an expression tree and is not meant to be used in 
 * a completed tree.
 * 
 * EXPR_ARRAY indicates that the subtree rooted at this ExprTreeNode encodes an array.
 * The structure of this sub-tree is left-leaning such that all right children are
 * sub-trees encoding array elements and left children are again EXPR_ARRAY nodes.  
 * The array encoding continues until one EXPR_ARRAY node has a NULL left node, 
 * indicating that the array ends. The first element of the array is the one at
 * the lowest level of the tree. The last element of the array is right-child of 
 * the second-top-most ELEM_ARRAY. The top-most ELEM_ARRAY node has a NULL right-child
 * to indicate the end of the array.
 */
enum ExprTreeType {
	EXPR_NONE = 0x0,
	EXPR_VAR,
	EXPR_INT,
	EXPR_MPZ,
	EXPR_MPQ,
	EXPR_ADD,
	EXPR_SUB,
	EXPR_MULT,
	EXPR_DIV,
	EXPR_EXP,
	EXPR_NEG,
	EXPR_SYM,
	EXPR_ARRAY,
	EXPR_DATA
};

/**
 * ExprTreeVal is the data element of a ExprTreeNode. 
 * the active union member is determined by the ExprTreeNode's type.
 * The data elements here have a one-to-one correspondence to values int he 
 * ExprTreeType enum. 
 * For example EXPR_VAR <-> var, EXPR_INT <-> i, EXPR_MPZ <-> z.
 *
 * The final element void* genericData is meant to be 1) a default NULL value
 * when the ExprTreeNode does not need any data itself (as is the case with an operator)
 * or in the future when some flexible, generic data needs to be stored.  
 */
union ExprTreeVal {
	std::string* var;
	long int i;
	mpz_class* z;
	mpq_class* q;
	Symbol* sym;
	void* genericData;
};

/**
 * ExprTreeNode is a single node in the bianry tree of an ExpressionTree.
 * Nodes have a recursive structure, where each node is the root of a sub-tree
 * with edges connecting *this node to the left and right sub-trees. 
 * 
 * Each tree has a type and a val. The type is used to specify how the val should 
 * interpreted. 
 *
 * See Also, The enum ExprTreeType, and the union ExprTreeVal.
 */
struct ExprTreeNode {
	ExprTreeNode* left;
	ExprTreeNode* right;
	ExprTreeType type;
	ExprTreeVal val;

	/**
	 * Construct a default, empty ExprTreeNode.
	 */
	ExprTreeNode();

	/**
	 * Construct an ExprTreeNode given its type, value, and left and right children.
	 * Value and children can be NULL.
	 */
	ExprTreeNode(ExprTreeType type, ExprTreeVal* val, ExprTreeNode* lNode, ExprTreeNode* rNode);

	/**
	 * Construct an ExprTreeNode of type EXPR_VAR with a copy of s as data.
	 */
	ExprTreeNode(const std::string& s);

	/**
	 * Construct an ExprTreeNode of type EXPR_INT with i as data.
	 */
	ExprTreeNode(long int i);

	/**
	 * Construct an ExprTreeNode of type EXPR_MPZ with a copy of z as data.
	 */
	ExprTreeNode(const mpz_class& z);

	/**
	 * Construct an ExprTreeNode of type EXPR_MPQ with a copy of q as data.
	 */	
	ExprTreeNode(const mpq_class& q);

 	/** 
 	 * Construct an ExpeTreeNode of type EXPR_SYM with a copy of s as data.
 	 */
	ExprTreeNode(const Symbol& s);

	/**
	 * ExprTreeNode destructor.
	 */
	~ExprTreeNode();

	/**
	 * Determine if *this or any of it's children is of type searchType.
	 * 
	 * returns true iff searchType found.
	 */ 
	bool findChildType(ExprTreeType searchType);

	/**
	 * Covnert *this to a string, depending on it's type and data.
	 * Only convert this particular node, not it's children. This allows
	 * for different tree traversals.
	 */
	std::string toString() const;


	/**
	 * Check to see if the expression rooted at *this encodes a variable
	 * or a variable to a power.
	 */
	bool isVar() const;

	/**
	 * Check to see if the expression rooted at *this encodes a constant. 
	 * That is, an integer, a rational number, etc.
	 */
	bool isConstant() const;

	/**
	 * Check to see if the expression rooted at *this encodes a monomial.
	 * 
	 * returns true iff *this encodes a monomial.
	 */
	bool isPolynomialTerm() const;

	/**
	 * Get a deep copy of the tree rooted at root.
	 *
	 * Node that if the node or any its children hold data
	 * as a void* then that data will NOT be deeply copied.
	 *
	 * returns the new root
	 */
	ExprTreeNode* deepCopy() const;

	///////////////////
	// Static Helpers
	///////////////////

	/**
	 * Creates a new ExprTreeNode and combines two sub-trees (lNode and rNode) 
	 * by adjoining them to the newly created node as their root. 
	 * The newly created node has type tType and val vVal (if vVal is not NULL). 
	 */
	static ExprTreeNode* combineExprTreeNodes(ExprTreeNode* lNode, ExprTreeNode* rNode, ExprTreeType tType, ExprTreeVal* vVal = NULL);
};

#endif
