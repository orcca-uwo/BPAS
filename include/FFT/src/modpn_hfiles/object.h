/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __object_h
#define __object_h 


// ====>>  During the summer, make these stuffs into C++  <<====


//#include "cexcept.h"
//define_exception_type(int);
//extern struct exception_context the_exception_context[1];



#ifndef plong
#define plong int32
#endif

typedef union operandUnion *operand;
typedef struct SL_Graph SLG;
typedef void *Pointer;
typedef long sfixn;

typedef struct dummy_struct DUMYO;
typedef struct sfixn_struct SFIXNO;
typedef struct variable_struct VARO;
typedef struct variablePow_struct VARPOWO;
typedef struct biPlus_struct BIPLUS;
typedef struct biSub_struct BISUB;
typedef struct biProd_struct BIPROD;


#define Ntypes 6


typedef
enum type {
  t_sfixn,
  t_var,
  t_varpow,
  t_biPlus,
  t_biSub,
  t_biProd,
}operandType;


#define	type_of(oper)	((operandType)(((operand)(oper))->DUMY.type))
#define	type_set(oper, t)	((((operand)(oper)  )->DUMY.type)=t)
#define	id_of(oper)	((operandType)(((operand)(oper))->DUMY.id))
#define	id_set(oper, id)	((((operand)(oper)  )->DUMY.id)=id)

#define is_var(oper)   (type_of(oper) == t_var)
#define is_sfixn(oper)  (type_of(oper) == t_sfixn)
#define is_varpow(oper)  (type_of(oper) == t_varpow)
#define is_biPlus(oper)  (type_of(oper) == t_biPlus)
#define is_biSub(oper)  (type_of(oper) == t_biSub)
#define is_biProd(oper)  (type_of(oper) == t_biProd)


#define sfixn_val(oper) ( ((oper)->SFIX).sfixnVal )
#define var_no(oper)  ( ((oper)->VAR).no )
#define varPow_e(oper)  ( ((oper)->VARPOW).e )
#define varPow_no(oper)  ( ((oper)->VARPOW).no )
#define biPlus_oper1(oper) ( ((oper)->BI_PLUS).oper1)
#define biPlus_oper2(oper) ( ((oper)->BI_PLUS).oper2)
#define biSub_oper1(oper) ( ((oper)->BI_SUB).oper1)
#define biSub_oper2(oper) ( ((oper)->BI_SUB).oper2)
#define biProd_oper1(oper) ( ((oper)->BI_PROD).oper1)
#define biProd_oper2(oper) ( ((oper)->BI_PROD).oper2)


// dirty. compare 1 to 2.
#define is_sameVar(oper1, oper2)  ( ((oper1)->VAR).no == ((oper2)->VAR).no ) 
// dirty. compare 1 to 2.
#define is_inSameVar(oper1, oper2)  ( ((oper1)->VARPOW).no == ((oper2)->VAR).no )


#define new_poly(oper)  oper=(operand)my_malloc(sizeof(POLYO)); \
			 type_set(oper, t_poly) 

#define new_poly_ini(oper, thepoly)  oper=(operand)my_malloc(sizeof(POLYO)); \
			         type_set(oper, t_poly); \
                                  poly_poly(oper) = thepoly 

#define new_sfixn(oper)  oper=(operand)my_malloc(sizeof(SFIXNO)); \
			 type_set(oper, t_sfixn) 

#define new_sfixn_ini(oper, v)  oper=(operand)my_malloc(sizeof(SFIXNO)); \
			         type_set(oper, t_sfixn); \
                                 (oper->SFIX).sfixnVal = v

#define	new_var(oper)   oper=(operand)my_malloc(sizeof(VARO));   \
	                 type_set(oper, t_var)

#define	new_var_ini(oper, newno)   oper=(operand)my_malloc(sizeof(VARO));   \
	                          type_set(oper, t_var);\
                                  (oper->VAR).no = newno

#define	new_varpow(oper)  oper=(operand)my_malloc(sizeof(VARPOWO)); \
                          type_set(oper, t_varpow)

#define	new_varpow_ini(oper, newe, newno)  oper=(operand)my_malloc(sizeof(VARPOWO)); \
                          type_set(oper, t_varpow); \
                          (oper->VARPOW).e = newe; \
                          (oper->VARPOW).no = newno

#define	new_biPlus(oper)  oper=(operand)my_malloc(sizeof(BIPLUS)); \
                           type_set(oper, t_biPlus)

#define	new_biPlus_ini(oper, newoper1, newoper2)  oper=(operand)my_malloc(sizeof(BIPLUS)); \
                               type_set(oper, t_biPlus); \
                               (oper->BI_PROD).oper1 = newoper1; \
                               (oper->BI_PROD).oper2 = newoper2                             

#define	new_biSub(oper)   oper=(operand)my_malloc(sizeof(BISUB));  \
                           type_set(oper, t_biSub)

#define	new_biSub_ini(oper, newoper1, newoper2)   oper=(operand)my_malloc(sizeof(BISUB));  \
                               type_set(oper, t_biSub); \
                               (oper->BI_PROD).oper1 = newoper1; \
                               (oper->BI_PROD).oper2 = newoper2

#define	new_biProd(oper)  oper=(operand)my_malloc(sizeof(BIPROD));  \
                            type_set(oper, t_biProd) 


#define	new_biProd_ini(oper, newoper1, newoper2)  oper=(operand)my_malloc(sizeof(BIPROD));  \
                            type_set(oper, t_biProd); \
                            (oper->BI_PROD).oper1 = newoper1; \
                            (oper->BI_PROD).oper2 = newoper2

/* type   -- data type.
   flag   -- bit-0 liveness.
             bit-1 tmpMark (derivative, ).
             bit-2~7 reserved.
   id     -- The ID of a object (a node in G.)
*/

#define FIRSTWORD unsigned char  type, flag; unsigned short int id  

#define tmpMarkMask 0x2
#define tmpMarkUnMask 0xfd

#define	is_TmpMarkOn(oper)	((((operand)(oper))->DUMY.flag) & tmpMarkMask)
#define	clear_TmpMark(oper)	((((operand)(oper))->DUMY.flag) &= tmpMarkUnMask)
#define	set_TmpMark(oper)	((((operand)(oper))->DUMY.flag) |= tmpMarkMask)

// dummy struct.

struct dummy_struct {
	FIRSTWORD;
};


// single fixnum.

struct sfixn_struct {
	FIRSTWORD;
	sfixn	sfixnVal;
};


// Variable.

struct variable_struct {
	FIRSTWORD;
        int32   no; // no=1 -> is x.
};



// VariablePow.

struct variablePow_struct {
	FIRSTWORD;
        int32   e;  // e is the exponent .
        int32   no; // no=1 -> is x.
};


// .

struct biPlus_struct {
	 FIRSTWORD;
         operand oper1;
         operand oper2; // when Uni-operation op2 is NULL.
};

// .

struct biSub_struct {
	 FIRSTWORD;
         operand oper1;
         operand oper2; // when Uni-operation op2 is NULL.
};


// .

struct biProd_struct {
	 FIRSTWORD;
         operand oper1;
         operand oper2; // when Uni-operation op2 is NULL.
};


// operand types.
union operandUnion {
        DUMYO   DUMY;
	SFIXNO  SFIX;
        VARO    VAR;
        VARPOWO VARPOW;
        BIPLUS  BI_PLUS;
        BISUB   BI_SUB;
        BIPROD  BI_PROD;
};


// A graph G for the Strait Line input.
// The strait line is a graph with following property.
//  0) encoded in the adjacency-list style. 
//  1) G is connected or unconnected.
//  2) A node's children must appear before this node in the nodes vector.
//  3) the sum of all subgraphs in G is the complete encoding for the give polynomial.
struct SL_Graph{
  int32 GN;          // Number of Nodes  G
  int32 GE;          // Number of Edges in G
  operand * Nodes; // a vector of nodes. Typical adjacency-list encoding. 
};




#define ROOT_G(slg)  ((slg->Nodes)[(slg->GN) - 1])








// operand types.
union operand2Union {
        DUMYO   DUMY;
	SFIXNO  SFIX;
        VARO    VAR;
        VARPOWO VARPOW;
        BIPLUS  BI_PLUS;
        BISUB   BI_SUB;
        BIPROD  BI_PROD;
};



#endif
