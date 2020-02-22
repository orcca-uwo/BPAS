/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __Types_h
#define __Types_h 


#define TRY64 1

#define WINDOWS 1

#define DEBUG 0
#define DEBUG1 0
#define DEBUG2 0
#define DEBUG3 0
#define DEBUG4 1
#define DEBUG5 0
#define DEBUG20 0
#define DEBUG21 0
#define DEBUG22 0
#define DEBUG23 0
#define DEBUG24 0
#define DEBUG25 0
#define DEBUG30 0
#define DEBUG31 0
#define DEBUG33 0
#define DEBUG35 0
#define DEBUG38 0
#define DEBUG39 0
#define NoOfCPU 4

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cexcept.h"
#include <signal.h>

//#ifdef WINDOWS
//#define my_malloc(A) my_aligned_malloc(A, 128)
//#define my_calloc(A, B) my_aligned_calloc(A, B, 128)
//#define my_free(A) _aligned_free(A);
//#endif



#ifdef _MSC_VER
#  define CDECL __cdecl
#else
#  define CDECL
#endif








define_exception_type(int);
extern struct exception_context the_exception_context[1];




#ifdef LINUXINTEL32
#undef WINDOWS
//typedef long  sfixn;
typedef long  sfixn;
typedef long long  longfixnum;
//typedef unsigned long  usfixn;
typedef unsigned long  usfixn;
typedef unsigned long long ulongfixnum;
#endif

#ifdef TRY64

#ifdef LINUXINTEL64
#undef WINDOWS
typedef __int64_t    sfixn;
typedef __int128_t   longfixnum;
typedef __uint64_t   usfixn;
typedef __uint128_t  ulongfixnum;
typedef __int32_t  int32;
typedef __uint32_t  uint32;
#endif


#else

#ifdef LINUXINTEL64
#undef WINDOWS
//typedef long  sfixn;
typedef int  sfixn;
typedef long  longfixnum;
//typedef unsigned long  usfixn;
typedef unsigned int  usfixn;
typedef unsigned long ulongfixnum;
#endif

#endif


#ifdef MAC32
#undef WINDOWS
//typedef long  sfixn;
typedef int32_t  sfixn;
typedef int64_t  longfixnum;

//typedef unsigned long  usfixn;
typedef uint32_t  usfixn;
typedef uint64_t ulongfixnum;
#endif

#ifdef SOLARIS64
#undef WINDOWS
//typedef long  sfixn;
typedef int32_t  sfixn;
typedef int64_t  longfixnum;

//typedef unsigned long  usfixn;
typedef uint32_t  usfixn;
typedef uint64_t ulongfixnum;
#endif


#ifdef PPC64
#undef WINDOWS
//typedef long  sfixn;
typedef int32_t  sfixn;
typedef int64_t  longfixnum;

//typedef unsigned long  usfixn;
typedef uint32_t  usfixn;
typedef uint64_t ulongfixnum;
#endif


#ifdef MAC64
#undef WINDOWS
//typedef long  sfixn;
typedef int32_t  sfixn;
typedef int64_t  longfixnum;

//typedef unsigned long  usfixn;
typedef uint32_t  usfixn;
typedef uint64_t ulongfixnum;
#endif


#ifdef WINDOWS
//typedef long  sfixn;
typedef _int32  sfixn; 
typedef _int64  longfixnum;
 
//typedef unsigned long  usfixn;
typedef unsigned _int32  usfixn;
typedef unsigned _int64 ulongfixnum;

#endif





 /**
 * preFFTRep: "cube" representation of multivariate dense  
 * @N: number of variables = number of dimensions in the cube
 * @defN: backup of N 
 * @bufSizs: 1-dim array of length N+1.
 *           bufSizs[i]+1 is the size of the buffer at dimension i.
 *           bufSizs[0] is not used.
 * @cuts:   1-dim array of length N+1. Gives the partial degrees.
 * @accum:  1-dim array of length N+1.
 *          accum[i] is the size of a coefficient w.r.t. the i-th variable.
 * @size:   the size of the coefficient buffer.
 * @defSize:  backup of size
 * @offset: shold be zero
 * @data: coefficient buffer
 * @tmpData: backup for the data pointer
 * @defData: backup for the data pointer
 * 
 * The elements of the cube are coefficients (not values).
 * 
 * 
 * Return value: 
 **/   
typedef struct preFFTRepST{
  sfixn N; 
  sfixn defN;
  sfixn * bufSizs; // bufSizes is a erroric name. It should be degrees.
                   // I.e. bufSizs is a vector keeps the partial degrees of 
                   // the give polynomias
  sfixn * cuts;  
  sfixn * accum;
  sfixn size;  
  sfixn datasize;
  sfixn defSize;
  sfixn offset;
  sfixn * data;
  sfixn * tmpData; 
  sfixn * defData; 
} preFFTRep;

#define N(Ptr)              (Ptr->N)
#define DFN(Ptr)            (Ptr->defN)
#define BUSZSI(Ptr, i)      ((Ptr->bufSizs)[i])
#define BUSZS(Ptr)          (Ptr->bufSizs)
#define BUSZSD(Ptr)         (Ptr.bufSizs)
#define BUSZSDI(Ptr, i)     ((Ptr.bufSizs)[i])
#define CUTSI(Ptr, i)       ((Ptr->cuts)[i])
#define CUTS(Ptr)           (Ptr->cuts)
#define CUTSD(Ptr)         (Ptr.cuts)
#define CUTSDI(Ptr, i)     ((Ptr.cuts)[i])
#define CUMI(Ptr, i)        ((Ptr->accum)[i])
#define CUM(Ptr)            (Ptr->accum)
#define CUMD(Ptr)            (Ptr.accum)
#define SIZ(Ptr)            (Ptr->size)
#define DFSIZ(Ptr)          (Ptr->defSize)
#define OFST(Ptr)           (Ptr->offset)
#define DATI(Ptr, i)        ((Ptr->data)[i])
#define DAT(Ptr)            (Ptr->data)
#define DATD(Ptr)            (Ptr.data)
#define TMPDAT(Ptr)         (Ptr->tmpData)
#define DEFDAT(Ptr)         (Ptr->defData)


 /**
 * KroFFTRep: 
 * @N: number of variables
 * @M: number of polynomias to be multiplied (must be 2 most of the time)
 * @es: 2^(es[i]) is the FFT size on dimension i 
 * @dims: dims[i]  FFT size on dimension i 
 * @size: size of the Kronecker encoding of any of the input polynomials
 * @Defsize: backup of size
 * @accum: 1-dim array of length N+1.
 *          accum[i] is the size of a coefficient w.r.t. 
 *          the i-th variable  of the product.
 * @datas:  array of pointers to the datas of the Kronecker encoding of the
 *          input polynomials. datas[0] and datas[1] point to the first
 *          and second data encodings.
 * @KN: FFT-size of the univariate product
 * @KE: 2^KE = KN
 * @KrootsPtr: powers of the primitive root of unity used for this produc
 * 
 * 
 * Return value: 
 **/   
typedef struct KroFFTRepST{
  sfixn N; 
  sfixn M; 
  sfixn * es; 
  sfixn * dims; 
  sfixn size; 
  sfixn Defsize; 
  sfixn * accum; 
  sfixn ** datas; 
  sfixn KN; 
  sfixn KE; 
  sfixn * KrootsPtr; 
} KroFFTRep;


//#define N(Ptr)              (Ptr->N)
#define M(Ptr)              (Ptr->M)
#define ES(Ptr)             (Ptr->es)
#define ESI(Ptr, i)         ((Ptr->es)[i])
#define DIMS(Ptr)           (Ptr->dims)
#define DIMSI(Ptr, i)       ((Ptr->dims)[i])
//#define SIZ(Ptr)            (Ptr->size)
//#define DFSIZ(Ptr)          (Ptr->defSize)
//#define CUMI(Ptr, i)        ((Ptr->accum)[i])
//#define CUM(Ptr)            (Ptr->accum)
#define DATS(Ptr)           (Ptr->datas)
#define DATSI(Ptr, i)       ((Ptr->datas)[i])
#define KN(Ptr)             (Ptr->KN)
#define KE(Ptr)             (Ptr->KE)
#define KROOTS(Ptr)         (Ptr->KrootsPtr)




 /**
 * KroTFTRep: 
 * @N: number of variables
 * @M: number of polynomias to be multiplied (must be 2 most of the time)
 * @es: 2^(es[i]) is the FFT size on dimension i 
 * @dims: dims[i]  FFT size on dimension i 
 * @ls: ls[i] TFT size on dimension i (= FFT size - number-of-leading-zeros) 
 * @size: size of the Kronecker encoding of any of the input polynomials
 * @Defsize: backup of size
 * @accum: 1-dim array of length N+1.
 *          accum[i] is the size of a coefficient w.r.t. 
 *          the i-th variable  of the product.
 * @datas:  array of pointers to the datas of the Kronecker encoding of the
 *          input polynomials. datas[0] and datas[1] point to the first
 *          and second data encodings.
 * @KN: FFT-size of the univariate product
 * @KE: 2^KE = KN
 * @KrootsPtr: powers of the primitive root of unity used for this produc
 * 
 * 
 * Return value: 
 **/   
typedef struct KroTFTRepST{
  sfixn N; 
  sfixn M; 
  sfixn * es; 
  sfixn * dims;
  sfixn * ls; 
  sfixn size; 
  sfixn Defsize; 
  sfixn * accum; 
  sfixn ** datas; 
  sfixn KN; 
  sfixn KE; 
  sfixn * KrootsPtr; 
} KroTFTRep;
#define LS(Ptr)           (Ptr->ls)
#define LSI(Ptr, i)       ((Ptr->ls)[i])


 /**
 * TriSetST: Data structure for zero-dim regular chains.
 *           Can be used for non zero-dim regular chain BUT this with care!
 *           Implictely we work with variables X1 < ... < XN.
 * @normalized: 1 means yes, 0 means no.
 *              Here normalized means reduced in the sense of Grobner bases
 *              and each initial is ONE.
 * @N:         the number of polynomials
 * @bounds:    bounds[i] + 1 i the degree of the i-th polynomial
 *             w.r.t. variable Xi.
 * @elems:     elems[i] is the polynomial whose main variable is Xi
 *             (if any, otherwise 0).
 * 
 * Return value: 
 **/   
typedef struct TriSetST{
  sfixn normalized;  // 1 means yes, 0 means no.
  sfixn N; 
  sfixn *bounds; 
  preFFTRep **elems; // elems[0] is useless.
}TriSet;
//#define N(Ptr)             (Ptr->N)
#define NMLZ(Ptr)             (Ptr->normalized)
#define BDS(Ptr)              (Ptr->bounds)
#define BDSI(Ptr, i)          ((Ptr->bounds)[i])
#define ELEM(Ptr)             (Ptr->elems)
#define ELEMI(Ptr, i)         ((Ptr->elems)[i])

// Pair of triangular sets, only used by lifting 
typedef struct RFuncTriSetST{
  TriSet *numeris;
  TriSet *denomis;
}RFuncTriSet;




 /**
 * TriRevInvSet: modular inverses of the elements of zero-dim reg chain
 *               to be used in the fast normal form (based on Newton iteration)
 * @N: number of polynomials = variables
 * @exist:  exist[i]=1  means the RevInv(T_i) is already computed.
 *          exist[0] is undefined
 * @NewtonLbounds: NewtonLbounds[i] is the smallest power of greater or equal
 *                 to the main degree of the  i-th quotient plus one
 * @NewtonSbounds: NewtonSbounds[i] is the degree of the  i-th quotient 
 *                 plus one.
 * @elems:         elems[i] is the inverse if the i-th polynomial in the
 *                 triangular set modulo the i-th variable and the lower
 *                 polynomials; moreover the coefficients w.r.t. 
 *                 the i-th variable are in reverse order.
 * 
 * Return value: 
 **/   
typedef struct TriRevInvSetST{
  sfixn N; // # of polys in this triangular set.
  sfixn * exist; // exist[0] is undefined. exist[i]=1 
                  // means the RevInv(T_i) is already computed.
  sfixn * NewtonLbounds; // power of 2 in Newton iteration.
  sfixn * NewtonSbounds; // real bound in Newton Iteration.
  preFFTRep ** elems; // keeps the inverses.
}TriRevInvSet;
//#define N(Ptr)              (Ptr->N)
#define EXSTI(Ptr, i)       ((Ptr->exist)[i])
#define EXST(Ptr)           (Ptr->exist)
#define NLB(Ptr)            (Ptr->NewtonLbounds)
#define NSB(Ptr)            (Ptr->NewtonSbounds)


// Pure Montgommery trick for integers
typedef struct MONTPRIMESTRUCT{
  sfixn  v;
  sfixn  Rpow;
  sfixn  Rmask;
  sfixn  Rinv;
  sfixn  V;
  sfixn  Vneg;
}  MONTP_GENE;

// prime= c.2^Npow+1.
// Used for the improved Montgommery trick for integers
typedef struct MONTPRIMEOPT2STRUCT{
  sfixn  P;
  sfixn  c;
  sfixn  Npow;
  sfixn  Rpow;
  sfixn  R_Npow;
  sfixn  Base_Rpow;// 32-Row on 32-bit machine.
  sfixn  Base_Npow;// 32-Row on 32-bit machine.
  sfixn  c_sft; // only for assembly subroutine.
  sfixn  c_pow; // 2^c_pow+1=c. for special case.
  sfixn  N2_Rpow;
  sfixn  Rmask;
  sfixn  R_Nmask;
  sfixn  Max_Root;
  sfixn  R2BRsft;
}  MONTP_OPT2_AS_GENE;



// obselete  
typedef struct MONTSMALLPRIMESTRUCT{
  sfixn  P;
  sfixn  c;
  sfixn  center;
  sfixn  N;
  sfixn  cNpow;
  sfixn  Npow;  // N = 2^Npow.
  sfixn  Nmask;
  sfixn  Ninv;
  sfixn  gate;
  sfixn  re;
  sfixn  center2;
  sfixn  reduMask;
  sfixn  * redu;
}  MONTP_SML;







// ====>>  During the summer, make these stuffs into C++  <<====


//#include "cexcept.h"
//define_exception_type(int);
//extern struct exception_context the_exception_context[1];

#ifndef plong
#define plong int32
#endif
// Add to fix the lifting allocation problem!
// Change-Code: lift-0
typedef union operandUnion operandObj;

typedef union operandUnion *operand;
typedef struct SL_Graph SLG;
typedef void *Pointer;

// Macros for the types of nodes in the DAG representation
// of polynomials as used in the lifting
typedef struct dummy_struct DUMYO;
typedef struct sfixn_struct SFIXNO;
typedef struct variable_struct VARO;
typedef struct variablePow_struct VARPOWO;
typedef struct biPlus_struct BIPLUSO;
typedef struct biSub_struct BISUBO;
typedef struct biProd_struct BIPRODO;
typedef struct polynomial_struct POLYO;
typedef struct pow_struct POWO;


#define Ntypes 8

// Type of a node in a DAG representing a polyomial
typedef
enum type {
  t_poly,   // 0
  t_sfixn,  // 1
  t_var,    // 2 (1..n  ~  x1 ~ xn)   
  t_varpow, // 3  var^e   has NOT been implemented or tested.
  t_biPlus, // 4
  t_biSub,  // 5
  t_biProd, // 6
  t_pow     // 7  operand^e           has NOT been implemented or tested.
}operandType;


// The following are macros for operating on the nodes of a DAG
#define	type_of(oper)	((operandType)(((operand)(oper))->DUMY.type))
#define	type_set(oper, t)	((((operand)(oper)  )->DUMY.type)=t)
#define	id_of(oper)	((operandType)(((operand)(oper))->DUMY.id))
#define	id_set(oper, theid)	((((operand)(oper)  )->DUMY.id)=theid)

#define is_poly(oper)   (type_of(oper) == t_poly)
#define is_var(oper)   (type_of(oper) == t_var)
#define is_sfixn(oper)  (type_of(oper) == t_sfixn)
#define is_varpow(oper)  (type_of(oper) == t_varpow)
#define is_biPlus(oper)  (type_of(oper) == t_biPlus)
#define is_biSub(oper)  (type_of(oper) == t_biSub)
#define is_biProd(oper)  (type_of(oper) == t_biProd)
#define is_pow(oper)  (type_of(oper) == t_pow)



#define poly_poly(oper) ( ((oper)->POLY).poly )
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
#define pow_e(oper)  ( ((oper)->VARPOW).e )
#define pow_base(oper)  ( ((oper)->VARPOW).base )


// dirty. compare 1 to 2.
#define is_sameVar(oper1, oper2)  ( ((oper1)->VAR).no == ((oper2)->VAR).no ) 
// dirty. compare 1 to 2.
#define is_inSameVar(oper1, oper2)  ( ((oper1)->VARPOW).no == ((oper2)->VAR).no )


#define new_poly(oper)  oper=(operand)my_calloc(1,sizeof(POLYO));	\
			 type_set(oper, t_poly) 

#define new_poly_ini(oper, thepoly)  oper=(operand)my_calloc(1,sizeof(POLYO)); \
			         type_set(oper, t_poly); \
                                 (oper->POLY).poly = thepoly


#define new_sfixn(oper)  oper=(operand)my_calloc(1,sizeof(SFIXNO));	\
			 type_set(oper, t_sfixn) 

#define new_sfixn_ini(oper, v)  oper=(operand)my_calloc(1,sizeof(SFIXNO)); \
			         type_set(oper, t_sfixn); \
                                 (oper->SFIX).sfixnVal = v

#define	new_var(oper)   oper=(operand)my_calloc(1,sizeof(VARO));   \
	                 type_set(oper, t_var)

#define	new_var_ini(oper, newno)   oper=(operand)my_calloc(1,sizeof(VARO)); \
	                          type_set(oper, t_var);\
                                  (oper->VAR).no = newno

#define	new_varpow(oper)  oper=(operand)my_calloc(1,sizeof(VARPOWO));	\
                          type_set(oper, t_varpow)

#define	new_varpow_ini(oper, newno, newe)  oper=(operand)my_calloc(1,sizeof(VARPOWO)); \
                          type_set(oper, t_varpow); \
                          (oper->VARPOW).e = newe; \
                          (oper->VARPOW).no = newno

#define	new_biPlus(oper)  oper=(operand)my_calloc(1,sizeof(BIPLUSO));	\
                           type_set(oper, t_biPlus)

#define	new_biPlus_ini(oper, newoper1, newoper2)  oper=(operand)my_calloc(1,sizeof(BIPLUSO)); \
                               type_set(oper, t_biPlus); \
                               (oper->BI_PROD).oper1 = newoper1; \
                               (oper->BI_PROD).oper2 = newoper2                             

#define	new_biSub(oper)   oper=(operand)my_calloc(1,sizeof(BISUBO));	\
                           type_set(oper, t_biSub)

#define	new_biSub_ini(oper, newoper1, newoper2)   oper=(operand)my_calloc(1,sizeof(BISUBO)); \
                               type_set(oper, t_biSub); \
                               (oper->BI_PROD).oper1 = newoper1; \
                               (oper->BI_PROD).oper2 = newoper2

#define	new_biProd(oper)  oper=(operand)my_calloc(1,sizeof(BIPRODO));	\
                            type_set(oper, t_biProd) 


#define	new_biProd_ini(oper, newoper1, newoper2)  oper=(operand)my_calloc(1,sizeof(BIPRODO)); \
                            type_set(oper, t_biProd); \
                            (oper->BI_PROD).oper1 = newoper1; \
                            (oper->BI_PROD).oper2 = newoper2

#define	new_pow(oper)  oper=(operand)my_calloc(1,sizeof(POWO));	\
                          type_set(oper, t_pow)

#define	new_pow_ini(oper, newe, newbase)  oper=(operand)my_calloc(1,sizeof(POWO)); \
                          type_set(oper, t_pow); \
                          (oper->VARPOW).e = newe; \
                          (oper->VARPOW).base = newbase



// FIRSTWORD is a data-structure storing the following information about
// each node in a DAG representing a polynomial.
/* type   -- data type.
   flag   -- bit-0 tmpMark0 (liveness).
             bit-1 tmpMark1 (derivative).
             bit-2 tmpMark2 (copy).
             bit-3 tmpMark3 (tmp dead mark).
             bit-4~7 reserved.
   id     -- The ID of a object (a node in G.)
*/
#define FIRSTWORD unsigned char type, flag; unsigned short int id  

#define tmpMarkMask 0x2
#define tmpMark2Mask 0x4
#define tmpMark3Mask 0x8
#define tmpMarkUnMask 0xfd
#define tmpMark2UnMask 0xfb
#define tmpMark3UnMask 0xf7




#define	is_TmpMark1On(oper)	((((operand)(oper))->DUMY.flag) & tmpMarkMask)
#define	is_TmpMark2On(oper)	((((operand)(oper))->DUMY.flag) & tmpMark2Mask)
#define	is_TmpMark3On(oper)	((((operand)(oper))->DUMY.flag) & tmpMark3Mask)
#define	clear_TmpMark1(oper)	((((operand)(oper))->DUMY.flag) &= tmpMarkUnMask)
#define	clear_TmpMark2(oper)	((((operand)(oper))->DUMY.flag) &= tmpMark2UnMask)
#define	clear_TmpMark3(oper)	((((operand)(oper))->DUMY.flag) &= tmpMark3UnMask)
#define	set_TmpMark1(oper)	((((operand)(oper))->DUMY.flag) |= tmpMarkMask)
#define	set_TmpMark2(oper)	((((operand)(oper))->DUMY.flag) |= tmpMark2Mask)
#define	set_TmpMark3(oper)	((((operand)(oper))->DUMY.flag) |= tmpMark3Mask)



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


// .

struct polynomial_struct {
	 FIRSTWORD;
         preFFTRep *poly;    
};




// pow.

struct pow_struct {
	FIRSTWORD;
        int32   e;  // e is the exponent .
        operand   base; // no=1 -> is x.
};

// Each of the DUMYO, .., POWO are macros
// for the node types defined above.
// operand types.
union operandUnion {
        DUMYO   DUMY;
        POLYO   POLY;
	SFIXNO  SFIX;
        VARO    VAR;
        VARPOWO VARPOW;
        BIPLUSO BI_PLUS;
        BISUBO   BI_SUB;
        BIPRODO  BI_PROD;
        POWO POW;
};


// A graph G for the Strait Line input.
// The strait line is a graph with following property.
//  0) encoded in the adjacency-list style. 
//  1) G is connected or unconnected.
//  2) A node's children must appear efore this node in the nodes vector.
//  3) the sum of all subgraphs in G is the complete encoding for the give polynomial.
struct SL_Graph{
  int32 GN;          // Number of Nodes  G
  int32 GE;          // Number of Edges in G
  operand * Nodes; // a vector of nodes. Typical adjacency-list encoding. 
                   // Nodes has size GN; Nodes[i] is the i-th nodes
                   // and Nodes[i] has ID i, for i=1...GN
};


// macros for operations on the DAG
#define GN(slg)  (slg->GN)
#define GE(slg)  (slg->GE)
#define ROOT_G(slg)  ((slg->Nodes)[(slg->GN) - 1])
#define setGN_G(slg, N)  ((slg->GN)=N)
#define NodeI(slg, i)  ((slg->Nodes)[i])

typedef struct polyMatrix_struct  POLYMATRIX;

// M by N matrix of polynomials
struct polyMatrix_struct{
  int32 M; // number of rows;
  int32 N; // number of columns;                  
  preFFTRep **entries;
};


// macros for matrix operations
#define ENTRYI_M(polyM, m, n) ((polyM->entries)[m*(polyM->N) + n])


typedef struct polyVector_SLG_struct  POLYVECTOR_SLG;

// Could be obselete ...
struct polyVector_SLG_struct{
  int32 M; // #;                  
  SLG **entries;
};

#define ENTRYI_V(polyV, i) ((polyV->entries)[i])
//#define M(polyV) (polyV->M)


typedef struct polyVector_struct  POLYVECTOR;

// Vector of polynomials
struct polyVector_struct{
  int32 M; // #;                  
  preFFTRep **entries;
};


typedef struct subProdTree_struct subProdTree;

// level_0 is useless
// level_1 is the bottom level.
// level_{h+1} is the root.
 /**
 * subProdTree: subproduct tree structure for degree one moduli polynomials
 * @h: height of the tree; the root is at level 0
 * @W: 1-dim array of size h+1;  W[n] is the node size of level n. 
 * @NoNodes: 1-dim array of size h+1; NoNodes[n] the # of nodes at level n.
 * @data: the actual data: coefficients for levels 0,1,... consecutively
 * @Bases: 1-dim array of size h+1; Bases[n] is the position of the first
 *         coefficient of the first polynomial at level n
 * 
 * Return value: 
 **/   
struct subProdTree_struct{
  sfixn h;
  sfixn *W; // node of size. W[n] is the node size of level n. 
  sfixn *NoNodes; // NoNodes[n] the # of nodes at level n.
  sfixn *data;
  sfixn *Bases;
};



typedef struct pts_trees_struct  PTS_TREE;

 /**
 * PTS_TREE: subprodtree structure for fast evaluation / interpolation
 * @no: number of variables = number of trees
 * @ptsPtr: array of arrays of values (= points)
 *           ptsPtr[i] gives the evaluation points in the i-th dim
 * @trees:  array of the corresponding subprodtrees
 * 
 * Return value: 
 **/   
struct pts_trees_struct{
  sfixn no; // no of trees / <set of points>.
  sfixn **ptsPtr; // has no+1 slots, first one is not used.
  subProdTree **trees; // has no+1 slots, first one is not used.
  sfixn times;
};




 /**
 * InterpRFRST: data-structure for rational function reconstruction 
 * @bound: degree of the modulo
 * @degNum:  actual degree of the numerator of the reconstructed fraction
 * @Num: numerator of the reconstructed fraction
 * @degDen: actual degree of the denominator of the reconstructed fraction
 * @Den: denominator of the reconstructed fraction
 * 
 * Return value: 
 **/  
typedef struct InterpRFRStruct{
  sfixn bound;
  sfixn degNum;
  sfixn *Num;
  sfixn degDen;
  sfixn *Den;
} InterpRFRST;


 /**
 * InterpRFRPreST:  data-structure for rational function interpolation
 *                  But more precisely for the modular computation of
 *                  iterated resultants modulo a 1-dimreg chain
 * @no: number of points
 * @tree: the corresponding subprodtree 
 * @points: the points 
 * @values: the images of the iterated resultant at the points 
 * @chains: the images of the input regular chain at the points 
 * @polys:   the images of the input polynomial at the points 
 * 
 * Return value: 
 **/   
typedef struct InterpRFRPreStruct{
  sfixn no;  // Number of points.
  subProdTree *tree;
  sfixn *points;
  sfixn *values;
  TriSet **chains;
  preFFTRep **polys; 
} InterpRFRPreST;


// Pair of a polynomial and a triangular set
struct RegularPairST{
  preFFTRep *poly;
  TriSet *ts;
};

typedef struct RegularPairST RegularPair;

// Pair of a polynomial list and a triangular set
struct RegularListPairST{
  sfixn no;
  preFFTRep **polyList;
  TriSet *ts;
};
typedef struct RegularListPairST  RegularListPair;



struct TaskPairST{
  int32 index;
  TriSet *ts;
};

typedef struct TaskPairST TaskPair;


// A node in a linked list /  queueu
struct LinearNodeST{
  void *element;
  struct LinearNodeST *next;
};

typedef struct LinearNodeST LinearNode;

struct LinkedQueueST{
  int32 count;
  LinearNode *front, *rear;
};

typedef struct LinkedQueueST LinkedQueue;

 /**
 * SCUBE: Data-structure for computing the subresultant chains
 *        of two multivariate polynomials using fast evaluation
 *        and interpolation.
 *        The least N-1 variables of input polynomials are evaluated 
 *        at sufficiently many points such that interpolation of their
 *        resultant can be performed.
 *        Essentially, this data-structure stores the subresultant chains
 *        of the evaluated input polynomials.
 * @w: length of the subresultant chain
 * @dim: If the input polynomials of N variables then N+1
 * @Sdims: 1-D array of length N+2
 *         Sdims[0] is not used
 *         Sdims[i] is the size of the i-th dimension of the SCube
 *         We have: Sdims[N] = Sdims[N+1] = w  
 * @Ssize: the total number slots in the SCube = the size of the SCube
 *         as a 1-dim array 
 * @Sdata: pointer to the actual data (= the multi-dim array)
 * @points_trees: pointer to the subproducttree data-structure used
 *                fast evaluation / interpolation 
 * @SLcs:  the interpolated initials of the successive regular subresultants
 * @SPolys: the  interpolated regular subresultants
 * @doneBfr: flag (see Maple code)
 * 
 * 
 * Return value: 
 **/   
struct ScubeST{
  sfixn w;  
  sfixn dim;
  sfixn *Sdims;
  sfixn Ssize;
  sfixn *Sdata;
  PTS_TREE *points_trees;
  preFFTRep **SLcs;
  preFFTRep **SPolys; 
  sfixn doneBfr;
};

typedef struct ScubeST SCUBE;








#endif
