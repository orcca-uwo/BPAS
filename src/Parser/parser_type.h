#ifndef PARSER_TYPE_H
#define PARSER_TYPE_H

#include "../../include/RationalNumberPolynomial/SMQP_Support-AA.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
	# define __BEGIN_DECLS extern "C" {
	# define __END_DECLS }
#else
	# define __BEGIN_DECLS /* empty */
	# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct _term{
    ratNum_t coef;
	degree_t *exp;
}term;

typedef struct _powervar{
	char *var;
	degree_t exp;
}powervar;

typedef struct _altarr_pack{
	AltArr_t* altarr_t_data;
	char** vars;
	int numVars;
}altarr_pack;

__END_DECLS

#endif
