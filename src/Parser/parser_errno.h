#ifndef PARSER_ERRNO_H
#define PARSER_ERRNO_H

#include <stdio.h>
#include <stdlib.h>

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

/* error type */
enum{
    PARSER_SUCCESS = 0, 
    PARSER_FAILURE = -1, 
    PARSER_GENERIC = 1, /*generic error value */
    PARSER_NOALLOC = 2  /*memory allocation failed, no space available */
};

/**
 * @brief error details printed with user defined reason
 * 
 * @param reason 
 * @param file 
 * @param fun 
 * @param line 
 */
void parser_error(const char* reason, const char* file, const char* fun, int line);

/**
 * @brief error details printed with description of the error type
 * 
 * @param parser_errno 
 * @param file 
 * @param fun 
 * @param line 
 */
void parser_error_with_errno_reason(const int parser_errno, const char* file, const char* fun, int line);


__END_DECLS

#endif