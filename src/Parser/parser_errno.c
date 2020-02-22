#include "parser_errno.h"

/**
 * @brief provided the error number decription of error type will be return .
 * 
 * @param parser_errno 
 * @return const char* 
 */
static const char* parser_errno_msg(const int parser_errno){
    switch(parser_errno){
        case PARSER_SUCCESS:
            return "success";
        case PARSER_FAILURE:
            return "failure";
        case PARSER_NOALLOC:
            return "memory allocation failed, no space available";
        case PARSER_GENERIC:
            return "generic error code";
        default:
            return "unknown error code";        
    }
}

/**
 * @brief Display parser error with user message or reason specfied.
 * 
 * @param reason 
 * @param file 
 * @param fun 
 * @param line 
 */
void parser_error(const char* reason, const char* file, const char* fun, int line){
    fprintf(stderr, "Parser error: [%s]: %s: @%d: %s\n", file, fun, line, reason);
    exit(PARSER_FAILURE);
}

/**
 * @brief Display parser error with default parser error number definition.
 * 
 * @param parser_errno 
 * @param file 
 * @param fun 
 * @param line 
 */
void parser_error_with_errno_reason(const int parser_errno, const char* file, const char* fun, int line){
    fprintf(stderr, "Parser error: [%s] : %s: @%d:  %s:\n", file, fun, line, parser_errno_msg(parser_errno));
    exit(PARSER_FAILURE);
}