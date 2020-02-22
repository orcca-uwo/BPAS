
find_package (BISON)
# message("BISON_EXECUTABLE= ${BISON_EXECUTABLE}")
if (NOT ${BISON_FOUND}) 
    message(WARNING "Could not find BISON on this machine... some code may not compile.")
endif()

find_package (FLEX)
# message("FLEX_EXECUTABLE= ${FLEX_EXECUTABLE}")
if (NOT ${FLEX_FOUND}) 
    message(WARNING "Could not find FLEX on this machine... some code may not compile.")
endif()


set(BPAS_FLEX_SRC "parser_lex.yy.c")
set(BPAS_BISON_SRC "parser_grammar.tab.c")
set(BPAS_BISON_HEADER "parser_grammar.tab.h")

set(BPAS_FLEX_LEX_FILE "parser_scanner.l")
set(BPAS_BISON_GRAMMAR_FILE "parser_grammar.y")

add_custom_command(
    COMMAND ${FLEX_EXECUTABLE} -o ${BPAS_FLEX_SRC} ${BPAS_FLEX_LEX_FILE}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_FLEX_LEX_FILE}
    OUTPUT 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_FLEX_SRC}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/Parser"
    COMMENT "Generating parser lexer code using flex..."
)
add_custom_target(GENERATE_LEX DEPENDS 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_FLEX_LEX_FILE}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_FLEX_SRC}
)

add_custom_command(
    COMMAND ${BISON_EXECUTABLE} -vd ${BPAS_BISON_GRAMMAR_FILE}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_BISON_GRAMMAR_FILE}
    OUTPUT 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_BISON_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_BISON_HEADER}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/Parser"
    COMMENT "Generating parser grammar code using bison..."
)
add_custom_target(GENERATE_GRAMMAR DEPENDS 
    GENERATE_LEX
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_BISON_GRAMMAR_FILE}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_BISON_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/${BPAS_BISON_HEADER}
)

#----------------------------------
# Ensure FFT headers are generated before we try to compile anything
add_dependencies(${BPAS_LIB_TARGET} GENERATE_LEX GENERATE_GRAMMAR)

target_sources(${BPAS_LIB_TARGET} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/parser_errno.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/parser_helper.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/parser_print.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/parser_lex.yy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/parser_grammar.tab.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Parser/bpas_parser.c
)
