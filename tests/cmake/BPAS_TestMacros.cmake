#---------------------------------------------------------
# Macros for easily generating tests.
# Assumes the following cmake variables are set:
#   BPAS_SANITY_TEST_TARGET - the target name of sanity tests
#   BPAS_VALIDATE_TEST_TARGET - the target name of validate tests
#   BPAS_SANITY_TEST_LIBS - the common libs to link against for sanity tests
#   BPAS_VALIDATE_TEST_LIBS - the common libs to link against for validate tests
#---------------------------------------------------------


#---------------------------------------------------------
# Add a sanity test target.
# _name The test name.
# _exename The exe name.
# ARGN :
#    FILES the source files for the test
#    ARGUMENTS Arguments for test executable
#    LINK_WITH link test executable with libraries
macro(BPAS_ADD_SANITY_TEST _name _exename)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs FILES ARGUMENTS LINK_WITH)
    cmake_parse_arguments(BPAS_ADD_SANITY_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT TARGET ${_exename})
        add_executable(${_exename} ${BPAS_ADD_SANITY_TEST_FILES})
        target_compile_options(${_exename} PRIVATE ${BPAS_SANITY_TEST_ARGS} )
        target_link_libraries(${_exename} ${BPAS_SANITY_TEST_LIBS} ${BPAS_ADD_SANITY_TEST_LINK_WITH})
        add_dependencies(${BPAS_SANITY_TEST_TARGET} ${_exename})
        add_dependencies(${BPAS_BUILD_TESTS_TARGET} ${_exename})
    endif()

    add_test(NAME "${BPAS_SANITY_TESTS_MARK}${_name}" COMMAND ${_exename} ${BPAS_ADD_SANITY_TEST_ARGUMENTS})

    if(BPAS_WITH_MAPLE)
        set_tests_properties("${BPAS_SANITY_TESTS_MARK}${_name}" PROPERTIES ENVIRONMENT "MAPLE=${MAPLE_HOME};LD_LIBRARY_PATH=${MAPLE_HOME}/bin.X86_64_LINUX:$ENV{LD_LIBRARY_PATH}")
    endif()

endmacro()


#---------------------------------------------------------
# Add a validate test target.
# _name The test name.
# _exename The exe name.
# ARGN :
#    FILES the source files for the test
#    ARGUMENTS Arguments for test executable
#    LINK_WITH link test executable with libraries
macro(BPAS_ADD_VALIDATE_TEST _name _exename)
    if("${MAPLE_HOME}" STREQUAL "")
	return()
    endif()
    set(options)
    set(oneValueArgs)
    set(multiValueArgs FILES ARGUMENTS LINK_WITH)
    cmake_parse_arguments(BPAS_ADD_VALIDATE_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT TARGET ${_exename})
        add_executable(${_exename} ${BPAS_ADD_VALIDATE_TEST_FILES})
        target_compile_options(${_exename} PRIVATE ${BPAS_VALIDATE_TEST_ARGS} )
        target_link_libraries(${_exename} ${BPAS_VALIDATE_TEST_LIBS} ${BPAS_ADD_VALIDATE_TEST_LINK_WITH})
        add_dependencies(${BPAS_VALIDATE_TEST_TARGET} ${_exename})
        add_dependencies(${BPAS_BUILD_TESTS_TARGET} ${_exename})
    endif()

    add_test(NAME "${BPAS_VALIDATE_TESTS_MARK}${_name}" COMMAND ${_exename} ${BPAS_ADD_VALIDATE_TEST_ARGUMENTS})

    # YES, Maple required that you set environment variables like this...
    # It is beyond annoying.
    set_tests_properties("${BPAS_VALIDATE_TESTS_MARK}${_name}" PROPERTIES ENVIRONMENT "MAPLE=${MAPLE_HOME};LD_LIBRARY_PATH=${MAPLE_HOME}/bin.X86_64_LINUX:$ENV{LD_LIBRARY_PATH}")
endmacro()
