
find_package (PythonInterp)
# message("PYTHON_EXECUTABLE= ${PYTHON_EXECUTABLE}")
if (NOT PYTHONINTERP_FOUND)
    message (FATAL_ERROR "Python not found on this machine... some code will not compile!")
endif()

set(BPAS_FFT_THRESHOLD "1024" CACHE STRING "Threshold used when building FFT")
mark_as_advanced(BPAS_FFT_THRESHOLD)


add_custom_command(
    COMMAND ${PYTHON_EXECUTABLE} generate_fft_iter_gen.py ${BPAS_FFT_THRESHOLD}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/generate_fft_iter_gen.py
    OUTPUT 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversal_${BPAS_FFT_THRESHOLD}.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter_${BPAS_FFT_THRESHOLD}.c
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/arraybitreversal_${BPAS_FFT_THRESHOLD}.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_iter_${BPAS_FFT_THRESHOLD}.h
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src"
    COMMENT "Generating code for general fft."
)
add_custom_target(GENERATE_FFT_ITER DEPENDS 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversal_${BPAS_FFT_THRESHOLD}.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter_${BPAS_FFT_THRESHOLD}.c
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/arraybitreversal_${BPAS_FFT_THRESHOLD}.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_iter_${BPAS_FFT_THRESHOLD}.h
)



add_custom_command(
    COMMAND ${PYTHON_EXECUTABLE} generate_fft.py ${BPAS_FFT_THRESHOLD}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/generate_fft.py
    OUTPUT 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/arraybitreversal.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_iter1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_iter2.h
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src"
    COMMENT "Generating code for fft."
)
add_custom_target(GENERATE_FFT DEPENDS 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/arraybitreversal.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_iter1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_iter2.h
)


add_custom_command(
    COMMAND ${PYTHON_EXECUTABLE} generate_fft_furer.py ${BPAS_FFT_THRESHOLD}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/generate_fft_furer.py
    OUTPUT
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversalfurer.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_furer1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_furer2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/arraybitreversalfurer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_furer1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_furer2.h
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src"
    COMMENT "Generating code for fft furer."
)
add_custom_target(GENERATE_FURER DEPENDS 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversalfurer.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_furer1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_furer2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/arraybitreversalfurer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_furer1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/fft_furer2.h
)

add_custom_command(
    COMMAND ${PYTHON_EXECUTABLE} generate_tft_relax.py ${BPAS_FFT_THRESHOLD}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/generate_fft_furer.py
    OUTPUT
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/tft_tree1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/tft_tree2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/tft_tree1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/tft_tree2.h
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src"
    COMMENT "Generating code for tft."
)
add_custom_target(GENERATE_TFT DEPENDS 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/tft_tree1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/tft_tree2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/tft_tree1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FFT/src/tft_tree2.h
)

#----------------------------------
# Ensure FFT headers are generated before we try to compile anything
add_dependencies(${BPAS_LIB_TARGET} GENERATE_FFT GENERATE_FURER GENERATE_TFT)

target_sources(${BPAS_LIB_TARGET} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/twocon_general.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/twocon_basic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_furer1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_furer2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/modpn_export.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/ser_general_routine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/ser_basic_routine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/PrimeField.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/general_routine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/basic_routine.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/pbpas_basic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/tft_tree1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/tft_tree2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/transpose.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/arraybitreversal_${BPAS_FFT_THRESHOLD}.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter_${BPAS_FFT_THRESHOLD}.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/FFT/src/fft_iter_main.c
)
