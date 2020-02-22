
if(NOT BPAS_BUILD_SERIAL)

target_sources(${BPAS_LIB_TARGET} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/Parallel/FunctionExecutorThread.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/Parallel/ExecutorThreadPool.cpp
)

endif()
