
cmake_minimum_required(VERSION 3.0)


project(NussinovFolding)
#project(distMat C CXX)

find_package(Threads REQUIRED)
add_executable(
    NussinovFolding 
    ../../source/nussinov_folding.cpp
    ../../source/nussinov_dp.cpp)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++1z")

if(THREADS_HAVE_PTHREAD_ARG)
    target_compile_options(NussinovFolding PUBLIC "-pthread")
endif()

if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(NussinovFolding "${CMAKE_THREAD_LIBS_INIT}")
endif()
    
