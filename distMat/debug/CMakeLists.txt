
cmake_minimum_required(VERSION 3.0)


project(distMat)
#project(distMat C CXX)

find_package(Threads REQUIRED)
add_executable(
    distMat 
    ../../source/distMat.cpp
    ../../source/externs.cpp)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++1z")

if(THREADS_HAVE_PTHREAD_ARG)
    target_compile_options(distMat PUBLIC "-pthread")
endif()

if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(distMat "${CMAKE_THREAD_LIBS_INIT}")
endif()
    
