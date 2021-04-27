cmake_minimum_required(VERSION 3.16)
project(CLion)


set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(GA_Study main.cpp main.h Population.cpp Population.h NatureSO.cpp NatureSO.h)

if (MSVC)
    target_compile_options(GA_Study PUBLIC "/utf-8")
endif ()
