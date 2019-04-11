cmake_minimum_required(VERSION 3.13)
project(MC)

include_directories(../stats/include/)
set(CMAKE_CXX_STANDARD 11)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/config.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_executable(MC main.cpp MonteCarlo.cpp MonteCarlo.hpp OptionInfo.cpp OptionInfo.hpp RNG.cpp RNG.hpp BlackScholes.cpp BlackScholes.hpp EulerMethod.hpp EulerMethod.cpp)