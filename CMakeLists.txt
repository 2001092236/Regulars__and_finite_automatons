cmake_minimum_required(VERSION 2.6)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Link runTests with what we want to test and the GTest and pthread library

include_directories(./)

add_executable(runTests ./tests.cpp ./NKA.cpp ./DFA_regular_converter.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
