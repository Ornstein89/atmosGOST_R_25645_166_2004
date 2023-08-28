cmake_minimum_required(VERSION 3.19)

project(test)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
  # URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(test_runner test_runner.cpp)

target_link_libraries(test_runner PRIVATE gtest_main)
target_link_libraries(test_runner PRIVATE atmosGOST_R_25645_166_2004)

#include(GoogleTest)
#gtest_discover_tests(test_runner)
add_test(NAME test COMMAND test_runner)
