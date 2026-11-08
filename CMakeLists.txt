cmake_minimum_required(VERSION 3.16)
project(phys CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(PHYS_USE_FLOAT "Use float instead of double" OFF)
option(PHYS_BUILD_TESTS "Build tests" ON)
option(PHYS_BUILD_DEMOS "Build demos" ON)

add_library(phys
    src/body.cpp
    src/broadphase.cpp
    src/collision.cpp
    src/raycast.cpp
    src/solver.cpp
    src/world.cpp)
target_include_directories(phys PUBLIC include)
if(PHYS_USE_FLOAT)
    target_compile_definitions(phys PUBLIC PHYS_USE_FLOAT)
endif()
if(MSVC)
    target_compile_options(phys PRIVATE /W4)
else()
    target_compile_options(phys PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(PHYS_BUILD_TESTS)
    enable_testing()
    add_executable(phys_tests tests/tests.cpp)
    target_link_libraries(phys_tests PRIVATE phys)
    add_test(NAME phys_tests COMMAND phys_tests)
endif()

if(PHYS_BUILD_DEMOS)
    add_executable(stack_demo demos/stack_demo.cpp)
    target_link_libraries(stack_demo PRIVATE phys)
endif()
