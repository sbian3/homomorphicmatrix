# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.16)

project(BenchTIFS VERSION 1.0 LANGUAGES CXX)

set(SEAL_BUILD_EXAMPLES ON)

# Import Microsoft SEAL
find_package(SEAL 3.6 REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(LT_SOURCE_FILES "")
set(LT_INCLUDES_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/seal_lt)
add_subdirectory(src)
add_library(seal_lt STATIC ${LT_SOURCE_FILES})
set_property(TARGET seal_lt PROPERTY CXX_STANDARD 17)
target_include_directories(seal_lt PUBLIC ${CMAKE_INSTALL_PREFIX}/include/SEAL-3.6)
target_include_directories(seal_lt PUBLIC src)
add_subdirectory(bench_hlt)
