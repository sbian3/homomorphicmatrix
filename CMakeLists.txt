# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.12)

project(BenchTIFS VERSION 3.6.4 LANGUAGES CXX)

# If not called from root CMakeLists.txt
if(NOT DEFINED SEAL_BUILD_EXAMPLES)
    set(SEAL_BUILD_EXAMPLES ON)

    # Import Microsoft SEAL
    find_package(SEAL 3.6.4 EXACT REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()

set(LT_SOURCE_FILES "")
set(LT_INCLUDES_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/seal_lt)
add_subdirectory(src)
add_library(seal_lt STATIC ${LT_SOURCE_FILES})
set_property(TARGET seal_lt PROPERTY CXX_STANDARD 17)
target_include_directories(seal_lt PUBLIC ${CMAKE_INSTALL_PREFIX}/include/SEAL-3.6)
target_include_directories(seal_lt PUBLIC src)
add_subdirectory(bench_hlt)
