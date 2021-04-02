# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.12)

project(BenchTIFS VERSION 3.6.2 LANGUAGES CXX)

# If not called from root CMakeLists.txt
if(NOT DEFINED SEAL_BUILD_EXAMPLES)
    set(SEAL_BUILD_EXAMPLES ON)

    # Import Microsoft SEAL
    find_package(SEAL 3.6.2 EXACT REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()

add_executable(direct_conv)

target_sources(direct_conv
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
    ${CMAKE_CURRENT_LIST_DIR}/direct_conv.cpp
    )

add_executable(general_lt)

target_sources(general_lt
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
    ${CMAKE_CURRENT_LIST_DIR}/general_matrix.cpp
    )

add_executable(packed_conv)

target_sources(packed_conv
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
    ${CMAKE_CURRENT_LIST_DIR}/packed_conv.cpp
    )

if(TARGET SEAL::seal)
    target_link_libraries(direct_conv PRIVATE SEAL::seal)
    target_link_libraries(general_lt PRIVATE SEAL::seal)
    target_link_libraries(packed_conv PRIVATE SEAL::seal)
elseif(TARGET SEAL::seal_shared)
    target_link_libraries(direct_conv PRIVATE SEAL::seal_shared)
    target_link_libraries(general_lt PRIVATE SEAL::seal_shared)
    target_link_libraries(packed_conv PRIVATE SEAL::seal_shared)
else()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif()
