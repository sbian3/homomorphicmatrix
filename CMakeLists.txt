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
target_include_directories(seal_lt PRIVATE /home/kazumasita/local/include/SEAL-3.6)

### direct_conv
add_executable(direct_conv)

target_sources(direct_conv
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
    ${CMAKE_CURRENT_LIST_DIR}/direct_conv.cpp
    )

### general_lt
add_executable(general_lt)

target_sources(general_lt
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
    ${CMAKE_CURRENT_LIST_DIR}/general_matrix.cpp
    )

### packed_conv
add_executable(packed_conv)

target_sources(packed_conv
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
    ${CMAKE_CURRENT_LIST_DIR}/packed_conv.cpp
    )

### swap_slot
add_executable(swap_slot)

target_sources(swap_slot
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
    ${CMAKE_CURRENT_LIST_DIR}/swap_slot.cpp
    )

if(TARGET SEAL::seal)
    target_link_libraries(direct_conv
                            PRIVATE SEAL::seal
                            PRIVATE seal_lt)
    target_link_libraries(general_lt
                            PRIVATE SEAL::seal
                            PRIVATE seal_lt)
    target_link_libraries(packed_conv
                            PRIVATE SEAL::seal
                            PRIVATE seal_lt)
    target_link_libraries(swap_slot
                            PRIVATE SEAL::seal
                            PRIVATE seal_lt)
elseif(TARGET SEAL::seal_shared)
    target_link_libraries(direct_conv PRIVATE SEAL::seal_shared)
    target_link_libraries(general_lt PRIVATE SEAL::seal_shared)
    target_link_libraries(packed_conv PRIVATE SEAL::seal_shared)
    target_link_libraries(swap_slot PRIVATE SEAL::seal_shared)
else()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif()
