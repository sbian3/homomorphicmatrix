# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.16)

project(BenchTIFS VERSION 1.0 LANGUAGES CXX)

set(SEAL_BUILD_EXAMPLES ON)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(LT_SOURCE_FILES "")
set(LT_INCLUDES_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/seal_lt)
add_subdirectory(src)

option(HLT_FETCH_THIRDPARTY "Download Thirdparty Library" ON)
if(HLT_FETCH_THIRDPARTY)
    include(FetchContent)
    mark_as_advanced(FETCHCONTENT_BASE_DIR)
    mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
    mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
    mark_as_advanced(FETCHCONTENT_QUIET)
    mark_as_advanced(SEAL_BUILD_DEPS)
    mark_as_advanced(SEAL_BUILD_TESTS)
    mark_as_advanced(SEAL_BUILD_BENCH)
    set(HLT_THIRDPARTY_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)
    message(${HLT_THIRDPARTY_DIR})
    set(FETCHCONTENT_BASE_DIR ${HLT_THIRDPARTY_DIR} CACHE STRING "" FORCE)
    set(THIRDPARTY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty)
    message(STATUS "SEAL: download ...")
    FetchContent_Declare(
        seal
        GIT_REPOSITORY https://github.com/microsoft/SEAL.git
        GIT_TAG 0b058d99b7f18a00e5ebb2b80caee593804b0500 # 3.6.4
        )
    FetchContent_GetProperties(seal)
    if(NOT seal_POPULATED)
        FetchContent_Populate(seal)
        add_subdirectory(
            ${seal_SOURCE_DIR}
            EXCLUDE_FROM_ALL
            )
    endif()

    message(STATUS "Microsoft GSL: download")
    FetchContent_Declare(
        msgsl
        GIT_REPOSITORY https://github.com/microsoft/GSL.git
        GIT_TAG        0f6dbc9e2915ef5c16830f3fa3565738de2a9230 # 3.1.0
        )
    FetchContent_GetProperties(msgsl)
    if(NOT msgsl_POPULATED)
        FetchContent_Populate(msgsl)

        set(GSL_CXX_STANDARD "14" CACHE STRING "" FORCE)
        set(GSL_TEST OFF CACHE BOOL "" FORCE)
        mark_as_advanced(GSL_CXX_STANDARD )
        mark_as_advanced(GSL_TEST)
        mark_as_advanced(FETCHCONTENT_SOURCE_DIR_MSGSL)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_MSGSL)

        add_subdirectory(
            ${msgsl_SOURCE_DIR}
            EXCLUDE_FROM_ALL)
    endif()
else()
    # You have already installed SEAL
    find_package(SEAL 3.6 REQUIRED)
endif()

add_library(seal_lt STATIC ${LT_SOURCE_FILES})
target_compile_features(seal_lt PUBLIC cxx_std_17)
target_include_directories(seal_lt PUBLIC src)

if(HLT_FETCH_THIRDPARTY)
    add_dependencies(seal_lt seal)
    target_include_directories(seal_lt PUBLIC $<BUILD_INTERFACE:${msgsl_SOURCE_DIR}/include>)
    target_include_directories(seal_lt PUBLIC $<BUILD_INTERFACE:${seal_SOURCE_DIR}/native/src>)
    target_include_directories(seal_lt PUBLIC $<BUILD_INTERFACE:$<TARGET_PROPERTY:seal,BINARY_DIR>/native/src>)
else()
    target_include_directories(seal_lt PUBLIC ${CMAKE_INSTALL_PREFIX}/include/SEAL-3.6)
endif()

add_subdirectory(bench_hlt)
option(HLT_BUILD_TEST "build test" OFF)
if(HLT_BUILD_TEST)
    add_subdirectory(test)
endif()
