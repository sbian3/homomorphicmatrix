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
elseif(TARGET SEAL::seal_shared)
    target_link_libraries(direct_conv PRIVATE SEAL::seal_shared)
    target_link_libraries(general_lt PRIVATE SEAL::seal_shared)
    target_link_libraries(packed_conv PRIVATE SEAL::seal_shared)
else()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif()
