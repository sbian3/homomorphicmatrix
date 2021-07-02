set(LT_SOURCE_FILES
    ${LT_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/convolution.cpp
    ${CMAKE_CURRENT_LIST_DIR}/packedconv.cpp
    ${CMAKE_CURRENT_LIST_DIR}/uintlinarith.cpp
    ${CMAKE_CURRENT_LIST_DIR}/common.cpp
    )

install(
    FILES
    ${CMAKE_CURRENT_LIST_DIR}/convolution.h
    ${CMAKE_CURRENT_LIST_DIR}/packedconv.h
    ${CMAKE_CURRENT_LIST_DIR}/uintlinarith.h
    ${CMAKE_CURRENT_LIST_DIR}/common.h
    ${CMAKE_CURRENT_LIST_DIR}/define_tifs.h
    DESTINATION
    ${LT_INCLUDES_INSTALL_DIR}/util
    )

set(LT_SOURCE_FILES ${LT_SOURCE_FILES} PARENT_SCOPE)
