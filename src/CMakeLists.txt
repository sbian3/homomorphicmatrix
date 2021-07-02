set(LT_SOURCE_FILES
    ${LT_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/decryptor_lt.cpp
    )

install(
    FILES
    ${CMAKE_CURRENT_LIST_DIR}/decryptor_lt.h
    DESTINATION
    ${LT_INCLUDES_INSTALL_DIR}
    )

add_subdirectory(util)
set(LT_SOURCE_FILES ${LT_SOURCE_FILES} PARENT_SCOPE)
