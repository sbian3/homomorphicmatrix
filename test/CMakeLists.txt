
find_package(GTest REQUIRED)
add_executable(hlt_test "")

target_sources(hlt_test
    PRIVATE
    sum.cpp
    ntt.cpp
    testmain.cpp
)

target_link_libraries(hlt_test PRIVATE SEAL::seal seal_lt GTest::gtest)
