
find_package(GTest REQUIRED)
add_executable(hlt_test "")

target_sources(hlt_test
    PRIVATE
    ntt.cpp
    matrix.cpp
    kernelInfo.cpp
    experiment.cpp
    testmain.cpp
)

target_link_libraries(hlt_test PRIVATE SEAL::seal seal_lt GTest::gtest)
