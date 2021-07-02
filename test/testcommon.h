
#include "gtest/gtest.h"
#include "seal/seal.h"
#include "seal/modulus.h"
#include "seal/util/ntt.h"
#include "seal/util/polyarithsmallmod.h"
#include "util/uintlinarith.h"
#include "util/packedconv.h"
#include "util/common.h"
#include <iostream>
#include <vector>
#include <cstdint>
#include <chrono>

using namespace std;
using namespace seal;
using namespace seal::util;

inline void ASSERT_ARR(ConstCoeffIter expected, ConstCoeffIter actual, uint64_t size){
    for(uint64_t i = 0; i < size;i++){
        ASSERT_EQ(expected[i], actual[i]);
    }
}

inline void ASSERT_MATRIX(vector<vector<uint64_t>>& expected, vector<vector<uint64_t>>& actual){
    for(uint64_t i = 0;i < expected.size();i++){
        for(uint64_t j = 0;j < expected[i].size();j++){
            ASSERT_EQ(expected[i][j], actual[i][j]);
        }
    }
}
