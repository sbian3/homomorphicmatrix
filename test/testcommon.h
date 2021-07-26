
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

inline vector<pair<uint64_t, uint64_t>> make_pair_vector(vector<uint64_t> values, vector<uint64_t> counts){
    assert(values.size() == counts.size());
    vector<pair<uint64_t, uint64_t>> pair_vector(values.size());
    for(uint64_t i = 0;i < values.size();i++){
        pair_vector[i] = make_pair(values[i], counts[i]);
    }
    return pair_vector;
}

inline void ASSERT_PAIR_VEC(vector<pair<uint64_t, uint64_t>> expected, vector<pair<uint64_t, uint64_t>> actual){
    ASSERT_EQ(expected.size(), actual.size());
    for(uint64_t i = 0;i < expected.size();i++){
        ASSERT_EQ(expected[i].first, actual[i].first);
        ASSERT_EQ(expected[i].second, actual[i].second);
    }
}

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
