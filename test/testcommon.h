
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

inline vector<uint64_t> get_nonzero_indexlist(vector<uint64_t> &input){
    if(input.empty()){
        cerr << "no data in input!" << endl;
    }
    vector<uint64_t> nonzero_index;
    for(uint64_t i = 0;i < input.size();i++){
        if(input[i] != 0){
            nonzero_index.push_back(i);
        }
    }
    return nonzero_index;
}


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

inline void ASSERT_ARR(ConstRNSIter expected, ConstRNSIter actual, uint64_t rns_size){
    auto coeff_size = expected.poly_modulus_degree();
    SEAL_ITERATE(iter(expected, actual), rns_size, [&](auto I){
            ASSERT_ARR(get<0>(I), get<1>(I), coeff_size);
            });
}

inline void ASSERT_MATRIX(vector<vector<uint64_t>>& expected, vector<vector<uint64_t>>& actual){
    for(uint64_t i = 0;i < expected.size();i++){
        for(uint64_t j = 0;j < expected[i].size();j++){
            ASSERT_EQ(expected[i][j], actual[i][j]);
        }
    }
}
