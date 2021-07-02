
#pragma once

#include "seal/util/common.h"
#include "seal/util/defines.h"
#include "seal/util/pointer.h"
#include "seal/util/iterator.h"
#include "seal/util/uintarithmod.h"
#include "seal/util/polyarithsmallmod.h"
#include "convolution.h"
#include "seal/ciphertext.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <numeric>
#include <string>
#include <sstream>
#include <random>
#include <vector>
#include <iostream>
#include <cassert>
#include <chrono>

using namespace std;
using namespace seal;

inline vector<uint64_t> sample_rn(uint64_t size, Modulus modulus){
    vector<uint64_t> ret(size);
    for(uint64_t i = 0;i < size;i++){
        ret[i] = rand() % modulus.value();
        if(ret[i] == 0){
            i--;
        }
    }
    return ret;
}

inline void sample_rn(util::CoeffIter array, uint64_t size, Modulus modulus){
    for(uint64_t i = 0;i < size;i++){
        array[i] = rand() % modulus.value();
        if(array[i] == 0){
            i--;
        }
    }
}

inline vector<vector<uint64_t>> sample_rn(uint64_t size_col, uint64_t size_row, Modulus modulus){
    vector<vector<uint64_t>> ret(size_col, vector<uint64_t>(size_row));
    for(uint64_t i = 0; i < size_col;i++){
        for(uint64_t j = 0; j < size_row;j++){
            ret[i][j] = rand() % modulus.value();
            if(ret[i][j]==0){
                j--;
            }
        }
    }
    return ret;
}
