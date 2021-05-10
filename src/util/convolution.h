
#pragma once

#include "seal/util/common.h"
#include "seal/util/defines.h"
#include "seal/util/pointer.h"
#include "seal/util/iterator.h"
#include "seal/util/uintarithmod.h"
#include "seal/ciphertext.h"
#include "seal/memorymanager.h"
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

namespace seal
{
    namespace util
    {
        ///////////////////
        // 1D Convolution
        ///////////////////
        void conv_negacyclic(vector<uint64_t> &kernel, CoeffIter encrypted, uint64_t poly_degree, const Modulus &modulus, CoeffIter result);

        inline void conv_negacyclic(vector<uint64_t> &kernel, RNSIter encrypted,uint64_t rns_count, ConstModulusIter mod_chain, RNSIter result){
            uint64_t coeff_degree = encrypted.poly_modulus_degree();
            SEAL_ITERATE(iter(encrypted, mod_chain, result),rns_count, [&](auto I){
                conv_negacyclic(kernel, get<0>(I), coeff_degree, get<1>(I), get<2>(I));
            });
        }

        inline void conv_negacyclic(vector<uint64_t> &kernel, Ciphertext &encrypted, ConstModulusIter mod_chain, Ciphertext &result){
            uint64_t poly_count = encrypted.size();
            uint64_t rns_count = encrypted.coeff_modulus_size();
            SEAL_ITERATE(iter(encrypted, result), poly_count, [&](auto I){
                conv_negacyclic(kernel, get<0>(I), rns_count, mod_chain, get<1>(I));
                    });
        }
        ////////////////////
        // 2D Convolution
        ////////////////////
        vector<uint64_t> transform_to_helix(vector<vector<uint64_t>> &input, uint64_t colsize_helix);

    }
}

