
#pragma once

#include "seal/seal.h"
#include "seal/util/common.h"
#include "seal/util/defines.h"
#include "seal/util/pointer.h"
#include "seal/util/iterator.h"
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

using namespace std;

namespace seal
{
    namespace util
    {
        // 
        // matrix initialization
        // 

        void init_matrix_identity(vector<vector<int64_t>>& matrix, uint64_t poly_modulus_degree, int64_t scale);

        void init_matrix_rotate(vector<vector<int64_t>>& matrix, uint64_t size, uint64_t left_rotate, int64_t scale);

        void init_matrix_with_coeff(vector<vector<int64_t>>& matrix, uint64_t size, ConstCoeffIter iter);

        void init_matrix_rand_mod(vector<vector<int64_t>>& matrix, uint64_t size, uint64_t mod);

        //
        // linear arithmetic
        //

        std::uint64_t inner_product_coeffmod(vector<int64_t> operand1, ConstCoeffIter operand2, std::size_t coeff_count, const Modulus &modulus);

        //
        // matrix and vector arithmetic
        //

        void matrix_dot_vector(vector<vector<int64_t>> matrix, ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, CoeffIter result);

        inline void matrix_dot_vector(vector<vector<int64_t>> matrix, uint64_t coeff_modulus_size, ConstRNSIter poly_rns, ConstModulusIter mod_chain, RNSIter result){
            // parameter validation
            // TODO: size check
            uint64_t coeff_count = poly_rns.poly_modulus_degree();
            SEAL_ITERATE(iter(poly_rns, mod_chain, result), coeff_modulus_size, [&](auto I){
                    matrix_dot_vector(matrix, get<0>(I), get<1>(I), coeff_count, get<2>(I));
                    });
        }


        // 
        // for decryptor
        //

        void secret_product_with_matrix(vector<vector<int64_t>> matrix,uint64_t coeff_degree, CoeffIter c, CoeffIter s, const Modulus& modulus, CoeffIter result);

        inline void secret_product_with_matrix_rns(vector<vector<int64_t>> matrix, uint64_t rns_count, RNSIter c, RNSIter s, ConstModulusIter mod_chain, RNSIter result){
            uint64_t coeff_degree = c.poly_modulus_degree();
            SEAL_ITERATE(iter(c, s, mod_chain, result), rns_count, [&](auto I){
                    secret_product_with_matrix(matrix, coeff_degree, get<0>(I), get<1>(I), get<2>(I), get<3>(I));
                    });
        }

        //
        // matrix and matrix arithmetic
        //

        void matrix_dot_product_mod(vector<vector<int64_t>> matrixL, vector<vector<int64_t>> matrixR, vector<vector<int64_t>>& result, uint64_t mod);

        //
        // print function
        // 

        void print_matrix(vector<vector<int64_t>>& matrix);


    } // namespace util
} // namespace seal
