
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

        void init_matrix_rotate(vector<vector<int64_t>>& matrix, uint64_t size, int64_t right_rotate, int64_t scale);

        void init_matrix_with_coeff(vector<vector<int64_t>>& matrix, uint64_t size, ConstCoeffIter iter);

        void init_matrix_rand_mod(vector<vector<int64_t>>& matrix, uint64_t size, uint64_t mod);

        //
        // linear arithmetic
        //

        std::uint64_t inner_product_coeffmod(vector<int64_t> operand1, seal::util::ConstCoeffIter operand2, std::size_t coeff_count, const Modulus &modulus);

        std::uint64_t inner_product_coeffmod(seal::util::ConstCoeffIter operand1, seal::util::ConstCoeffIter operand2, std::size_t coeff_count, const Modulus &modulus);

        void matrix_dot_product_mod(vector<vector<int64_t>> matrixL, vector<vector<int64_t>> matrixR, vector<vector<int64_t>>& result, uint64_t mod);

        void matrix_dot_vector(seal::util::ConstRNSIter matrix, seal::util::ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, seal::util::CoeffIter result);

        void matrix_dot_vector(vector<vector<int64_t>> matrix, seal::util::ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, seal::util::CoeffIter result);

        void matrix_dot_vector(seal::util::ConstRNSIter matrix, uint64_t coeff_modulus_size, uint64_t coeff_count, seal::util::ConstRNSIter poly_rns, util::ConstModulusIter mod_chain, seal::util::RNSIter result);


        void matrix_dot_vector(vector<vector<int64_t>> matrix, uint64_t coeff_modulus_size, uint64_t coeff_count, seal::util::ConstRNSIter poly_rns, util::ConstModulusIter mod_chain, seal::util::RNSIter result);

        //
        // print function
        // 

        void print_matrix(vector<vector<int64_t>>& matrix);

        void print_iter(util::CoeffIter operand1, uint64_t coeff_count);

        inline void print_iter(util::RNSIter operand1, uint64_t rns_count, uint64_t coeff_count){
            SEAL_ITERATE(operand1, rns_count, [&](auto I){
                    print_iter(I, coeff_count);
                    cout << "end of RNS...." << endl;
                    });
        }
    } // namespace util
} // namespace seal
