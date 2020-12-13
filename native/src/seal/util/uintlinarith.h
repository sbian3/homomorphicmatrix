
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
#include <chrono>

using namespace std;

namespace seal
{
    namespace util
    {
        // 
        // matrix initialization
        // 

        void init_matrix_identity(vector<vector<uint64_t>>& matrix, uint64_t poly_modulus_degree, uint64_t scale);

        void init_matrix_rotate(vector<vector<uint64_t>>& matrix, uint64_t size, uint64_t left_rotate, uint64_t scale,const Modulus &modulus);

        void init_matrix_with_coeff(vector<vector<uint64_t>>& matrix, uint64_t size, ConstCoeffIter iter,const Modulus &modulus);
        void init_matrix_with_coeff(vector<vector<uint64_t>>& matrix, uint64_t size_matrix, ConstCoeffIter iter, uint64_t size_kernel, const Modulus &modulus);
        void init_matrix_rand_mod(vector<vector<uint64_t>>& matrix, uint64_t size, uint64_t mod);
        void init_matrix_rotate_partial(vector<vector<uint64_t>> &matrix, uint64_t size_kernel, uint64_t left_rotate, uint64_t start_col, uint64_t start_row, uint64_t scale, const Modulus &modulus);
        void init_matrix_rotate_partial(vector<vector<uint64_t>> &matrix, vector<uint64_t> kernel, uint64_t start_col, uint64_t start_row, const Modulus &modulus);
        void init_matrix_diagonal(vector<vector<uint64_t>> &matrix, uint64_t size, uint64_t scalar, uint64_t right_rotate);
        void copy_matrix(vector<vector<uint64_t>> &dest, vector<vector<uint64_t>> src, uint64_t start_col, uint64_t start_row);
        void init_matrix_2dconv(vector<vector<uint64_t>> &matrix, uint64_t input_size, vector<vector<uint64_t>> kernel);

        //
        // linear arithmetic
        //

        std::uint64_t inner_product_coeffmod(ConstCoeffIter operand1, ConstCoeffIter operand2, std::size_t coeff_count, const Modulus &modulus);

        //
        // matrix and vector arithmetic
        //

        inline void matrix_dot_vector(ConstRNSIter matrix, ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, CoeffIter result){
            // TODO: parameter validation

            SEAL_ITERATE(iter(matrix, result), coeff_count, [&](auto I){
                    get<1>(I) = inner_product_coeffmod(get<0>(I), poly_vector, coeff_count, modulus);
                    });
        }

        inline void matrix_dot_vector(ConstRNSIter matrix, uint64_t coeff_modulus_size, ConstRNSIter poly_rns, ConstModulusIter mod_chain, RNSIter result){
            // parameter validation
            // TODO: size check
            uint64_t coeff_count = poly_rns.poly_modulus_degree();
            SEAL_ITERATE(iter(poly_rns, mod_chain, result), coeff_modulus_size, [&](auto I){
                    matrix_dot_vector(matrix, get<0>(I), get<1>(I), coeff_count, get<2>(I));
                    });
        }

        inline void matrix_dot_vector(vector<vector<uint64_t>> matrix, ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, CoeffIter result){
            // TODO: parameter validation

            SEAL_ITERATE(iter(matrix, result), coeff_count, [&](auto I){
                    get<1>(I) = inner_product_coeffmod(get<0>(I), poly_vector, coeff_count, modulus);
                    });
        }

        inline void matrix_dot_vector(vector<vector<uint64_t>> matrix, uint64_t coeff_modulus_size, ConstRNSIter poly_rns, ConstModulusIter mod_chain, RNSIter result){
            uint64_t coeff_count = poly_rns.poly_modulus_degree();
            SEAL_ITERATE(iter(poly_rns, mod_chain, result), coeff_modulus_size, [&](auto I){
                    matrix_dot_vector(matrix, get<0>(I), get<1>(I), coeff_count, get<2>(I));
                    });
        }

        //
        // matrix and matrix arithmetic
        //


        void matrix_dot_product_mod(vector<vector<uint64_t>> matrixL, vector<vector<uint64_t>> matrixR, vector<vector<uint64_t>>& result,const Modulus &modulus);

        void matrix_dot_product_mod_t(vector<vector<uint64_t>> matrixL, vector<vector<uint64_t>> matrixtR, vector<vector<uint64_t>>& result, Modulus &modulus);

        //
        // Convolution
        //

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

        // 
        // for decryptor
        //

        inline void generate_c1(vector<vector<uint64_t>> matrix, uint64_t coeff_degree, CoeffIter c, const Modulus &modulus, vector<vector<uint64_t>> &new_c1){
            vector<vector<uint64_t>> A(coeff_degree, vector<uint64_t>(coeff_degree));
            init_matrix_with_coeff(A, coeff_degree, c, modulus);
            matrix_dot_product_mod(matrix, A, new_c1, modulus);
        }

        inline void generate_c1_conv(vector<uint64_t> &kernel, uint64_t coeff_degree, CoeffIter c, const Modulus &modulus, vector<vector<uint64_t>> &new_c1){
            vector<uint64_t> conved_c1(coeff_degree);
            conv_negacyclic(kernel, c, coeff_degree, modulus, conved_c1.data());
            init_matrix_with_coeff(new_c1, coeff_degree, conved_c1.data(), modulus);
        }

        inline void secret_product_with_matrix_rns(vector<vector<uint64_t>> matrix, uint64_t rns_count, RNSIter c, RNSIter s, ConstModulusIter mod_chain, RNSIter result){
            uint64_t coeff_degree = c.poly_modulus_degree();
            SEAL_ITERATE(iter(c, s, mod_chain, result), rns_count, [&](auto I){
                    auto time_start = chrono::high_resolution_clock::now();
                    vector<vector<uint64_t>> new_c1(coeff_degree, vector<uint64_t>(coeff_degree));
                    generate_c1(matrix, coeff_degree, get<0>(I), get<2>(I), new_c1);
                    auto time_c1 = chrono::high_resolution_clock::now();
                    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_c1 - time_start);
                    cout << "c1 generation: " << time_diff.count() << "ms" << endl;
                    matrix_dot_vector(new_c1, get<1>(I), get<2>(I), coeff_degree, get<3>(I));
                    });
        }

        inline void secret_pruduct_with_kernel(vector<uint64_t> &kernel, uint64_t rns_count, RNSIter c, RNSIter s, ConstModulusIter mod_chain, RNSIter result){
            uint64_t coeff_degree = c.poly_modulus_degree();
            SEAL_ITERATE(iter(c, s, mod_chain, result), rns_count, [&](auto I){
                    auto time_start = chrono::high_resolution_clock::now();
                    vector<vector<uint64_t>> new_c1(coeff_degree, vector<uint64_t>(coeff_degree));
                    generate_c1_conv(kernel, coeff_degree, get<0>(I), get<2>(I), new_c1);
                    auto time_c1 = chrono::high_resolution_clock::now();
                    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_c1 - time_start);
                    cout << "c1 generation: " << time_diff.count() << "ms" << endl;
                    matrix_dot_vector(new_c1, get<1>(I), get<2>(I), coeff_degree, get<3>(I));
                    });
        }

        //
        // print function
        // 

        void print_matrix(vector<vector<uint64_t>>& matrix);
        void print_iter(CoeffIter operand1, uint64_t coeff_count);

        inline void print_iter(RNSIter operand1, uint64_t rns_count){
            uint64_t coeff_count = operand1.poly_modulus_degree();
            SEAL_ITERATE(operand1, rns_count, [&](auto I){
                    print_iter(I, coeff_count);
                    cout << "end of RNS...." << endl;
                    });
        }

        inline void print_iter(PolyIter operand1, uint64_t poly_count){
            uint64_t rns_count = operand1.coeff_modulus_size();
            SEAL_ITERATE(operand1, poly_count, [&](auto I){
                    print_iter(I, rns_count);
                    });
        }

    } // namespace util
} // namespace seal
