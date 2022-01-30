
#pragma once

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
#include "uintlinarith.h"
#include "seal/util/uintarithmod.h"
#include "seal/memorymanager.h"
#include "util/uintlinarith.h"

using namespace std;

namespace seal
{
    namespace util
    {

        vector<uint64_t> create_diagonal_scalars(const vector<uint64_t> &kernel, const uint64_t colsize, const uint64_t rowsize, const Modulus &modulus, vector<uint64_t> &diagonal_list);

        class KernelInfo{
            public:
                uint64_t input_size;
                uint64_t kernel_size;
                uint64_t block_size;
                // toeplitz part of multiplied matrix
                vector<uint64_t> toeplitz_diagonal_scalars;
                KernelInfo(){

                }
                KernelInfo(uint64_t input_s, uint64_t block_s, uint64_t st_c,uint64_t  st_r,uint64_t si_c,uint64_t si_r,vector<uint64_t> dat, Modulus mod):
                    input_size(input_s), block_size(block_s), start_col(st_c), start_row(st_r), size_col(si_c), size_row(si_r), data(dat), modulus(mod){
                        kernel_size = dat.size();
                        diagonal_scalars = vector<uint64_t>(size_col + size_row - 1);
                        index = create_diagonal_scalars(data, size_col, size_row, modulus, diagonal_scalars);
                        //util::print_vector(diagonal_list, diagonal_list.size());
                    }
                vector<uint64_t> diagonal_scalars;
                vector<uint64_t> index;

                uint64_t get_colsize(){
                    return size_col;
                }

                uint64_t get_rowsize(){
                    return size_row;
                }

                uint64_t get_startcol(){
                    return start_col;
                }

                uint64_t get_startrow(){
                    return start_row;
                }

                void getParamsforSubmatrix(uint64_t &submat_startcol, uint64_t &submat_colsize){
                    submat_startcol = start_row;
                    submat_colsize = size_row;
                }

                // for matrix-vector multiplication routine
                // return vector<pair<index, value>>
                vector<pair<uint64_t, uint64_t>> make_rowpair(Modulus mod){
                    vector<pair<uint64_t, uint64_t>> ret(data.size());
                    ret[0] = make_pair(start_row, data[0]);
                    // negate seconds data and later
                    for(uint64_t i = 1;i< data.size();i++){
                        ret[i] = make_pair(start_row + size_row - i, negate_uint_mod(data[i], mod));
                    }
                    return ret;
                }

                void pair_nextcol(vector<pair<uint64_t, uint64_t>> &pair_kernel, Modulus mod){
                    for(uint64_t i = 0;i < pair_kernel.size();i++){
                        uint64_t next_index = pair_kernel[i].first + 1;
                        uint64_t next_value = pair_kernel[i].second;
                        if(next_index > start_row + size_row){
                            cerr << "pair_nextcol: invalid index " << endl;
                            return;
                        }
                        if(next_index == start_row + size_row){
                            next_index = start_row;
                            next_value = negate_uint_mod(next_value, mod);
                        }
                        pair_kernel[i] = make_pair(next_index, next_value);
                    }
                }

                // pair(value, valuelen)
                void get_toeplitz(vector<vector<pair<uint64_t, uint64_t>>>& pairs, uint64_t poly_degree){
                    uint64_t toeplitz_len = input_size + poly_degree - 1;
                    vector<uint64_t> toeplitz_tmp(toeplitz_len);
                    for(uint64_t i = 0;i < toeplitz_len;i++){
                        if(pairs[i].empty()){
                            throw std::out_of_range("no data in pair");
                        }
                        for(uint64_t j = pairs[i].size()-1;j >=0;j--){
                            if(pairs[i][j].first != 0 && pairs[i][j].second != 0){
                                toeplitz_tmp[i] = pairs[i][j].first;
                                break;
                            }
                        }
                    }
                    toeplitz_diagonal_scalars = toeplitz_tmp;
                }

                void print(){
                    cout << "start col: " << start_col << endl;
                    cout << "start row: " << start_row << endl;
                    cout << "diagonal: ";
                    for(uint64_t k = 0;k < diagonal_scalars.size();k++){
                        cout << diagonal_scalars[k] << " ";
                    }
                    cout << endl;
                    cout << "index: ";
                    for(uint64_t k = 0;k < index.size();k++){
                        cout << index[k] << " ";
                    }
                    cout << endl;
                }

            private:
                uint64_t size_col;
                uint64_t size_row;
                uint64_t start_col;
                uint64_t start_row;
                Modulus modulus;
                vector<uint64_t> data;
        };

        vector<uint64_t> create_diagonal_from_submatrix(CoeffIter a, uint64_t poly_degree, uint64_t start_col, uint64_t colsize, Modulus &modulus);

        //vector<uint64_t> matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus);

        inline void matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus, vector<vector<pair<uint64_t, uint64_t>>> &diagonalpairlist, uint64_t index_pairlist){
            // assert list_R is larger than kernel
            assert(kernel_L_indexes.size() <= list_R.size());
            //assert(colsize_R <= kernel_L.size());
        
#if HLT_DEBUG_PRINT == 1
            cout << "---------calculating a diagonal-----" << endl;
            cout << "kernel: " << endl;
            print_vector(kernel_L, kernel_L.size());
            cout << "list_R:"  << endl;
            print_vector(list_R, list_R.size());
            cout << "offset: " << offset << endl;
#endif
            // calculate element wise product
            // optimize complexity remembering kernel nonzero elements
            //uint64_t wise_prod_len = kernel_L.size() <= list_R.size()? kernel_L.size(): list_R.size();
#if HLT_DEBUG_TIME == 1
            auto diagonal_begin = chrono::high_resolution_clock::now();
#endif
#if HLT_DEBUG_TIME == 1
            auto mul_start = chrono::high_resolution_clock::now();
#endif
#if HLT_DEBUG_TIME == 1
            auto mul_end = chrono::high_resolution_clock::now();
#endif
            uint64_t wise_prod_len;
            uint64_t innerp_size = colsize_R;
            uint64_t first_right_edge;
            bool end_in_firstinner = false;
            if(offset >= 0){
                if(offset + kernel_L.size() > list_R.size()){
                    wise_prod_len = list_R.size() - offset;
                }else{
                    wise_prod_len = kernel_L.size(); 
                }
            }else{
                wise_prod_len = kernel_L.size() + offset;
            }
            if (wise_prod_len >= innerp_size) {
              first_right_edge = innerp_size;
            } else {
              cout << "end in first inner" << endl;
              first_right_edge = wise_prod_len;
              end_in_firstinner = true;
            }
            if(offset < 0){
                first_right_edge -= offset;
            }
            uint64_t partial_sum = 0;
            uint64_t index_iterator_right = 0;
            uint64_t index_iterator_left = 0;
            // calculate first inner prod
            for(uint64_t i = 0; i < kernel_L_indexes.size() && kernel_L_indexes[i] < first_right_edge;i++){
                uint64_t prod;
                if(offset >= 0){
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                }else{
                    // boarder check and mul
                    if(kernel_L_indexes[i] <= -offset-1){
                      index_iterator_right++;
                      index_iterator_left++;
                      continue;
                    }
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                }
                partial_sum = util::add_uint_mod(partial_sum, prod, modulus);
                index_iterator_right++;
            }
#if HLT_DEBUG_TIME == 1
            auto innerp_end = chrono::high_resolution_clock::now();
#endif
#if HLT_DEBUG_PRINT == 1
            cout << "index_iterator_right: " << index_iterator_right << endl;
            cout << "wise_prod_len: " << wise_prod_len << endl;;
            //for(uint64_t i = 0; i < wise_prod_index.size();i++){
            //    cout << "index " << wise_prod_index[i] << ": " << wise_prod[wise_prod_index[i]] << endl;;
            //}
            cout << "prod_times: " << innerp_size << endl;
            //cout << "end index of nonzero wise_prod: " << wise_prod_index[wise_prod_index.size()- 1] << endl;
#endif

            // slide window
            // we need O(1) to calc a next diagonal element
            uint64_t jump_right;
            uint64_t jump_left;
            bool is_right_edge = false;
            bool is_left_edge = false;
            uint64_t pair_num = 0;
            uint64_t i;
            if(offset < 0){
                i = -offset;
                wise_prod_len += i;
            }else{
                i = 0;
            }
            for (; i + innerp_size - 1 < wise_prod_len; i++) {
                // how many index should we jump to calc next partial_sum?
                uint64_t window_left = i;
                uint64_t window_right = i + innerp_size - 1;
                jump_right =
                    kernel_L_indexes[index_iterator_right] - window_right;
                jump_left = kernel_L_indexes[index_iterator_left] - window_left + 1;
                if (is_right_edge) {
                    jump_right = wise_prod_len - window_right;
                }
                if(is_left_edge){
                    jump_left = wise_prod_len - window_left;
                }
                uint64_t jump_len = min(jump_right, jump_left);
#if HLT_DEBUG_PRINT == 1
                cout << "--loop: i=" << i << "--" << endl;
                cout << "iterator_of_index: (" << index_iterator_left << ", "
                    << index_iterator_right << ") " << endl;
                cout << "kernel indexes: ("
                    << kernel_L_indexes[index_iterator_left] << ", "
                    << kernel_L_indexes[index_iterator_right] << ") " << endl;
                cout << "windows: (" << i << ", " << i + innerp_size - 1 << ")"
                    << endl;
                cout << "jump_right: " << jump_right << endl;
                cout << "jump_left: " << jump_left << endl;
                cout << "jump_len: " << jump_len << endl;
                cout << "partial_sum: " << partial_sum << endl;
                cout << "make pair: " << partial_sum << ", " << jump_len << endl;
                cout << endl;
#endif
                // make pair
                auto pair = make_pair(partial_sum, jump_len);
                diagonalpairlist[index_pairlist].push_back(pair);
                pair_num++;

                // update partial sum
                if (jump_len == jump_right) {
                    uint64_t list_R_coeff = list_R[kernel_L_indexes[index_iterator_right] + offset];
                    auto coeff_prod = util::multiply_uint_mod(
                            kernel_L[kernel_L_indexes[index_iterator_right]],
                            list_R_coeff, modulus);
                    partial_sum =
                        util::add_uint_mod(partial_sum, coeff_prod, modulus);
                    if (index_iterator_right == kernel_L_indexes.size() - 1) {
#if HLT_DEBUG_PRINT == 1
                        cout << "index_right is on edge!" << endl;
#endif
                        is_right_edge = true;
                    } else {
                        index_iterator_right++;
                    }
                }
                if (jump_len == jump_left) {
                    uint64_t list_R_coeff = list_R[kernel_L_indexes[index_iterator_left] + offset];
                    auto coeff_prod = util::multiply_uint_mod(
                            kernel_L[kernel_L_indexes[index_iterator_left]],
                            list_R_coeff, modulus);
                    partial_sum =
                        util::sub_uint_mod(partial_sum, coeff_prod, modulus);
                    if (index_iterator_left == kernel_L_indexes.size() - 1) {
#if HLT_DEBUG_PRINT == 1
                        cout << "index_left is on edge!" << endl;
#endif
                        is_left_edge = true;
                    } else {
                        index_iterator_left++;
                    }
                }
                i = i + jump_len - 1;
            }
            //cout << "pair_num: " << pair_num << endl;
#if HLT_DEBUG_TIME == 1
            auto slide_end   = chrono::high_resolution_clock::now();
            auto begin_diff  = chrono::duration_cast<chrono::nanoseconds>(mul_start - diagonal_begin);
            auto mul_diff    = chrono::duration_cast<chrono::nanoseconds>(mul_end - mul_start);
            auto innerp_diff = chrono::duration_cast<chrono::nanoseconds>(innerp_end - mul_end);
            auto slide_diff  = chrono::duration_cast<chrono::nanoseconds>(slide_end - innerp_end);
            cout <<  "begin: " << begin_diff.count() << " mul : " << mul_diff.count() << " innerp: " << innerp_diff.count() << " slide: " << slide_diff.count() << " sum: " << begin_diff.count() + mul_diff.count() + innerp_diff.count() + slide_diff.count() << endl;
#endif
        }

        void diagonallist_to_matrix(vector<vector<uint64_t>> &diagonallist, uint64_t start_col, uint64_t start_row, uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result);

        void diagonallist_to_matrix(vector<vector<pair<uint64_t, uint64_t>>> &diagonallist, uint64_t start_col, uint64_t start_row , uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result);

        vector<vector<uint64_t>> scalars_to_diagonallist(vector<uint64_t> scalars, uint64_t colsize, uint64_t rowsize);
        inline uint64_t get_blocksize(uint64_t input_dim, uint64_t kernel_dim, uint64_t padding){
            return input_dim + kernel_dim - 1 + padding;
        }

        vector<uint64_t> pack_input(const vector<vector<uint64_t>> &input, vector<KernelInfo> &kernel_infos, uint64_t poly_size);

        vector<KernelInfo> pack_kernel(vector<vector<uint64_t>> &kernels, vector<vector<uint64_t>> &inputs, Modulus modulus, uint64_t poly_degree);

        void pack_kernel_to_matrix(vector<KernelInfo> kernelinfos, vector<vector<uint64_t>> &matrix);

        void matrix_dot_matrix_toeplitz_mod(vector<KernelInfo> &kernel_infos, CoeffIter c1, uint64_t poly_degree, vector<vector<uint64_t>> &result, Modulus &modulus);


        
        uint64_t kernel_innerprod(vector<pair<uint64_t, uint64_t>> rowinfo, CoeffIter coeff_vec, Modulus modulus);

        inline void kernel_matrix_dot_vector(vector<KernelInfo> kernel_infos, CoeffIter coeff_vec, Modulus modulus, CoeffIter result){
            uint64_t result_index = 0;
            for(uint64_t i = 0;i < kernel_infos.size();i++){
                KernelInfo kinfo = kernel_infos[i];
                vector<pair<uint64_t, uint64_t>> krow = kinfo.make_rowpair(modulus);
                for(uint64_t j = 0;j < kinfo.get_colsize();j++){
                    // inner_prod
                    result[result_index] = kernel_innerprod(krow, coeff_vec, modulus);
                    // next column
                    kinfo.pair_nextcol(krow, modulus);
                    result_index++;
                }
            }
        }

        inline void kernel_matrix_dot_vector(vector<KernelInfo> kernel_infos, uint64_t coeff_modulus_size, RNSIter poly_rns, ConstModulusIter mod_chain, RNSIter result){
            SEAL_ITERATE(iter(poly_rns, mod_chain, result), coeff_modulus_size, [&](auto I){
                    kernel_matrix_dot_vector(kernel_infos, get<0>(I), get<1>(I), get<2>(I));
                    });
        }

        inline void make_packedconv_matrixproduct(vector<KernelInfo> &kernel_infos, Ciphertext &encrypted, uint64_t poly_degree, vector<vector<uint64_t>> &result, Modulus modulus){
            PolyIter cipher_poly(encrypted);
            cipher_poly++;
            util::matrix_dot_matrix_toeplitz_mod(kernel_infos, **cipher_poly, poly_degree, result, modulus);
        }

        void packedconv_matrix_dot_vector(vector<vector<uint64_t>> &matrix_multed, vector<KernelInfo> kernelinfos, CoeffIter vector_iter, uint64_t vector_size, CoeffIter destination, const Modulus &modulus, MemoryPoolHandle pool_);

        inline void packedconv_matrix_dot_vector(vector<vector<uint64_t>> &matrix_multed, vector<KernelInfo> kernelinfos, RNSIter vector_rns, uint64_t rns_len, RNSIter destination, ConstModulusIter mod_chain, MemoryPoolHandle pool_){
            auto poly_degree = vector_rns.poly_modulus_degree();
            SEAL_ITERATE(iter(vector_rns, destination, mod_chain), rns_len, [&](auto I){
                    packedconv_matrix_dot_vector(matrix_multed, kernelinfos, get<0>(I), poly_degree, get<1>(I), get<2>(I), pool_);
                    });
        }

    }
}
