
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

using namespace std;

namespace seal
{
    namespace util
    {

        vector<uint64_t> create_diagonal_list(vector<uint64_t> kernel, uint64_t colsize, uint64_t rowsize, Modulus &modulus, vector<uint64_t> &diagonal_list);

        class KernelInfo{
            public:
                uint64_t input_size;
                uint64_t block_size;
                vector<uint64_t> toeplitz;
                KernelInfo(){

                }
                KernelInfo(uint64_t input_s, uint64_t block_s, uint64_t st_c,uint64_t  st_r,uint64_t si_c,uint64_t si_r,vector<uint64_t> dat, Modulus mod):
                    input_size(input_s), block_size(block_s), start_col(st_c), start_row(st_r), size_col(si_c), size_row(si_r), data(dat), modulus(mod){
                        diagonal_list = vector<uint64_t>(size_col + size_row - 1);
                        index = create_diagonal_list(data, size_col, size_row, modulus, diagonal_list);
                        //util::print_vector(diagonal_list, diagonal_list.size());
                    }
                vector<uint64_t> diagonal_list;
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

                // 行列ベクトル積計算に用いる
                // pair<index, value>のベクトルを返す
                vector<pair<uint64_t, uint64_t>> make_rowpair(Modulus mod){
                    vector<pair<uint64_t, uint64_t>> ret(data.size());
                    // 最初は通常の値
                    ret[0] = make_pair(start_row, data[0]);
                    // 2つ目以降はnegateして右端へ
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

                void get_toeplitz(vector<vector<pair<uint64_t, uint64_t>>>& pairs, uint64_t poly_degree){
                    uint64_t toeplitz_len = input_size + poly_degree - 1;
                    vector<uint64_t> toeplitz_tmp(toeplitz_len);
                    for(uint64_t i = 0;i < toeplitz_len;i++){
                        toeplitz_tmp[i] = pairs[i].back().first;
                    }
                    toeplitz = toeplitz_tmp;
                }

                void print(){
                    cout << "start col: " << start_col << endl;
                    cout << "start row: " << start_row << endl;
                    cout << "diagonal: ";
                    for(uint64_t k = 0;k < diagonal_list.size();k++){
                        cout << diagonal_list[k] << " ";
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

        vector<uint64_t> matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus);

        void matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus, vector<pair<uint64_t, uint64_t>> &diagonalpairlist);

        void diagonallist_to_matrix(vector<vector<uint64_t>> &diagonallist, uint64_t start_col, uint64_t start_row, uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result);

        void diagonallist_to_matrix(vector<vector<pair<uint64_t, uint64_t>>> &diagonallist, uint64_t start_col, uint64_t start_row , uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result);

        vector<vector<uint64_t>> scalars_to_diagonallist(vector<uint64_t> scalars, uint64_t colsize, uint64_t rowsize);
        inline uint64_t get_blocksize(uint64_t input_dim, uint64_t kernel_dim, uint64_t padding){
            return input_dim + kernel_dim - 1 + padding;
        }

        vector<uint64_t> pack_input(const vector<vector<uint64_t>> input, uint64_t block_size, uint64_t poly_size);

        vector<KernelInfo> pack_kernel(vector<vector<uint64_t>> kernels, uint64_t input_size, Modulus modulus);

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

    }
}
