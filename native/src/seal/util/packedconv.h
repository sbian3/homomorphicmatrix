
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
#include "seal/util/uintlinarith.h"
#include "seal/util/uintarithmod.h"

using namespace std;

namespace seal
{
    namespace util
    {

        class KernelInfo{
            public:
                KernelInfo(){

                }
                KernelInfo(uint64_t st_c,uint64_t  st_r,uint64_t si_c,uint64_t si_r,vector<uint64_t> data, Modulus mod):
                    start_col(st_c), start_row(st_r), size_col(si_c), size_row(si_r), data(data), modulus(mod){
                        diagonal_list = vector<uint64_t>(size_col + size_row - 1);
                        index = create_diagonal_list(data, size_col, size_row, modulus, diagonal_list);
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
                        if(next_index < 0 || next_index > start_row + size_row){
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

        vector<uint64_t> create_diagonal_list(vector<uint64_t> kernel, uint64_t colsize, uint64_t rowsize, Modulus &modulus, vector<uint64_t> &diagonal_list);
        vector<uint64_t> create_diagonal_from_submatrix(CoeffIter a, uint64_t poly_degree, uint64_t start_col, uint64_t colsize, Modulus &modulus);
        vector<uint64_t> matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> kernel_L, vector<uint64_t> kernel_L_indexes, vector<uint64_t> list_R, Modulus & modulus);
        void diagonallist_to_matrix(vector<vector<uint64_t>> diagonallist, uint64_t start_col, uint64_t start_row, uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result);
        vector<vector<uint64_t>> scalars_to_diagonallist(vector<uint64_t> scalars, uint64_t colsize, uint64_t rowsize);
    }
}
