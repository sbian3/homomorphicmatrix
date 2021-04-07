#include "packedconv.h"

namespace seal
{
    namespace util
    {
        ////////////////////////////////////////////
        //
        // diagonal functions for packed convolution
        //
        ////////////////////////////////////////////

        // result list is reversed(index is also)
        vector<uint64_t> create_diagonal_list(vector<uint64_t> kernel, uint64_t colsize, uint64_t rowsize, Modulus &modulus, vector<uint64_t> &diagonal_list){
            vector<uint64_t> indexes;
            diagonal_list[colsize-1] = kernel[0];
            for(uint64_t i = 1;i< kernel.size();i++){
                diagonal_list[colsize-i-1] = kernel[i];
                diagonal_list[colsize + rowsize - 1 - i] = util::negate_uint_mod(kernel[i], modulus);
            }
            reverse(diagonal_list.begin(), diagonal_list.end());
            indexes.reserve(kernel.size());
            for(uint64_t i = 0;i < diagonal_list.size();i++){
                if(diagonal_list[i] != 0) indexes.push_back(i);
            }
            return indexes;
        }

        vector<uint64_t> create_diagonal_from_submatrix(CoeffIter a, uint64_t poly_degree, uint64_t start_col, uint64_t colsize, Modulus &modulus){
            vector<uint64_t> ret;
            uint64_t n = poly_degree;
            int64_t index = start_col + colsize-1;
            ret.reserve(n + colsize);
            // 縦へ移動
            for(;index > start_col;index--){
                ret.push_back(a[index]);
            }
            // 0へ移動
            uint64_t counter = 0;
            while(index >= 0){
                ret.push_back(a[index]);
                index--;
                counter++;
            }
            for(uint64_t i = 0;i< n-counter;i++){
                ret.push_back(util::negate_uint_mod(a[n-1-i], modulus));
            }
            return ret;
        }

        // assume kernel_L is reversed( indexes is also )
        // offsetは0が中心．
        vector<uint64_t> matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> kernel_L, vector<uint64_t> kernel_L_indexes, vector<uint64_t> list_R, Modulus & modulus){
            // RはL以上を想定
            assert(kernel_L_indexes.size() <= list_R.size());
        
            // 要素積を計算
            // kernelの非ゼロ要素の場所を覚えることで計算時間を短縮
            //uint64_t wise_prod_len = kernel_L.size() <= list_R.size()? kernel_L.size(): list_R.size();
            uint64_t wise_prod_len;
            if(offset >= 0){
                if(offset + kernel_L.size() > list_R.size()){
                    wise_prod_len = list_R.size() - offset;
                }else{
                    wise_prod_len = kernel_L.size(); 
                }
            }else{
                wise_prod_len = kernel_L.size() + offset;
            }
            vector<uint64_t> wise_prod(wise_prod_len);
            for(uint64_t i = 0;i < kernel_L_indexes.size();i++){
                uint64_t prod;
                if(offset >= 0){
                    // boarder check and mul
                    if(kernel_L_indexes[i]+offset >= list_R.size() || kernel_L_indexes[i] >= kernel_L.size()) break;
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                }else{
                    // boarder check and mul
                    if(kernel_L_indexes[i] <= -offset-1) continue;
                    if(kernel_L_indexes[i] >= kernel_L.size() || kernel_L_indexes[i]+offset >= list_R.size()) break;
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                }
                if(offset < 0){
                    wise_prod[kernel_L_indexes[i]+offset] = prod;
                }else{
                    wise_prod[kernel_L_indexes[i]] = prod;
                }
            }
            vector<uint64_t> diagonal;
            uint64_t partial_sum = 0;
            uint64_t prod_times = colsize_R;
            // 最初の内積を計算
            for(uint64_t i = 0;i < prod_times;i++){
                partial_sum = util::add_uint_mod(partial_sum, wise_prod[i], modulus);
            }
            diagonal.reserve(rowsize_R);
            diagonal.push_back(partial_sum);
            // 次の対角成分はO(1)で計算
            for(uint64_t i = 0;i < wise_prod.size() - prod_times;i++){
                partial_sum = util::add_uint_mod(partial_sum, wise_prod[i+prod_times], modulus);
                partial_sum = util::sub_uint_mod(partial_sum, wise_prod[i], modulus);
                diagonal.push_back(partial_sum);
            }
            return diagonal;
        }
        
        // 対角成分ベクトルを行列に書き込む
        // diagonallist: 左端からの対角成分ベクトル
        // resultは長方形を仮定
        void diagonallist_to_matrix(vector<vector<uint64_t>> diagonallist, uint64_t start_col, uint64_t start_row, uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result){
            assert(start_col + colsize <= result.size());
            assert(start_row + rowsize <= result[0].size()); 
            assert(diagonallist.size() == colsize + rowsize - 1);
            // (0, 0)までの対角成分
            for(uint64_t i = 0;i < colsize-1;i++){
                for(uint64_t j = 0;j < diagonallist[i].size();j++){
                    if(start_col + colsize - 1 -i + j >= result.size() || start_row + j >= result[0].size()){
                        cerr << "Error: diagonallist_to_matrix: out of range" << endl;
                        return;
                    }
                    result[start_col + colsize - 1 -i + j][start_row + j] = diagonallist[i][j];
                }
            }
            // (0, 0)から右の成分
            for(uint64_t i=colsize-1;i<diagonallist.size();i++){
                for(uint64_t j = 0;j < diagonallist[i].size();j++){
                    if(start_col + j >= result.size() || start_row + (i-colsize+1) + j >= result[0].size()){
                        cerr << "Error: diagonallist_to_matrix: out of range" << endl;
                        return;
                    }
                    result[start_col + j][start_row + (i-colsize+1) + j] = diagonallist[i][j];
                }
            }
        }

        vector<vector<uint64_t>> scalars_to_diagonallist(vector<uint64_t> scalars, uint64_t colsize, uint64_t rowsize){
            vector<vector<uint64_t>> diagonals(colsize + rowsize - 1);
            // 小さいほうの長さがmax
            uint64_t maxlength = colsize >= rowsize? rowsize: colsize;
            uint64_t keeplength = colsize >= rowsize? colsize-rowsize: rowsize-colsize;
            // 縦(最上も含む)
            uint64_t i = 0;
            for(;i < colsize;i++){
                for(uint64_t j = 0;j < i+1;j++){
                    diagonals[i].push_back(scalars[i]);
                }
            }
            // キープ
            for(;i < colsize + keeplength;i++){
                for(uint64_t j = 0;j < maxlength;j++){
                    diagonals[i].push_back(scalars[i]);
                }
            }
            // しぼむ
            while(maxlength > 0){
                maxlength--;
                for(uint64_t j = 0;j < maxlength;j++){
                    diagonals[i].push_back(scalars[i]);
                }
                i++;
            }
            return diagonals;
        }

    }
}
