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


        vector<uint64_t> pack_input(const vector<vector<uint64_t>> input, uint64_t block_size, uint64_t poly_size){
            assert(input.size() * block_size <= poly_size);
            vector<uint64_t> packed_input(poly_size);
            for(uint64_t i = 0;i < input.size();i++){
                uint64_t input_start = block_size * i;
                for(uint64_t j = 0;j < input[i].size();j++){
                    packed_input[input_start + j] = input[i][j];
                }
            }
            return packed_input;
        }

        // kernel: 初期化したkernel vectorのvector. 長さはpack_num個
        // block_size: 1ブロックの大きさ．
        // return: KernelInfoのベクトルを生成して返す
        vector<KernelInfo> pack_kernel(vector<vector<uint64_t>> kernels, uint64_t block_size, Modulus modulus){
            uint64_t packing_num = kernels.size();
            vector<KernelInfo> kernel_info(packing_num);
            uint64_t start_col = 0;
            uint64_t start_row = 0;
            for(uint64_t i = 0;i < packing_num;i++){
                KernelInfo kinfo(start_col, start_row, block_size, block_size, kernels[i], modulus);
                kernel_info[i] = kinfo;
                start_col += block_size;
                start_row += block_size;
            }
            return kernel_info;
        }

        void pack_kernel_to_matrix(vector<KernelInfo> kernelinfos, vector<vector<uint64_t>> &matrix){
            for(uint64_t i = 0;i < kernelinfos.size();i++){
                KernelInfo kinfo = kernelinfos[i];
                // kernel scalar is reversed in default
                // so, scalar must be reversed again
                vector<uint64_t> k_scalar = kinfo.diagonal_list;
                reverse(k_scalar.begin(), k_scalar.end());
                vector<vector<uint64_t>> diagonals = util::scalars_to_diagonallist(k_scalar, kinfo.get_colsize(), kinfo.get_rowsize());
                util::diagonallist_to_matrix(diagonals, kinfo.get_startcol(), kinfo.get_startrow(), kinfo.get_colsize(), kinfo.get_rowsize(),matrix);
            }
        }

        // 行列積結果を対角成分のみから計算する．
        void matrix_dot_matrix_toeplitz_mod(vector<KernelInfo> kernel_infos, CoeffIter c1, uint64_t poly_degree, vector<vector<uint64_t>> &result, Modulus &modulus){
            // for each block
            for(uint64_t i = 0;i < kernel_infos.size();i++){
                // get diagonal lists for kernel
                KernelInfo kinfo = kernel_infos[i];
                vector<uint64_t> kernel_diagonal_list = kinfo.diagonal_list;
                vector<uint64_t> kernel_index = kinfo.index;
                uint64_t colsize_K = kinfo.get_colsize();

                // diagonal list of c1
                uint64_t submat_startcol,submat_startrow, submat_colsize, submat_rowsize;
                submat_rowsize = poly_degree;
                submat_startrow = 0;
                kinfo.getParamsforSubmatrix(submat_startcol, submat_colsize);
                vector<uint64_t> diagonal_c1 = util::create_diagonal_from_submatrix(c1, poly_degree , submat_startcol, submat_colsize, modulus);

                // calc diagonal of product
                vector<vector<uint64_t>> matrix_product_diagonals(colsize_K + submat_rowsize - 1);
                uint64_t index = 0;
                int64_t k = static_cast<int64_t>(colsize_K);
                k = -k+1;
                for(;k<static_cast<int64_t>(submat_rowsize);k++){
                    vector<uint64_t> diagonal_vec;
                    diagonal_vec = util::matrix_product_diagonal(k, submat_colsize, submat_rowsize, kernel_diagonal_list, kernel_index, diagonal_c1, modulus);
                    matrix_product_diagonals[index] = diagonal_vec;
                    index++;
                }

                // write diagonals to result matrix
                util::diagonallist_to_matrix(matrix_product_diagonals, submat_startcol, submat_startrow, colsize_K, submat_rowsize, result);
            }
        }

        uint64_t kernel_innerprod(vector<pair<uint64_t, uint64_t>> rowinfo, CoeffIter coeff_vec, Modulus modulus){
            uint64_t sum = 0;
            for(uint64_t i = 0;i < rowinfo.size();i++){
                pair<uint64_t, uint64_t> row = rowinfo[i];
                sum = multiply_add_uint_mod(row.second, coeff_vec[row.first],sum, modulus);
            }
            return sum;
        }


    }
}