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
        vector<uint64_t> matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus){
            bool print_debug = false;
            bool print_time = false;
            // assert list_R is larger than kernel
            assert(kernel_L_indexes.size() <= list_R.size());
        
            // calculate element wise product
            // optimize complexity remembering kernel nonzero elements
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
            //cout << "wise_prod_len: " << wise_prod_len << endl;
            vector<uint64_t> wise_prod(wise_prod_len);
            // non-zero index
            vector<uint64_t> wise_prod_index;
            auto mul_start = chrono::high_resolution_clock::now();
            for(uint64_t i = 0;i < kernel_L_indexes.size();i++){
                uint64_t prod;
                if(offset >= 0){
                    // boarder check and mul
                    if(kernel_L_indexes[i]+offset >= list_R.size() || kernel_L_indexes[i] >= kernel_L.size()) break;
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                    wise_prod[kernel_L_indexes[i]] = prod;
                    wise_prod_index.push_back(kernel_L_indexes[i]);
                }else{
                    // boarder check and mul
                    if(kernel_L_indexes[i] <= -offset-1) continue;
                    if(kernel_L_indexes[i] >= kernel_L.size() || kernel_L_indexes[i]+offset >= list_R.size()) break;
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                    wise_prod[kernel_L_indexes[i]+offset] = prod;
                    wise_prod_index.push_back(kernel_L_indexes[i] + offset);
                }
            }
            auto mul_end = chrono::high_resolution_clock::now();
            uint64_t prod_times = colsize_R;
            uint64_t partial_sum = 0;
            uint64_t index_of_index = 0;
            // calculate first inner prod
            for(uint64_t i = 0; i < wise_prod_index.size() && wise_prod_index[i] < prod_times;i++){
                partial_sum = util::add_uint_mod(partial_sum, wise_prod[wise_prod_index[i]], modulus);
                index_of_index++;
            }
            auto innerp_end = chrono::high_resolution_clock::now();
            if(print_debug){
                cout << "---------calculating a diagonal-----" << endl;
                cout << "index_of_index: " << index_of_index << endl;
                cout << "wise_prod.size: " << wise_prod.size() << endl;;
                for(uint64_t i = 0; i < wise_prod_index.size();i++){
                    cout << "index " << wise_prod_index[i] << ": " << wise_prod[wise_prod_index[i]] << endl;;
                }
                cout << "prod_times: " << prod_times << endl;
                cout << "end index of nonzero wise_prod: " << wise_prod_index[wise_prod_index.size()- 1] << endl;
            }
            // all elements end in first inner prod(constant vector)
            if(index_of_index  == wise_prod_index.size() && wise_prod_index[wise_prod_index.size()-1] < prod_times){
                if(print_debug){
                    cout << "end in first inner prod: return" << endl;
                }
                vector<uint64_t> diagonal(wise_prod_len - prod_times + 1, partial_sum);
                return diagonal;
            }
            vector<uint64_t> diagonal(wise_prod_len - prod_times + 1);
            diagonal[0] = partial_sum;
            // we need O(1) to calc a next diagonal element
            // TODO: skip constant slide
            for(uint64_t i = 0;i < wise_prod.size() - prod_times;i++){
                //cout << "i: " << i << endl;
                //cout << "index of index: " << index_of_index << endl;
                if(partial_sum == 0){
                    if(index_of_index == wise_prod_index.size() - 1){
                        //cout << "index of index: out" << endl;
                        break;
                    }
                    if(wise_prod_index[index_of_index] - i > prod_times ){
                        i = wise_prod_index[index_of_index] - prod_times;
                        //cout << "jumped" << endl;
                    }
                }
                // move window to right
                if(wise_prod[i+prod_times] != 0){
                    index_of_index++;
                }
                partial_sum = util::add_uint_mod(partial_sum, wise_prod[i+prod_times], modulus);
                partial_sum = util::sub_uint_mod(partial_sum, wise_prod[i], modulus);
                //cout << "partial sum: " << partial_sum << endl;
                diagonal[i+1] = partial_sum;
            }
            auto slide_end = chrono::high_resolution_clock::now();
            auto mul_diff = chrono::duration_cast<chrono::nanoseconds>(mul_end - mul_start);
            auto innerp_diff = chrono::duration_cast<chrono::nanoseconds>(innerp_end - mul_end);
            auto slide_diff = chrono::duration_cast<chrono::nanoseconds>(slide_end - innerp_end);
            cout << "mul : " << mul_diff.count() << " innerp: " << innerp_diff.count() << " slide: " << slide_diff.count() << " sum: " << mul_diff.count() + innerp_diff.count() + slide_diff.count() << endl;
            //cout << "------expected-----" << endl;
            //// expected loop
            //for(uint64_t i = 0;i < wise_prod.size() - prod_times;i++){
            //    cout << "i: " << i << endl;
            //    // move window to right
            //    partial_sum_expect = util::add_uint_mod(partial_sum_expect, wise_prod[i+prod_times], modulus);
            //    partial_sum_expect = util::sub_uint_mod(partial_sum_expect, wise_prod[i], modulus);
            //    cout << "partial sum: " << partial_sum_expect << endl;
            //    diagonal[i+1] = partial_sum_expect;
            //}
            //cout << endl;
            if(print_debug){
                util::print_vector(diagonal, diagonal.size());
            }
            return diagonal;
        }
        
        void matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus, vector<pair<uint64_t, uint64_t>> &diagonalpairlist){
            // assert list_R is larger than kernel
            assert(kernel_L_indexes.size() <= list_R.size());
        
            // calculate element wise product
            // optimize complexity remembering kernel nonzero elements
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
            // non-zero index list of wise_prod
            vector<uint64_t> wise_prod_index;
#if TIFS_DEBUG_TIME == 1
            auto mul_start = chrono::high_resolution_clock::now();
#endif
            for(uint64_t i = 0;i < kernel_L_indexes.size();i++){
                uint64_t prod;
                if(offset >= 0){
                    // boarder check and mul
                    if(kernel_L_indexes[i]+offset >= list_R.size() || kernel_L_indexes[i] >= kernel_L.size()) break;
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                    wise_prod[kernel_L_indexes[i]] = prod;
                    wise_prod_index.push_back(kernel_L_indexes[i]);
                }else{
                    // boarder check and mul
                    if(kernel_L_indexes[i] <= -offset-1) continue;
                    if(kernel_L_indexes[i] >= kernel_L.size() || kernel_L_indexes[i]+offset >= list_R.size()) break;
                    prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
                    wise_prod[kernel_L_indexes[i]+offset] = prod;
                    wise_prod_index.push_back(kernel_L_indexes[i] + offset);
                }
            }
#if TIFS_DEBUG_TIME
            auto mul_end = chrono::high_resolution_clock::now();
#endif
            uint64_t prod_times = colsize_R;
            uint64_t partial_sum = 0;
            uint64_t index_of_index_right = 0;
            // calculate first inner prod
            for(uint64_t i = 0; i < wise_prod_index.size() && wise_prod_index[i] < prod_times;i++){
                partial_sum = util::add_uint_mod(partial_sum, wise_prod[wise_prod_index[i]], modulus);
                index_of_index_right++;
            }
#if TIFS_DEBUG_TIME
            auto innerp_end = chrono::high_resolution_clock::now();
#endif
#if TIFS_DEBUG
            cout << "---------calculating a diagonal-----" << endl;
            cout << "index_of_index: " << index_of_index_right << endl;
            cout << "wise_prod.size: " << wise_prod.size() << endl;;
            for(uint64_t i = 0; i < wise_prod_index.size();i++){
                cout << "index " << wise_prod_index[i] << ": " << wise_prod[wise_prod_index[i]] << endl;;
            }
            cout << "prod_times: " << prod_times << endl;
            cout << "end index of nonzero wise_prod: " << wise_prod_index[wise_prod_index.size()- 1] << endl;
#endif
            uint64_t return_len = wise_prod_len - prod_times + 1;
            // all elements end in first inner prod(constant vector)
            if(index_of_index_right  == wise_prod_index.size() && wise_prod_index[wise_prod_index.size()-1] < prod_times){
#if TIFS_DEBUG
                cout << "end in first inner prod: return" << endl;
#endif
                diagonalpairlist.push_back(make_pair(partial_sum, return_len));
                return;
            }
            uint64_t index_of_index_left = 0;
            // slide window
            // we need O(1) to calc a next diagonal element
            uint64_t jump_upper;
            uint64_t jump_bottom;
            bool nonzeroleft_over = false;
            uint64_t pair_num = 0;
            for(uint64_t i = 0;i < wise_prod.size() - prod_times+1;i++){
                // how many index should we jump to calc next partial_sum?
                jump_upper = wise_prod_index[index_of_index_right] - (i+prod_times-1);
                jump_bottom = wise_prod_index[index_of_index_left] - i + 1;
                if(nonzeroleft_over){
                    jump_upper = wise_prod.size() - (i + prod_times - 1);
                }
                uint64_t jump_len = min(jump_upper, jump_bottom);
#if TIFS_DEBUG == 1
                cout << "--loop: i=" << i << "--" << endl;
                cout << "index_of_index: (" << index_of_index_left << ", " << index_of_index_right << ") " << endl;
                cout << "windows: (" << i << ", " << i + prod_times - 1 << ")" << endl;
                cout << "jump_upper: " << jump_upper << endl;
                cout << "jump_bottom: " << jump_bottom << endl;
                cout << "jump_len: " << jump_len << endl;
                cout << "partial_sum: " << partial_sum << endl;
                cout << "make pair: " << partial_sum << ", " << jump_len << endl;
#endif
                auto pair = make_pair(partial_sum, jump_len);
                pair_num++;
                diagonalpairlist.push_back(pair);
                if(jump_len == jump_upper){
                    partial_sum = util::add_uint_mod(partial_sum, wise_prod[wise_prod_index[index_of_index_right]], modulus);
                    if(index_of_index_right == wise_prod_index.size()-1){
#if TIFS_DEBUG == 1
                        cout << "index_right is edge!" << endl;
#endif
                        nonzeroleft_over = true;
                    }else{
                        index_of_index_right++;
                    }
                }
                if(jump_len == jump_bottom){
                    partial_sum = util::sub_uint_mod(partial_sum, wise_prod[wise_prod_index[index_of_index_left]], modulus);
                    index_of_index_left++;
                }
                i = i + jump_len - 1;
            }
            //cout << "pair_num: " << pair_num << endl;
#if TIFS_DEBUG_TIME
            auto slide_end = chrono::high_resolution_clock::now();
            auto mul_diff = chrono::duration_cast<chrono::nanoseconds>(mul_end - mul_start);
            auto innerp_diff = chrono::duration_cast<chrono::nanoseconds>(innerp_end - mul_end);
            auto slide_diff = chrono::duration_cast<chrono::nanoseconds>(slide_end - innerp_end);
            cout << "mul : " << mul_diff.count() << " innerp: " << innerp_diff.count() << " slide: " << slide_diff.count() << " sum: " << mul_diff.count() + innerp_diff.count() + slide_diff.count() << endl;
#endif
        }

        // 対角成分ベクトルを行列に書き込む
        // diagonallist: 左端からの対角成分ベクトル
        // resultは長方形を仮定
        void diagonallist_to_matrix(vector<vector<uint64_t>> &diagonallist, uint64_t start_col, uint64_t start_row, uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result){
            assert(start_col + colsize <= result.size());
            assert(start_row + rowsize <= result[0].size()); 
            assert(diagonallist.size() == colsize + rowsize - 1);
            // from (start_col + colsize - 1, start_row ) to (start_col, start_row)
            for(uint64_t i = 0;i < colsize-1;i++){
                for(uint64_t j = 0;j < diagonallist[i].size();j++){
                    if(start_col + colsize - 1 -i + j >= result.size() || start_row + j >= result[0].size()){
                        cerr << "Error: diagonallist_to_matrix: out of range" << endl;
                        return;
                    }
                    result[start_col + colsize - 1 -i + j][start_row + j] = diagonallist[i][j];
                }
            }
            // from (start_col, start_row) to (start_col, start_row + rowsize - 1)
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

        void diagonallist_to_matrix(vector<vector<pair<uint64_t, uint64_t>>> &diagonallist, uint64_t start_col, uint64_t start_row , uint64_t colsize, uint64_t rowsize, vector<vector<uint64_t>> &result){
            assert(start_col + colsize <= result.size());
            assert(start_row + rowsize <= result[0].size()); 
            assert(diagonallist.size() == colsize + rowsize - 1);
            // from (start_col + colsize - 1, start_row ) to (start_col, start_row)
            for(uint64_t i = 0;i < colsize-1;i++){
                // iterate all pairs in vector[i]
                uint64_t index_resultdiagonal = 0;
                for(uint64_t j = 0;j < diagonallist[i].size();j++){
                    // pair to result matrix
                    auto pair = diagonallist[i][j];
                    auto pair_value = pair.first;
                    auto pair_valuelen = pair.second;
                    for(uint64_t k = 0;k < pair_valuelen;k++){
                        if(start_col + colsize - 1 -i + index_resultdiagonal >= result.size() || start_row + index_resultdiagonal >= result[0].size()){
                            cerr << "Error: diagonallist_to_matrix: out of range" << endl;
                            return;
                        }
                        result[start_col + colsize - 1 -i + index_resultdiagonal][start_row + index_resultdiagonal] = pair_value;
                        index_resultdiagonal++;
                    }
                }
            }
            // from (start_col, start_row) to (start_col, start_row + rowsize - 1)
            for(uint64_t i=colsize-1;i<diagonallist.size();i++){
                uint64_t index_resultdiagonal = 0;
                for(uint64_t j = 0;j < diagonallist[i].size();j++){
                    // pair to result matrix
                    auto pair = diagonallist[i][j];
                    auto pair_value = pair.first;
                    auto pair_valuelen = pair.second;
                    for(uint64_t k = 0;k < pair_valuelen;k++){
                        if(start_col + index_resultdiagonal >= result.size() || start_row + (i-colsize+1)+ index_resultdiagonal >= result[0].size()){
                            cerr << "Error: diagonallist_to_matrix: out of range" << endl;
                            return;
                        }
                        result[start_col + index_resultdiagonal][start_row + (i-colsize+1)+ index_resultdiagonal] = pair_value;
                        index_resultdiagonal++;
                    }
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

        // matrix dot matrix product using toeplitz algorithm
        void matrix_dot_matrix_toeplitz_mod(vector<KernelInfo> &kernel_infos, CoeffIter c1, uint64_t poly_degree, vector<vector<uint64_t>> &result, Modulus &modulus){
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
                vector<vector<pair<uint64_t, uint64_t>>> matrix_product_diagonals(colsize_K + submat_rowsize - 1);
                uint64_t index = 0;
                int64_t k = static_cast<int64_t>(colsize_K);
                k = -k+1;
                for(;k<static_cast<int64_t>(submat_rowsize);k++){
                    vector<pair<uint64_t, uint64_t>> diagonal_pairs;
#if TIFS_DEBUG_TIME == 1
                    auto diagonal_start = chrono::high_resolution_clock::now();
#endif
                    util::matrix_product_diagonal(k, submat_colsize, submat_rowsize, kernel_diagonal_list, kernel_index, diagonal_c1, modulus, diagonal_pairs);
#if TIFS_DEBUG_TIME == 1
                    auto diagonal_end = chrono::high_resolution_clock::now();
                    auto diagonal_diff = chrono::duration_cast<chrono::nanoseconds>(diagonal_end - diagonal_start);
                    cout << "calc one diagonal vector: " << diagonal_diff.count() << endl;
#endif
                    matrix_product_diagonals[index] = diagonal_pairs;
                    index++;
                }
                //cout << "calc diagonals: " << lt_diff.count() << " ms" <<endl;

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
