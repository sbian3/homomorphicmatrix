#include "packedconv.h"
#include "define_tifs.h"

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
        vector<uint64_t> create_diagonal_scalars(const vector<uint64_t> &kernel, const uint64_t colsize, const uint64_t rowsize, const Modulus &modulus, vector<uint64_t> &diagonal_list){
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
            // get elements from lower left
            for(;index > start_col;index--){
                ret.push_back(a[index]);
            }
            // go to diagonal element
            uint64_t counter = 0;
            while(index >= 0){
                ret.push_back(a[index]);
                index--;
                counter++;
            }
            // go to upper right
            for(uint64_t i = 0;i< n-counter;i++){
                ret.push_back(util::negate_uint_mod(a[n-1-i], modulus));
            }
            return ret;
        }

        // WARN: currently not used
        // assume kernel_L is reversed( indexes is also )
        // offset: 0 is center of diagonal
        //vector<uint64_t> matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus){
        //    bool print_debug = false;
        //    // assert list_R is larger than kernel
        //    assert(kernel_L_indexes.size() <= list_R.size());
        //
        //    // calculate element wise product
        //    // optimize complexity remembering kernel nonzero elements
        //    //uint64_t wise_prod_len = kernel_L.size() <= list_R.size()? kernel_L.size(): list_R.size();
        //    uint64_t wise_prod_len;
        //    if(offset >= 0){
        //        if(offset + kernel_L.size() > list_R.size()){
        //            wise_prod_len = list_R.size() - offset;
        //        }else{
        //            wise_prod_len = kernel_L.size(); 
        //        }
        //    }else{
        //        wise_prod_len = kernel_L.size() + offset;
        //    }
        //    //cout << "wise_prod_len: " << wise_prod_len << endl;
        //    vector<uint64_t> wise_prod(wise_prod_len);
        //    // non-zero index
        //    vector<uint64_t> wise_prod_index;
        //    auto mul_start = chrono::high_resolution_clock::now();
        //    for(uint64_t i = 0;i < kernel_L_indexes.size();i++){
        //        uint64_t prod;
        //        if(offset >= 0){
        //            // boarder check and mul
        //            if(kernel_L_indexes[i]+offset >= list_R.size() || kernel_L_indexes[i] >= kernel_L.size()) break;
        //            prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
        //            wise_prod[kernel_L_indexes[i]] = prod;
        //            wise_prod_index.push_back(kernel_L_indexes[i]);
        //        }else{
        //            // boarder check and mul
        //            if(kernel_L_indexes[i] <= -offset-1) continue;
        //            if(kernel_L_indexes[i] >= kernel_L.size() || kernel_L_indexes[i]+offset >= list_R.size()) break;
        //            prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
        //            wise_prod[kernel_L_indexes[i]+offset] = prod;
        //            wise_prod_index.push_back(kernel_L_indexes[i] + offset);
        //        }
        //    }
        //    auto mul_end = chrono::high_resolution_clock::now();
        //    uint64_t prod_times = colsize_R;
        //    uint64_t partial_sum = 0;
        //    uint64_t index_of_index = 0;
        //    // calculate first inner prod
        //    for(uint64_t i = 0; i < wise_prod_index.size() && wise_prod_index[i] < prod_times;i++){
        //        partial_sum = util::add_uint_mod(partial_sum, wise_prod[wise_prod_index[i]], modulus);
        //        index_of_index++;
        //    }
        //    auto innerp_end = chrono::high_resolution_clock::now();
        //    if(print_debug){
        //        cout << "---------calculating a diagonal-----" << endl;
        //        cout << "index_of_index: " << index_of_index << endl;
        //        cout << "wise_prod.size: " << wise_prod.size() << endl;;
        //        for(uint64_t i = 0; i < wise_prod_index.size();i++){
        //            cout << "index " << wise_prod_index[i] << ": " << wise_prod[wise_prod_index[i]] << endl;;
        //        }
        //        cout << "prod_times: " << prod_times << endl;
        //        cout << "end index of nonzero wise_prod: " << wise_prod_index[wise_prod_index.size()- 1] << endl;
        //    }
        //    // all elements end in first inner prod(constant vector)
        //    if(index_of_index  == wise_prod_index.size() && wise_prod_index[wise_prod_index.size()-1] < prod_times){
        //        if(print_debug){
        //            cout << "end in first inner prod: return" << endl;
        //        }
        //        vector<uint64_t> diagonal(wise_prod_len - prod_times + 1, partial_sum);
        //        return diagonal;
        //    }
        //    vector<uint64_t> diagonal(wise_prod_len - prod_times + 1);
        //    diagonal[0] = partial_sum;
        //    // we need O(1) to calc a next diagonal element
        //    // TODO: skip constant slide
        //    for(uint64_t i = 0;i < wise_prod.size() - prod_times;i++){
        //        //cout << "i: " << i << endl;
        //        //cout << "index of index: " << index_of_index << endl;
        //        if(partial_sum == 0){
        //            if(index_of_index == wise_prod_index.size() - 1){
        //                //cout << "index of index: out" << endl;
        //                break;
        //            }
        //            if(wise_prod_index[index_of_index] - i > prod_times ){
        //                i = wise_prod_index[index_of_index] - prod_times;
        //                //cout << "jumped" << endl;
        //            }
        //        }
        //        // move window to right
        //        if(wise_prod[i+prod_times] != 0){
        //            index_of_index++;
        //        }
        //        partial_sum = util::add_uint_mod(partial_sum, wise_prod[i+prod_times], modulus);
        //        partial_sum = util::sub_uint_mod(partial_sum, wise_prod[i], modulus);
        //        //cout << "partial sum: " << partial_sum << endl;
        //        diagonal[i+1] = partial_sum;
        //    }
        //    auto slide_end = chrono::high_resolution_clock::now();
        //    auto mul_diff = chrono::duration_cast<chrono::nanoseconds>(mul_end - mul_start);
        //    auto innerp_diff = chrono::duration_cast<chrono::nanoseconds>(innerp_end - mul_end);
        //    auto slide_diff = chrono::duration_cast<chrono::nanoseconds>(slide_end - innerp_end);
        //    cout << "mul : " << mul_diff.count() << " innerp: " << innerp_diff.count() << " slide: " << slide_diff.count() << " sum: " << mul_diff.count() + innerp_diff.count() + slide_diff.count() << endl;
        //    //cout << "------expected-----" << endl;
        //    //// expected loop
        //    //for(uint64_t i = 0;i < wise_prod.size() - prod_times;i++){
        //    //    cout << "i: " << i << endl;
        //    //    // move window to right
        //    //    partial_sum_expect = util::add_uint_mod(partial_sum_expect, wise_prod[i+prod_times], modulus);
        //    //    partial_sum_expect = util::sub_uint_mod(partial_sum_expect, wise_prod[i], modulus);
        //    //    cout << "partial sum: " << partial_sum_expect << endl;
        //    //    diagonal[i+1] = partial_sum_expect;
        //    //}
        //    //cout << endl;
        //    if(print_debug){
        //        util::print_vector(diagonal, diagonal.size());
        //    }
        //    return diagonal;
        //}
        
        // kernel_L: diagona_scalars
        // kernel_L_indexes: index of nonezero element in kernel_L
        // case1: offset > 0
        //           | -----kernel-------|
        // |-------list_R------------|
        //           |-----product---|
        //
        // case2: offset < 0
        // |---kernel---------|
        //       |-----list_R--------|
        //       |---product--|
        //void matrix_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> &kernel_L, vector<uint64_t> &kernel_L_indexes, vector<uint64_t> &list_R, Modulus & modulus, vector<pair<uint64_t, uint64_t>> &diagonalpairlist){
        //    // assert list_R is larger than kernel
        //    assert(kernel_L_indexes.size() <= list_R.size());
        //    //assert(colsize_R <= kernel_L.size());
        //
#if HLT_//DEBUG_PRINT == 1
        //    cout << "---------calculating a diagonal-----" << endl;
        //    cout << "kernel: " << endl;
        //    print_vector(kernel_L, kernel_L.size());
        //    cout << "list_R:"  << endl;
        //    print_vector(list_R, list_R.size());
        //    cout << "offset: " << offset << endl;
#endif
        //    // calculate element wise product
        //    // optimize complexity remembering kernel nonzero elements
        //    //uint64_t wise_prod_len = kernel_L.size() <= list_R.size()? kernel_L.size(): list_R.size();
#if HLT_//DEBUG_TIME == 1
        //    auto diagonal_begin = chrono::high_resolution_clock::now();
#endif
#if HLT_//DEBUG_TIME == 1
        //    auto mul_start = chrono::high_resolution_clock::now();
#endif
#if HLT_//DEBUG_TIME == 1
        //    auto mul_end = chrono::high_resolution_clock::now();
#endif
        //    uint64_t wise_prod_len;
        //    uint64_t innerp_size = colsize_R;
        //    uint64_t first_right_edge;
        //    bool end_in_firstinner = false;
        //    if(offset >= 0){
        //        if(offset + kernel_L.size() > list_R.size()){
        //            wise_prod_len = list_R.size() - offset;
        //        }else{
        //            wise_prod_len = kernel_L.size(); 
        //        }
        //    }else{
        //        wise_prod_len = kernel_L.size() + offset;
        //    }
        //    if (wise_prod_len >= innerp_size) {
        //      first_right_edge = innerp_size;
        //    } else {
        //      cout << "end in first inner" << endl;
        //      first_right_edge = wise_prod_len;
        //      end_in_firstinner = true;
        //    }
        //    if(offset < 0){
        //        first_right_edge -= offset;
        //    }
        //    uint64_t partial_sum = 0;
        //    uint64_t index_iterator_right = 0;
        //    uint64_t index_iterator_left = 0;
        //    // calculate first inner prod
        //    for(uint64_t i = 0; i < kernel_L_indexes.size() && kernel_L_indexes[i] < first_right_edge;i++){
        //        uint64_t prod;
        //        if(offset >= 0){
        //            prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
        //        }else{
        //            // boarder check and mul
        //            if(kernel_L_indexes[i] <= -offset-1){
        //              index_iterator_right++;
        //              index_iterator_left++;
        //              continue;
        //            }
        //            prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
        //        }
        //        partial_sum = util::add_uint_mod(partial_sum, prod, modulus);
        //        index_iterator_right++;
        //    }
#if HLT_//DEBUG_TIME == 1
        //    auto innerp_end = chrono::high_resolution_clock::now();
#endif
#if HLT_//DEBUG_PRINT == 1
        //    cout << "index_iterator_right: " << index_iterator_right << endl;
        //    cout << "wise_prod_len: " << wise_prod_len << endl;;
        //    //for(uint64_t i = 0; i < wise_prod_index.size();i++){
        //    //    cout << "index " << wise_prod_index[i] << ": " << wise_prod[wise_prod_index[i]] << endl;;
        //    //}
        //    cout << "prod_times: " << innerp_size << endl;
        //    //cout << "end index of nonzero wise_prod: " << wise_prod_index[wise_prod_index.size()- 1] << endl;
#endif

        //    // slide window
        //    // we need O(1) to calc a next diagonal element
        //    uint64_t jump_right;
        //    uint64_t jump_left;
        //    bool is_right_edge = false;
        //    bool is_left_edge = false;
        //    uint64_t pair_num = 0;
        //    uint64_t i;
        //    if(offset < 0){
        //        i = -offset;
        //        wise_prod_len += i;
        //    }else{
        //        i = 0;
        //    }
        //    for (; i + innerp_size - 1 < wise_prod_len; i++) {
        //        // how many index should we jump to calc next partial_sum?
        //        uint64_t window_left = i;
        //        uint64_t window_right = i + innerp_size - 1;
        //        jump_right =
        //            kernel_L_indexes[index_iterator_right] - window_right;
        //        jump_left = kernel_L_indexes[index_iterator_left] - window_left + 1;
        //        if (is_right_edge) {
        //            jump_right = wise_prod_len - window_right;
        //        }
        //        if(is_left_edge){
        //            jump_left = wise_prod_len - window_left;
        //        }
        //        uint64_t jump_len = min(jump_right, jump_left);
#if HLT_//DEBUG_PRINT == 1
        //        cout << "--loop: i=" << i << "--" << endl;
        //        cout << "iterator_of_index: (" << index_iterator_left << ", "
        //            << index_iterator_right << ") " << endl;
        //        cout << "kernel indexes: ("
        //            << kernel_L_indexes[index_iterator_left] << ", "
        //            << kernel_L_indexes[index_iterator_right] << ") " << endl;
        //        cout << "windows: (" << i << ", " << i + innerp_size - 1 << ")"
        //            << endl;
        //        cout << "jump_right: " << jump_right << endl;
        //        cout << "jump_left: " << jump_left << endl;
        //        cout << "jump_len: " << jump_len << endl;
        //        cout << "partial_sum: " << partial_sum << endl;
        //        cout << "make pair: " << partial_sum << ", " << jump_len << endl;
        //        cout << endl;
#endif
        //        // make pair
        //        auto pair = make_pair(partial_sum, jump_len);
        //        diagonalpairlist.push_back(pair);
        //        pair_num++;

        //        // update partial sum
        //        if (jump_len == jump_right) {
        //            uint64_t list_R_coeff = list_R[kernel_L_indexes[index_iterator_right] + offset];
        //            auto coeff_prod = util::multiply_uint_mod(
        //                    kernel_L[kernel_L_indexes[index_iterator_right]],
        //                    list_R_coeff, modulus);
        //            partial_sum =
        //                util::add_uint_mod(partial_sum, coeff_prod, modulus);
        //            if (index_iterator_right == kernel_L_indexes.size() - 1) {
#if HLT_//DEBUG_PRINT == 1
        //                cout << "index_right is on edge!" << endl;
#endif
        //                is_right_edge = true;
        //            } else {
        //                index_iterator_right++;
        //            }
        //        }
        //        if (jump_len == jump_left) {
        //            uint64_t list_R_coeff = list_R[kernel_L_indexes[index_iterator_left] + offset];
        //            auto coeff_prod = util::multiply_uint_mod(
        //                    kernel_L[kernel_L_indexes[index_iterator_left]],
        //                    list_R_coeff, modulus);
        //            partial_sum =
        //                util::sub_uint_mod(partial_sum, coeff_prod, modulus);
        //            if (index_iterator_left == kernel_L_indexes.size() - 1) {
#if HLT_//DEBUG_PRINT == 1
        //                cout << "index_left is on edge!" << endl;
#endif
        //                is_left_edge = true;
        //            } else {
        //                index_iterator_left++;
        //            }
        //        }
        //        i = i + jump_len - 1;
        //    }
#if HLT_//DEBUG_TIME == 1
        //    auto slide_end   = chrono::high_resolution_clock::now();
        //    auto begin_diff  = chrono::duration_cast<chrono::nanoseconds>(mul_start - diagonal_begin);
        //    auto mul_diff    = chrono::duration_cast<chrono::nanoseconds>(mul_end - mul_start);
        //    auto innerp_diff = chrono::duration_cast<chrono::nanoseconds>(innerp_end - mul_end);
        //    auto slide_diff  = chrono::duration_cast<chrono::nanoseconds>(slide_end - innerp_end);
        //    cout <<  "begin: " << begin_diff.count() << " mul : " << mul_diff.count() << " innerp: " << innerp_diff.count() << " slide: " << slide_diff.count() << " sum: " << begin_diff.count() + mul_diff.count() + innerp_diff.count() + slide_diff.count() << endl;
#endif
        //}

        // convert diagonal value vectors to matrix
        // diagonallist: diagonal element vectors(first element is lower left)
        // assume: all rowsize of result is same
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
                    if(pair_value == 0 && pair_valuelen == 0) break;
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
                    if(pair_value == 0 && pair_valuelen == 0) break;
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

        // convert scalar list to diagonal vectors
        vector<vector<uint64_t>> scalars_to_diagonallist(vector<uint64_t> scalars, uint64_t colsize, uint64_t rowsize){
            vector<vector<uint64_t>> diagonals(colsize + rowsize - 1);
            uint64_t maxlength = colsize >= rowsize? rowsize: colsize;
            uint64_t keeplength = colsize >= rowsize? colsize-rowsize: rowsize-colsize;
            // lower left to uppwer left
            uint64_t i = 0;
            for(;i < colsize;i++){
                for(uint64_t j = 0;j < i+1;j++){
                    diagonals[i].push_back(scalars[i]);
                }
            }
            // move right
            for(;i < colsize + keeplength;i++){
                for(uint64_t j = 0;j < maxlength;j++){
                    diagonals[i].push_back(scalars[i]);
                }
            }
            // diagonal size decreases
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

        // kernel: vectors of filter weights.
        // block_size: size of 1d convolution result
        vector<KernelInfo> pack_kernel(vector<vector<uint64_t>> kernels, uint64_t input_size, Modulus modulus){
            uint64_t packing_num = kernels.size();
            vector<KernelInfo> kernel_info(packing_num);
            uint64_t start_col = 0;
            uint64_t start_row = 0;
            for(uint64_t i = 0;i < packing_num;i++){
                uint64_t block_size = get_blocksize(input_size, kernels[i].size(), 0 );
                KernelInfo kinfo(input_size, block_size, start_col, start_row, block_size, block_size, kernels[i], modulus);
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
                vector<uint64_t> k_scalar = kinfo.diagonal_scalars;
                reverse(k_scalar.begin(), k_scalar.end());
                vector<vector<uint64_t>> diagonals = util::scalars_to_diagonallist(k_scalar, kinfo.get_colsize(), kinfo.get_rowsize());
                util::diagonallist_to_matrix(diagonals, kinfo.get_startcol(), kinfo.get_startrow(), kinfo.get_colsize(), kinfo.get_rowsize(),matrix);
            }
        }

        // matrix dot matrix product using toeplitz algorithm
        void matrix_dot_matrix_toeplitz_mod(vector<KernelInfo> &kernel_infos, CoeffIter c1, uint64_t poly_degree, vector<vector<uint64_t>> &result,Modulus &modulus){
            // for each block
            for(uint64_t i = 0;i < kernel_infos.size();i++){
#if HLT_DEBUG_TIME == 1
                auto first_settings_begin = chrono::high_resolution_clock::now();
#endif
                // get diagonal lists for kernel
                KernelInfo kinfo = kernel_infos[i];
                vector<uint64_t> kernel_index = kinfo.index;
                uint64_t colsize_K = kinfo.get_colsize();

                // diagonal list of c1
                uint64_t submat_startcol,submat_startrow, submat_colsize, submat_rowsize;
                submat_rowsize = poly_degree;
                submat_startrow = 0;
                kinfo.getParamsforSubmatrix(submat_startcol, submat_colsize);
                vector<uint64_t> diagonal_c1 = util::create_diagonal_from_submatrix(c1, poly_degree , submat_startcol, submat_colsize, modulus);
                //vector<uint64_t> diagonal_scalars = kernel_infos[i].diagonal_scalars;

                // calc diagonal of product
                // pair(value, value_len)
                uint64_t diagonal_vector_size = colsize_K + submat_rowsize - 1;
                uint64_t diagonal_element_size = kernel_infos[i].kernel_size;
                vector<vector<pair<uint64_t, uint64_t>>> matrix_product_diagonals(diagonal_vector_size);
                //cout << "diagonal element_size: " << diagonal_element_size << endl;
                for(int i = 0;i < diagonal_vector_size;i++){
                    matrix_product_diagonals[i].reserve(diagonal_element_size);
                }
                //vector<vector<uint64_t>> matrix_product_diagonals(colsize_K + submat_rowsize - 1);
                uint64_t index = 0;
                int64_t k = static_cast<int64_t>(colsize_K);
                k = -k+1;
#if HLT_DEBUG_TIME == 1
                auto diagonal_all_start = chrono::high_resolution_clock::now();
#endif
#if HLT_DEBUG_TIME == DEBUG_TIME_SUM_DIAGONAL
                vector<uint64_t> times_diagonal_calc(submat_rowsize + colsize_K+1);
                vector<uint64_t> times_prepare_vector(submat_rowsize + colsize_K+1);
                vector<uint64_t> times_save_vector(submat_rowsize + colsize_K+1);
                int iter_diagonal = 0;
#endif
                for(;k<static_cast<int64_t>(submat_rowsize);k++){
#if HLT_DEBUG_TIME == DEBUG_TIME_SUM_DIAGONAL
                    //cout << "alloc " << kernel_index.size() << "pairs" << endl;
                    auto prepare_vector = chrono::high_resolution_clock::now();
#endif
                    //matrix_product_diagonals[index].reserve(diagonal_element_size);
#if HLT_DEBUG_TIME == DEBUG_TIME_SUM_DIAGONAL
                    auto diagonal_start = chrono::high_resolution_clock::now();
#endif
                    //vector<uint64_t> diagonal_pairs;
                    util::matrix_product_diagonal(k, submat_colsize, submat_rowsize, kernel_infos[i].diagonal_scalars, kernel_index, diagonal_c1, modulus, matrix_product_diagonals, index);
                    // digaonal_pairs = util::matrix_product_diagonal(k, submat_colsize, submat_rowsize, kernel_diagonal_list, kernel_index, diagonal_c1, modulus, diagonal_pairs);
                    //cout << "kernel_index_size: " << kernel_index.size() << endl;
                    //cout << "diagonal_pairs,size = " << diagonal_pairs.size() << endl;
#if HLT_DEBUG_TIME == DEBUG_TIME_SUM_DIAGONAL
                    auto save_start = chrono::high_resolution_clock::now();
#endif
                    //matrix_product_diagonals[index] = diagonal_pairs;
                    index++;
#if HLT_DEBUG_TIME == DEBUG_TIME_SUM_DIAGONAL
                    auto diagonal_end = chrono::high_resolution_clock::now();

                    auto diagonal_diff = chrono::duration_cast<chrono::nanoseconds>(diagonal_end - diagonal_start);
                    auto alloc_diff = chrono::duration_cast<chrono::nanoseconds>(diagonal_start- prepare_vector);
                    auto save_diff = chrono::duration_cast<chrono::nanoseconds>(diagonal_end- save_start);
                    times_diagonal_calc[iter_diagonal] = diagonal_diff.count();
                    times_prepare_vector[iter_diagonal] = alloc_diff.count();
                    times_save_vector[iter_diagonal] = save_diff.count();
                    iter_diagonal++;
                    //cout << "calc one diagonal vector: " << diagonal_diff.count() << "ns"  << endl;
#endif
                }
#if HLT_DEBUG_PRINT == 1
                print_pair_vectors(matrix_product_diagonals);
#endif

#if HLT_DEBUG_TIME == 1
                auto get_toeplitz_start = chrono::high_resolution_clock::now();
#endif
                kernel_infos[i].get_toeplitz(matrix_product_diagonals, poly_degree);
#if HLT_DEBUG_TIME == 1
                // write diagonals to result matrix
                auto write_matrix_start = chrono::high_resolution_clock::now();
#endif
                util::diagonallist_to_matrix(matrix_product_diagonals, submat_startcol, submat_startrow, colsize_K, submat_rowsize, result);
#if HLT_DEBUG_PRINT
                util::print_matrix(result, 0, 0, 20, 20);
#endif
#if HLT_DEBUG_TIME == 1
                auto first_settings_diff = chrono::duration_cast<chrono::microseconds>(diagonal_all_start - first_settings_begin);
                auto write_matrix_end = chrono::high_resolution_clock::now();
                auto write_diff = chrono::duration_cast<chrono::microseconds>(write_matrix_end - write_matrix_start);
                auto diagonal_diff = chrono::duration_cast<chrono::microseconds>(get_toeplitz_start - diagonal_all_start);
                auto toeplitz_diff = chrono::duration_cast<chrono::microseconds>(write_matrix_start - get_toeplitz_start);
                cout << "first settings: " << first_settings_diff.count() << "us" << endl;
                cout << "all diagonal: " << diagonal_diff.count() << "us" << endl;
                cout << "get toeplitz: " << toeplitz_diff.count() << "us" << endl;
                cout << "write matrix: " << write_diff.count()    << "us" << endl;
#endif
#if HLT_DEBUG_TIME == DEBUG_TIME_SUM_DIAGONAL
                uint64_t sum_diagonal_time = 0;
                uint64_t sum_prepare_time = 0;
                uint64_t sum_save_time = 0;
                for(int i = 0;i < times_diagonal_calc.size();i++){
                    //cout << "diagonal time[" << i << "]: " << times_diagonal_calc[i] << "ns" << endl;
                    sum_diagonal_time += times_diagonal_calc[i];
                }
                for(int i = 0;i < times_prepare_vector.size();i++){
                    //cout << "diagonal time[" << i << "]: " << times_diagonal_calc[i] << "ns" << endl;
                    sum_prepare_time += times_prepare_vector[i];
                }
                for(int i = 0;i < times_save_vector.size();i++){
                    sum_save_time    += times_save_vector[i];
                }
                cout << "sum diagonal time: " << sum_diagonal_time/1000 << "us" << endl;
                cout << "sum prepare time: "  << sum_prepare_time/1000 << "us" << endl;
                cout << "sum save time: "  << sum_save_time/1000 << "us" << endl;
#endif
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

        void packedconv_matrix_dot_vector(vector<vector<uint64_t>> &matrix_multed, vector<KernelInfo> kernelinfos, CoeffIter vector_iter, uint64_t poly_degree, CoeffIter destination, const Modulus &modulus, MemoryPoolHandle pool_){
            uint64_t offset_begin = 0;
            for(uint64_t i = 0;i < kernelinfos.size();i++){
                uint64_t kernel_range = kernelinfos[i].kernel_size-1;
                uint64_t block_size = kernelinfos[i].block_size;
                matrix_dot_vector(matrix_multed, vector_iter, poly_degree, offset_begin, kernel_range, destination, modulus);
                CoeffIter dest_toeplitz_begin = destination + offset_begin + kernel_range;
                toeplitz_dot_vector(kernelinfos[i].toeplitz_diagonal_scalars, vector_iter, kernelinfos[i].input_size, poly_degree, modulus, dest_toeplitz_begin, pool_);
                offset_begin += block_size;
            }
        }

        }
    }
