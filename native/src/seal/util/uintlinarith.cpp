
#include "seal/util/uintlinarith.h"
#include <chrono>


namespace seal
{
    namespace util
    {

        //
        // matrix initialization
        //

        void init_matrix_identity(vector<vector<uint64_t>>& matrix, uint64_t poly_modulus_degree, uint64_t scale){
            for(auto i = 0U;i < poly_modulus_degree;i++){
                for(auto j = 0U;j < poly_modulus_degree;j++){
                    if(i == j)
                        matrix[i][j] = scale;
                    else
                        matrix[i][j] = 0;
                }
            }
        }

        void init_matrix_rotate(vector<vector<uint64_t>>& matrix, uint64_t size, uint64_t left_rotate, uint64_t scale,const Modulus &modulus){
            for(auto i = 0U;i < size;i++){
                uint64_t tmp_scale = scale;
                uint64_t col_index = i + left_rotate;
                if(col_index >= size){
                    col_index = col_index % size;
                    tmp_scale = negate_uint_mod(tmp_scale, modulus);
                }
                matrix[col_index][i] = tmp_scale;
            }
        }

        void init_matrix_rotate_partial(vector<vector<uint64_t>> &matrix, uint64_t size_kernel, uint64_t left_rotate, uint64_t start_col, uint64_t start_row, uint64_t scale, const Modulus &modulus){
            // assert
            // 1.start_row + size_kernel <= matrix.row.size
            // 2.start_col + size_kernel <= matrix.col.size
            for(auto i = 0U;i < size_kernel;i++){
                uint64_t tmp_scale = scale;
                uint64_t col_index = start_col + i + left_rotate;
                if(col_index >= start_col + size_kernel){
                    col_index = col_index - size_kernel;
                    tmp_scale = negate_uint_mod(tmp_scale, modulus);
                }
                matrix[col_index][start_row + i] = tmp_scale;
            }
        }

        void init_matrix_rotate_partial(vector<vector<uint64_t>> &matrix, vector<uint64_t> kernel, uint64_t start_col, uint64_t start_row, const Modulus &modulus){
            uint64_t kernel_size = kernel.size();
            for(auto i = 0U;i < kernel_size;i++){
                init_matrix_rotate_partial(matrix, kernel_size, i, start_col, start_row, kernel[i], modulus);
            }
        }

        void init_matrix_with_coeff(vector<vector<uint64_t>>& matrix, uint64_t size, ConstCoeffIter iter, const Modulus &modulus){
            for(uint64_t i = 0;i < size;i++){
                init_matrix_rotate(matrix, size, i, iter[i], modulus);
            }
        }

        void init_matrix_with_coeff(vector<vector<uint64_t>>& matrix, uint64_t size_matrix, ConstCoeffIter iter, uint64_t size_kernel, const Modulus &modulus){
            for(uint64_t i = 0;i < size_kernel;i++){
                init_matrix_rotate(matrix, size_matrix, i, iter[i], modulus);
            }
        }

        void init_matrix_rand_mod(vector<vector<uint64_t>>& matrix, uint64_t size, uint64_t mod){
            random_device rnd;
            for(auto i = 0U;i < size;i++){
                for(auto j = 0U;j < size;j++){
                    matrix[i][j] = rnd() % mod;
                }
            }
        }

        void init_matrix_diagonal(vector<vector<uint64_t>> &matrix, uint64_t size, uint64_t scalar, uint64_t right_rotate){
            // assert that matrix is big enough
            for(auto i = 0U;i < size;i++){
                matrix[i][i + right_rotate] = scalar;
            }
        }

        void copy_matrix(vector<vector<uint64_t>> &dest, vector<vector<uint64_t>> src, uint64_t start_col, uint64_t start_row){
            for(auto i = 0U;i < src.size();i++){
                for(auto j = 0U;j < src[i].size();j++){
                    dest[start_col + i][start_row + j] = src[i][j];
                }
            }
        }

        void init_matrix_2dconv(vector<vector<uint64_t>> &matrix, uint64_t input_size, vector<vector<uint64_t>> kernel){
            // assert that....
            // matrix is square
            // kernel is square
            uint64_t matrix_size = input_size;
            uint64_t kernel_size = kernel.size();
            uint64_t dest_size = matrix_size - kernel_size + 1;
            for(auto i = 0U;i < kernel_size;i++){
               vector<uint64_t> kernel_row = kernel[i]; 
               // first, generate and init block matrix
               vector<vector<uint64_t>> block(dest_size, vector<uint64_t>(matrix_size));
               for(auto j = 0U;j < kernel_size;j++){
                    init_matrix_diagonal(block, dest_size, kernel_row[j], j);
               }
               // copy it in diagonal order
               for(auto j = 0U;j < dest_size;j++){
                    copy_matrix(matrix, block, dest_size * j,matrix_size * (j + i));
               }
            }
        }

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
        
        // 対角成分ベクトルを乗算結果の行列に書き込む
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
        ///////////////////////////
        //
        // linear arithmetic
        //
        ///////////////////////////

        std::uint64_t inner_product_coeffmod(ConstCoeffIter operand1, ConstCoeffIter operand2, std::size_t coeff_count, const Modulus &modulus){
            //parameter validation
            const uint64_t modulus_value = modulus.value();
            const uint64_t const_ratio_0 = modulus.const_ratio()[0];
            const uint64_t const_ratio_1 = modulus.const_ratio()[1];

            uint64_t result = 0;
            SEAL_ITERATE(iter(operand1, operand2), coeff_count, [&](auto I){
                    // WARNING: 1 code line below cannot calculate as I thought.
                    // result += get<0>(I) * get<1>(I) % modulus_value;
                    //
                    //
                    //
                    // COPIED FROM  dyadic_product_coeffmod
                    //
                    //
                    //

                    // Reduces z using base 2^64 Barrett reduction
                    unsigned long long z[2], tmp1, tmp2[2], tmp3,tmp4, carry;
                    util::multiply_uint64(get<0>(I), get<1>(I), z);

                    // Multiply input and const_ratio
                    // Round 1
                    util::multiply_uint64_hw64(z[0], const_ratio_0, &carry);
                    util::multiply_uint64(z[0], const_ratio_1, tmp2);
                    tmp3 = tmp2[1] + util::add_uint64(tmp2[0], carry, &tmp1);

                    // Round 2
                    util::multiply_uint64(z[1], const_ratio_0, tmp2);
                    carry = tmp2[1] + util::add_uint64(tmp1, tmp2[0], &tmp1);

                    // This is all we care about
                    tmp1 = z[1] * const_ratio_1 + tmp3 + carry;

                    // Barrett subtraction
                    tmp3 = z[0] - tmp1 * modulus_value;

                    // Claim: One more subtraction is enough
                    tmp4 =
                        tmp3 - (modulus_value & static_cast<uint64_t>(-static_cast<int64_t>(tmp3 >= modulus_value)));
                    result = util::add_uint_mod(result, tmp4, modulus);
            });
            return result;
        }

        ///////////////////////////////
        // 
        // matrix and matrix arithmetic
        //
        ///////////////////////////////

        void matrix_dot_product_mod(vector<vector<uint64_t>> matrixL, vector<vector<uint64_t>> matrixR, vector<vector<uint64_t>>& result,const Modulus &modulus){
            auto time_start = chrono::high_resolution_clock::now();
            assert(matrixL[0].size() == matrixR.size());
            for(auto i = 0U;i < matrixL.size();i++){
                for(auto j = 0U;j < matrixR[0].size();j++){
                    uint64_t tmp_sum = 0;
                    if(matrixL[i].size() != matrixR.size()){
                        cout << "Error: Left and Right matrix size does not match!!" << endl;
                        return;
                    }
                    for(auto k = 0U;k < matrixR.size();k++){
                        tmp_sum = multiply_add_uint_mod(matrixL[i][k], matrixR[k][j], tmp_sum, modulus);
                    }
                    result[i][j] = tmp_sum;
                }
            }
            auto time_end = chrono::high_resolution_clock::now();
            auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
            cout << "matrix dot product: " << time_diff.count() << "milliseconds" << endl;
        }

        void matrix_dot_product_mod_t(vector<vector<uint64_t>> matrixL, vector<vector<uint64_t>> matrixtR, vector<vector<uint64_t>>& result, Modulus &modulus){
            auto time_start = chrono::high_resolution_clock::now();
            assert(matrixL[0].size() == matrixtR.size());
            for(auto i = 0U;i < matrixL.size();i++){
                for(auto j = 0U;j < matrixtR.size();j++){
                    uint64_t tmp_sum = 0;
                    for(auto k = 0U;k < matrixtR[j].size();k++){
                       tmp_sum = multiply_add_uint_mod(matrixL[i][k], matrixtR[j][k], tmp_sum, modulus);
                    }
                    result[i][j] = tmp_sum;
                }
            }
            auto time_end = chrono::high_resolution_clock::now();
            auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
            cout << "matrix dot product: " << time_diff.count() << "milliseconds" << endl;
        }
        
        void matrix_dot_matrix_toeplitz_mod(vector<vector<uint64_t>> matrixL, vector<vector<uint64_t>> matrixR, vector<vector<uint64_t>>& result,const Modulus &modulus){
        }

        ///////////////////////
        //
        // Convolution
        //
        //////////////////////

        void conv_negacyclic(vector<uint64_t> &kernel, CoeffIter encrypted, uint64_t poly_degree, const Modulus &modulus, CoeffIter result){
            uint64_t size_poly = poly_degree;
            // make sure result is 0 initialized
            for(uint64_t i = 0;i < poly_degree;i++){
                result[i] = 0;
            }
            // direct convolution (kernel * coeffiter)
            for(uint64_t i = 0U;i < kernel.size();i++){
                for(uint64_t j = 0U;j < size_poly;j++){
                    uint64_t tmp_prod = util::multiply_uint_mod(kernel[i], encrypted[j], modulus);
                    uint64_t index = i+j;
                    if((index) >= size_poly){
                        tmp_prod = util::negate_uint_mod(tmp_prod, modulus);
                        index -= size_poly;
                    }
                    result[index] = util::add_uint_mod(result[index], tmp_prod, modulus);
                }
            }
        }

        ////////////////////////////
        // 
        // print function
        //
        ///////////////////////////

        void print_iter(CoeffIter operand1, uint64_t coeff_count){
            SEAL_ITERATE(operand1, coeff_count, [&](auto I){
                    cout << I << " ";
                    });
            cout << endl;
        }

        void print_matrix(vector<vector<uint64_t>>& matrix){
            for(auto i = 0U;i < matrix.size();i++){
                auto row = matrix[i];
                for( auto j = 0U;j < row.size();j++ ){
                    cout << row.at(j) << " ";
                }
                cout << endl;
            } 
        }

    } // namespace util
} // namespace seal

