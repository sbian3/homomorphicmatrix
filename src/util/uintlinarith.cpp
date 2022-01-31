#include "uintlinarith.h"
#include <chrono>
#include <iomanip>


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

        void init_matrix_negate(vector<vector<uint64_t>>& matrix, uint64_t size, uint64_t left_rotate, uint64_t scale,const Modulus &modulus){
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

        void init_matrix_circ(vector<vector<uint64_t>>& matrix, uint64_t size, ConstCoeffIter iter, const Modulus &modulus){
            for(uint64_t i = 0;i < size;i++){
                init_matrix_negate(matrix, size, i, iter[i], modulus);
            }
        }

        void init_matrix_circ(vector<vector<uint64_t>>& matrix, uint64_t size_matrix, ConstCoeffIter iter, uint64_t size_kernel, const Modulus &modulus){
            for(uint64_t i = 0;i < size_kernel;i++){
                init_matrix_negate(matrix, size_matrix, i, iter[i], modulus);
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

        vector<vector<uint64_t>> toeplitz_to_matrix(vector<uint64_t> toeplitz, uint64_t toeplitz_rowsize, uint64_t toeplitz_colsize){
            vector<vector<uint64_t>> matrix(toeplitz_rowsize, vector<uint64_t>(toeplitz_colsize));
            assert(toeplitz.size() == toeplitz_rowsize + toeplitz_colsize - 1);
            uint64_t max_diagonal_len = toeplitz_rowsize;
            uint64_t maxdiagonal_count;
            if(toeplitz_rowsize > toeplitz_colsize){
                maxdiagonal_count = toeplitz_rowsize - toeplitz_colsize + 1;
            }else{
                maxdiagonal_count = toeplitz_colsize - toeplitz_rowsize + 1;
            }
            // start from lower left
            uint64_t counter = 1;
            for(;counter < max_diagonal_len;counter++){
                for(uint64_t i = 0;i < counter;i++){
                    matrix[toeplitz_rowsize+i-counter][i] = toeplitz[counter-1];
                }
            }
            // diagonal len is max
            for(uint64_t i = 0;i < maxdiagonal_count;i++){
                for(uint64_t j = 0;j < max_diagonal_len;j++){
                    matrix[j][i+j] = toeplitz[max_diagonal_len + i- 1];
                }
            }
            // end to upper right
            counter--;
            for(uint64_t k = 0;counter > 0;k++,counter--){
                for(uint64_t i = 0;i < counter;i++){
                    matrix[i][maxdiagonal_count+i+k] = toeplitz[max_diagonal_len + maxdiagonal_count + k - 1];
                }
            }
            return matrix;
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
        // matrix and vector arithmetic
        //
        ///////////////////////////////
        void toeplitz_to_circ(vector<uint64_t> &toeplitz, uint64_t toeplitz_rowsize, uint64_t toeplitz_colsize, CoeffIter circ, Modulus modulus){
            assert(toeplitz.size() == toeplitz_rowsize + toeplitz_colsize - 1);
            //assert(circ.size() == toeplitz_colsize * 2);
            uint64_t circ_size = toeplitz_colsize * 2;
            uint64_t iter_toep = toeplitz_rowsize - 1;
            for(uint64_t i = 0;i < toeplitz_rowsize;i++){
                circ[i] = toeplitz[iter_toep];
                iter_toep--;
            }
            iter_toep = toeplitz_rowsize;
            for(uint64_t i = circ_size - 1;i > toeplitz_colsize;i--){
                circ[i] = negate_uint_mod(toeplitz[iter_toep], modulus);
                iter_toep++;
            }
        }

        void toeplitz_dot_vector(vector<uint64_t> &toeplitz, CoeffIter right_vec_coeff, uint64_t toeplitz_rowsize, uint64_t toeplitz_colsize, const Modulus &modulus, CoeffIter destination, MemoryPoolHandle pool_){
            uint64_t right_vec_coeff_size = toeplitz_colsize;
            uint64_t circ_size = get_bigger_poweroftwo(toeplitz_colsize) * 2;
            uint64_t coeff_count_power = get_power_of_two(circ_size);
#if HLT_DEBUG_PRINT == DEBUG_PRINT_DEC
#endif
#if HLT_DEBUG_TIME == DEBUG_TIME_DEC_NTT
            auto prepare_begin = chrono::high_resolution_clock::now();
#endif
            Pointer<NTTTables> ntt_tables = allocate<NTTTables>(pool_, coeff_count_power, modulus, pool_);
#if HLT_DEBUG_TIME == DEBUG_TIME_DEC_NTT
            auto prepare_end = chrono::high_resolution_clock::now();
#endif
            SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(right_vec, circ_size, pool_);
            util::set_poly(right_vec_coeff, right_vec_coeff_size, 1, right_vec);
            SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(circ, circ_size, pool_);
            SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(dest_tmp, circ_size, pool_);
            toeplitz_to_circ(toeplitz, toeplitz_rowsize, toeplitz_colsize, circ, modulus);
            //assert(circ.size() == right_vec.size());

#if HLT_DEBUG_TIME == DEBUG_TIME_DEC_NTT
            auto ntt_begin = chrono::high_resolution_clock::now();
#endif
            // multiply circ and right_vec
            ntt_negacyclic_harvey(circ, *ntt_tables);
            ntt_negacyclic_harvey(right_vec, *ntt_tables);
            dyadic_product_coeffmod(circ, right_vec, circ_size, modulus, dest_tmp);
            inverse_ntt_negacyclic_harvey(dest_tmp, *ntt_tables);
            util::set_poly(dest_tmp, toeplitz_rowsize, 1, destination);
#if HLT_DEBUG_TIME == DEBUG_TIME_DEC_NTT
            auto ntt_end = chrono::high_resolution_clock::now();
            auto ntt_time = chrono::duration_cast<chrono::microseconds>(ntt_end - ntt_begin);
            auto prepare_time = chrono::duration_cast<chrono::microseconds>(prepare_end - prepare_begin);
            cout << "prepare time: " << prepare_time.count() << "us" << endl;
            cout << "ntt time: " << ntt_time.count() << "us" << endl;
#endif
        }

        ///////////////////////////////
        // 
        // matrix and matrix arithmetic
        //
        ///////////////////////////////

        void matrix_dot_matrix_mod(vector<vector<uint64_t>> &matrixL, vector<vector<uint64_t>>& matrixR, vector<vector<uint64_t>>& result,const Modulus &modulus){
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

        void matrix_dot_matrix_mod_t(vector<vector<uint64_t>> &matrixL, vector<vector<uint64_t>> &matrixtR, vector<vector<uint64_t>>& result, Modulus &modulus){
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

        ///////////////////////
        //
        // Convolution
        //
        //////////////////////


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
                    cout << " " << row.at(j);
                }
                cout << endl;
            } 
        }

    } // namespace util
} // namespace seal

