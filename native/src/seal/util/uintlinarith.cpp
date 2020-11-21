
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

        void init_matrix_with_coeff(vector<vector<uint64_t>>& matrix, uint64_t size, ConstCoeffIter iter, const Modulus &modulus){
            for(uint64_t i = 0;i < size;i++){
                init_matrix_rotate(matrix, size, i, iter[i], modulus);
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

        //
        // linear arithmetic
        //

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

        // 
        // matrix and matrix arithmetic
        //

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

        // 
        // print function
        //

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

