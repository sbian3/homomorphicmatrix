
#include "seal/util/linarith.h"


namespace seal
{
    namespace util
    {

        void init_matrix_identity(vector<vector<int64_t>>& matrix, uint64_t poly_modulus_degree, int64_t scale){
            for(auto i = 0U;i < poly_modulus_degree;i++){
                for(auto j = 0U;j < poly_modulus_degree;j++){
                    if(i == j)
                        matrix[i][j] = scale;
                    else
                        matrix[i][j] = 0;
                }
            }
        }

        void init_matrix_rotate(vector<vector<int64_t>>& matrix, uint64_t size, uint64_t left_rotate, int64_t scale){
            for(auto i = 0U;i < size;i++){
                int64_t tmp_scale = scale;
                uint64_t col_index = i + left_rotate;
                if(col_index >= size){
                    col_index = col_index % size;
                    tmp_scale = -tmp_scale;
                }
                matrix[col_index][i] = tmp_scale;
            }
        }

        void init_matrix_with_coeff(vector<vector<int64_t>>& matrix, uint64_t size, ConstCoeffIter iter){
            for(uint64_t i = 0;i < size;i++){
                init_matrix_rotate(matrix, size, i, (int64_t)iter[i]);
            }
        }

        void init_matrix_rand_mod(vector<vector<int64_t>>& matrix, uint64_t size, uint64_t mod){
            random_device rnd;
            for(auto i = 0U;i < size;i++){
                for(auto j = 0U;j < size;j++){
                    matrix[i][j] = rnd() % mod;
                }
            }
        }

        std::uint64_t inner_product_coeffmod(vector<int64_t> operand1, ConstCoeffIter operand2, std::size_t coeff_count, const Modulus &modulus){
            const uint64_t modulus_value = modulus.value();
            int64_t result = 0;
            for(size_t i = 0;i < coeff_count;i++){
                result += operand1[i] * static_cast<int64_t>(operand2[i]); 
                result = result % static_cast<int64_t>(modulus_value);
            }
            if(result < 0){
                result += modulus_value;
            }
            return static_cast<uint64_t>(result);
        }

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
        void matrix_dot_product_mod(vector<vector<int64_t>> matrixL, vector<vector<int64_t>> matrixR, vector<vector<int64_t>>& result, uint64_t mod){
            assert(matrixL[0].size() == matrixR.size());
            for(auto i = 0U;i < matrixL.size();i++){
                for(auto j = 0U;j < matrixR[0].size();j++){
                    int64_t tmp_sum = 0;
                    for(auto k = 0U;k < matrixR.size();k++){
                        tmp_sum += matrixL[i][k] * matrixR[k][j]; 
                        tmp_sum = tmp_sum % static_cast<int64_t>(mod);
                    }
                    result[i][j] = tmp_sum;
                }
            }
        }


        void matrix_dot_vector(ConstRNSIter matrix, ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, CoeffIter result){
            // TODO: parameter validation

            SEAL_ITERATE(iter(matrix, result), coeff_count, [&](auto I){
                    get<1>(I) = inner_product_coeffmod(get<0>(I), poly_vector, coeff_count, modulus);
                    });
        }

        void matrix_dot_vector(vector<vector<int64_t>> matrix, ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, CoeffIter result){
            // TODO: parameter validation

            for(uint64_t i = 0;i < coeff_count;i++){
                result[i] = inner_product_coeffmod(matrix[i], poly_vector, coeff_count, modulus);
            }
        }

        void matrix_dot_vector(ConstRNSIter matrix, uint64_t coeff_modulus_size, ConstRNSIter poly_rns, ConstModulusIter mod_chain, RNSIter result){
            // parameter validation
            // TODO: size check
            uint64_t coeff_count = poly_rns.poly_modulus_degree();
            SEAL_ITERATE(iter(poly_rns, mod_chain, result), coeff_modulus_size, [&](auto I){
                    matrix_dot_vector(matrix, get<0>(I), get<1>(I), coeff_count, get<2>(I));
                    });
        }

        void matrix_dot_vector(vector<vector<int64_t>> matrix, uint64_t coeff_modulus_size, ConstRNSIter poly_rns, ConstModulusIter mod_chain, RNSIter result){
            // parameter validation
            // TODO: size check
            uint64_t coeff_count = poly_rns.poly_modulus_degree();
            SEAL_ITERATE(iter(poly_rns, mod_chain, result), coeff_modulus_size, [&](auto I){
                    matrix_dot_vector(matrix, get<0>(I), get<1>(I), coeff_count, get<2>(I));
                    });
        }

        void secret_product_with_matrix(vector<vector<int64_t>> matrix,uint64_t coeff_degree, CoeffIter c, CoeffIter s, const Modulus& modulus, CoeffIter result){

            uint64_t mod_value = modulus.value();
            vector<vector<int64_t>> A(coeff_degree, vector<int64_t>(coeff_degree));
            vector<vector<int64_t>> B(coeff_degree, vector<int64_t>(coeff_degree));
            init_matrix_with_coeff(A, coeff_degree, c);
            matrix_dot_product_mod(matrix, A, B, mod_value);
            matrix_dot_vector(B, s, modulus, coeff_degree,result );
        }

        void secret_product_with_matrix_rns(vector<vector<int64_t>> matrix, uint64_t rns_count, RNSIter c, RNSIter s, ConstModulusIter mod_chain, RNSIter result){
            uint64_t coeff_degree = c.poly_modulus_degree();
            SEAL_ITERATE(iter(c, s, mod_chain, result), rns_count, [&](auto I){
                    secret_product_with_matrix(matrix, coeff_degree, get<0>(I), get<1>(I), get<2>(I), get<3>(I));
                    });
        }

        //
        // print functions
        //

        void print_iter(CoeffIter operand1, uint64_t coeff_count){
            SEAL_ITERATE(operand1, coeff_count, [&](auto I){
                    cout << I << " ";
                    });
            cout << endl;
        }


        void print_matrix(vector<vector<int64_t>>& matrix){
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

