
#include "seal/util/linarith.h"
#include <chrono>


namespace seal
{
    namespace util
    {

        // 
        // matrix initialization
        //

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


        void matrix_dot_product_mod(vector<vector<int64_t>> matrixL, vector<vector<int64_t>> matrixR, vector<vector<int64_t>>& result, uint64_t mod){
            auto time_start = chrono::high_resolution_clock::now();
            assert(matrixL[0].size() == matrixR.size());
            for(auto i = 0U;i < matrixL.size();i++){
                for(auto j = 0U;j < matrixR[0].size();j++){
                    int64_t tmp_sum = 0;
                    if(matrixL[i].size() != matrixR.size()){
                        cout << "Error: Left and Right matrix size does not match!!" << endl;
                        return;
                    }
                    for(auto k = 0U;k < matrixR.size();k++){
                        tmp_sum += matrixL[i][k] * matrixR[k][j]; 
                        tmp_sum = tmp_sum % static_cast<int64_t>(mod);
                    }
                    result[i][j] = tmp_sum;
                }
            }
            auto time_end = chrono::high_resolution_clock::now();
            auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
            cout << "matrix dot product: " << time_diff.count() << "milliseconds" << endl;
        }


        void matrix_dot_vector(vector<vector<int64_t>> matrix, ConstCoeffIter poly_vector, const Modulus& modulus, uint64_t coeff_count, CoeffIter result){
            // TODO: parameter validation

            for(uint64_t i = 0;i < coeff_count;i++){
                result[i] = inner_product_coeffmod(matrix[i], poly_vector, coeff_count, modulus);
            }
        }

        //
        // functions for decryptor
        //

        void secret_product_with_matrix(vector<vector<int64_t>> matrix,uint64_t coeff_degree, CoeffIter c, CoeffIter s, const Modulus& modulus, CoeffIter result){

            uint64_t mod_value = modulus.value();
            vector<vector<int64_t>> A(coeff_degree, vector<int64_t>(coeff_degree));
            vector<vector<int64_t>> B(coeff_degree, vector<int64_t>(coeff_degree));
            init_matrix_with_coeff(A, coeff_degree, c);
            matrix_dot_product_mod(matrix, A, B, mod_value);
            matrix_dot_vector(B, s, modulus, coeff_degree,result );
        }


        //
        // print functions
        //

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

