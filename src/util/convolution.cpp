#include "convolution.h"
#include "define_tifs.h"


namespace seal
{
    namespace util
    {

        ///////////////////
        // 1D Convolution
        ///////////////////

        // direct convolution (kernel * coeffiter)
        void conv_negacyclic(vector<uint64_t> &kernel, CoeffIter encrypted, uint64_t poly_degree, const Modulus &modulus, CoeffIter result){
            uint64_t size_poly = poly_degree;
            for(uint64_t j = 0U;j < size_poly;j++){
                for(uint64_t i = 0U;i < kernel.size();i++){
                    uint64_t index = i+j;
                    uint64_t tmp_prod = kernel[i];
                    if((index) >= size_poly){
                        tmp_prod = negate_uint_mod(tmp_prod, modulus);
                        index -= size_poly;
                    }
                    result[index] = multiply_add_uint_mod(tmp_prod, encrypted[j], result[index], modulus);
                }
            }
        }

        ////////////////////
        // 2D Convolution
        ////////////////////
        vector<uint64_t> transform_to_helix(vector<vector<uint64_t>> &input, uint64_t colsize_helix){
            uint64_t inputmat_colsize = input.size();
            // assert every input[i] is a same length vector
            uint64_t inputmat_rowsize = input[0].size();
            // concatenate last zeroes
            vector<uint64_t> helix(colsize_helix * (inputmat_rowsize-1) + inputmat_colsize);

            for(uint64_t i = 0;i < inputmat_rowsize;i++){
                for(uint64_t j = 0;j < inputmat_colsize;j++){
                    helix[colsize_helix * i + j] = input[j][i];
                }
            }
            return helix;
        }

    }
}
