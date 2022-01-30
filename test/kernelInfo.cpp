#include "testcommon.h"

TEST(KernelInfotest, init){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6}, {1, 2, 3, 7} };
    vector<vector<uint64_t>> inputs  = { {4, 3, 2,}, {1, 2, 4}, {5, 3, 6, 7} };
    uint64_t pack_num = kernels.size();
    Modulus modulus(11);
    uint64_t poly_degree = 2048;
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, inputs, modulus, poly_degree);
    uint64_t start_col = 0;
    uint64_t start_row = 0;
    for(uint64_t i = 0; i < pack_num;i++){
        ASSERT_EQ(kernels[i].size(), kernel_infos[i].kernel_size);
        ASSERT_EQ(inputs[i].size(), kernel_infos[i].input_size);
        ASSERT_EQ(kernels[i].size()+inputs[i].size()-1, kernel_infos[i].block_size);
        ASSERT_EQ(start_col, kernel_infos[i].get_startcol());
        ASSERT_EQ(start_row, kernel_infos[i].get_startrow());
        start_col += kernels[i].size()+inputs[i].size()-1;
        start_row += kernels[i].size()+inputs[i].size()-1;
    }
}

TEST(KernelInfotest ,pack_input){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6}};
    vector<vector<uint64_t>> inputs  = { {4, 3, 2,}, {1, 2, 4}};
    Modulus modulus(11);
    uint64_t poly_degree = 16;
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, inputs, modulus, poly_degree);
    auto actual = pack_input(inputs, kernel_infos, poly_degree);
    vector<uint64_t> expected = {4, 3, 2, 0, 0, 1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0};
    ASSERT_ARR(actual.data(), expected.data(), poly_degree);
}

TEST(KernelInfotest, get_toeplitz){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6} };
    vector<vector<uint64_t>> inputs  = { {5, 1, 3}, {1, 6, 2} };
    Modulus modulus(7);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, inputs, modulus, 1024);
    uint64_t poly_degree = 3;
    vector<vector<pair<uint64_t, uint64_t>>> pairs = 
       {{make_pair(1, 1), make_pair(2, 2), make_pair(0, 0)},
        {make_pair(3, 3), make_pair(4, 1)},
        {make_pair(2, 2), make_pair(5, 2), make_pair(0, 0)},
        {make_pair(6, 2), make_pair(6, 8)},
        {make_pair(4, 6), make_pair(2, 2)},
        {make_pair(2, 1)},
        {make_pair(4, 9)},
        {make_pair(4, 4)},
        {make_pair(2, 3)},
        {make_pair(4, 9)}
    };
    kernel_infos[0].get_toeplitz(pairs, poly_degree);
    vector<uint64_t> expect = {2, 4, 5, 6, 2};
    for(uint64_t i = 0;i < expect.size();i++){
        ASSERT_EQ(expect[i], kernel_infos[0].toeplitz_diagonal_scalars[i]);
    }
}
