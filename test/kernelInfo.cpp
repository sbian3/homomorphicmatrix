#include "testcommon.h"

TEST(KernelInfotest, init){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6} };
    uint64_t pack_num = kernels.size();
    uint64_t input_size = 32 - kernels[0].size()+ 1;
    Modulus modulus(7);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, input_size, modulus);
    for(uint64_t i = 0; i < pack_num;i++){
        ASSERT_EQ(kernels[i].size(), kernel_infos[i].kernel_size);
        ASSERT_EQ(input_size, kernel_infos[i].input_size);
        ASSERT_EQ(32, kernel_infos[i].block_size);
    }
}

TEST(KernelInfotest, get_toeplitz){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6} };
    uint64_t input_size = 3;
    Modulus modulus(7);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, input_size, modulus);
    uint64_t poly_degree = 3;
    vector<vector<pair<uint64_t, uint64_t>>> pairs = 
    { { make_pair(1, 1), make_pair(2, 2) },
        { make_pair(3, 3), make_pair(4, 1) },
        {make_pair(2,2), make_pair(5, 2)},
        {make_pair(6, 2), make_pair(6,8)},
        {make_pair(4, 6), make_pair(2,2)},
        {make_pair(2,1)},
        {make_pair(4, 9)},
        {make_pair(4,4)},
        {make_pair(2,3)},
        {make_pair(4, 9)}
    };
    kernel_infos[0].get_toeplitz(pairs, poly_degree);
    vector<uint64_t> expect = {2, 4, 5, 6, 2};
    for(uint64_t i = 0;i < expect.size();i++){
        ASSERT_EQ(expect[i], kernel_infos[0].toeplitz_diagonal_scalars[i]);
    }
}
