#include "testcommon.h"

TEST(Packing, kernel_dot_c1){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
    vector<uint64_t> c1 = {1, 4, 2, 5, 1, 3, 6, 2, 10, 3, 11, 9, 2, 7, 2, 5, 7, 1, 2, 4, 8, 1, 2};
    uint64_t poly_degree = c1.size();
    vector<vector<uint64_t>> kernels = { {3, 1, 2}, {3, 1, 2}};
    uint64_t input_dim = 7;
    uint64_t block_size = get_blocksize(input_dim, kernels[0].size(), 0);
    Modulus modulus(13);
    vector<KernelInfo> kernel_info = pack_kernel(kernels, block_size, modulus);
    vector<vector<uint64_t>> result(poly_degree, vector<uint64_t>(poly_degree));
    matrix_dot_matrix_toeplitz_mod(kernel_info, c1.data(), poly_degree, result, modulus);
    //util::print_matrix(result);
    //cout << "negacyclic: " << endl;
    //SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(nega_result, poly_degree, pool_);
    //conv_negacyclic(kernels[0], c1.data(), poly_degree, modulus, nega_result);
    //util::print_iter(nega_result, poly_degree);
}

