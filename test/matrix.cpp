#include "testcommon.h"

TEST(Arithtest, inner_product){
    vector<uint64_t> left_vec = { 1, 2, 3, 4 };
    vector<uint64_t> right_vec = { 2, 3, 4, 5 };
    auto modulus = Modulus(0x7e00001ULL);
    uint64_t actual = inner_product_coeffmod(left_vec.data(), right_vec.data(), left_vec.size(), modulus);
    ASSERT_EQ(40, actual);
    right_vec = { 2, 3, 4, negate_uint_mod(6, modulus) };
    actual = inner_product_coeffmod(left_vec.data(), right_vec.data(), left_vec.size(), modulus);
    ASSERT_EQ(negate_uint_mod(4, modulus), actual);
}

TEST(Matrixtest, init_circulant_matrix){
    uint64_t poly_degree = 10;
    Modulus modulus(23);

    vector<vector<uint64_t>> actual(poly_degree, vector<uint64_t>(poly_degree));
    vector<uint64_t> coeff(poly_degree);
    for(uint64_t i = 0;i < coeff.size();i++){
        coeff[i] = i+1;
    }

    init_matrix_circ(actual, poly_degree, coeff.data(), modulus);
    vector<vector<uint64_t>> expected = 
    { 
        { 1, 13, 14, 15, 16, 17, 18, 19, 20, 21},
        { 2,  1, 13, 14, 15, 16, 17, 18, 19, 20},
        { 3,  2,  1, 13, 14, 15, 16, 17, 18, 19},
        { 4,  3,  2,  1, 13, 14, 15, 16, 17, 18},
        { 5,  4,  3,  2,  1, 13, 14, 15, 16, 17},
        { 6,  5,  4,  3,  2,  1, 13, 14, 15, 16},
        { 7,  6,  5,  4,  3,  2,  1, 13, 14, 15},
        { 8,  7,  6,  5,  4,  3,  2,  1, 13, 14},
        { 9,  8,  7,  6,  5,  4,  3,  2,  1, 13},
        {10,  9,  8,  7,  6,  5,  4,  3,  2,  1},
    };
    ASSERT_MATRIX(expected, actual);
}

TEST(Matrixtest, toeplitz_to_matrix){
    vector<uint64_t> toeplitz = {0, 1, 2, 3, 4, 5, 6, 7};
    uint64_t toeplitz_rowsize = 4;
    uint64_t toeplitz_colsize = 5;
    auto actual = toeplitz_to_matrix(toeplitz, toeplitz_rowsize, toeplitz_colsize);
    vector<vector<uint64_t>> expected = 
    { 
        {3, 4, 5, 6, 7},
        {2, 3, 4, 5, 6},
        {1, 2, 3, 4, 5},
        {0, 1, 2, 3, 4}
    };
    ASSERT_MATRIX(expected, actual);
    toeplitz = {0, 1, 2, 3, 4, 5, 6};
    toeplitz_rowsize = 4;
    toeplitz_colsize = 4;
    expected = 
    { 
        {3, 4, 5, 6},
        {2, 3, 4, 5},
        {1, 2, 3, 4},
        {0, 1, 2, 3}
    };
    ASSERT_MATRIX(expected, actual);
}

TEST(Matrixtest, matrix_dot_vector){
    MemoryPoolHandle pool_ = MemoryPoolHandle::Global();
    uint64_t array_size = 10;
    uint64_t coeff_degree = array_size;
    auto modulus = Modulus(0x7e00001ULL);
    vector<std::uint64_t> arr(array_size);
    vector<vector<uint64_t>> matrix(coeff_degree, vector<uint64_t>(coeff_degree));

    for(uint64_t i = 0;i < array_size;i++){
        arr[i] = i+1;
    }
    util::init_matrix_identity(matrix, coeff_degree, 2);
    matrix[0][1] = 1;

    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(actual, coeff_degree, pool_);
    util::matrix_dot_vector(matrix, coeff_degree, arr.data(), modulus, coeff_degree, actual);
    vector<uint64_t> expected = { 4, 4, 6, 8, 10, 12, 14, 16, 18, 20};
    ASSERT_ARR(expected, actual, expected.size());
}

TEST(Matrixtest, matrix_dot_vector_rnd){
    MemoryPoolHandle pool_ = MemoryPoolHandle::Global();
    auto modulus = Modulus(0x7e00001ULL);

    uint64_t vector_len = 4;
    vector<vector<uint64_t>> matrix = 
    {
        { 2, 3, negate_uint_mod(4, modulus) },
        {11, 8, 7},
        { 1, 3, 9},
        { 4, negate_uint_mod(2, modulus), 3}
    };
    vector<uint64_t> arr = {3, 7, 5};

    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(actual, vector_len, pool_);
    util::matrix_dot_vector(matrix, vector_len, arr.data(), modulus, arr.size(), actual);
    vector<uint64_t> expected = {7, 124, 69, 13};
    ASSERT_ARR(expected, actual, vector_len);
}

TEST(Matrixtest, matrix_dot_vector_time){
    MemoryPoolHandle pool_ = MemoryPoolHandle::Global();
    uint64_t coeff_degree = 1024;
    uint64_t matrix_valid_rowsize = 16;
    auto modulus = Modulus(0x7e00001ULL);
    vector<std::uint64_t> arr(matrix_valid_rowsize);
    vector<vector<uint64_t>> matrix(coeff_degree, vector<uint64_t>(coeff_degree));

    for(uint64_t i = 0;i < arr.size();i++){
        arr[i] = i+1;
    }
    util::init_matrix_identity(matrix, coeff_degree, 2);
    matrix[0][1] = 1;

    SEAL_ALLOCATE_GET_COEFF_ITER(poly_vector, coeff_degree, pool_);
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(result, coeff_degree, pool_);
    util::set_poly(arr.data(), coeff_degree, 1, poly_vector);
    auto time_s = chrono::high_resolution_clock::now();
    util::matrix_dot_vector(matrix, matrix_valid_rowsize, poly_vector, modulus, coeff_degree, result);
    auto time_e = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_e - time_s);
    cout << "time: " << time_diff.count() << " us" << endl;
}

TEST(Matrixtest, packedconv_matrix_dot_vector){
    MemoryPoolHandle pool_ = MemoryPoolHandle::Global();
    uint64_t poly_degree = 16; // should be power of two
    auto modulus = Modulus(0x7e00001ULL);
    vector<vector<uint64_t>> kernels = { {3, 1, 2}};
    uint64_t input_dim = 7;
    uint64_t block_size = get_blocksize(input_dim, kernels[0].size(), 0);
    vector<uint64_t> c1(poly_degree);
    sample_rn(c1.data(), c1.size(), modulus);
    vector<KernelInfo> kernel_info = pack_kernel(kernels, block_size, modulus);
    vector<vector<uint64_t>> conved_C1(poly_degree, vector<uint64_t>(poly_degree));
    matrix_dot_matrix_toeplitz_mod(kernel_info, c1.data(), poly_degree, conved_C1, modulus);
    vector<uint64_t> right_vec(poly_degree);
    sample_rn(right_vec.data(), right_vec.size(), modulus);
    ASSERT_EQ(poly_degree, right_vec.size());

    // actual
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(actual, poly_degree, pool_);
    util::matrix_dot_vector(conved_C1, kernel_info[0].kernel_size-1, right_vec.data(), modulus, poly_degree, actual);
    //print_iter(actual, poly_degree);
    CoeffIter actual_toeplitz = actual + kernel_info[0].kernel_size-1;
    util::toeplitz_dot_vector(kernel_info[0].toeplitz_diagonal_scalars, right_vec.data(), block_size, poly_degree, modulus, actual_toeplitz, pool_);

    // expected
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(expected, poly_degree, pool_);
    util::matrix_dot_vector(conved_C1, poly_degree, right_vec.data(), modulus, poly_degree, expected);

    ASSERT_ARR(expected, actual, poly_degree);
}

TEST(Matrixtest, matrix_product_diagonal){
    auto modulus = Modulus(0x7e00001ULL);
    vector<uint64_t> kernel_L = {1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 1, 9, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<uint64_t> kernel_L_indexes = get_nonzero_indexlist(kernel_L);
    vector<uint64_t> list_R = { 2, 3, 1, 6, 3, 2, 6, 4, 12, 2, 8, 7,9, 2, 12, 7, 5, 2, 1, 13, 10, 4, 2, 4, 3, 1, 2, 8, 6, 12, 7, 3, 5, 11, 2, 6, 4};
    uint64_t colsize_R = 4;
    uint64_t rowsize_R = 7;
    vector<pair<uint64_t, uint64_t>> diagonalpairlist;
    diagonalpairlist.reserve(kernel_L_indexes.size());
    // offset = 0
    // expected: (11, 1), (9, 1), (3, 1), (0, 4), (8, 1), (71, 1), (161, 2), (153, 1), (90, 1), (0, 9)
    int64_t offset = 0;
    matrix_product_diagonal(offset, colsize_R, rowsize_R, kernel_L, kernel_L_indexes, list_R, modulus, diagonalpairlist);
    vector<uint64_t> values = { 11, 9, 3, 0, 8, 71, 161, 153, 90, 0};
    vector<uint64_t> counts = {1, 1, 1, 4, 1, 1, 2, 1, 1, 9};
    vector<pair<uint64_t, uint64_t>> expected = make_pair_vector(values, counts);
    ASSERT_PAIR_VEC(expected, diagonalpairlist);
    // offset < 0
    // expected: (6, 1), (0, 4), (12, 1), (30, 1), (110, 2), (98, 1), (80, 1), (0, 7)
    diagonalpairlist.clear();
    offset = -2;
    matrix_product_diagonal(offset, colsize_R, rowsize_R, kernel_L, kernel_L_indexes, list_R, modulus, diagonalpairlist);
    values = {6, 0, 12, 30, 110, 98, 80, 0};
    counts = {1, 4, 1, 1, 2, 1, 1, 7};
    expected = make_pair_vector(values, counts);
    ASSERT_PAIR_VEC(expected, diagonalpairlist);
    diagonalpairlist.clear();

    // offset > 0
    // expected: (22, 1), (21, 1), (9, 1), (0, 4), (9, 1), (27, 1), (147, 2), (138, 1), (120, 1), (0, 9)
    offset = 2;
    matrix_product_diagonal(offset, colsize_R, rowsize_R, kernel_L, kernel_L_indexes, list_R, modulus, diagonalpairlist);
    values = {22, 21, 9, 0, 9, 27, 147, 138, 120, 0};
    counts = {1, 1, 1, 4, 1, 1, 2, 1, 1, 9};
    expected = make_pair_vector(values, counts);
    ASSERT_PAIR_VEC(expected, diagonalpairlist);
    diagonalpairlist.clear();

    //kernel_L = {1, 2, 3, 0, 0, 0};
    //kernel_L_indexes = get_nonzero_indexlist(kernel_L);
    offset = -10;
    matrix_product_diagonal(offset, 15, rowsize_R, kernel_L, kernel_L_indexes, list_R, modulus, diagonalpairlist);
    expected = {make_pair(39, 1)};
    ASSERT_PAIR_VEC(expected, diagonalpairlist);
    diagonalpairlist.clear();

    offset = 10;
    matrix_product_diagonal(offset, 15, rowsize_R, kernel_L, kernel_L_indexes, list_R, modulus, diagonalpairlist);
    print_pair_vector(diagonalpairlist);
}
