#include "testcommon.h"

TEST(Matrixtest, init_circulant_matrix){
    uint64_t poly_degree = 10;
    Modulus modulus(23);

    vector<vector<uint64_t>> actual(poly_degree, vector<uint64_t>(poly_degree));
    vector<uint64_t> coeff(poly_degree);
    for(uint64_t i = 0;i < coeff.size();i++){
        coeff[i] = i+1;
    }

    init_matrix_with_coeff(actual, poly_degree, coeff.data(), modulus);
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

    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(result, coeff_degree, pool_);
    util::matrix_dot_vector(matrix, coeff_degree, arr.data(), modulus, coeff_degree, result);
    util::print_iter(result, coeff_degree);
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
