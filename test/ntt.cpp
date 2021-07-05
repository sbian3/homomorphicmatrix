#include "testcommon.h"


TEST(NTTtest, multiply_poly){
    uint64_t poly_degree = 8;
    uint64_t bigger_power_of_two = get_bigger_poweroftwo(poly_degree);
    uint64_t coeff_count_power = get_power_of_two(bigger_power_of_two);
    MemoryPoolHandle pool_ = MemoryPoolHandle::Global();
    auto modulus = Modulus(0x7e00001);
    NTTTables ntt_tables_(coeff_count_power, modulus, pool_);

    vector<uint64_t> left_vec(bigger_power_of_two);
    vector<uint64_t> right_vec(bigger_power_of_two);
    vector<uint64_t> result(bigger_power_of_two);
    left_vec[0] = 5;
    left_vec[1] = 2;
    right_vec[0] = 1;
    right_vec[1] = 3;
    ntt_negacyclic_harvey(left_vec.data(), ntt_tables_);
    ntt_negacyclic_harvey(right_vec.data(), ntt_tables_);
    dyadic_product_coeffmod(left_vec.data(), right_vec.data(), bigger_power_of_two, modulus, result);
    inverse_ntt_negacyclic_harvey(result.data(), ntt_tables_);
    ASSERT_EQ(5, result[0]);
    ASSERT_EQ(17, result[1]);
    ASSERT_EQ(6, result[2]);
}

TEST(NTTtest, multiply_poly_negacyclic){
    uint64_t poly_degree = 2;
    uint64_t coeff_count_power = get_power_of_two(poly_degree);
    MemoryPoolHandle pool_ = MemoryPoolHandle::Global();
    auto modulus = Modulus(0x7e00001);
    NTTTables ntt_tables_(coeff_count_power, modulus, pool_);

    vector<uint64_t> left_vec(poly_degree);
    vector<uint64_t> right_vec(poly_degree);
    vector<uint64_t> result(poly_degree);
    left_vec[0] = 5;
    left_vec[1] = 2;
    right_vec[0] = 1;
    right_vec[1] = 3;
    ntt_negacyclic_harvey(left_vec.data(), ntt_tables_);
    ntt_negacyclic_harvey(right_vec.data(), ntt_tables_);
    dyadic_product_coeffmod(left_vec.data(), right_vec.data(), poly_degree, modulus, result);
    inverse_ntt_negacyclic_harvey(result.data(), ntt_tables_);
    ASSERT_EQ(negate_uint_mod(1, modulus), result[0]);
    ASSERT_EQ(17, result[1]);
}


TEST(NTTtest, poweroftwo){
    auto actual = get_bigger_poweroftwo(63);
    ASSERT_EQ(actual, uint64_t(64));
    actual = get_bigger_poweroftwo(32);
    ASSERT_EQ(actual, uint64_t(32));
}

TEST(NTTtest, toeplitz_to_circ){
    uint64_t toep_rowsize = 5;
    uint64_t toep_colsize = 8;
    uint64_t toep_size = toep_rowsize + toep_colsize - 1;
    Modulus modulus(23);

    vector<uint64_t> toep(toep_size);
    vector<uint64_t> actual(toep_colsize * 2);
    for(uint64_t i = 0;i < toep_size;i++){
        toep[i] = i+1;
    }
    toeplitz_to_circ(toep, toep_rowsize, toep_colsize, actual, modulus);

    vector<uint64_t> expected = { 5, 4, 3, 2, 1, 0, 0, 0, 0, 11, 12, 13, 14, 15, 16, 17 };
    ASSERT_ARR(expected.data(), actual.data(), expected.size());
}

TEST(NTTtest, toeplitz_dot_vector){
    // params
    MemoryPoolHandle pool_ = MemoryPoolHandle::Global();
    uint64_t poly_degree = 1024;
    uint64_t toeplitz_rowsize = 16;
    uint64_t toeplitz_colsize = poly_degree;
    auto modulus = Modulus(0x7e00001ULL);

    vector<uint64_t> toeplitz(toeplitz_rowsize + toeplitz_colsize - 1);
    vector<uint64_t> right_vec(poly_degree);
    sample_rn(toeplitz.data(), toeplitz.size(), modulus);
    sample_rn(right_vec.data(), right_vec.size(), modulus);
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(actual, toeplitz_rowsize, pool_);

    // actual
    auto time_s = chrono::high_resolution_clock::now();
    toeplitz_dot_vector(toeplitz, right_vec.data(), toeplitz_rowsize, toeplitz_colsize, modulus, actual, pool_);
    auto time_e = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_e - time_s);
    cout << "time: " << time_diff.count() << " us" << endl;

    // expect
    auto toeplitz_matrix = toeplitz_to_matrix(toeplitz, toeplitz_rowsize, toeplitz_colsize);
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(expected, toeplitz_rowsize, pool_);
    matrix_dot_vector(toeplitz_matrix, toeplitz_rowsize, right_vec.data(), modulus, toeplitz_colsize, expected);

    ASSERT_ARR(expected, actual, toeplitz_rowsize);
}

