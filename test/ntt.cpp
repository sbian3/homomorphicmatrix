#include "testcommon.h"


TEST(NTTtest, multiply_poly){
    uint64_t poly_degree = 1024;
    uint64_t coeff_count_power = get_power_of_two(poly_degree);
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    print_vector(result, 4);
    ASSERT_EQ(5, result[0]);
}

TEST(NTTtest, multiply_poly_negacyclic){
    uint64_t poly_degree = 2;
    uint64_t coeff_count_power = get_power_of_two(poly_degree);
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    //print_vector(result, 2);
    ASSERT_EQ(negate_uint_mod(1, modulus), result[0]);
    ASSERT_EQ(17, result[1]);
}

TEST(NTTtest, toeplitz_vector){
    uint64_t poly_degree = 10;
    Modulus modulus(11);
    vector<vector<uint64_t>> matrix(poly_degree, vector<uint64_t>(poly_degree));
    vector<uint64_t> coeff(poly_degree);
    for(uint64_t i = 0;i < coeff.size();i++){
        coeff[i] = i+1;
    }
    init_matrix_with_coeff(matrix, poly_degree, coeff.data(), modulus);
    print_matrix(matrix);
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
    vector<uint64_t> toep(toep_size);
    for(uint64_t i = 0;i < toep_size;i++){
        toep[i] = i+1;
    }
    vector<uint64_t> circ(toep_colsize * 2);
    Modulus modulus(23);
    toeplitz_to_circ(toep, toep_rowsize, toep_colsize, circ, modulus);
    vector<uint64_t> expect = { 5, 4, 3, 2, 1, 0, 0, 0, 0, 11, 12, 13, 14, 15, 16, 17 };
    for(uint64_t i = 0;i < circ.size();i++){
        ASSERT_EQ(expect[i], circ[i]);
    }
}

TEST(NTTtest, toeplitz_vector_mult){
    uint64_t poly_degree = 8;
    uint64_t toeplitz_rowsize = 5;
    uint64_t toeplitz_colsize = poly_degree;
    vector<uint64_t> toeplitz(toeplitz_rowsize + toeplitz_colsize - 1);
    vector<uint64_t> right_vec(poly_degree);
    vector<uint64_t> result(poly_degree);
    auto modulus = Modulus(0x7e00001);

    for(uint64_t i = 0;i < toeplitz.size();i++){
        toeplitz[i] = i + 1;
    }
    for(uint64_t i = 0;i < right_vec.size();i++){
        right_vec[i] = i + 2;
    }
    toeplitz_dot_vector(toeplitz, right_vec.data(), toeplitz_rowsize, toeplitz_colsize, modulus, result.data()); 
    print_vector(result, result.size());
}
