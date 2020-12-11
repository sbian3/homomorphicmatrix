#include "testutil.h"
#include "seal/util/uintlinarith.h"

void test_matrix_initialization(){
    uint64_t size=3;
    vector<vector<int64_t>> matrix(size, vector<int64_t>(size));
    //util::init_matrix_rotate(matrix, size, 1, 1);
    util::print_matrix(matrix);
    //util::init_matrix_rotate(matrix, size, -1, 1);
    util::print_matrix(matrix);
    util::init_matrix_identity(matrix, size, 2);
    util::print_matrix(matrix);
    //init_matrix_identity_rnd(matrix, size);
    util::print_matrix(matrix);
}

void test_matrix_init_partial(){
    uint64_t size=20;
    vector<vector<uint64_t>> matrix(size, vector<uint64_t>(size));
    vector<uint64_t> kernel = {1,2,3};
    uint64_t modulus = 7;
    util::init_matrix_rotate_partial(matrix, kernel, 10, 10, modulus);
    util::print_matrix(matrix);
}

void test_matconv(){
    uint64_t size=3;
    vector<vector<int64_t>> matrix(size, vector<int64_t>(size));
    uint64_t rotate = 0;
    util::init_matrix_rotate(matrix, size, rotate, 1);
    rotate++;
    util::init_matrix_rotate(matrix, size, rotate, 2);
    rotate++;
    util::init_matrix_rotate(matrix, size, rotate, 3);
    util::print_matrix(matrix);
}

void test_print_iter(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    // how many "c"?
    uint64_t poly_count = 4;
    // polynomial degree
    uint64_t mod_degree = 10;
    // rns count
    uint64_t coeff_mod_size = 3;
    uint64_t array_size = poly_count * mod_degree * coeff_mod_size;

    vector<std::uint64_t> arr(array_size);

    for(uint64_t i = 0;i < array_size;i++){
        arr[i] = i+1;
    }

    SEAL_ALLOCATE_GET_POLY_ITER(poly_iter, poly_count, mod_degree, coeff_mod_size, pool_);
    util::set_poly_array(arr.data(), poly_count, mod_degree, coeff_mod_size, poly_iter);
    poly_iter++;
    SEAL_ITERATE(poly_iter, 1, [&](auto I){
            util::print_iter(I, coeff_mod_size);
              cout << "end of c ...." << endl;
            });
}

void test_innerprod(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 10;
    // polynomial degree
    uint64_t coeff_degree = array_size;
    Modulus modulus(100);

    vector<std::uint64_t> arr(array_size);
    vector<std::uint64_t> arr2(array_size);

    uint64_t expected_result = 0;
    for(uint64_t i = 0;i < array_size;i++){
        arr[i] = i+1;
        arr2[i] = i+2;
        expected_result += arr[i] * arr2[i] % modulus.value();
    }
    expected_result %= modulus.value();

    SEAL_ALLOCATE_GET_COEFF_ITER(iter1, coeff_degree, pool_);
    SEAL_ALLOCATE_GET_COEFF_ITER(iter2, coeff_degree, pool_);
    util::set_poly(arr.data(), coeff_degree, 1, iter1);
    util::set_poly(arr2.data(), coeff_degree, 1, iter2);
    util::print_iter(iter1, coeff_degree);
    util::print_iter(iter2, coeff_degree);
    uint64_t result = inner_product_coeffmod(iter1, iter2, coeff_degree, modulus);
    cout << "innerprod result: " << result << endl;
    cout << "expected result: " << expected_result << endl;
}

void test_conviter(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 10;
    // polynomial degree
    uint64_t coeff_degree = array_size;
    Modulus modulus(100);

    vector<std::uint64_t> arr(array_size);
    vector<std::uint64_t> arr2(array_size);

    uint64_t expected_result = 0;
    for(uint64_t i = 0;i < array_size;i++){
        arr[i] = i+1;
        arr2[i] = i+2;
        expected_result += arr[i] * arr2[i] % modulus.value();
    }
    expected_result %= modulus.value();

    SEAL_ALLOCATE_GET_COEFF_ITER(iter1, coeff_degree, pool_);
    SEAL_ALLOCATE_GET_COEFF_ITER(iter2, coeff_degree, pool_);
    util::set_poly(arr2.data(), coeff_degree, 1, iter2);
    util::print_iter(iter2, coeff_degree);
    util::conv_negacyclic(arr, iter2, coeff_degree, modulus, iter1);
    util::print_iter(iter1, coeff_degree);
}

void test_innerprod_vector(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 10;
    // polynomial degree
    uint64_t coeff_degree = array_size;
    Modulus modulus(100);

    vector<std::int64_t> arr(array_size);
    vector<std::uint64_t> arr2(array_size);

    int64_t expected_result = 0;
    for(size_t i = 0;i < array_size;i++){
        arr[i] = -i+1;
        arr2[i] = i+2;
        expected_result += (arr[i] * static_cast<int64_t>(arr2[i]));
        expected_result %= modulus.value();
    }

    SEAL_ALLOCATE_GET_COEFF_ITER(iter2, coeff_degree, pool_);
    util::set_poly(arr2.data(), coeff_degree, 1, iter2);
    util::print_iter(iter2, coeff_degree);
    uint64_t result = inner_product_coeffmod(arr, iter2, coeff_degree, modulus);
    cout << "innerprod result: " << result << endl;
    cout << "expected result: " << expected_result << endl;
}

void test_init_matrix(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 10;
    // polynomial degree
    uint64_t coeff_degree = array_size;
    vector<std::uint64_t> arr(array_size);
    vector<vector<int64_t>> matrix(coeff_degree, vector<int64_t>(coeff_degree));

    for(size_t i = 0;i < array_size;i++){
        arr[i] = i+1;
    }
    SEAL_ALLOCATE_GET_COEFF_ITER(coeff_iter, coeff_degree, pool_);
    util::set_poly(arr.data(), coeff_degree, 1, coeff_iter);
    util::init_matrix_with_coeff(matrix, coeff_degree, coeff_iter);
    util::print_matrix(matrix);
}

void test_init_matrix_uint(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 10;
    Modulus modulus(15);
    // polynomial degree
    uint64_t coeff_degree = array_size;
    vector<std::uint64_t> arr(array_size);
    vector<vector<uint64_t>> matrix(coeff_degree, vector<uint64_t>(coeff_degree));

    for(size_t i = 0;i < array_size;i++){
        arr[i] = i+1;
    }
    SEAL_ALLOCATE_GET_COEFF_ITER(coeff_iter, coeff_degree, pool_);
    util::set_poly(arr.data(), coeff_degree, 1, coeff_iter);
    util::init_matrix_with_coeff(matrix, coeff_degree, coeff_iter, modulus);
    util::print_matrix(matrix);
}

void test_init_matrix_uint_by_kernel(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 3;
    Modulus modulus(15);
    // polynomial degree
    uint64_t coeff_degree = 10;
    vector<std::uint64_t> arr(array_size);
    vector<vector<uint64_t>> matrix(coeff_degree, vector<uint64_t>(coeff_degree));

    for(size_t i = 0;i < array_size;i++){
        arr[i] = i+1;
    }
    SEAL_ALLOCATE_GET_COEFF_ITER(coeff_iter, coeff_degree, pool_);
    util::set_poly(arr.data(), coeff_degree, 1, coeff_iter);
    util::init_matrix_with_coeff(matrix, coeff_degree, arr.data(), array_size, modulus);
    util::print_matrix(matrix);
}


void test_matrix_conversion_with_coeffiter(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 10;
    uint64_t coeff_degree = array_size;
    Modulus modulus(100);
    vector<std::uint64_t> arr(array_size);
    vector<vector<int64_t>> matrix(coeff_degree, vector<int64_t>(coeff_degree));

    for(uint64_t i = 0;i < array_size;i++){
        arr[i] = i+1;
    }
    util::init_matrix_identity(matrix, coeff_degree, 2);
    matrix[0][1] = 1;

    SEAL_ALLOCATE_GET_COEFF_ITER(poly_vector, coeff_degree, pool_);
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(result, coeff_degree, pool_);
    util::set_poly(arr.data(), coeff_degree, 1, poly_vector);
    util::matrix_dot_vector(matrix, poly_vector, modulus, coeff_degree, result);
    util::print_iter(result, coeff_degree);
}

void test_secret_product(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 3;
    uint64_t coeff_degree = array_size;
    Modulus modulus(100);
    vector<std::uint64_t> arr(array_size);
    arr[0] = 1; arr[1] = 0; arr[2] = 1;
    vector<uint64_t> s_arr = {static_cast<uint64_t>(-1), 0, 1};
    vector<vector<int64_t>> matrix(coeff_degree, vector<int64_t>(coeff_degree));
    util::init_matrix_identity(matrix, coeff_degree, 1);
    util::CoeffIter c = util::CoeffIter(arr);
    print_iter(c, coeff_degree);
    util::CoeffIter s = util::CoeffIter(s_arr);
    print_iter(s, coeff_degree);
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(result, coeff_degree, pool_);
    util::secret_product_with_matrix(matrix, coeff_degree, c, s, modulus, result);
    util::print_iter(result, coeff_degree);
}

void test_matrix_conversion_with_rnsiter(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    // modulus chain
    vector<Modulus> mod_chain = {10, 15, 100};
    // polynomial degree
    uint64_t coeff_degree = 10;
    // rns count
    uint64_t coeff_mod_size = mod_chain.size();

    uint64_t array_size = coeff_degree * coeff_mod_size;
    vector<std::uint64_t> arr(array_size);
    vector<vector<int64_t>> matrix(coeff_degree, vector<int64_t>(coeff_degree));

    for(uint64_t i = 0;i < coeff_mod_size;i++){
        for(uint64_t j=0;j < coeff_degree;j++){
            arr[i* coeff_degree+j] = j;
        }
    }
    util::init_matrix_identity(matrix, coeff_degree, 2);
    matrix[0][1] = 1;
    SEAL_ALLOCATE_GET_RNS_ITER(rns_iter, coeff_degree, coeff_mod_size, pool_);
    util::set_poly(arr.data(),  coeff_degree, coeff_mod_size, rns_iter);
    cout << "print rns iter first" << endl;
    print_iter(rns_iter, coeff_mod_size);
    SEAL_ALLOCATE_ZERO_GET_RNS_ITER(result, coeff_degree, coeff_mod_size, pool_);
    util::matrix_dot_vector(matrix, coeff_mod_size,  rns_iter, mod_chain, result);
    cout << "print converted rns iter " << endl;
    print_iter(result, coeff_mod_size);
}

void test_matrix_dot_product(){
    uint64_t coeff_degree = 1000;
    uint64_t modulus = 10000;
    vector<vector<int64_t>> matrix(coeff_degree, vector<int64_t>(coeff_degree));
    vector<vector<int64_t>> matrix2(coeff_degree, vector<int64_t>(coeff_degree));
    vector<vector<int64_t>> result(coeff_degree, vector<int64_t>(coeff_degree));
    util::init_matrix_rand_mod(matrix, coeff_degree, modulus);
    cout << "first matrix: " << endl;
    //util::print_matrix(matrix);
    util::init_matrix_rand_mod(matrix2, coeff_degree, modulus);
    cout << "second matrix: " << endl;
    //util::print_matrix(matrix2);
    cout << "calculating matrix dot product..." << endl;
    util::matrix_dot_product_mod(matrix, matrix2, result, modulus);
    //util::print_matrix(result);
}

void test_matrix_dot_product_uint(uint64_t coeff_degree){
    Modulus modulus = 10000;
    vector<vector<uint64_t>> matrix(coeff_degree, vector<uint64_t>(coeff_degree));
    vector<vector<uint64_t>> matrix2(coeff_degree, vector<uint64_t>(coeff_degree));
    vector<vector<uint64_t>> result(coeff_degree, vector<uint64_t>(coeff_degree));
    util::init_matrix_rand_mod(matrix, coeff_degree, modulus.value());
    cout << "first matrix: " << endl;
    //util::print_matrix(matrix);
    util::init_matrix_rand_mod(matrix2, coeff_degree, modulus.value());
    cout << "second matrix: " << endl;
    //util::print_matrix(matrix2);
    cout << "calculating matrix dot product..." << endl;
    util::matrix_dot_product_mod(matrix, matrix2, result, modulus);
    //util::print_matrix(result);
}

void test_matrix_dot_product_uint_t(uint64_t coeff_degree){
    Modulus modulus = 10000;
    vector<vector<uint64_t>> matrix(coeff_degree, vector<uint64_t>(coeff_degree));
    vector<vector<uint64_t>> matrix2(coeff_degree, vector<uint64_t>(coeff_degree));
    vector<vector<uint64_t>> result(coeff_degree, vector<uint64_t>(coeff_degree));
    util::init_matrix_rand_mod(matrix, coeff_degree, modulus.value());
    cout << "first matrix: " << endl;
    //util::print_matrix(matrix);
    util::init_matrix_rand_mod(matrix2, coeff_degree, modulus.value());
    cout << "second matrix: " << endl;
    //util::print_matrix(matrix2);
    cout << "calculating matrix dot product..." << endl;
    util::matrix_dot_product_mod_t(matrix, matrix2, result, modulus);
    //util::print_matrix(result);
}

void matrix_product_benchmark(){
    uint64_t coeff_degree = 100;
    cout << "coeff degree: " << coeff_degree << endl;
    test_matrix_dot_product_uint_t(coeff_degree);
    coeff_degree *= 2;
    cout << "coeff degree: " << coeff_degree << endl;
    test_matrix_dot_product_uint_t(coeff_degree);
    coeff_degree *= 2;
    cout << "coeff degree: " << coeff_degree << endl;
    test_matrix_dot_product_uint_t(coeff_degree);
    coeff_degree *= 2;
    cout << "coeff degree: " << coeff_degree << endl;
    test_matrix_dot_product_uint_t(coeff_degree);
}

void test_inverse(){
    uint64_t operand = 5;
    uint64_t mod = 11;
    uint64_t result =  util::negate_uint_mod(operand, mod);
    cout << "-" << operand << "mod" << mod <<  "= " << result;
}

void test_util_dot_product_mod(){
    constexpr uint64_t count = 101;
    uint64_t arr1[count];
    uint64_t arr2[count];
    uint64_t expected = 0;
    Modulus modulus(22);
    for(auto i = 0U;i < count;i++){
        arr1[i] = i+1;
        arr2[i] = i+2;
        expected += (i+1) * (i+2);
        expected = expected % modulus.value();
    }
    uint64_t result = util::dot_product_mod(arr1, arr2, count, modulus);
    cout << "expected result: " << expected << endl;
    cout << "result: " << result << endl;
}

void benchmark_singlefunction(){
    uint64_t count;
    cout << "input count: ";
    cin >> count;
    vector<uint64_t> arr1(count);
    vector<uint64_t> arr2(count);
    uint64_t expected = 0;
    uint64_t result = 0;
    Modulus modulus(22);
    for(auto i = 0U;i < count;i++){
        arr1[i] = i+1;
        arr2[i] = i+2;
        expected += (i+1) * (i+2);
        expected = expected % modulus.value();
    }
    auto time_start = chrono::high_resolution_clock::now();
    for(auto i = 0U;i < count;i++){
        result = util::multiply_add_uint_mod(arr1[i], arr2[i], result, modulus);
    }
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
    cout << "expect: " << expected << endl;
    cout << "result: " << result << endl;
    cout << "time: " << time_diff.count() << "ms" << endl;
}

//void test_conv_nega(){
//    vector<uint64_t> kernel = {1,1,2};
//    vector<uint64_t> input = {1, 2, 3, 4};
//    vector<uint64_t> result(20);
//    util::conv_negacyclic(kernel, input, 9, result);
//    for(uint64_t i = 0U;i < result.size();i++){
//        cout << result[i];
//        if(i != result.size() -1){
//            cout << " ";
//        }
//    }
//    cout << endl;
//}