#include "testprimitives.h"
#include "testutil.h"
#include "seal/util/uintlinarith.h"
#include "seal/util/packedconv.h"

using namespace seal::util;

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
    uint64_t size=5;
    uint64_t conv_size = 5;
    Modulus modulus(11);
    vector<vector<int64_t>> matrix(size, vector<int64_t>(size));
    vector<uint64_t> a(conv_size);
    for(uint64_t i = 0;i<conv_size;i++){
        a[i] = rand() % modulus.value();
    }
    for(uint64_t i = 0;i<conv_size;i++){
    }
    uint64_t rotate = 0;
    util::init_matrix_rotate(matrix, size, rotate, 1);
    rotate++;
    util::init_matrix_rotate(matrix, size, rotate, 2);
    rotate++;
    util::init_matrix_rotate(matrix, size, rotate, 3);
    util::print_matrix(matrix);
}

void test_print_iter(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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

void test_init_matrix_coeff(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
    uint64_t array_size = 10;
    Modulus modulus(15);
    cout << "modulus: " << modulus.value() << endl;
    // polynomial degree
    uint64_t coeff_degree = array_size;
    vector<std::uint64_t> arr(array_size);
    vector<vector<uint64_t>> matrix(coeff_degree, vector<uint64_t>(coeff_degree));

    for(size_t i = 0;i < array_size;i++){
        arr[i] = i;
    }
    SEAL_ALLOCATE_GET_COEFF_ITER(coeff_iter, coeff_degree, pool_);
    util::set_poly(arr.data(), coeff_degree, 1, coeff_iter);
    util::init_matrix_with_coeff(matrix, coeff_degree, coeff_iter, modulus);
    util::print_matrix(matrix);
}

void test_init_matrix_uint_by_kernel(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
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

void test_dianonal_kernel(){
    vector<uint64_t> kernel = {1, 2, 3};
    uint64_t matrix_size = 28;
    uint64_t dest_size = matrix_size - kernel.size() + 1;
    vector<vector<uint64_t>> block(dest_size, vector<uint64_t>(matrix_size));
    for(auto i = 0U; i < kernel.size();i++){
        util::init_matrix_diagonal(block, dest_size, kernel[i], i);
    }
    util::print_matrix(block);
}

void test_init_matrix_2dconv(){
    //vector<vector<uint64_t>> kernel = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    vector<vector<uint64_t>> kernel = {{1, 2}, {4, 5}};
    uint64_t input_size = 3;
    uint64_t matrix_size = input_size * input_size;
    uint64_t dest_size = input_size - kernel.size() + 1;
    uint64_t matrix_colsize = dest_size * dest_size;
    vector<vector<uint64_t>> matrix(matrix_colsize, vector<uint64_t>(matrix_size));
    util::init_matrix_2dconv(matrix, input_size, kernel);
    util::print_matrix(matrix);
}

void test_create_diagonal(){
    vector<uint64_t> kernel = {1,2,3};
    uint64_t colsize = 5;
    uint64_t rowsize = 5;
    Modulus modulus(7);
    vector<uint64_t> diagonal(colsize + rowsize - 1);
    vector<uint64_t> indexes = create_diagonal_list(kernel, colsize, rowsize, modulus, diagonal);
    cout << "diagonal: ";
    for(uint64_t i = 0;i < diagonal.size();i++){
        cout << diagonal[i] << " ";
    }
    cout << endl;

    cout << "index: ";
    for(uint64_t i = 0;i < indexes.size();i++){
        cout << indexes[i] << " ";
    }
    cout << endl;
}

void test_diagonal_from_submat(){
    uint64_t start_col = 4;
    uint64_t colsize = 3;
    uint64_t a_size = 10;
    Modulus modulus(15);
    vector<uint64_t> a(a_size);
    for(auto i = 0L;i < a.size();i++){
        a[i] = i;
    }
    vector<uint64_t> diagonal;
    diagonal = create_diagonal_from_submatrix(a, a_size, start_col, colsize, modulus);
    cout << "start_col: " << start_col << " " << "colsize: " << colsize << endl;
    cout << "diagonal list of submatrix: ";
    for(auto i =0L;i < diagonal.size();i++){
        cout << diagonal[i] << " ";
    }
    cout << endl;
}

void test_prod_diagonal(){
    vector<uint64_t> kernel = {1, 2, 4};
    uint64_t colsize_K = 4;
    uint64_t rowsize_K = 5;
    uint64_t colsize_R = 5;
    uint64_t rowsize_R = 6;
    vector<uint64_t> diagonallist_R = {3, 2, 1, 5, 3, 8, 3, 4, 7, 1};
    Modulus modulus(11);
    vector<uint64_t> kernel_diagonallist(colsize_K + rowsize_K - 1);
    vector<uint64_t> indexes = create_diagonal_list(kernel, colsize_K, rowsize_K, modulus, kernel_diagonallist);
    int64_t k = static_cast<int64_t>(colsize_K);
    k  = -k + 1;
    vector<vector<uint64_t>> matrix_product_diagonals(colsize_K + rowsize_R-1);
    uint64_t index = 0;
    for(;k < static_cast<int64_t>(rowsize_R);k++){
        vector<uint64_t> diagonal_vec;
        diagonal_vec = util::matrix_product_diagonal(k, colsize_R, rowsize_R, kernel_diagonallist, indexes, diagonallist_R, modulus);
        cout << "offset: " << k << endl;
        util::print_iter(diagonal_vec, diagonal_vec.size());
        matrix_product_diagonals[index] = diagonal_vec;
        index++;
    }
    cout << "result matrix" << endl;
    vector<vector<uint64_t>> result(colsize_K , vector<uint64_t>(rowsize_R));
    util::diagonallist_to_matrix(matrix_product_diagonals, 0, 0, colsize_K, rowsize_R, result);
    util::print_matrix(result);
}

/////////////////
// tests for kernelinfo
/////////////////
void test_init_kernelinfo(){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6} };
    uint64_t block_size = 32;
    Modulus modulus(7);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, block_size, modulus);
    for(uint64_t i = 0;i < kernel_infos.size();i++){
        kernel_infos[i].print();
    }
}

void test_kernelinfo_to_rowpair(){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6} };
    uint64_t block_size = 32;
    Modulus modulus(7);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, block_size, modulus);
    for(uint64_t i = 0;i < kernel_infos.size();i++){
        KernelInfo kinfo = kernel_infos[i];
        cout << "----kernel info[ " << i << " ]----" << endl;
        // for each row
        vector<pair<uint64_t, uint64_t>> rowvec = kinfo.make_rowpair(modulus);
        for(uint64_t j = 0;j < block_size;j++){
            for(uint64_t k = 0;k < rowvec.size();k++){
                cout << "(index, value) = "  << rowvec[k].first << ", " << rowvec[k].second << endl;
            }
            cout << "next row" << endl;
            kinfo.pair_nextcol(rowvec, modulus);
        }
    }
}

void test_pack_kernel_to_matrix(){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6} };
    uint64_t block_size = 6;
    uint64_t matrix_size = 18;
    Modulus modulus(7);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, block_size, modulus);
    vector<vector<uint64_t>> matrix(matrix_size, vector<uint64_t>(matrix_size));
    pack_kernel_to_matrix(kernel_infos, matrix);
    util::print_matrix(matrix);
}


void test_kernel_dot_c1(){
    vector<uint64_t> c1 = {1, 4, 2, 5, 1, 3};
    uint64_t poly_degree = c1.size();
    vector<vector<uint64_t>> kernels = { {3, 1, 2}, {2, 3, 5} };
    uint64_t block_size = 3;
    Modulus modulus(7);
    vector<KernelInfo> kernel_info = pack_kernel(kernels, block_size, modulus);
    vector<vector<uint64_t>> result(poly_degree, vector<uint64_t>(poly_degree));
    matrix_dot_matrix_toeplitz_mod(kernel_info, c1, poly_degree, result, modulus);
    util::print_matrix(result);
}

void test_pack_input(){
    vector<vector<uint64_t>> inputs = { {1 , 2, 3}, {4, 5, 6} };
    uint64_t block_size = 6;
    uint64_t poly_size = 20;
    vector<uint64_t> plain_vec = pack_input(inputs, block_size, poly_size);
    for(uint64_t i = 0;i < plain_vec.size();i++){
       cout << plain_vec[i] << " "; 
    }
    cout << endl;
}

void vector_to_plain(){
    print_example_banner("Packed Convolution Benchmark");
    vector<uint64_t> coeff = {1, 2, 3,4, 10};
    Plaintext plain(coeff);
    //print_plain(plain, coeff.size());
}
