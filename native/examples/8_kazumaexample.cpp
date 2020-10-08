#include "seal/util/defines.h"
#include "seal/util/linarith.h"
#include "examples.h"
#include <cassert>

using namespace std;
using namespace seal;

void convert(Plaintext& src, Plaintext& dst, vector<vector<int64_t>>& matrix, uint64_t poly_degree, uint64_t mod_plain){
    assert(matrix.size() == poly_degree);
    assert(matrix[0].size() == poly_degree);
    assert(dst.coeff_count() == poly_degree);
    uint64_t sum ;
    for(auto i = 0U;i < poly_degree;i++){
        sum = 0;
        for( auto j = 0U;j < src.coeff_count();j++ ){
            sum += src[j] * matrix.at(i).at(j) % mod_plain;
        }
        dst[i] = sum % mod_plain;
    }
}

void convert(Ciphertext& src, Ciphertext& dst, vector<vector<int64_t>>& matrix, uint64_t poly_degree,vector<Modulus> mod_chain){
    assert(matrix.size() == poly_degree);
    assert(matrix[0].size() == poly_degree);
    uint64_t sum;
    for(auto size = 0U;size < src.size();size++){
        for(auto coeff_m = 0U; coeff_m < src.coeff_modulus_size();coeff_m++){
            for(auto i = 0U;i < poly_degree;i++){
                sum = 0;
                for(auto j = 0U;j < poly_degree;j++){
                    sum += (src[(size * src.coeff_modulus_size() + coeff_m ) * poly_degree + j] * matrix.at(i).at(j)) % mod_chain[coeff_m].value();
                }
                dst[(size * src.coeff_modulus_size() + coeff_m ) * poly_degree + i] = sum % mod_chain[coeff_m].value();
            }

        }
    }
}

void print_plain_coefficients(Plaintext plain){
    // print all coefficients
    cout << "Plain: [";
    for(int i = 0;i < plain.coeff_count();i++){
        cout << plain[i] << " ";
    }
    cout << "]" << endl;
}

void print_plain(Plaintext plain, uint64_t num){
    if(num > plain.coeff_count()){
        cout << "print_plain: warn: num is bigger than plain coefficient size!" << endl;
        num = plain.coeff_count();
    }
    cout << "Plaintext Coefficients(first " << num << " elements): [";
    for(auto i = 0U;i < num;i++){
        cout << plain[i] << " ";
    }
    cout << "]" << endl;
}

void print_cipher_coefficients(Ciphertext cipher, uint64_t poly_modulus){
    cout << "Cipher Coefficients: [";
    for(auto i = 0U;i < cipher.coeff_modulus_size() * cipher.size() * poly_modulus;i++){
        cout << cipher[i] << " ";
    }
    cout << "]" << endl;
}

void print_cipher(Ciphertext cipher, uint64_t num){
    cout << "Ciphertext Coefficients(first " << num << " elements): [";
    for(auto i = 0U;i < num;i++){
        cout << cipher[i] << " ";
    }
    cout << "]" << endl;
}


bool compare_cipher(Ciphertext cipher1, Ciphertext cipher2, uint64_t poly_modulus){
    uint64_t cipher1_count = cipher1.coeff_modulus_size() * cipher1.size() * poly_modulus;
    uint64_t cipher2_count = cipher2.coeff_modulus_size() * cipher2.size() * poly_modulus;
    if(cipher1_count != cipher2_count){
        cout << "coefficient number does not match: " << cipher1_count << " " << cipher2_count << endl;
        return false;
    }
    for(auto i = 0U;i < cipher1_count;i++){
        if(cipher1[i] != cipher2[i]){
            cout << i << "-th coefficient does not match!" << endl;
            cout << "cipher1: " << cipher1[i] << " cipher2: " << cipher2[i] << endl;
            return false;
        }
    }
    return true;
}

bool compare_cipher_all(Ciphertext cipher1, Ciphertext cipher2, uint64_t poly_modulus){
    uint64_t cipher1_count = cipher1.coeff_modulus_size() * cipher1.size() * poly_modulus;
    uint64_t cipher2_count = cipher2.coeff_modulus_size() * cipher2.size() * poly_modulus;
    if(cipher1_count != cipher2_count){
        cout << "coefficient number does not match: " << cipher1_count << " " << cipher2_count << endl;
        return false;
    }
    bool flag = true;
    for(auto i = 0U;i < cipher1_count;i++){
        if(cipher1[i] != cipher2[i]){
            cout << i << "-th coefficient does not match: " << "(" << cipher1[i] << ", " << cipher2[i] << ")" <<endl;
            flag = false;
        }
    }
    return flag;
}

void print_cipher_info(Ciphertext& c, Decryptor& decryptor){
    cout << "Size: " << c.size() << endl;
    cout << "Coeff modulus size: " << c.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = c.size() * c.poly_modulus_degree() * c.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(c) << " bits" << endl;
    cout << "is NTT form? " << c.is_ntt_form() << endl;
}


void init_matrix_identity_rnd(vector<vector<int64_t>>& matrix, uint64_t poly_modulus_degree){
    int64_t scale = 2;
    for(auto i = 0U;i < poly_modulus_degree;i++){
        for(auto j = 0U;j < poly_modulus_degree;j++){
            if(i == j)
                matrix[i][j] = scale;
            else
                matrix[i][j] = 0;
        }
    }
    matrix[0][0] = 1;
}




void test_conversion(){
    uint64_t size=3;
    vector<vector<int64_t>> matrix(size, vector<int64_t>(size));
    util::init_matrix_rotate(matrix, size, 1, 1);
    util::print_matrix(matrix);
    util::init_matrix_rotate(matrix, size, -1, 1);
    util::print_matrix(matrix);
    util::init_matrix_identity(matrix, size, 2);
    util::print_matrix(matrix);
    init_matrix_identity_rnd(matrix, size);
    util::print_matrix(matrix);
}

void test_matconv(){
    uint64_t size=3;
    vector<vector<int64_t>> matrix(size, vector<int64_t>(size));
    util::init_matrix_rotate(matrix, size, 0, 1);
    util::init_matrix_rotate(matrix, size, -1, 2);
    util::init_matrix_rotate(matrix, size, -2, 3);
    util::print_matrix(matrix);
}





void test_print_iter(){
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW, true);
    uint64_t array_size = 120;
    // how many "c"?
    uint64_t poly_count = 4;
    // polynomial degree
    uint64_t mod_degree = 10;
    // rns count
    uint64_t coeff_mod_size = 3;

    vector<std::uint64_t> arr(array_size);

    for(uint64_t i = 0;i < array_size;i++){
        arr[i] = i+1;
    }

    SEAL_ALLOCATE_GET_POLY_ITER(poly_iter, poly_count, mod_degree, coeff_mod_size, pool_);
    util::set_poly_array(arr.data(), poly_count, mod_degree, coeff_mod_size, poly_iter);
    SEAL_ITERATE(poly_iter, poly_count, [&](auto I){
            util::print_iter(I, coeff_mod_size, mod_degree);
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
    print_iter(rns_iter, coeff_mod_size, coeff_degree);
    SEAL_ALLOCATE_ZERO_GET_RNS_ITER(result, coeff_degree, coeff_mod_size, pool_);
    util::matrix_dot_vector(matrix, coeff_mod_size, coeff_degree, rns_iter, mod_chain, result);
    cout << "print converted rns iter " << endl;
    print_iter(result, coeff_mod_size, coeff_degree);
}

void test_matrix_dot_product(){
    uint64_t coeff_degree = 3;
    uint64_t modulus = 10;
    vector<vector<int64_t>> matrix(coeff_degree, vector<int64_t>(coeff_degree));
    vector<vector<int64_t>> matrix2(coeff_degree, vector<int64_t>(coeff_degree));
    vector<vector<int64_t>> result(coeff_degree, vector<int64_t>(coeff_degree));
    util::init_matrix_rand_mod(matrix, coeff_degree, modulus);
    cout << "first matrix: " << endl;
    util::print_matrix(matrix);
    util::init_matrix_identity(matrix2, coeff_degree, 2);
    cout << "second matrix: " << endl;
    util::print_matrix(matrix2);
    util::matrix_dot_product_mod(matrix, matrix2, result, modulus);
    util::print_matrix(result);
}

void test_matrix_conversion(){
    print_example_banner("matrix_conversion");

    // parameter setting
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = 4096;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    //uint64_t modulo_cipher = 0xfffffffd8001;
    //uint64_t modulo_cipher = 0xfffffffd8001;
    //0x3ffffffff040001
    //vector<Modulus> mod_chain = {Modulus(modulo_cipher)};
    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 7;
    parms.set_plain_modulus(plaintext_modulus);
    auto context = SEALContext::Create(parms);
    print_line(__LINE__);
    cout << "Set encryption parameters and print" << endl;
    print_parameters(context);
    cout << "Parameter validation: " << context->parameter_error_message() << endl;

    cout << endl;

    // generate encryption helper
    cout << "keygen step" << endl;
    KeyGenerator keygen(context);
    cout << "pubkey " << endl;
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // generate plaintext x
    print_line(__LINE__);
    string x = "6x^2 + 1x^1 + 2";
    cout << "Input plaintext: ";
    Plaintext x_plain(x);
    cout << "Express x = " + x + " as a plaintext polynomial " + x_plain.to_string() + "." << endl;
    cout << "Coeff count: " << x_plain.coeff_count() << endl;
    cout << "is NTT form? " << x_plain.is_ntt_form() << endl;
    print_plain(x_plain, 10);

    // generate transform matrix
    vector<vector<int64_t>> matrix(poly_modulus_degree, vector<int64_t>(poly_modulus_degree));
    init_matrix_identity_rnd(matrix, poly_modulus_degree);
    //print_matrix(matrix);

    // convert plaintext by matrix
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    convert(x_plain, copied_plain, matrix, poly_modulus_degree, plaintext_modulus);
    cout << "Copied and converted plaintext" << endl;
    //print_plain_coefficients(copied_plain);

    // encrypt x
    Ciphertext x_encrypted;
    cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;
    cout << "is NTT form? " << x_encrypted.is_ntt_form() << endl;
    print_cipher(x_encrypted, 5);
    cout << "4095-th coeff is " << x_encrypted[4095] << endl;

    // convert ciphertext by matrix
    cout << "----convert ciphertext by matrix---" << endl;
    Ciphertext x_enc_converted = Ciphertext(x_encrypted);
    convert(x_encrypted, x_enc_converted, matrix, poly_modulus_degree, mod_chain);
    print_cipher(x_enc_converted, 5);
    cout << "ciphertext info---->" << endl;
    print_cipher_info(x_enc_converted, decryptor);

    // multiply ciphertext by constant using seal function
    cout << "generate constant multiplied ciphertext by using seal" << endl;
    Ciphertext x_constant;
    Plaintext constant("2");
    evaluator.multiply_plain(x_encrypted, constant, x_constant);
    // constant addition
    cout << "add const 1" << endl;
    Plaintext addconst("1");
    Ciphertext addconst_enc;
    encryptor.encrypt(addconst, addconst_enc);
    print_cipher(addconst_enc, 5);
    evaluator.add_plain_inplace(x_constant, addconst);
    print_cipher(x_constant, 5);
    print_cipher_info(x_constant, decryptor);

    // compare cipher
    //cout << "Compare: " << compare_cipher_all(x_enc_converted, x_constant, poly_modulus_degree) << endl;

    // decryption normal x
    print_line(__LINE__);
    Plaintext x_decrypted;
    cout << "    + decryption of x_encrypted: ";
    decryptor.decrypt(x_encrypted, x_decrypted);
    cout << x_decrypted.to_string() << " ...... Correct?" << endl;

    // decrypt constant multiplied x
    cout << "seal generated constant multiplied x" << endl;
    Plaintext x_constant_dec;
    decryptor.decrypt(x_constant, x_constant_dec);
    print_plain(x_constant_dec, 10);

    // decryption converted x
    Plaintext x_converted_decrypted;
    cout << "decryption of x_copied: " << endl;;
    decryptor.decrypt(x_enc_converted, x_converted_decrypted);

    // compare converted plain and decryption of x_converted
    cout << "Converted plain: " << endl;
    print_plain(copied_plain, 10);
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_converted_decrypted, 20);
}

void example_kazuma(){
    //matrix_conversion();
    //test_conversion();
    //test_matconv();
    //iterator_toy();
    //test_innerprod_vector();
    //test_matrix_conversion_with_rnsiter();
    //test_init_matrix();
    //test_matrix_dot_product();
    test_print_iter();
}
