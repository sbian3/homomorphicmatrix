#include "testprimitives.h"
#include "testutil.h"
#include "seal/util/uintlinarith.h"

using namespace seal::util;

void test_scalars_to_diagonallist(){
    vector<uint64_t> scalars = { 3, 2, 1, 5, 3, 8, 3, 4, 7, 1 };
    uint64_t colsize = 5;
    uint64_t rowsize = 6;
    vector<vector<uint64_t>> diagonallists = util::scalars_to_diagonallist(scalars, colsize, rowsize);
    util::print_matrix(diagonallists);
}

inline vector<uint64_t> sample_rn(uint64_t size, Modulus modulus){
    vector<uint64_t> ret(size);
    for(uint64_t i = 0;i < size;i++){
        ret[i] = rand() % modulus.value();
    }
    return ret;
}

void print_plain(Plaintext plain, uint64_t num){
    if(num > plain.coeff_count()){
        std::cout << "print_plain: warn: num is bigger than plain coefficient size!" << std::endl;
        num = plain.coeff_count();
    }
    std::cout << "Plaintext Coefficients(first " << num << " elements): [";
    for(auto i = 0U;i < num;i++){
        std::cout << plain[i] << " ";
    }
    std::cout << "]" << std::endl;
}

void set_zero(util::RNSIter rns, uint64_t rns_count){
    SEAL_ITERATE(rns, rns_count, [&](auto I){
            
            });

}

void set_zero_cipher(Ciphertext &encrypted,uint64_t poly_modulus_degree, uint64_t rns_count, uint64_t c_size, uint64_t input_dim){
    util::CoeffIter coeff_iter(encrypted.data(0));
    coeff_iter = coeff_iter + (input_dim);
    util::set_zero_uint(poly_modulus_degree - input_dim, coeff_iter);
}

void test_delete_cipher(uint64_t input_dim, uint64_t kernel_dim, uint64_t poly_modulus_degree){
    print_example_banner("Direct convolution of ciphertext Benchmark");

    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);

    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    //vector<Modulus> mod_chain =  {Modulus(0xffffff00000001)};
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 1032193;
    parms.set_plain_modulus(plaintext_modulus);
    SEALContext context(parms);
    print_line(__LINE__);
    cout << "Set encryption parameters and print" << endl;
    print_parameters(context);
    cout << "Parameter validation: " << context.parameter_error_message() << endl;

    cout << endl;

    // generate encryption helper
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // sample plaintext x
    //cout << "Input plaintext: ";
    Plaintext x_plain(sample_rn(input_dim, Modulus(9)));
    //cout << "plaintext polynomial " + x_plain.to_string() + "." << endl;
    //cout << "Coeff count: " << x_plain.coeff_count() << endl;
    print_plain(x_plain, 10);

    // encrypt x
    Ciphertext x_encrypted;
    //cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    //uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    //cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;
    set_zero_cipher(x_encrypted, poly_modulus_degree, 1, 2, input_dim);
    util::PolyIter x_iter(x_encrypted);
    util::print_iter(x_iter, 2);

    // convolve encrypted x
    cout << "decryption of x_tranformed: " << endl;
    Plaintext x_decrypted;
    decryptor.decrypt(x_encrypted, x_decrypted);

    print_plain(x_decrypted, 20);
}

void test_kernel_matrix_dot_vector(){
    vector<vector<uint64_t>> kernels = { {1 , 0, 3}, {4, 2, 1} };
    uint64_t block_size = 10;
    uint64_t matrix_size = 300;
    Modulus modulus(1023);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, block_size, modulus);
    vector<vector<uint64_t>> kernel_matrix(matrix_size, vector<uint64_t>(matrix_size));
    pack_kernel_to_matrix(kernel_infos, kernel_matrix);
    cout << "kernel matrix: " << endl;
    print_matrix(kernel_matrix);
    MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);
    SEAL_ALLOCATE_GET_COEFF_ITER(vec_right, matrix_size, pool_);
    SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(result_actual, matrix_size, pool_);
    vector<uint64_t> result_expected(matrix_size);

    // right vector generation
    for(uint64_t i = 0;i < matrix_size;i++){
        vec_right[i] = rand() % modulus.value();
    }
    cout << "vector (right)" << endl;
    print_iter(vec_right, matrix_size);

    // expect result
    util::matrix_dot_vector(kernel_matrix, vec_right, modulus, matrix_size, result_expected);
    cout << "expect" << endl;
    print_iter(result_expected, matrix_size);

    // actual result
    util::kernel_matrix_dot_vector(kernel_infos, vec_right, modulus, result_actual);
    cout << "actual" << endl;
    print_iter(result_actual, matrix_size);
}

void test_print_vector(){
    vector<uint64_t> vec = {1, 3, 2, 4, 5, 2, 10, 2};
    util::print_vector(vec, 3);
    util::print_vector(vec, 30);
}

// 行列操作など基礎的な関数のテストや実験のためのもの
int main(){
    cout << "test_primitives" << endl;
    //test_scalars_to_diagonallist();
    //test_delete_cipher(5, 0, 2048);
    //test_kernel_matrix_dot_vector();
    test_print_vector();
}
