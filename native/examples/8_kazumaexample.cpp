#include "examples.h"
#include <cassert>
#include "simd.h"

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
    for(auto i = 0U;i < plain.coeff_count();i++){
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

void test_bfv_matrix(){
    print_example_banner("matrix_conversion");

    // parameter setting
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree;
    cout << "poly_modulus_degree: ";
    cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

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
    print_plain(x_plain, 10);

    // generate transform matrix
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    util::init_matrix_identity(matrix, poly_modulus_degree, 2);
    matrix[0][1] = 1;
    //print_matrix(matrix);

    // convert plaintext by matrix
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    //convert(x_plain, copied_plain, matrix, poly_modulus_degree, plaintext_modulus);
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

    // decryption normal x
    print_line(__LINE__);
    Plaintext x_decrypted;
    cout << "decryption of x_encrypted: ";
    auto time_start = chrono::high_resolution_clock::now();
    decryptor.decrypt_bfv_with_matrix(x_encrypted, x_decrypted, matrix);
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);

    // compare converted plain and decryption of x_converted
    //cout << "Converted plain: " << endl;
    //print_plain(copied_plain, 10);
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_decrypted, 20);
}

void test_lin_conv_packing(){
    print_example_banner("matrix_conversion");

    // parameter setting
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree;
    cout << "poly_modulus_degree: ";
    cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

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
    string x = "3x^514 + 2x^513 + 1x^512 + 6x^2 + 1x^1 + 2";
    cout << "Input plaintext: ";
    Plaintext x_plain(x);
    //x_plain[poly_modulus_degree/2] = 1;
    //x_plain[poly_modulus_degree/2 + 1] = 2;
    //x_plain[poly_modulus_degree/2 + 2] = 3;

    cout << "Express x = " + x + " as a plaintext polynomial " + x_plain.to_string() + "." << endl;
    cout << "Coeff count: " << x_plain.coeff_count() << endl;
    print_plain(x_plain, 10);

    // generate transform matrix using kernel
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    vector<uint64_t> kernel1 = {1,2,3};
    vector<uint64_t> kernel2 = {2,3,4};
    util::init_matrix_rotate_partial(matrix, kernel1, 0, 0, plaintext_modulus);
    util::init_matrix_rotate_partial(matrix, kernel2, poly_modulus_degree/2, poly_modulus_degree/2, plaintext_modulus);

    // encrypt x
    Ciphertext x_encrypted;
    cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // decryption normal x
    print_line(__LINE__);
    Plaintext x_decrypted;
    cout << "decryption of x_encrypted: ";
    auto time_start = chrono::high_resolution_clock::now();
    decryptor.decrypt_bfv_with_matrix(x_encrypted, x_decrypted, matrix);
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);

    // compare converted plain and decryption of x_converted
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_decrypted, 516);
}

void test_bfv_conv(){
    print_example_banner("matrix_conversion");

    // parameter setting
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree;
    cout << "poly_modulus_degree: ";
    cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

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
    print_plain(x_plain, 10);

    // generate transform kernel(filter)
    vector<uint64_t> kernel = {1, 1, 2, 3};
    cout << "kernel(conv filter): ";
    for(auto i = 0U;i < kernel.size();i++){
        cout << kernel[i] << " ";
    }
    cout << endl;

    // convert plaintext by matrix
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    //convert(x_plain, copied_plain, matrix, poly_modulus_degree, plaintext_modulus);
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

    // decryption normal x
    print_line(__LINE__);
    Plaintext x_decrypted;
    cout << "decryption of x_encrypted: ";
    auto time_start = chrono::high_resolution_clock::now();
    decryptor.decrypt_bfv_with_kernel(x_encrypted, x_decrypted, kernel);
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);

    // compare converted plain and decryption of x_converted
    //cout << "Converted plain: " << endl;
    //print_plain(copied_plain, 10);
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_decrypted, 20);
}

void test_ntt_conv(){
    print_example_banner("matrix_conversion");

    // parameter setting
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree;
    cout << "poly_modulus_degree: ";
    cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

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
    print_plain(x_plain, 10);

    // generate transform kernel(filter)
    string kernel_st = "3x^3 + 2x^2 + 1x^1 + 1";
    Plaintext kernel(kernel_st);
    print_plain(kernel, 10);

    // convert plaintext by matrix
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    cout << "Copied and converted plaintext" << endl;

    // encrypt x
    Ciphertext x_encrypted;
    cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // conv plaintext = multiply plain function
    auto time_start = chrono::high_resolution_clock::now();
    evaluator.multiply_plain_inplace(x_encrypted, kernel);
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "NTT Conv: " << time_diff.count() << "us"<<endl;

    // decryption normal x
    print_line(__LINE__);
    Plaintext x_decrypted;
    cout << "decryption of x_tranformed: " << endl;
    decryptor.decrypt(x_encrypted, x_decrypted);
    print_plain(x_decrypted, 20);
}

void test_conv_cipher_direct(){
    print_example_banner("strait convolution of ciphertext");

    // parameter setting
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree;
    cout << "poly_modulus_degree: ";
    cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 1032193;
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
    print_plain(x_plain, 10);

    // convert plaintext by matrix
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    cout << "Copied and converted plaintext" << endl;

    // encrypt x
    Ciphertext x_encrypted;
    cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // convolve encrypted x
    vector<uint64_t> kernel = {1,1,2,3};
    Ciphertext conved_x(x_encrypted);
    for(uint64_t i = 0;i < cipher_coeffsize;i++){
        conved_x[i] = 0;
    }
    auto time_start = chrono::high_resolution_clock::now();
    util::conv_negacyclic(kernel, x_encrypted, mod_chain, conved_x);
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "straight convolution: " << time_diff.count() << "us" << endl;

    cout << "decryption of x_tranformed: " << endl;
    Plaintext x_conved_decrypted;
    decryptor.decrypt(conved_x, x_conved_decrypted);
    print_plain(x_conved_decrypted, 20);
}

void example_kazuma(){
    // matrix_conversion();
    //test_conversion();
    //test_matconv();
    //iterator_toy();
    //test_innerprod_vector();
    //test_matrix_conversion_with_rnsiter();
    //test_init_matrix();
    //test_init_matrix_uint();
    //test_init_matrix_uint_by_kernel();
    //test_matrix_dot_product();
    //test_print_iter();
    //test_bfv_matrix();
    //test_secret_product();
    //test_inverse();
    //test_util_dot_product_mod(); 
    //benchmark_singlefunction();
    //test_batch_encoder();
    //test_batch_matrix();
    //test_conv_nega();
    //test_conv_cipher_direct();
    //test_batch_convolution();
    //test_conviter();
    //test_bfv_conv();
    //test_ntt_conv();
    //test_matrix_init_partial();
    test_lin_conv_packing();
}

int main(){
    cout << "kazuma example!" << endl;
}
