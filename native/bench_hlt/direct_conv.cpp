#include "bench.h"

using namespace std;
using namespace seal;


void conv_cipher_direct(vector<uint64_t> input, vector<uint64_t> kernel, uint64_t poly_modulus_degree, vector<uint64_t> &decrypted, int64_t &latency_lt, int64_t &latency_dec){

    bool print_data = false;
    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    uint64_t output_dim = input.size() + kernel.size() - 1;
    if(output_dim > poly_modulus_degree){
        throw invalid_argument("polynomial degree is too small");
    }

    //vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    vector<Modulus> mod_chain =  select_modchain(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 1024;
    parms.set_plain_modulus(plaintext_modulus);
    SEALContext context(parms);
    if(print_data){
        print_parameters(context);
        cout << "Parameter validation: " << context.parameter_error_message() << endl;
        cout << endl;
    }

    // generate encryption helper
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // plaintext x
    Plaintext x_plain(input);
    if(print_data){
        print_plain(x_plain, input.size());
    }

    // encrypt x
    Ciphertext x_encrypted;
    encryptor.encrypt(x_plain, x_encrypted);
    if(print_data){
        cout << "----Encrypt x_plain to x_encrypted.----" << endl;
        cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
        cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;
    }

    // convolve encrypted x
    if(print_data){
        cout << "print kernel Coefficients: ";
        print_iter(kernel, kernel.size());
    }
    Ciphertext conved_x(x_encrypted);
    util::set_zero_uint(conved_x.size() * conved_x.poly_modulus_degree() * conved_x.coeff_modulus_size(), conved_x.data());
    auto lt_start = chrono::high_resolution_clock::now();
    util::conv_negacyclic(kernel, x_encrypted, mod_chain, conved_x);
    auto lt_end = chrono::high_resolution_clock::now();

    Plaintext x_conved_decrypted;
    auto dec_start = chrono::high_resolution_clock::now();
    decryptor.decrypt(conved_x, x_conved_decrypted);
    auto dec_end = chrono::high_resolution_clock::now();

    // time result
    auto lt_diff = chrono::duration_cast<chrono::microseconds>(lt_end - lt_start);
    auto dec_diff = chrono::duration_cast<chrono::microseconds>(dec_end - dec_start);
    latency_lt = lt_diff.count();
    latency_dec = dec_diff.count();

    if(print_data){
        cout << TIME_LABEL_LT << lt_diff.count() << US << endl;
        cout << TIME_LABEL_DEC << dec_diff.count() << US << endl;
        print_plain(x_conved_decrypted, input.size() + kernel.size() - 1);
    }
    decrypted.assign(x_conved_decrypted.data(), x_conved_decrypted.data() + output_dim);
}

void test_directconv(){
    // assert plaintext_modulus is 1024
    uint64_t poly_degree = 1024;
    vector<uint64_t> input = {3, 1, 6, 1, 3, 0, 6, 1, 2, 2};
    vector<uint64_t> kernel = {3, 4, 3, 4};
    int64_t latency_lt, latency_dec;
    vector<uint64_t> decrypted(input.size() + kernel.size() - 1);
    conv_cipher_direct(input, kernel, poly_degree, decrypted, latency_lt, latency_dec);
    vector<uint64_t> expect = { 9, 15, 31, 42, 35, 39, 31, 39, 28, 41, 18, 14, 8 };
    for(uint64_t i = 0;i < expect.size();i++){
        if(expect[i] != decrypted[i]){
            cerr << "test failed: " << i << "th number" << endl;
            cerr << "expected: " << expect[i] << endl;
            cerr << "actual: " << decrypted[i] << endl;
            return;
        }
    }
    cout << "test passed!!" << endl;
}

int main(int argc, char* argv[]){
    uint64_t bench_times = 10;
    uint64_t input_dim, kernel_dim, poly_degree;
    int ret = check_args(argc, argv, input_dim, kernel_dim, poly_degree);
    if(ret == 0){
        return 1;
    }
    //test_directconv();
    print_example_banner("Direct convolution of ciphertext Benchmark");
    cout << "input_dim: " << input_dim << endl;
    cout << "kernel_dim: " << kernel_dim << endl;
    cout << "poly_degree: " << poly_degree << endl;
    cout << "benchmark for " << bench_times << " times" << endl;

    vector<uint64_t> decrypted(input_dim + kernel_dim - 1);
    uint64_t latency_lt_sum = 0;
    uint64_t latency_dec_sum = 0;

    for(uint64_t i = 0;i < bench_times;i++){
        Modulus sample_mod(7);
        vector<uint64_t> input = sample_rn(input_dim, sample_mod);
        vector<uint64_t> kernel = sample_rn(kernel_dim, sample_mod);
        int64_t latency_lt, latency_dec;
        conv_cipher_direct(input, kernel, poly_degree, decrypted, latency_lt, latency_dec);
        latency_lt_sum += static_cast<uint64_t>(latency_lt);
        latency_dec_sum += static_cast<uint64_t>(latency_dec);
    }
    double latency_lt = latency_lt_sum/static_cast<double>(bench_times);
    double latency_dec = latency_dec_sum/static_cast<double>(bench_times);
    cout << TIME_LABEL_LT << latency_lt << US << endl;
    cout << TIME_LABEL_DEC << latency_dec << US << endl;
}
