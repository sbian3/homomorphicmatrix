#include "bench.h"

void print_input_kernel(vector<vector<uint64_t>> input, vector<vector<uint64_t>> kernel){
    cout << "inputs: " << endl;
    util::print_matrix(input);
    cout << "kernels: " << endl;
    util::print_matrix(kernel);
}

// Benchmark
void bench_packed_conv(vector<vector<uint64_t>> input, vector<vector<uint64_t>> kernel, uint64_t pack_num, uint64_t poly_modulus_degree, vector<uint64_t> &decrypted, int64_t &latency_lt, int64_t &latency_dec, bool print_data){
    if(print_data){
        print_input_kernel(input, kernel);
    }

    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    vector<Modulus> mod_chain =  select_modchain(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 7;
    parms.set_plain_modulus(plaintext_modulus);
    SEALContext context(parms);
    if(print_data){
        print_parameters(context);
        cout << "Set encryption parameters and print" << endl;
        cout << "Parameter validation: " << context.parameter_error_message() << endl;
        cout << endl;
    }

    // generate encryption helper
    auto keygen_start = chrono::high_resolution_clock::now();
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor_LT decryptor(context, secret_key);
    auto keygen_end = chrono::high_resolution_clock::now();

    // pack kernels
    vector<KernelInfo> kernelinfos = pack_kernel(kernel, input, parms.coeff_modulus()[0], poly_modulus_degree);
    vector<uint64_t> packed_input = pack_input(input, kernelinfos, poly_modulus_degree);
    Plaintext x_plain(packed_input);
    uint64_t result_len_packed = kernelinfos.back().get_startcol() + kernelinfos.back().get_colsize();


    // generate transform matrix
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    pack_kernel_to_matrix(kernelinfos, matrix);

    // encrypt x
    auto enc_start = chrono::high_resolution_clock::now();
    Ciphertext x_encrypted;
    encryptor.encrypt(x_plain, x_encrypted);
    auto enc_end = chrono::high_resolution_clock::now();
    if(print_data){
        cout << "----Encrypt x_plain to x_encrypted.----" << endl;
        //cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;
    }

    // lt
    Ciphertext x_enc_lin(x_encrypted);
    util::set_zero_poly(poly_modulus_degree, x_encrypted.coeff_modulus_size(), x_enc_lin.data());
    vector<vector<uint64_t>> matrix_conved(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    //cout << "decryption of x_encrypted: ";
    auto lt_start = chrono::high_resolution_clock::now();
    make_packedconv_matrixproduct(kernelinfos, x_encrypted, poly_modulus_degree, matrix_conved, parms.coeff_modulus()[0]);
    auto lt_half = chrono::high_resolution_clock::now();
    decryptor.lt_packedconv(x_encrypted, kernelinfos, x_enc_lin);
    auto lt_end = chrono::high_resolution_clock::now();

    // decrypt
    Plaintext x_decrypted;
    auto dec_start = chrono::high_resolution_clock::now();
    decryptor.decrypt_bfv_lt_toeplitz(kernelinfos, x_enc_lin, matrix_conved, poly_modulus_degree, x_decrypted);
    auto dec_end = chrono::high_resolution_clock::now();

    // time result
    auto keygen_diff = chrono::duration_cast<chrono::microseconds>(keygen_end - keygen_start);
    auto enc_diff = chrono::duration_cast<chrono::microseconds>(enc_end - enc_start);
    auto lt_diff = chrono::duration_cast<chrono::microseconds>(lt_end - lt_start);
    auto lt_c0 = chrono::duration_cast<chrono::microseconds>(lt_end - lt_half);
    auto dec_diff = chrono::duration_cast<chrono::microseconds>(dec_end - dec_start);
    latency_lt = lt_diff.count();
    latency_dec = dec_diff.count();
    cout << "keygen: " << keygen_diff.count() << US << endl;
    cout << "encrypt: " << enc_diff.count() << US << endl;
    cout << "c0 matrix_vector product: " << lt_c0.count() << US << endl;
    cout << "Decryption: " << latency_dec << US << endl;
    if(print_data){
        cout << TIME_LABEL_LT << lt_diff.count() << US << endl;
        cout << TIME_LABEL_DEC << dec_diff.count() << US << endl;
    }

    // plaintext check
    if(print_data){
        cout << "----Decryption---- " << endl;
        print_plain(x_decrypted, result_len_packed);
    }

    // put decrypted result to vector
    decrypted.assign(x_decrypted.data(), x_decrypted.data() + result_len_packed);
}

bool pass_test_packedconv(){
    cout << "packedconv test" << endl;
    uint64_t pack_num = 2;
    uint64_t poly_degree = 1024;
    vector<vector<uint64_t>> input = { {1, 4, 2}, {5, 1, 3}};
    vector<vector<uint64_t>> kernel = { {3, 2, 1}, {3, 2, 5}};
    vector<uint64_t> decrypted(10);
    int64_t time_lt, time_dec;
    bench_packed_conv(input, kernel, pack_num, poly_degree , decrypted, time_lt, time_dec, false);

    // decrypted shold be     [3, 0, 1, 1, 2, 1, 6, 1, 4, 1]
    vector<uint64_t> expect = {3, 0, 1, 1, 2, 1, 6, 1, 4, 1};
    for(uint64_t i = 0;i < expect.size();i++){
        if(expect[i] != decrypted[i]){
            print_iter(decrypted, decrypted.size());
            cerr << "test failed: " << i << "th number" << endl;
            cerr << "expected: " << expect[i] << endl;
            cerr << "actual: " << decrypted[i] << endl;
            return false;
        }
    }
    cout << "test passed!!" << endl;
    return true;
}

bool pass_test_packedconv_1pack(){
    cout << "packedconv test" << endl;
    uint64_t pack_num = 1;
    uint64_t poly_degree = 1024;
    vector<vector<uint64_t>> input = { {1, 4, 2}};
    vector<vector<uint64_t>> kernel = { {3, 2, 1}};
    vector<uint64_t> decrypted(10);
    int64_t time_lt, time_dec;
    bench_packed_conv(input, kernel, pack_num, poly_degree , decrypted, time_lt, time_dec, false);

    // decrypted shold be [3, 0, 1, 1, 2]
    vector<uint64_t> expect = {3, 0, 1, 1, 2};
    for(uint64_t i = 0;i < expect.size();i++){
        if(expect[i] != decrypted[i]){
            cerr << "test failed: " << i << "th number" << endl;
            cerr << "expected: " << expect[i] << endl;
            cerr << "actual: " << decrypted[i] << endl;
            return false;
        }
    }
    cout << "test passed!!" << endl;
    return true;
}
//////////////////
// main
/////////////////
int main(int argc, char* argv[]){
    //test_toeplitz_dot_vector();
    uint64_t bench_times = 100;
    if(!pass_test_packedconv()){
        cerr << "test failed!" << endl;
        return 1;
    }
    uint64_t input_dim, kernel_dim, pack_num, poly_degree;
    if(argc != 5){
        cerr << "please input two numbers.argc: " << argc << endl;
        return 1;
    }
    input_dim = stoull(argv[1]);
    kernel_dim = stoull(argv[2]);
    poly_degree = stoull(argv[3]);
    pack_num = stoull(argv[4]);

    print_example_banner("Homomorphic packed convolution Benchmark");
    // sample input and kernel data

    uint64_t output_dim = input_dim + kernel_dim - 1;
    vector<uint64_t> decrypted(output_dim);

    uint64_t latency_lt_sum = 0;
    uint64_t latency_dec_sum = 0;

    // benchmark impl
    for(uint64_t i = 0;i < bench_times;i++){
        Modulus sample_mod(7);
        vector<vector<uint64_t>> input = sample_rn(pack_num, input_dim, sample_mod);
        vector<vector<uint64_t>> kernel = sample_rn(pack_num, kernel_dim, sample_mod);
        int64_t latency_lt, latency_dec;
        bench_packed_conv(input, kernel, pack_num, poly_degree, decrypted, latency_lt, latency_dec, false);
        latency_lt_sum += static_cast<uint64_t>(latency_lt);
        latency_dec_sum += static_cast<uint64_t>(latency_dec);
    }

    double latency_lt = latency_lt_sum/static_cast<double>(bench_times);
    double latency_dec = latency_dec_sum/static_cast<double>(bench_times);
    cout << TIME_LABEL_LT << latency_lt << US << endl;
    cout << TIME_LABEL_DEC << latency_dec << US << endl;
}
