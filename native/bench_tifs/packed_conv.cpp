#include "bench.h"
#include "seal/util/packedconv.h"

void make_packedconv_matrixproduct(vector<KernelInfo> kernel_infos, Ciphertext &encrypted, uint64_t poly_degree, vector<vector<uint64_t>> &result, Modulus modulus){
    PolyIter cipher_poly(encrypted);
    cipher_poly++;
    util::matrix_dot_matrix_toeplitz_mod(kernel_infos, **cipher_poly, poly_degree, result, modulus);
}


void print_input_kernel(vector<vector<uint64_t>> input, vector<vector<uint64_t>> kernel){
    cout << "inputs: " << endl;
    util::print_matrix(input);
    cout << "kernels: " << endl;
    util::print_matrix(kernel);
}

// Benchmark
void bench_packed_conv(vector<vector<uint64_t>> input, vector<vector<uint64_t>> kernel, uint64_t pack_num, uint64_t poly_modulus_degree, vector<uint64_t> &decrypted){
    print_example_banner("Homomorphic packed convolution Benchmark");
    bool print_data = true;

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

    // generate plaintext x
    uint64_t block_size = get_blocksize(input[0].size(), kernel[0].size(), 0);
    if(block_size * pack_num > poly_modulus_degree){
        throw invalid_argument("polynomial degree is too small");
    }
    vector<uint64_t> packed_input = pack_input(input, block_size, poly_modulus_degree);
    Plaintext x_plain(packed_input);
    //print_plain(x_plain, block_size * pack_num);

    // pack kernels
    vector<KernelInfo> kernelinfos = pack_kernel(kernel, block_size, parms.coeff_modulus()[0]);

    // generate transform matrix
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    pack_kernel_to_matrix(kernelinfos, matrix);

    // encrypt x
    Ciphertext x_encrypted;
    cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // lt
    Ciphertext x_enc_lin(x_encrypted);
    util::set_zero_poly(poly_modulus_degree, x_encrypted.coeff_modulus_size(), x_enc_lin.data());
    vector<vector<uint64_t>> matrix_conved(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    //cout << "decryption of x_encrypted: ";
    auto lt_start = chrono::high_resolution_clock::now();
    make_packedconv_matrixproduct(kernelinfos, x_encrypted, poly_modulus_degree, matrix_conved, parms.coeff_modulus()[0]);
    decryptor.lt_packedconv(x_encrypted, kernelinfos, x_enc_lin);
    auto lt_end = chrono::high_resolution_clock::now();

    // decrypt
    Plaintext x_decrypted;
    decryptor.decrypt_bfv_lt(x_enc_lin, matrix_conved, x_decrypted);
    auto dec_start = chrono::high_resolution_clock::now();
    decryptor.decrypt_bfv_lt(x_enc_lin, matrix_conved, x_decrypted);
    auto dec_end = chrono::high_resolution_clock::now();

    // time result
    auto lt_diff = chrono::duration_cast<chrono::milliseconds>(lt_end - lt_start);
    auto dec_diff = chrono::duration_cast<chrono::milliseconds>(dec_end - dec_start);
    cout << TIME_LABEL_LT << lt_diff.count() << "ms" << endl;
    cout << TIME_LABEL_DEC << dec_diff.count() << "ms" << endl;

    // plaintext check
    cout << "----Decryption---- " << endl;
    if(print_data){
        print_plain(x_decrypted, block_size * pack_num);
    }

    // put decrypted result to vector
    decrypted.assign(x_decrypted.data(), x_decrypted.data() + block_size * pack_num);
}

void test_packedconv(){
    cout << "packedconv test" << endl;
    uint64_t pack_num = 2;
    uint64_t poly_degree = 1024;
    vector<vector<uint64_t>> input = { {1, 4, 2}, {5, 1, 3} };
    vector<vector<uint64_t>> kernel = { {3, 2, 1}, {3, 2, 5} };
    vector<uint64_t> decrypted(10);
    bench_packed_conv(input, kernel, pack_num, poly_degree , decrypted);

    // decrypted shold be [3, 0, 1, 1, 2, 1, 6, 1, 4, 1]
    vector<uint64_t> expect = {3, 0, 1, 1, 2, 1, 6, 1, 4, 1};
    for(uint64_t i = 0;i < decrypted.size();i++){
        if(expect[i] != decrypted[i]){
            cerr << "test failed: " << i << "th number" << endl;
            cerr << "expected: " << expect[i] << endl;
            cerr << "actual: " << decrypted[i] << endl;
            return;
        }
    }
    cout << "test passed!!" << endl;
}


//////////////////
// main
/////////////////
int main(int argc, char* argv[]){
    uint64_t input_dim, kernel_dim, pack_num, poly_degree;
    if(argc != 5){
        cerr << "please input two numbers.argc: " << argc << endl;
        return 1;
    }
    input_dim = stoull(argv[1]);
    kernel_dim = stoull(argv[2]);
    poly_degree = stoull(argv[3]);
    pack_num = stoull(argv[4]);

    uint64_t output_dim = input_dim + kernel_dim - 1;
    // sample input and kernel data
    Modulus sample_mod(7);
    vector<vector<uint64_t>> input = sample_rn(pack_num, input_dim, sample_mod);
    vector<vector<uint64_t>> kernel = sample_rn(pack_num, kernel_dim, sample_mod);
    vector<uint64_t> decrypted(output_dim);

    // benchmark impl
    //bench_packed_conv(input, kernel, pack_num, poly_degree, decrypted);
    test_packedconv();
}
