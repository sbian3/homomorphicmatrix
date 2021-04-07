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
void bench_packed_conv(uint64_t input_dim, uint64_t kernel_dim, uint64_t pack_num, uint64_t poly_modulus_degree){
    print_example_banner("Homomorphic packed convolution Benchmark");

    // sample input and kernel
    Modulus sample_mod(7);
    vector<vector<uint64_t>> input = sample_rn(pack_num, input_dim, sample_mod);
    vector<vector<uint64_t>> kernel = sample_rn(pack_num, kernel_dim, sample_mod);
    //print_input_kernel(input, kernel);


    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    //size_t poly_modulus_degree;
    //cout << "poly_modulus_degree: ";
    //cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    //vector<Modulus> mod_chain =  {Modulus(0xffffff00000001)};
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 7;
    parms.set_plain_modulus(plaintext_modulus);
    cout << "Plaintext modulus: " << plaintext_modulus << endl;
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

    // generate plaintext x
    uint64_t block_size = get_blocksize(input_dim, kernel_dim, 0);
    if(block_size * pack_num > poly_modulus_degree){
        throw invalid_argument("polynomial degree is too small");
    }
    vector<uint64_t> packed_input = pack_input(input, block_size, poly_modulus_degree);
    //cout << "Input plaintext: ";
    Plaintext x_plain(packed_input);
    //cout << "Express x = as a plaintext polynomial " + x_plain.to_string() + "." << endl;
    //cout << "Coeff count: " << x_plain.coeff_count() << endl;
    //print_plain(x_plain, block_size * pack_num);

    // pack kernels
    vector<KernelInfo> kernelinfos = pack_kernel(kernel, block_size, plaintext_modulus);

    // generate transform matrix
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    pack_kernel_to_matrix(kernelinfos, matrix);

    // encrypt x
    Ciphertext x_encrypted;
    cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    //cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    //uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    //cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // lt
    Ciphertext x_enc_lin(x_encrypted);
    util::set_zero_poly(poly_modulus_degree, x_encrypted.coeff_modulus_size(), x_enc_lin.data());
    vector<vector<uint64_t>> matrix_conved(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    cout << "decryption of x_encrypted: ";
    auto lt_start = chrono::high_resolution_clock::now();
    make_packedconv_matrixproduct(kernelinfos, x_encrypted, poly_modulus_degree, matrix_conved, parms.coeff_modulus()[0]);
    decryptor.linear_trans(x_encrypted, matrix, x_enc_lin);
    auto lt_end = chrono::high_resolution_clock::now();

    // decrypt
    Plaintext x_decrypted;
    auto dec_start = chrono::high_resolution_clock::now();
    decryptor.decrypt_bfv_lt(x_enc_lin, matrix_conved, x_decrypted);
    auto dec_end = chrono::high_resolution_clock::now();

    // time result
    auto lt_diff = chrono::duration_cast<chrono::milliseconds>(lt_end - lt_start);
    auto dec_diff = chrono::duration_cast<chrono::milliseconds>(dec_end - dec_start);
    cout << "Linear transformation: " << lt_diff.count() << "ms" << endl;
    cout << "Decryption: " << dec_diff.count() << "ms" << endl;

    // plaintext check
    cout << "decryption of x_tranformed: " << endl;
    //print_plain(x_decrypted, block_size * pack_num);
}



/////////////////
// tests
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
    print_plain(plain, coeff.size());
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
    //bench_packed_conv(input_dim, kernel_dim, pack_num, poly_degree);
    test_kernelinfo_to_rowpair();
}
