#include "bench.h"
#include "seal/util/packedconv.h"

void make_circ_matrix(Ciphertext &encrypted, uint64_t poly_modulus_degree, Modulus modulus, vector<vector<uint64_t>> &circ){
    PolyIter polyiter(encrypted);
    polyiter++;
    util::init_matrix_with_coeff(circ, poly_modulus_degree, **polyiter, modulus);
}


void swap_cipher(Ciphertext &encrypted, uint64_t i , uint64_t j){
    PolyIter polyiter(encrypted);
    CoeffIter coeffiter = **polyiter;
    uint64_t tmp = coeffiter[i];
    coeffiter[i] = coeffiter[j];
    coeffiter[j] = tmp;
}

void swap_matrix(vector<vector<uint64_t>> &matrix, uint64_t i , uint64_t j){
    vector<uint64_t> tmp = matrix[i];
    matrix[i] = matrix[j];
    matrix[j] = tmp;
}

// Benchmark
void bench_swap_slot(vector<uint64_t> &input, uint64_t poly_modulus_degree, vector<uint64_t> &decrypted, int64_t &latency_lt, int64_t &latency_dec, bool print_data){

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
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // generate plaintext x
    Plaintext x_plain(input);


    // encrypt x
    Ciphertext x_encrypted;
    encryptor.encrypt(x_plain, x_encrypted);
    if(print_data){
        cout << "----Encrypt x_plain to x_encrypted.----" << endl;
        cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;
    }

    // generate transform matrix
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    make_circ_matrix(x_encrypted, poly_modulus_degree, mod_chain[0] , matrix);

    // lt
    //Ciphertext x_enc_lin(x_encrypted);
    //util::set_zero_poly(poly_modulus_degree, x_encrypted.coeff_modulus_size(), x_enc_lin.data());
    //vector<vector<uint64_t>> matrix_conved(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    //cout << "decryption of x_encrypted: ";
    uint64_t swap_i = 0;
    uint64_t swap_j = 2;
    auto lt_start = chrono::high_resolution_clock::now();
    swap_cipher(x_encrypted, swap_i, swap_j);
    auto lt_half = chrono::high_resolution_clock::now();
    swap_matrix(matrix, swap_i, swap_j);
    auto lt_end = chrono::high_resolution_clock::now();

    // decrypt
    Plaintext x_decrypted;
    auto dec_start = chrono::high_resolution_clock::now();
    decryptor.decrypt_bfv_lt(x_encrypted, matrix, poly_modulus_degree, x_decrypted);
    auto dec_end = chrono::high_resolution_clock::now();

    // time result
    auto lt_diff = chrono::duration_cast<chrono::microseconds>(lt_end - lt_start);
    auto lt_c0 = chrono::duration_cast<chrono::microseconds>(lt_end - lt_half);
    auto dec_diff = chrono::duration_cast<chrono::microseconds>(dec_end - dec_start);
    latency_lt = lt_diff.count();
    latency_dec = dec_diff.count();
    //cout << "c0 matrix_vector product: " << lt_c0.count() << US << endl;
    if(print_data){
        cout << TIME_LABEL_LT << lt_diff.count() << US << endl;
        cout << TIME_LABEL_DEC << dec_diff.count() << US << endl;
    }

    // plaintext check

    // put decrypted result to vector
    decrypted.assign(x_decrypted.data(), x_decrypted.data() + input.size());
}

bool pass_test_packedconv(){
    cout << "packedconv test" << endl;
    uint64_t poly_degree = 1024;
    vector<uint64_t> input = {1, 4, 2, 5};
    vector<uint64_t> decrypted(10);
    int64_t time_lt, time_dec;
    bench_swap_slot(input,poly_degree , decrypted, time_lt, time_dec, false);
    print_vector(decrypted, 4);

    vector<uint64_t> expect = {2, 4, 1, 5};
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
    uint64_t bench_times = 10;
    if(!pass_test_packedconv()){
        cerr << "test failed!" << endl;
        return 1;
    }
    uint64_t input_dim, poly_degree;
    if(argc != 3){
        cerr << "please input two numbers.argc: " << argc << endl;
        return 1;
    }
    input_dim = stoull(argv[1]);
    poly_degree = stoull(argv[2]);

    print_example_banner("Homomorphic packed convolution Benchmark");
    // sample input and kernel data

    vector<uint64_t> decrypted(input_dim);

    uint64_t latency_lt_sum = 0;
    uint64_t latency_dec_sum = 0;

    // benchmark impl
    for(uint64_t i = 0;i < bench_times;i++){
        Modulus sample_mod(7);
        vector<uint64_t> input = sample_rn(input_dim, sample_mod);
        int64_t latency_lt, latency_dec;
        bench_swap_slot(input, poly_degree, decrypted, latency_lt, latency_dec, false);
        latency_lt_sum += static_cast<uint64_t>(latency_lt);
        latency_dec_sum += static_cast<uint64_t>(latency_dec);
    }

    double latency_lt = latency_lt_sum/static_cast<double>(bench_times);
    double latency_dec = latency_dec_sum/static_cast<double>(bench_times);
    cout << TIME_LABEL_LT << latency_lt << US << endl;
    cout << TIME_LABEL_DEC << latency_dec << US << endl;
}
