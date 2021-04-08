#include "bench.h"

using namespace std;
using namespace seal;


void test_conv_cipher_direct(uint64_t input_dim, uint64_t kernel_dim, uint64_t poly_modulus_degree){
    print_example_banner("Direct convolution of ciphertext Benchmark");

    bool print_arr = true;
    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    //size_t poly_modulus_degree;
    //cout << "poly_modulus_degree: ";
    //cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    cout << "input_dim: " << input_dim << endl;
    cout << "kernel_dim: " << kernel_dim << endl;

    //vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    vector<Modulus> mod_chain =  select_modchain(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 1032193;
    parms.set_plain_modulus(plaintext_modulus);
    SEALContext context(parms);
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
    Plaintext x_plain(sample_rn(input_dim, Modulus(7)));
    if(print_arr){
        print_plain(x_plain, input_dim);
    }
    //cout << "plaintext polynomial " + x_plain.to_string() + "." << endl;
    //cout << "Coeff count: " << x_plain.coeff_count() << endl;
    //print_plain(x_plain, 10);

    // encrypt x
    Ciphertext x_encrypted;
    //cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    //uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    //cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // convolve encrypted x
    vector<uint64_t> kernel = sample_rn(kernel_dim, Modulus(7));
    if(print_arr){
        print_iter(kernel, kernel_dim);
    }
    Ciphertext conved_x(x_encrypted);
    util::set_zero_uint(conved_x.size() * conved_x.poly_modulus_degree() * conved_x.coeff_modulus_size(), conved_x.data());
    auto lt_start = chrono::high_resolution_clock::now();
    util::conv_negacyclic(kernel, x_encrypted, mod_chain, conved_x);
    auto lt_end = chrono::high_resolution_clock::now();

    cout << "decryption of x_tranformed: " << endl;
    Plaintext x_conved_decrypted;
    auto dec_start = chrono::high_resolution_clock::now();
    decryptor.decrypt(conved_x, x_conved_decrypted);
    auto dec_end = chrono::high_resolution_clock::now();

    // time result
    auto lt_diff = chrono::duration_cast<chrono::microseconds>(lt_end - lt_start);
    auto dec_diff = chrono::duration_cast<chrono::microseconds>(dec_end - dec_start);
    cout << "Linear transformation: " << lt_diff.count() << "us" << endl;
    cout << "Decryption: " << dec_diff.count() << "us" << endl;

    if(print_arr){
        print_plain(x_conved_decrypted, input_dim + kernel_dim - 1);
    }
}

int main(int argc, char* argv[]){
    uint64_t input_dim, kernel_dim, poly_degree;
    int ret = check_args(argc, argv, input_dim, kernel_dim, poly_degree);
    if(ret == 0){
        return 1;
    }
    test_conv_cipher_direct(input_dim, kernel_dim, poly_degree);
}
