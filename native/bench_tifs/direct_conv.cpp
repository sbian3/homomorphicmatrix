#include "bench.h"

using namespace std;
using namespace seal;


void test_conv_cipher_direct(uint64_t input_dim, uint64_t kernel_dim){
    print_example_banner("Direct convolution of ciphertext Benchmark");

    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree;
    cout << "poly_modulus_degree: ";
    cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
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
    cout << "keygen step" << endl;
    KeyGenerator keygen(context);
    cout << "pubkey " << endl;
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // generate plaintext x
    print_line(__LINE__);
    cout << "Input plaintext: ";
    Plaintext x_plain(sample_rn(input_dim, plaintext_modulus));
    cout << "plaintext polynomial " + x_plain.to_string() + "." << endl;
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
    vector<uint64_t> kernel = sample_rn(kernel_dim, plaintext_modulus);
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

int main(int argc, char* argv[]){
    uint64_t input_dim, kernel_dim;
    int ret = check_args(argc, argv, input_dim, kernel_dim);
    if(ret == 0){
        return 1;
    }
    test_conv_cipher_direct(input_dim, kernel_dim);
}
