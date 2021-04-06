#include "bench.h"

using namespace std;
using namespace seal;

void make_convedmatrix(vector<vector<uint64_t>> matrix, uint64_t poly_modulus_degree, Ciphertext &encrypted, Modulus modulus, vector<vector<uint64_t>> &matrix_conved){
    PolyIter enc_poly(encrypted);
    enc_poly++;
    util::matrix_dot_convedcoeff(matrix, poly_modulus_degree, **enc_poly, modulus, matrix_conved);
}

void test_bfv_matrix(){
    print_example_banner("Homomorphic general Linear Transformation Benchmark");

    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree;
    cout << "poly_modulus_degree: ";
    cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 7;
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

    // lt
    cout << "decryption of x_encrypted: ";
    auto time_start = chrono::high_resolution_clock::now();
    Plaintext x_decrypted;
    decryptor.decrypt_bfv_with_matrix(x_encrypted, x_decrypted, matrix);
    // decrypt
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);

    vector<vector<uint64_t>> matrix_conved(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    cout << "decryption of x_encrypted: ";
    make_convedmatrix(matrix, poly_modulus_degree, x_encrypted, parms.coeff_modulus()[0], matrix_conved);
    print_matrix(matrix_conved);
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_decrypted, 20);
}

void test_decrypt_separate(uint64_t input_dim, uint64_t kernel_dim, uint64_t poly_modulus_degree){
    print_example_banner("Homomorphic general Linear Transformation Benchmark");

    if(poly_modulus_degree < input_dim){
        throw invalid_argument("input_dim is too large");
    }

    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);

    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 7;
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
    vector<uint64_t> input = sample_rn(input_dim, Modulus(7));
    //string x = "6x^2 + 1x^1 + 2";
    //cout << "Input plaintext: ";
    Plaintext x_plain(input);
    //cout << "Express x = " + x + " as a plaintext polynomial " + x_plain.to_string() + "." << endl;
    //cout << "Coeff count: " << x_plain.coeff_count() << endl;
    print_plain(x_plain, 10);

    // generate transform matrix
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    util::init_matrix_identity(matrix, poly_modulus_degree, 2);
    matrix[0][1] = 1;

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
    //uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    //cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // lt
    Ciphertext x_enc_lin(x_encrypted);
    util::set_zero_poly(poly_modulus_degree, x_encrypted.coeff_modulus_size(), x_enc_lin.data());
    vector<vector<uint64_t>> matrix_conved(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    cout << "decryption of x_encrypted: ";
    auto lt_start = chrono::high_resolution_clock::now();
    make_convedmatrix(matrix, poly_modulus_degree, x_encrypted, parms.coeff_modulus()[0], matrix_conved);
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

    // compare converted plain and decryption of x_converted
    //cout << "Converted plain: " << endl;
    //print_plain(copied_plain, 10);
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_decrypted, 20);
}

int main(int argc, char* argv[]){
    uint64_t input_dim, kernel_dim,poly_modulus_degree;
    int ret = check_args(argc, argv, input_dim, kernel_dim, poly_modulus_degree);
    if(ret == 0){
        return 1;
    }
    //test_bfv_matrix();
    test_decrypt_separate(input_dim, kernel_dim, poly_modulus_degree);
}
