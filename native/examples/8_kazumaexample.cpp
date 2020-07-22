#include "seal/util/defines.h"
#include "examples.h"
#include <cassert>

using namespace std;
using namespace seal;

void convert(Plaintext& src, Plaintext& dst, vector<vector<uint64_t>>& matrix, uint64_t poly_degree){
    assert(matrix.size() == poly_degree);
    assert(matrix[0].size() == poly_degree);
    assert(dst.coeff_count() == poly_degree);
    uint64_t sum ;
    for(int i = 0;i < poly_degree;i++){
        sum = 0;
        for( int j = 0;j < src.coeff_count();j++ ){
            sum += src[j] * matrix.at(i).at(j);
        }
        dst[i] = sum;
    }
}

void print_plain_coefficients(Plaintext plain){
    // print all coefficients
    cout << "Plain: [";
    for(int i = 0;i < plain.coeff_count();i++){
        cout << plain[i] << " ";
    }
    cout << "]" << endl;
}

void matrix_conversion(){
    print_example_banner("matrix_conversion");
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = 4096;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    vector<Modulus> mod_chain = {Modulus(0x3ffffffdf0001)};
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 1009;
    parms.set_plain_modulus(plaintext_modulus);
    auto context = SEALContext::Create(parms);
    print_line(__LINE__);
    cout << "Set encryption parameters and print" << endl;
    print_parameters(context);
    cout << "Parameter validation (success): " << context->parameter_error_message() << endl;

    cout << endl;

    cout << "keygen step" << endl;
    KeyGenerator keygen(context);
    cout << "pubkey " << endl;
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    print_line(__LINE__);
    /* just encrypt and decrypt */
    string x = "5x^2 + 1x^1 + 2";
    cout << "Input plaintext: ";
    Plaintext x_plain(x);
    cout << "Express x = " + x + " as a plaintext polynomial " + x_plain.to_string() + "." << endl;
    cout << "Coeff count: " << x_plain.coeff_count() << endl;
    print_plain_coefficients(x_plain);
    // generate random vector
    vector<vector<uint64_t>> matrix(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree));
    std::random_device rnd;
    for(int i = 0;i < poly_modulus_degree;i++){
        for(int j = 0;j < poly_modulus_degree;j++){
            matrix[i][j] = rnd() % plaintext_modulus;
        }
    }
    // create copy and resize
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    convert(x_plain, copied_plain, matrix, poly_modulus_degree);
    cout << "Copied and converted plaintext" << endl;
    print_plain_coefficients(copied_plain);
    print_line(__LINE__);
    // encryption
    Ciphertext x_encrypted;
    cout << "Encrypt x_plain to x_encrypted." << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    // print ciphertext coefficients
    for(uint64_t i = 0;i < 10;i++){
        cout << "Cipher[" << i << "]: " << x_encrypted[i] << endl;
    }
    cout << "    + size of freshly encrypted x: " << x_encrypted.size() << endl;
    cout << "    + noise budget in x_sq_plus_one: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;
    // decryption
    Plaintext x_decrypted;
    cout << "    + decryption of x_encrypted: ";
    decryptor.decrypt(x_encrypted, x_decrypted);
    cout << x_decrypted.to_string() << " ...... Correct?" << endl;

}

void integer_encoder(){
    print_example_banner("Kazuma example: integer encoder");
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = 4096;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    //parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    vector<Modulus> mod_chain = {Modulus(1073153)};
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 1024;
    parms.set_plain_modulus(plaintext_modulus);
    auto context = SEALContext::Create(parms);
    print_line(__LINE__);
    cout << "Set encryption parameters and print" << endl;
    print_parameters(context);
    cout << "Parameter validation (success): " << context->parameter_error_message() << endl;

    cout << endl;

    // setup environment
    cout << "keygen step" << endl;
    KeyGenerator keygen(context);
    cout << "pubkey " << endl;
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    IntegerEncoder encoder_int(context);
    print_line(__LINE__);
    // just encrypt and decrypt
    cout << "Input plain value: ";
    int plainvalue;
    cin >> plainvalue;
    Plaintext x_plain("5x^1 + 1");
    //Plaintext x_plain = encoder_int.encode(plainvalue);
    cout << "Express x = " << plainvalue <<" as a plaintext polynomial " << x_plain.to_string() + "." << endl;
    cout << "Coeff count: " << x_plain.coeff_count() << endl;
    // print all coefficients
    for(int i = 0;i < x_plain.coeff_count();i++){
        cout << "Plain[" << i << "]: " << x_plain[i] << endl;
    }
    print_line(__LINE__);
    Ciphertext x_encrypted;
    cout << "Encrypt x_plain to x_encrypted." << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    // print ciphertext coefficients
    for(uint64_t i = 0;i < 10;i++){
        cout << "Cipher[" << i << "]: " << x_encrypted[i] << endl;
    }
    cout << "    + size of freshly encrypted x: " << x_encrypted.size() << endl;
    Plaintext x_decrypted;
    cout << "    + decryption of x_encrypted: ";
    decryptor.decrypt(x_encrypted, x_decrypted);
    //cout << encoder_int.decode_int32(x_decrypted) << endl;
    cout << "Decrypted as a plaintext polynomial " << x_decrypted.to_string() + "." << endl;

    /* compute (3x^2 + 1) and decrypt*/
    print_line(__LINE__);
    cout << "Compute x_sq_plus_one (3x^2+1)." << endl;
    Ciphertext x_sq_plus_one;
    evaluator.square(x_encrypted, x_sq_plus_one);
    Plaintext three("3");
    evaluator.multiply_plain_inplace(x_sq_plus_one, three);
    Plaintext plain_one("1");
    evaluator.add_plain_inplace(x_sq_plus_one, plain_one);
    cout << "Coeff modulus size: " << x_sq_plus_one.coeff_modulus_size() << endl;
    cout << "    + size of x_sq_plus_one: " << x_sq_plus_one.size() << endl;
    cout << "    + noise budget in x_sq_plus_one: " << decryptor.invariant_noise_budget(x_sq_plus_one) << " bits"
         << endl;
    Plaintext decrypted_result;
    cout << "    + decryption of x_sq_plus_one: ";
    decryptor.decrypt(x_sq_plus_one, decrypted_result);
    cout << "0x" << decrypted_result.to_string() << " ...... Correct." << endl;

}

void example_kazuma1(){
    print_example_banner("Kazuma example");
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = 4096;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    //vector<Modulus> mod_chain = {Modulus(0x7e00001)};
    //    //vector<Modulus> mod_chain = {Modulus(0x3ffffac001)};
    vector<Modulus> mod_chain = {Modulus(0x3ffffffdf0001)};
    //vector<Modulus> mod_chain1 = CoeffModulus::BFVDefault(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 1009;
    parms.set_plain_modulus(plaintext_modulus);
    auto context = SEALContext::Create(parms);
    print_line(__LINE__);
    cout << "Set encryption parameters and print" << endl;
    print_parameters(context);
    cout << "Parameter validation (success): " << context->parameter_error_message() << endl;

    cout << endl;

    cout << "keygen step" << endl;
    KeyGenerator keygen(context);
    cout << "pubkey " << endl;
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    print_line(__LINE__);
    /* just encrypt and decrypt */
    string x = "5x^2 + 2";
    cout << "Input plaintext: ";
    Plaintext x_plain(x);
    cout << "Express x = " + x + " as a plaintext polynomial " + x_plain.to_string() + "." << endl;
    cout << "Coeff count: " << x_plain.coeff_count() << endl;
    // print all coefficients
    for(int i = 0;i < x_plain.coeff_count();i++){
        cout << "Plain[" << i << "]: " << x_plain[i] << endl;
    }
    print_line(__LINE__);
    // encryption
    Ciphertext x_encrypted;
    cout << "Encrypt x_plain to x_encrypted." << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    // print ciphertext coefficients
    for(uint64_t i = 0;i < 10;i++){
        cout << "Cipher[" << i << "]: " << x_encrypted[i] << endl;
    }
    cout << "    + size of freshly encrypted x: " << x_encrypted.size() << endl;
    cout << "    + noise budget in x_sq_plus_one: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;
    // decryption
    Plaintext x_decrypted;
    cout << "    + decryption of x_encrypted: ";
    decryptor.decrypt(x_encrypted, x_decrypted);
    cout << x_decrypted.to_string() << " ...... Correct?" << endl;

    /* compute (3x^2 + 1) and decrypt*/
    print_line(__LINE__);
    cout << "Compute x_sq_plus_one (3x^2+1)." << endl;
    Ciphertext x_sq_plus_one;
    evaluator.square(x_encrypted, x_sq_plus_one);
    Plaintext three("3");
    evaluator.multiply_plain_inplace(x_sq_plus_one, three);
    Plaintext plain_one("1");
    evaluator.add_plain_inplace(x_sq_plus_one, plain_one);
    cout << "Coeff modulus size: " << x_sq_plus_one.coeff_modulus_size() << endl;
    cout << "    + size of x_sq_plus_one: " << x_sq_plus_one.size() << endl;
    cout << "    + noise budget in x_sq_plus_one: " << decryptor.invariant_noise_budget(x_sq_plus_one) << " bits" << endl;
    Plaintext decrypted_result;
    cout << "    + decryption of x_sq_plus_one: ";
    decryptor.decrypt(x_sq_plus_one, decrypted_result);
    cout << "0x" << decrypted_result.to_string() << " ...... Correct." << endl;

    /* compute x^3 */
    //cout << "Computer x^3" << endl;
    //Ciphertext x_squared;
    //evaluator.square(x_encrypted, x_squared);
    //Ciphertext x_boxed;
    //evaluator.multiply(x_squared, x_encrypted, x_boxed);
    //cout << "Coeff modulus size: " << x_boxed.coeff_modulus_size() << endl;
    //cout << " + size of x^3: " << x_boxed.size() << endl;
    //cout << " + noise budget in x^3 :" << decryptor.invariant_noise_budget(x_boxed) << endl;
    //cout << " + decrypt x^3: ";
    //Plaintext box_decrypted;
    //decryptor.decrypt(x_boxed, box_decrypted);
    //cout << "0x" << box_decrypted.to_string() <<  endl;
}

void example_kazuma(){
    //integer_encoder();
    //simple_encryption();
    matrix_conversion();
}
