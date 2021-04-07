#include "testprimitives.h"
#include "testutil.h"
#include "seal/util/uintlinarith.h"

void test_scalars_to_diagonallist(){
    vector<uint64_t> scalars = { 3, 2, 1, 5, 3, 8, 3, 4, 7, 1 };
    uint64_t colsize = 5;
    uint64_t rowsize = 6;
    vector<vector<uint64_t>> diagonallists = util::scalars_to_diagonallist(scalars, colsize, rowsize);
    util::print_matrix(diagonallists);
}

inline vector<uint64_t> sample_rn(uint64_t size, Modulus modulus){
    vector<uint64_t> ret(size);
    for(uint64_t i = 0;i < size;i++){
        ret[i] = rand() % modulus.value();
    }
    return ret;
}

void print_plain(Plaintext plain, uint64_t num){
    if(num > plain.coeff_count()){
        std::cout << "print_plain: warn: num is bigger than plain coefficient size!" << std::endl;
        num = plain.coeff_count();
    }
    std::cout << "Plaintext Coefficients(first " << num << " elements): [";
    for(auto i = 0U;i < num;i++){
        std::cout << plain[i] << " ";
    }
    std::cout << "]" << std::endl;
}

void test_delete_cipher(uint64_t input_dim, uint64_t kernel_dim, uint64_t poly_modulus_degree){
    print_example_banner("Direct convolution of ciphertext Benchmark");

    // parameter setting
    EncryptionParameters parms(scheme_type::bfv);
    //size_t poly_modulus_degree;
    //cout << "poly_modulus_degree: ";
    //cin >> poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    //vector<Modulus> mod_chain =  {Modulus(0xffffff00000001)};
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
    //cout << "plaintext polynomial " + x_plain.to_string() + "." << endl;
    //cout << "Coeff count: " << x_plain.coeff_count() << endl;
    print_plain(x_plain, 10);

    // convert plaintext by matrix
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    cout << "Copied and converted plaintext" << endl;

    // encrypt x
    Ciphertext x_encrypted;
    //cout << "----Encrypt x_plain to x_encrypted.----" << endl;
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "Coeff modulus size: " << x_encrypted.coeff_modulus_size() << endl;
    //uint64_t cipher_coeffsize = x_encrypted.size() * x_encrypted.poly_modulus_degree() * x_encrypted.coeff_modulus_size();
    //cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(x_encrypted) << " bits" << endl;

    // convolve encrypted x
    cout << "decryption of x_tranformed: " << endl;
    Plaintext x_decrypted;
    decryptor.decrypt(x_encrypted, x_decrypted);

    print_plain(x_decrypted, 20);
}

// 行列操作など基礎的な関数のテストのためのもの
int main(){
    cout << "test_primitives" << endl;
    //test_scalars_to_diagonallist();
    test_delete_cipher(5, 0, 2048);
}
