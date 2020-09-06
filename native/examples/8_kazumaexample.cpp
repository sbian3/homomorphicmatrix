#include "seal/util/defines.h"
#include "examples.h"
#include <cassert>

using namespace std;
using namespace seal;

void convert(Plaintext& src, Plaintext& dst, vector<vector<int64_t>>& matrix, uint64_t poly_degree, uint64_t mod_plain){
    assert(matrix.size() == poly_degree);
    assert(matrix[0].size() == poly_degree);
    assert(dst.coeff_count() == poly_degree);
    uint64_t sum ;
    for(auto i = 0U;i < poly_degree;i++){
        sum = 0;
        for( auto j = 0U;j < src.coeff_count();j++ ){
            sum += src[j] * matrix.at(i).at(j) % mod_plain;
        }
        dst[i] = sum % mod_plain;
    }
}

void convert(Ciphertext& src, Ciphertext& dst, vector<vector<int64_t>>& matrix, uint64_t poly_degree,vector<Modulus> mod_chain){
    assert(matrix.size() == poly_degree);
    assert(matrix[0].size() == poly_degree);
    uint64_t sum;
    for(auto size = 0U;size < src.size();size++){
        for(auto coeff_m = 0U; coeff_m < src.coeff_modulus_size();coeff_m++){
            for(auto i = 0U;i < poly_degree;i++){
                sum = 0;
                for(auto j = 0U;j < poly_degree;j++){
                    sum += (src[(size * src.coeff_modulus_size() + coeff_m ) * poly_degree + j] * matrix.at(i).at(j)) % mod_chain[coeff_m].value();
                }
                dst[(size * src.coeff_modulus_size() + coeff_m ) * poly_degree + i] = sum % mod_chain[coeff_m].value();
            }

        }
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

void print_plain(Plaintext plain, uint64_t num){
    if(num > plain.coeff_count()){
        cout << "print_plain: warn: num is bigger than plain coefficient size!" << endl;
        num = plain.coeff_count();
    }
    cout << "Plaintext Coefficients(first " << num << " elements): [";
    for(auto i = 0U;i < num;i++){
        cout << plain[i] << " ";
    }
    cout << "]" << endl;
}

void print_cipher_coefficients(Ciphertext cipher, uint64_t poly_modulus){
    cout << "Cipher Coefficients: [";
    for(auto i = 0U;i < cipher.coeff_modulus_size() * cipher.size() * poly_modulus;i++){
        cout << cipher[i] << " ";
    }
    cout << "]" << endl;
}

void print_cipher(Ciphertext cipher, uint64_t num){
    cout << "Ciphertext Coefficients(first " << num << " elements): [";
    for(auto i = 0U;i < num;i++){
        cout << cipher[i] << " ";
    }
    cout << "]" << endl;
}

void print_matrix(vector<vector<int64_t>>& matrix){
    for(auto i = 0U;i < matrix.size();i++){
        auto row = matrix[i];
        for( auto j = 0U;j < row.size();j++ ){
            cout << row.at(j) << " ";
        }
        cout << endl;
    } 
}

bool compare_cipher(Ciphertext cipher1, Ciphertext cipher2, uint64_t poly_modulus){
    uint64_t cipher1_count = cipher1.coeff_modulus_size() * cipher1.size() * poly_modulus;
    uint64_t cipher2_count = cipher2.coeff_modulus_size() * cipher2.size() * poly_modulus;
    if(cipher1_count != cipher2_count){
        cout << "coefficient number does not match: " << cipher1_count << " " << cipher2_count << endl;
        return false;
    }
    for(auto i = 0U;i < cipher1_count;i++){
        if(cipher1[i] != cipher2[i]){
            cout << i << "-th coefficient does not match!" << endl;
            cout << "cipher1: " << cipher1[i] << " cipher2: " << cipher2[i] << endl;
            return false;
        }
    }
    return true;
}

bool compare_cipher_all(Ciphertext cipher1, Ciphertext cipher2, uint64_t poly_modulus){
    uint64_t cipher1_count = cipher1.coeff_modulus_size() * cipher1.size() * poly_modulus;
    uint64_t cipher2_count = cipher2.coeff_modulus_size() * cipher2.size() * poly_modulus;
    if(cipher1_count != cipher2_count){
        cout << "coefficient number does not match: " << cipher1_count << " " << cipher2_count << endl;
        return false;
    }
    bool flag = true;
    for(auto i = 0U;i < cipher1_count;i++){
        if(cipher1[i] != cipher2[i]){
            cout << i << "-th coefficient does not match: " << "(" << cipher1[i] << ", " << cipher2[i] << ")" <<endl;
            flag = false;
        }
    }
    return flag;
}

void print_cipher_info(Ciphertext& c, Decryptor& decryptor){
    cout << "Size: " << c.size() << endl;
    cout << "Coeff modulus size: " << c.coeff_modulus_size() << endl;
    uint64_t cipher_coeffsize = c.size() * c.poly_modulus_degree() * c.coeff_modulus_size();
    cout << "Coeff size: " << cipher_coeffsize << endl;
    cout << "noise budget in ciphertext: " << decryptor.invariant_noise_budget(c) << " bits" << endl;
    cout << "is NTT form? " << c.is_ntt_form() << endl;
}

void init_matrix_identity(vector<vector<int64_t>>& matrix, uint64_t poly_modulus_degree){
    int64_t scale = 2;
    for(auto i = 0U;i < poly_modulus_degree;i++){
        for(auto j = 0U;j < poly_modulus_degree;j++){
            //matrix[i][j] = rnd() % 5;
            //if(matrix[i][j] == 4 || matrix[i][j] == 3 || matrix[i][j] == 2)
            //    matrix[i][j] = 0;
            if(i == j)
                matrix[i][j] = scale;
            else
                matrix[i][j] = 0;
        }
    }
    //matrix[0][1] = 1;
}

void init_matrix_identity_rnd(vector<vector<int64_t>>& matrix, uint64_t poly_modulus_degree){
    int64_t scale = 2;
    for(auto i = 0U;i < poly_modulus_degree;i++){
        for(auto j = 0U;j < poly_modulus_degree;j++){
            if(i == j)
                matrix[i][j] = scale;
            else
                matrix[i][j] = 0;
        }
    }
    matrix[0][0] = 1;
}

void init_matrix_rotate(vector<vector<int64_t>>& matrix, uint64_t size, int64_t right_rotate, int64_t scale){
    for(auto i = 0U;i < size;i++){
        for(auto j = 0U;j < size;j++){
            int64_t ii = i + right_rotate;
            bool reverse = false;
            if(ii < 0){
                ii+= size;
                reverse = true;
            }
            if(ii >= size){
                reverse = true;
            }
            if(j == ii% size){
                if(reverse)
                    matrix[i][j] = scale * -1;
                else
                    matrix[i][j] = scale;
            }
        }
    }
}

void init_matrix_rand(vector<vector<uint64_t>>& matrix, uint64_t size, uint64_t mod){
    std::random_device rnd;
    for(auto i = 0U;i < size;i++){
        for(auto j = 0U;j < size;j++){
            matrix[i][j] = rnd() % mod;
        }
    }
}

void test_conversion(){
    uint64_t size=3;
    vector<vector<int64_t>> matrix(size, vector<int64_t>(size));
    init_matrix_rotate(matrix, size, 1, 1);
    print_matrix(matrix);
    init_matrix_rotate(matrix, size, -1, 1);
    print_matrix(matrix);
    init_matrix_identity(matrix, size);
    print_matrix(matrix);
    init_matrix_identity_rnd(matrix, size);
    print_matrix(matrix);
}

void test_matconv(){
    uint64_t size=3;
    vector<vector<int64_t>> matrix(size, vector<int64_t>(size));
    init_matrix_rotate(matrix, size, 0, 1);
    init_matrix_rotate(matrix, size, -1, 2);
    init_matrix_rotate(matrix, size, -2, 3);
    print_matrix(matrix);
}

void matrix_conversion(){
    print_example_banner("matrix_conversion");

    // parameter setting
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = 4096;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    //uint64_t modulo_cipher = 0xfffffffd8001;
    //uint64_t modulo_cipher = 0xfffffffd8001;
    //0x3ffffffff040001
    //vector<Modulus> mod_chain = {Modulus(modulo_cipher)};
    vector<Modulus> mod_chain = CoeffModulus::BFVDefault(poly_modulus_degree);
    parms.set_coeff_modulus(mod_chain);
    uint64_t plaintext_modulus = 7;
    parms.set_plain_modulus(plaintext_modulus);
    auto context = SEALContext::Create(parms);
    print_line(__LINE__);
    cout << "Set encryption parameters and print" << endl;
    print_parameters(context);
    cout << "Parameter validation: " << context->parameter_error_message() << endl;

    cout << endl;

    // generate encryption helper
    cout << "keygen step" << endl;
    KeyGenerator keygen(context);
    cout << "pubkey " << endl;
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
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
    cout << "is NTT form? " << x_plain.is_ntt_form() << endl;
    print_plain(x_plain, 10);

    // generate transform matrix
    vector<vector<int64_t>> matrix(poly_modulus_degree, vector<int64_t>(poly_modulus_degree));
    init_matrix_identity_rnd(matrix, poly_modulus_degree);
    //print_matrix(matrix);

    // convert plaintext by matrix
    Plaintext copied_plain = Plaintext(x_plain);
    copied_plain.resize(poly_modulus_degree);
    convert(x_plain, copied_plain, matrix, poly_modulus_degree, plaintext_modulus);
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
    cout << "is NTT form? " << x_encrypted.is_ntt_form() << endl;
    print_cipher(x_encrypted, 5);
    cout << "4095-th coeff is " << x_encrypted[4095] << endl;

    // convert ciphertext by matrix
    cout << "----convert ciphertext by matrix---" << endl;
    Ciphertext x_enc_converted = Ciphertext(x_encrypted);
    convert(x_encrypted, x_enc_converted, matrix, poly_modulus_degree, mod_chain);
    print_cipher(x_enc_converted, 5);
    cout << "ciphertext info---->" << endl;
    print_cipher_info(x_enc_converted, decryptor);

    // multiply ciphertext by constant using seal function
    cout << "generate constant multiplied ciphertext by using seal" << endl;
    Ciphertext x_constant;
    Plaintext constant("2");
    evaluator.multiply_plain(x_encrypted, constant, x_constant);
    // constant addition
    cout << "add const 1" << endl;
    Plaintext addconst("1");
    Ciphertext addconst_enc;
    encryptor.encrypt(addconst, addconst_enc);
    print_cipher(addconst_enc, 5);
    evaluator.add_plain_inplace(x_constant, addconst);
    print_cipher(x_constant, 5);
    print_cipher_info(x_constant, decryptor);

    // compare cipher
    //cout << "Compare: " << compare_cipher_all(x_enc_converted, x_constant, poly_modulus_degree) << endl;

    // decryption normal x
    print_line(__LINE__);
    Plaintext x_decrypted;
    cout << "    + decryption of x_encrypted: ";
    decryptor.decrypt(x_encrypted, x_decrypted);
    cout << x_decrypted.to_string() << " ...... Correct?" << endl;

    // decrypt constant multiplied x
    cout << "seal generated constant multiplied x" << endl;
    Plaintext x_constant_dec;
    decryptor.decrypt(x_constant, x_constant_dec);
    print_plain(x_constant_dec, 10);

    // decryption converted x
    Plaintext x_converted_decrypted;
    cout << "decryption of x_copied: " << endl;;
    decryptor.decrypt(x_enc_converted, x_converted_decrypted);

    // compare converted plain and decryption of x_converted
    cout << "Converted plain: " << endl;
    print_plain(copied_plain, 10);
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_converted_decrypted, 20);
}

void example_kazuma(){
    //integer_encoder();
    //simple_encryption();
    //matrix_conversion();
    //test_conversion();
    test_matconv();
}
