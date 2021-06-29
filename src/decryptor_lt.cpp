#include "decryptor_lt.h"
#include "seal/util/polyarithsmallmod.h"


Decryptor_LT::Decryptor_LT(const SEALContext &context, const SecretKey &secret_key): context_(context){
        // Verify parameters
        if (!context_.parameters_set())
        {
            throw invalid_argument("encryption parameters are not set correctly");
        }
        if (!is_valid_for(secret_key, context_))
        {
            throw invalid_argument("secret key is not valid for encryption parameters");
        }

        auto &parms = context_.key_context_data()->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();

        // Set the secret_key_array to have size 1 (first power of secret)
        // and copy over data
        secret_key_array_ = allocate_poly(coeff_count, coeff_modulus_size, pool_);
        set_poly(secret_key.data().data(), coeff_count, coeff_modulus_size, secret_key_array_.get());
        secret_key_array_size_ = 1;
}

// Decrypt linear transformed ciphertext
void Decryptor_LT::decrypt_bfv_lt(Ciphertext &encrypted, std::vector<std::vector<uint64_t>> &matrix_conved, uint64_t colsize, Plaintext &destination){
    // Verify that encrypted is valid.
    if (!is_valid_for(encrypted, context_))
    {
        throw invalid_argument("encrypted is not valid for encryption parameters");
    }

    // Additionally check that ciphertext doesn't have trivial size
    if (encrypted.size() < SEAL_CIPHERTEXT_SIZE_MIN)
    {
        throw invalid_argument("encrypted is empty");
    }

    //
    // copied from bfv_decrypt
    // 
    if (encrypted.is_ntt_form())
    {
        throw invalid_argument("encrypted cannot be in NTT form");
    }

    auto &context_data = *context_.get_context_data(encrypted.parms_id());
    auto &parms = context_data.parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t coeff_modulus_size = coeff_modulus.size();

    // Firstly find c_0 + c_1 *s + ... + c_{count-1} * s^{count-1} mod q
    // This is equal to Delta m + v where ||v|| < Delta/2.
    // Add Delta / 2 and now we have something which is Delta * (m + epsilon) where epsilon < 1
    // Therefore, we can (integer) divide by Delta and the answer will round down to m.

    // Make a temp destination for all the arithmetic mod qi before calling FastBConverse
    SEAL_ALLOCATE_ZERO_GET_RNS_ITER(tmp_dest_modq, coeff_count, coeff_modulus_size, pool_);

    // original dot_product function.
    //dot_product_with_matrix(encrypted, tmp_dest_modq, matrix, pool_);
    dot_product_with_secret_lt(encrypted, matrix_conved, colsize, tmp_dest_modq, pool_);

    // Allocate a full size destination to write to
    destination.parms_id() = parms_id_zero;
    destination.resize(coeff_count);

    // Divide scaling variant using BEHZ FullRNS techniques
    context_data.rns_tool()->decrypt_scale_and_round(tmp_dest_modq, destination.data(), pool_);

    // How many non-zero coefficients do we really have in the result?
    size_t plain_coeff_count = get_significant_uint64_count_uint(destination.data(), coeff_count);

    // Resize destination to appropriate size
    destination.resize(max(plain_coeff_count, size_t(1)));
}

// for linear transformation
void Decryptor_LT::dot_product_with_secret_lt(Ciphertext &encrypted, std::vector<std::vector<uint64_t>> matrix_conved, uint64_t colsize, util::RNSIter destination, MemoryPoolHandle pool){
    auto &context_data = *context_.get_context_data(encrypted.parms_id());
    auto &parms = context_data.parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t coeff_modulus_size = coeff_modulus.size();
    size_t key_coeff_modulus_size = context_.key_context_data()->parms().coeff_modulus().size();
    size_t encrypted_size = encrypted.size();
    //auto is_ntt_form = encrypted.is_ntt_form();

    auto ntt_tables = context_data.small_ntt_tables();

    if(encrypted_size != 2){
        throw invalid_argument("encrypted_size must be 2");
    }

    //compute_secret_key_array(encrypted_size - 1);
    SEAL_ALLOCATE_GET_RNS_ITER(secret_key_array, coeff_count, coeff_modulus_size, pool);
    set_poly(secret_key_array_.get(), coeff_count, coeff_modulus_size, secret_key_array);
    // transform secret key array into non-NTT form
    inverse_ntt_negacyclic_harvey(secret_key_array, coeff_modulus_size, ntt_tables);

    PolyIter cipher_polyiter(encrypted);
    set_poly(cipher_polyiter, coeff_count, coeff_modulus_size, destination);
    cipher_polyiter++;
    SEAL_ALLOCATE_ZERO_GET_RNS_ITER(C1_s, coeff_count,coeff_modulus_size, pool);
    auto matrix_s = chrono::high_resolution_clock::now();
    SEAL_ITERATE(iter(*cipher_polyiter, secret_key_array, coeff_modulus, C1_s), coeff_modulus_size, [&](auto I){
            matrix_dot_vector(matrix_conved, colsize, get<1>(I), get<2>(I), coeff_count, get<3>(I));
            });
    auto matrix_e = chrono::high_resolution_clock::now();
    auto matrix_diff = chrono::duration_cast<chrono::microseconds>(matrix_e - matrix_s);
    cout << "decrypt: matrix_dot_vector: " << matrix_diff.count() << " us" << endl;
    // add c0 and c1
    add_poly_coeffmod(destination, C1_s, coeff_modulus_size, coeff_modulus, destination);
}

void Decryptor_LT::linear_trans(Ciphertext &encrypted, vector<std::vector<uint64_t>> lt_matrix, Ciphertext &lt_cipher){
    auto &context_data = *context_.get_context_data(encrypted.parms_id());
    auto &parms = context_data.parms();
    auto &coeff_modulus = parms.coeff_modulus();
    auto coeff_degree = parms.poly_modulus_degree();
    size_t coeff_modulus_size = coeff_modulus.size();

    PolyIter cipher_polyiter(encrypted);
    PolyIter lt_cipher_polyiter(lt_cipher);
    util::matrix_dot_vector(lt_matrix, coeff_degree, coeff_modulus_size, *cipher_polyiter, coeff_modulus, *lt_cipher_polyiter);
}

// WIP
void Decryptor_LT::lt_packedconv(Ciphertext &encrypted, std::vector<util::KernelInfo> kernel_infos, Ciphertext &lt_cipher){
    auto &context_data = *context_.get_context_data(encrypted.parms_id());
    auto &parms = context_data.parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    PolyIter cipher_polyiter(encrypted);
    PolyIter lt_cipher_polyiter(lt_cipher);
    util::kernel_matrix_dot_vector(kernel_infos, **cipher_polyiter, coeff_modulus[0], **lt_cipher_polyiter);
    //util::matrix_dot_vector(lt_matrix, coeff_modulus_size, *cipher_polyiter, coeff_modulus, *lt_cipher_polyiter);
}
