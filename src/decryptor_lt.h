#pragma once

#include "seal/ciphertext.h"
#include "seal/context.h"
#include "seal/encryptionparams.h"
#include "seal/memorymanager.h"
#include "seal/modulus.h"
#include "seal/plaintext.h"
#include "seal/randomgen.h"
#include "seal/secretkey.h"
#include "seal/util/defines.h"
#include "seal/util/iterator.h"
#include "seal/util/locks.h"
#include "seal/util/ntt.h"
#include "seal/util/rns.h"
#include "util/packedconv.h"
#include <memory>
#include <vector>

using namespace seal;
using namespace seal::util;

class Decryptor_LT
{
    public:

        Decryptor_LT(const SEALContext &context, const SecretKey &secret_key);

        void linear_trans(Ciphertext &encrypted, std::vector<std::vector<uint64_t>> lt_matrix, Ciphertext &lt_cipher);

        void lt_packedconv(Ciphertext &encrypted, std::vector<util::KernelInfo> kernel_infos, Ciphertext &lt_cipher);

        void decrypt_bfv_lt(Ciphertext &encrypted, std::vector<std::vector<uint64_t>> &matrix_conved, uint64_t colsize, Plaintext &destination);

        void decrypt_bfv_lt_toeplitz(Ciphertext &encrypted, std::vector<std::vector<uint64_t>> &matrix_conved, uint64_t colsize, Plaintext &destination);

    private:
        void dot_product_with_secret_lt(Ciphertext &encrypted, std::vector<std::vector<uint64_t>> &matrix_conved, uint64_t colsize, util::RNSIter destination, MemoryPoolHandle pool);

        void dot_product_with_secret_lt_toeplitz(Ciphertext &encrypted, std::vector<std::vector<uint64_t>> matrix_conved, uint64_t colsize, util::RNSIter destination, MemoryPoolHandle pool);

        MemoryPoolHandle pool_ = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);

        SEALContext context_;

        std::size_t secret_key_array_size_ = 0;

        util::Pointer<std::uint64_t> secret_key_array_;
};
