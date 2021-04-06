#include "bench.h"

class KernelInfo{
    public:
        KernelInfo(){

        }
        KernelInfo(uint64_t st_c,uint64_t  st_r,uint64_t si_c,uint64_t si_r,vector<uint64_t> data, Modulus mod):
            start_col(st_c), start_row(st_r), size_col(si_c), size_row(si_r), data(data), modulus(mod){
                diagonal_list = vector<uint64_t>(size_col + size_row - 1);
                index = util::create_diagonal_list(data, size_col, size_row, modulus, diagonal_list);
        }
        vector<uint64_t> diagonal_list;
        vector<uint64_t> index;

        uint64_t get_colsize(){
            return size_col;
        }

        uint64_t get_rowsize(){
            return size_row;
        }

        uint64_t get_startcol(){
            return start_col;
        }

        uint64_t get_startrow(){
            return start_row;
        }

        void getParamsforSubmatrix(uint64_t &submat_startcol, uint64_t &submat_colsize){
            submat_startcol = start_row;
            submat_colsize = size_row;
        }

        void print(){
            cout << "start col: " << start_col << endl;
            cout << "start row: " << start_row << endl;
            cout << "diagonal: ";
            for(uint64_t k = 0;k < diagonal_list.size();k++){
                cout << diagonal_list[k] << " ";
            }
            cout << endl;
            cout << "index: ";
            for(uint64_t k = 0;k < index.size();k++){
                cout << index[k] << " ";
            }
            cout << endl;
        }

    private:
        uint64_t size_col;
        uint64_t size_row;
        uint64_t start_col;
        uint64_t start_row;
        Modulus modulus;
        vector<uint64_t> data;
};

inline uint64_t get_blocksize(uint64_t input_dim, uint64_t kernel_dim, uint64_t padding){
    return input_dim + kernel_dim - 1 + padding;
}

vector<uint64_t> pack_input(const vector<vector<uint64_t>> input, uint64_t block_size, uint64_t poly_size){
    assert(input.size() * block_size <= poly_size);
    vector<uint64_t> packed_input(poly_size);
    for(uint64_t i = 0;i < input.size();i++){
        uint64_t input_start = block_size * i;
        for(uint64_t j = 0;j < input[i].size();j++){
            packed_input[input_start + j] = input[i][j];
        }
    }
    return packed_input;
}

// kernel: 初期化したkernel vectorのvector. 長さはpack_num個
// block_size: 1ブロックの大きさ．
// return: KernelInfoのベクトルを生成して返す
vector<KernelInfo> pack_kernel(vector<vector<uint64_t>> kernels, uint64_t block_size, Modulus modulus){
    uint64_t packing_num = kernels.size();
    vector<KernelInfo> kernel_info(packing_num);
    uint64_t start_col = 0;
    uint64_t start_row = 0;
    for(uint64_t i = 0;i < packing_num;i++){
        KernelInfo kinfo(start_col, start_row, block_size, block_size, kernels[i], modulus);
        kernel_info[i] = kinfo;
        start_col += block_size;
        start_row += block_size;
    }
    return kernel_info;
}

void pack_kernel_to_matrix(vector<KernelInfo> kernelinfos, vector<vector<uint64_t>> &matrix){
    for(uint64_t i = 0;i < kernelinfos.size();i++){
        KernelInfo kinfo = kernelinfos[i];
        // kernel scalar is reversed in default
        // so, scalar must be reversed again
        vector<uint64_t> k_scalar = kinfo.diagonal_list;
        reverse(k_scalar.begin(), k_scalar.end());
        vector<vector<uint64_t>> diagonals = util::scalars_to_diagonallist(k_scalar, kinfo.get_colsize(), kinfo.get_rowsize());
        util::diagonallist_to_matrix(diagonals, kinfo.get_startcol(), kinfo.get_startrow(), kinfo.get_colsize(), kinfo.get_rowsize(),matrix);
    }
}

// 行列積結果を対角成分のみから計算する．
void matrix_dot_matrix_toeplitz_mod(vector<KernelInfo> kernel_infos, CoeffIter c1, uint64_t poly_degree, vector<vector<uint64_t>> &result, Modulus &modulus){
    // for each block
    for(uint64_t i = 0;i < kernel_infos.size();i++){
        // get diagonal lists for kernel
        KernelInfo kinfo = kernel_infos[i];
        vector<uint64_t> kernel_diagonal_list = kinfo.diagonal_list;
        vector<uint64_t> kernel_index = kinfo.index;
        uint64_t colsize_K = kinfo.get_colsize();

        // diagonal list of c1
        uint64_t submat_startcol,submat_startrow, submat_colsize, submat_rowsize;
        submat_rowsize = poly_degree;
        submat_startrow = 0;
        kinfo.getParamsforSubmatrix(submat_startcol, submat_colsize);
        vector<uint64_t> diagonal_c1 = util::create_diagonal_from_submatrix(c1, poly_degree , submat_startcol, submat_colsize, modulus);

        // calc diagonal of product
        vector<vector<uint64_t>> matrix_product_diagonals(colsize_K + submat_rowsize - 1);
        uint64_t index = 0;
        int64_t k = static_cast<int64_t>(colsize_K);
        k = -k+1;
        for(k;k<static_cast<int64_t>(submat_rowsize);k++){
            vector<uint64_t> diagonal_vec;
            diagonal_vec = util::matrix_product_diagonal(k, submat_colsize, submat_rowsize, kernel_diagonal_list, kernel_index, diagonal_c1, modulus);
            matrix_product_diagonals[index] = diagonal_vec;
            index++;
        }
        
        // write diagonals to result matrix
        util::diagonallist_to_matrix(matrix_product_diagonals, submat_startcol, submat_startrow, colsize_K, submat_rowsize, result);
    }
}

void make_packedconv_matrixproduct(vector<KernelInfo> kernel_infos, Ciphertext &encrypted, uint64_t poly_degree, vector<vector<uint64_t>> &result, Modulus modulus){
    PolyIter cipher_poly(encrypted);
    cipher_poly++;
    matrix_dot_matrix_toeplitz_mod(kernel_infos, **cipher_poly, poly_degree, result, modulus);
}

void print_input_kernel(vector<vector<uint64_t>> input, vector<vector<uint64_t>> kernel){
    cout << "inputs: " << endl;
    util::print_matrix(input);
    cout << "kernels: " << endl;
    util::print_matrix(kernel);
}


void bench_packed_conv(uint64_t input_dim, uint64_t kernel_dim, uint64_t pack_num){
    print_example_banner("Homomorphic packed convolution Benchmark");

    // sample input and kernel
    Modulus sample_mod(7);
    vector<vector<uint64_t>> input = sample_rn(pack_num, input_dim, sample_mod);
    vector<vector<uint64_t>> kernel = sample_rn(pack_num, kernel_dim, sample_mod);
    print_input_kernel(input, kernel);


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
    cout << "Input plaintext: ";
    Plaintext x_plain(packed_input);
    cout << "Express x = as a plaintext polynomial " + x_plain.to_string() + "." << endl;
    cout << "Coeff count: " << x_plain.coeff_count() << endl;
    print_plain(x_plain, 10);

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
    auto time_start = chrono::high_resolution_clock::now();
    make_packedconv_matrixproduct(kernelinfos, x_encrypted, poly_modulus_degree, matrix_conved, parms.coeff_modulus()[0]);
    decryptor.linear_trans(x_encrypted, matrix, x_enc_lin);
    //print_iter(PolyIter(x_encrypted), 2);
    //print_iter(PolyIter(x_enc_lin), 2);
    // decrypt
    Plaintext x_decrypted;
    decryptor.decrypt_bfv_lt(x_enc_lin, matrix_conved, x_decrypted);
    auto time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);

    // compare converted plain and decryption of x_converted
    //cout << "Converted plain: " << endl;
    //print_plain(copied_plain, 10);
    cout << "decryption of x_tranformed: " << endl;
    print_plain(x_decrypted, block_size * pack_num);
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
    uint64_t input_dim, kernel_dim, pack_num;
    if(argc != 4){
        cerr << "please input two numbers (input_dim, kernel_dim, pack_num)" << endl;
        return 1;
    }
    input_dim = stoull(argv[1]);
    kernel_dim = stoull(argv[2]);
    pack_num = stoull(argv[3]);
    bench_packed_conv(input_dim, kernel_dim, pack_num);
}
