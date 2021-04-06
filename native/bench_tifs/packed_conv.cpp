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
int main(){
    test_kernel_dot_c1();
}
