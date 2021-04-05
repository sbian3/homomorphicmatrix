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
        uint64_t start_col;
        uint64_t start_row;
    private:
        uint64_t size_col;
        uint64_t size_row;
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

void test_init_kernelinfo(){
    vector<vector<uint64_t>> kernels = { {1 , 2, 3}, {4, 5, 6} };
    uint64_t block_size = 32;
    Modulus modulus(7);
    vector<KernelInfo> kernel_infos = pack_kernel(kernels, block_size, modulus);
    for(uint64_t i = 0;i < kernel_infos.size();i++){
        cout << "start col: " << kernel_infos[i].start_col << endl;
        cout << "start row: " << kernel_infos[i].start_row << endl;
        cout << "diagonal: ";
        vector<uint64_t> diagonal = kernel_infos[i].diagonal_list;
        for(uint64_t k = 0;k < diagonal.size();k++){
            cout << diagonal[k] << " ";
        }
        cout << endl;
        cout << "index: ";
        vector<uint64_t> index = kernel_infos[i].index;
        for(uint64_t k = 0;k < index.size();k++){
            cout << index[k] << " ";
        }
        cout << endl;
    }
}

int main(){
    print_example_banner("Packed Convolution Benchmark");
    vector<uint64_t> coeff = {1, 2, 3,4, 10};
    Plaintext plain(coeff);
    print_plain(plain, coeff.size());
    test_init_kernelinfo();
}
