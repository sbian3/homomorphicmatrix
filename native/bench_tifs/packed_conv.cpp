#include "bench.h"

uint64_t get_blocksize(uint64_t input_dim, uint64_t kernel_dim){

}

vector<uint64_t> pack_input(vector<vector<uint64_t>> input, uint64_t block_size, uint64_t poly_size){

}

vector<uint64_t> pack_kernel(vector<vector<uint64_t>> kernel, uint64_t block_size, uint64_t poly_size){

}

int main(){
    print_example_banner("Packed Convolution Benchmark");
    vector<uint64_t> coeff = {1, 2, 3,4, 10};
    Plaintext plain(coeff);
    print_plain(plain, coeff.size());
}
