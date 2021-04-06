#include "testprimitives.h"
#include "testutil.h"
#include "seal/util/uintlinarith.h"



void test_prod_diagonal(){
    vector<uint64_t> kernel = {1, 2, 4};
    uint64_t colsize_K = 4;
    uint64_t rowsize_K = 5;
    uint64_t colsize_R = 5;
    uint64_t rowsize_R = 6;
    vector<uint64_t> diagonallist_R = {3, 2, 1, 5, 3, 8, 3, 4, 7, 1};
    Modulus modulus(11);
    vector<uint64_t> kernel_diagonallist(colsize_K + rowsize_K - 1);
    vector<uint64_t> indexes = util::create_diagonal_list(kernel, colsize_K, rowsize_K, modulus, kernel_diagonallist);
    int64_t k = static_cast<int64_t>(colsize_K);
    k  = -k + 1;
    vector<vector<uint64_t>> matrix_product_diagonals(colsize_K + rowsize_R-1);
    uint64_t index = 0;
    for(;k < static_cast<int64_t>(rowsize_R);k++){
        vector<uint64_t> diagonal_vec;
        diagonal_vec = util::matrix_product_diagonal(k, colsize_R, rowsize_R, kernel_diagonallist, indexes, diagonallist_R, modulus);
        cout << "offset: " << k << endl;
        util::print_iter(diagonal_vec, diagonal_vec.size());
        matrix_product_diagonals[index] = diagonal_vec;
        index++;
    }
    cout << "result matrix" << endl;
    vector<vector<uint64_t>> result(colsize_K , vector<uint64_t>(rowsize_R));
    util::diagonallist_to_matrix(matrix_product_diagonals, 0, 0, colsize_K, rowsize_R, result);
    util::print_matrix(result);
}

int main(){
    cout << "test_primitives" << endl;
    test_diagonal_from_submat();
    test_prod_diagonal();
}
