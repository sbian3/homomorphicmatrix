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

// 行列操作など基礎的な関数のテストのためのもの
int main(){
    cout << "test_primitives" << endl;
    test_scalars_to_diagonallist();
}
