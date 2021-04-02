#include "testprimitives.h"
#include "testutil.h"
#include "seal/util/uintlinarith.h"



void test_diagonal_partial(){
    uint64_t start_col = 4;
    uint64_t colsize = 3;
    uint64_t a_size = 10;
    Modulus modulus(15);
    vector<uint64_t> a(a_size);
    for(auto i = 0L;i < a.size();i++){
        a[i] = i;
    }
    vector<uint64_t> diagonal = util::create_diagonal_partial(a, start_col, colsize, modulus);
    cout << "start_col: " << start_col << " " << "colsize: " << colsize << endl;
    cout << "diagonal list of submatrix: ";
    for(auto i =0L;i < diagonal.size();i++){
        cout << diagonal[i] << " ";
    }
    cout << endl;
}

int main(){
    cout << "test_primitives" << endl;
    test_init_matrix_uint();
    test_diagonal_partial();
}
