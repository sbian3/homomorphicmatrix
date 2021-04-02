#include "testprimitives.h"
#include "testutil.h"
#include "seal/util/uintlinarith.h"


// assume kernel_L is reversed( indexes is also )
vector<uint64_t> calc_product_diagonal(int64_t offset, uint64_t colsize_R, uint64_t rowsize_R, vector<uint64_t> kernel_L, vector<uint64_t> kernel_L_indexes, vector<uint64_t> list_R, Modulus & modulus){
    // RはL以上を想定
    assert(kernel_L_indexes.size() <= list_R.size());

    // 要素積を計算
    // kernelの非ゼロ要素の場所を覚えることで計算時間を短縮
    //uint64_t wise_prod_len = kernel_L.size() <= list_R.size()? kernel_L.size(): list_R.size();
    uint64_t wise_prod_len;
    if(offset >= 0){
        if(offset + kernel_L.size() > list_R.size()){
            wise_prod_len = list_R.size() - offset;
        }else{
            wise_prod_len = kernel_L.size(); 
        }
    }else{
        wise_prod_len = kernel_L.size() + offset;
    }
    vector<uint64_t> wise_prod(wise_prod_len);
    for(uint64_t i = 0;i < kernel_L_indexes.size();i++){
        uint64_t prod;
        if(offset >= 0){
            if(kernel_L_indexes[i]+offset >= list_R.size() || kernel_L_indexes[i] >= kernel_L.size()) break;
            prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
        }else{
            if(kernel_L_indexes[i] <= -offset-1) continue;
            if(kernel_L_indexes[i]-offset >= kernel_L.size() || kernel_L_indexes[i] >= list_R.size()) break;
            prod = util::multiply_uint_mod(kernel_L[kernel_L_indexes[i]], list_R[kernel_L_indexes[i]+offset], modulus);
        }
        if(offset < 0){
            wise_prod[kernel_L_indexes[i]+offset] = prod;
        }else{
            wise_prod[kernel_L_indexes[i]] = prod;
        }
    }
    vector<uint64_t> diagonal;
    uint64_t partial_sum = 0;
    uint64_t prod_times = rowsize_R;
    // 最初の内積を計算
    for(uint64_t i = 0;i < prod_times;i++){
        partial_sum = util::add_uint_mod(partial_sum, wise_prod[i], modulus);
    }
    diagonal.reserve(rowsize_R);
    diagonal.push_back(partial_sum);
    // 次の対角成分はO(1)で計算
    for(uint64_t i = 0;i < wise_prod.size() - prod_times;i++){
        partial_sum = util::add_uint_mod(partial_sum, wise_prod[i+prod_times], modulus);
        partial_sum = util::sub_uint_mod(partial_sum, wise_prod[i], modulus);
        diagonal.push_back(partial_sum);
    }
    return diagonal;
}

void test_prod_diagonal(){
    vector<uint64_t> kernel = {1, 2, 4};
    uint64_t colsize_K = 5;
    uint64_t rowsize_K = 5;
    uint64_t colsize_R = 5;
    uint64_t rowsize_R = 5;
    vector<uint64_t> list_R = {3, 2, 1, 5, 3, 8, 3, 4, 7};
    Modulus modulus(11);
    vector<uint64_t> kernel_diagonallist(colsize_K + rowsize_K - 1);
    vector<uint64_t> indexes = util::create_diagonal_list(kernel, colsize_K, rowsize_K, modulus, kernel_diagonallist);
    vector<uint64_t> diagonal_0 = calc_product_diagonal(-3, colsize_R, rowsize_R, kernel_diagonallist, indexes, list_R, modulus);
    for(auto i = 0;i < diagonal_0.size();i++){
        cout << diagonal_0[i] << " ";
    }
    cout << endl;
}

int main(){
    cout << "test_primitives" << endl;
    test_prod_diagonal();
}
