#include "bench.h"
#include <chrono>

using namespace std;
using namespace seal;

inline std::uint64_t multiply_uint_mod2k(
        std::uint64_t operand1, std::uint64_t operand2, const Modulus &modulus)
{
    unsigned long long z[2];
    util::multiply_uint64(operand1, operand2, z);
    return z[0] & (modulus.value() - 1);
}

inline std::uint64_t uint64_mod2k(
        std::uint64_t operand1, const Modulus &modulus)
{
    return operand1 & (modulus.value() - 1);
}

void test_mod2k(){
    uint64_t operand1 = 122;
    uint64_t operand2 = 7;
    Modulus modulus(8);
    uint64_t operand1_mod = uint64_mod2k(operand1, modulus);
    uint64_t mult_mod = multiply_uint_mod2k(operand1, operand2, modulus);
    cout << "operand1: " << operand1 << endl;
    cout << "operand2: " << operand2 << endl;
    cout << operand1 << " mod " << modulus.value() << "= " << operand1_mod << endl;
    cout << operand1 << " * " << operand2 << " mod " << modulus.value() << "= " << mult_mod << endl;
}

void multiply_mod_loop(std::uint64_t operand1,std::uint64_t operand2,const Modulus &modulus, uint64_t loop_cnt){
}

void bench_mod(){
    uint64_t operand1 = 122;
    uint64_t operand2 = 7;
    Modulus modulus(8);
    
    uint64_t mult_mod;
    uint64_t loop_cnt = 2048 * 25;
    auto mod_2k_begin = chrono::high_resolution_clock::now();
    for(uint64_t i = 0;i < loop_cnt;i++){
        mult_mod= multiply_uint_mod2k(operand1, operand2, modulus);
    }
    auto mod_2k_end = chrono::high_resolution_clock::now();

    auto mod_barrett_begin = chrono::high_resolution_clock::now();
    for(uint64_t i = 0;i < loop_cnt;i++){
        mult_mod= multiply_uint_mod(operand1, operand2, modulus);
    }
    auto mod_barrett_end = chrono::high_resolution_clock::now();
    
    auto barrett_time = chrono::duration_cast<chrono::microseconds>(mod_barrett_end - mod_barrett_begin);
    auto mod_2k_time = chrono::duration_cast<chrono::microseconds>(mod_2k_end - mod_2k_begin);
    cout << "mod_barrett: " << barrett_time.count() << US << endl;
    cout << "mod 2k:      " << mod_2k_time.count()  << US << endl;
}


int main(int argc, char* argv[]){
    bench_mod();
}
