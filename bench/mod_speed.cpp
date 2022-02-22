#include "bench.h"
#include <chrono>

using namespace std;
using namespace seal;
typedef chrono::nanoseconds unit_mult_bench;
typedef chrono::nanoseconds  unit_mod_bench;

inline std::uint64_t mod_barrett64(
        std::uint64_t operand1, const Modulus &modulus)
{
    return util::barrett_reduce_64(operand1, modulus);
}

inline std::uint64_t mod2k_uint64(
        std::uint64_t operand1, const Modulus &modulus)
{
    return operand1 & (modulus.value() - 1);
}

inline std::uint64_t multiply_sample(
        std::uint64_t operand1, std::uint64_t operand2, const Modulus &modulus)
{
    unsigned long long z[2];
    util::multiply_uint64(operand1, operand2, z);
    return z[0];
}


inline std::uint64_t multiply_uint_mod2k(
        std::uint64_t operand1, std::uint64_t operand2, const Modulus &modulus)
{
    unsigned long long z[2];
    util::multiply_uint64(operand1, operand2, z);
    return mod2k_uint64(z[0], modulus);
}


void test_mod2k(){
    uint64_t operand1 = 122;
    uint64_t operand2 = 7;
    Modulus modulus(8);
    uint64_t operand1_mod = mod2k_uint64(operand1, modulus);
    uint64_t mult_mod = multiply_uint_mod2k(operand1, operand2, modulus);
    cout << "operand1: " << operand1 << endl;
    cout << "operand2: " << operand2 << endl;
    cout << operand1 << " mod " << modulus.value() << "= " << operand1_mod << endl;
    cout << operand1 << " * " << operand2 << " mod " << modulus.value() << "= " << mult_mod << endl;
}

void multiply_mod_loop(std::uint64_t operand1,std::uint64_t operand2,const Modulus &modulus, uint64_t loop_cnt){
}

void bench_mod(uint64_t loop_cnt){
    vector<uint64_t> operand1(loop_cnt);
    for(int i = 0;i < loop_cnt;i++){
        operand1[i] = i +1;
    }
    uint64_t operand2 = 7;
    Modulus modulus(16);
    
    uint64_t mult_mod;
    auto mod_2k_begin = chrono::high_resolution_clock::now();
    for(uint64_t i = 0;i < loop_cnt;i++){
        mult_mod+= mod2k_uint64(operand1[i], modulus);
    }
    auto mod_2k_end = chrono::high_resolution_clock::now();
    cout << mult_mod << endl;
    mult_mod = 0;

    auto mod_barrett_begin = chrono::high_resolution_clock::now();
    for(uint64_t i = 0;i < loop_cnt;i++){
        mult_mod+= mod_barrett64(operand1[i], modulus);
    }
    auto mod_barrett_end = chrono::high_resolution_clock::now();
    cout << mult_mod << endl;
    mult_mod = 0;

    auto mult_mod_2k_begin = chrono::high_resolution_clock::now();
    for(uint64_t i = 0;i < loop_cnt;i++){
        mult_mod+= multiply_uint_mod2k(operand1[i], operand2, modulus);
    }
    auto mult_mod_2k_end = chrono::high_resolution_clock::now();
    cout << mult_mod << endl;
    mult_mod = 0;

    auto mult_mod_barrett_begin = chrono::high_resolution_clock::now();
    for(uint64_t i = 0;i < loop_cnt;i++){
        mult_mod+= multiply_uint_mod(operand1[i], operand2, modulus);
    }
    auto mult_mod_barrett_end = chrono::high_resolution_clock::now();
    cout << mult_mod << endl;
    
    auto barrett_time      = chrono::duration_cast<unit_mod_bench> (mod_barrett_end      - mod_barrett_begin);
    auto mod_2k_time       = chrono::duration_cast<unit_mod_bench> (mod_2k_end           - mod_2k_begin);
    auto mult_barrett_time = chrono::duration_cast<unit_mult_bench>(mult_mod_barrett_end - mult_mod_barrett_begin);
    auto mult_mod_2k_time  = chrono::duration_cast<unit_mult_bench>(mult_mod_2k_end      - mult_mod_2k_begin);
    cout << "mod_barrett(64bit) " << barrett_time.count()/static_cast<double>(loop_cnt) << NS << endl;
    cout << "mod 2k:            " << mod_2k_time.count()/static_cast<double>(loop_cnt) << NS << endl;
    cout << "mult_mod_barrett:  " << mult_barrett_time.count()/static_cast<double>(loop_cnt) << NS << endl;
    cout << "mult_mod 2k:       " << mult_mod_2k_time.count()/static_cast<double>(loop_cnt)  << NS << endl;
}


int main(int argc, char* argv[]){
    if(argc < 2){
        cout << "Usage: ./exe loop_count" << endl;
        return 0;
    }
    uint64_t loop_cnt = atoll(argv[1]);
    bench_mod(loop_cnt);
}
