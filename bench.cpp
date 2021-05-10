#include "bench.h"


void print_plain(Plaintext plain, uint64_t num){
    if(num > plain.coeff_count()){
        std::cout << "print_plain: warn: num is bigger than plain coefficient size!" << std::endl;
        num = plain.coeff_count();
    }
    std::cout << "Plaintext Coefficients(first " << num << " elements): [";
    for(auto i = 0U;i < num;i++){
        std::cout << plain[i] << " ";
    }
    std::cout << "]" << std::endl;
}
