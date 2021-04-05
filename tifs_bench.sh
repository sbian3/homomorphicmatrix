#!/bin/bash

input_n=(4 8 16 32 64)
input_dims=()
kernel_dim=4
err=tifs_result/err.log

# calculate input dimensions (n^2)
for n in ${input_n[@]}
do
    echo $n
    nn=$(($n ** 2))
    input_dims=("${input_dims[@]}" $nn)
done

# assume executable file is build/bin/direct_conv
function direct_conv(){
    mkdir -p tifs_result/
    for input_dim in ${input_dims[@]}
    do
        echo "Direct Convolution: (input, kernel) = ($input_dim, $kernel_dim)"
        output_path="tifs_result/direct_conv$input_dim.txt"
        env build/bin/direct_conv $input_dim $kernel_dim 2>$err # > $output_path
    done
}

direct_conv
