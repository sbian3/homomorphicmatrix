#!/bin/bash

input_n=(4 8 16 32 64 128)
input_dims=()
kernel_dim=4
poly_degrees=(1024 2048 4096 8192 16384 32768)
poly_degree=32768
pack_num=2
err=tifs_result/err.log

# calculate input dimensions (n^2)
for n in ${input_n[@]}
do
    nn=$(($n ** 2))
    input_dims=("${input_dims[@]}" $nn)
done

# assume executable file is build/bin/direct_conv
function direct_conv(){
    mkdir -p tifs_result/direct_conv
    for input_dim in ${input_dims[@]}
    do
        echo "Direct Convolution: (input, kernel) = ($input_dim, $kernel_dim)"
        output_path="tifs_result/direct_conv/direct$input_dim.txt"
        env build/bin/direct_conv $input_dim $kernel_dim $poly_degree 2>$err  > $output_path
    done
}

function packed_conv(){
    mkdir -p tifs_result/packed_conv
    for input_dim in ${input_dims[@]}
    do
        echo "Packed Convolution: (input, kernel, pack_num) = ($input_dim, $kernel_dim, $pack_num)"
        echo "Polynomial Length: $poly_degree"
        output_path="tifs_result/packed_conv/pack$input_dim.txt"
        env build/bin/packed_conv $input_dim $kernel_dim $poly_degree $pack_num 2>$err > $output_path
    done
}

function general_lt(){
    mkdir -p tifs_result/general_lt
    # In general_lt, poly_degree 4096 or longer is too long to measure
    poly_degree=1024
    for input_dim in ${input_dims[@]}
    do
        echo "General Linear Transformation: ()"
        output_path="tifs_result/general_lt/general_lt$input_dim.txt"
        env build/bin/general_lt $input_dim $kernel_dim $poly_degree 2>$err > $output_path
    done
}

#direct_conv
#packed_conv
general_lt
