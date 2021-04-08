#!/bin/bash

input_n=(4 8 16 32 64 128)
input_dims=()
kernel_dim=4
poly_degrees=(1024 2048 4096 8192 16384 32768)
poly_degree=16384
pack_num=2
inputdim_multipoly=16
err=tifs_result/err.log

# calculate input dimensions (n^2)
for n in ${input_n[@]}
do
    nn=$(($n ** 2))
    input_dims=("${input_dims[@]}" $nn)
done

# assume executable file is build/bin/direct_conv
function directconv_multiinput(){
    mkdir -p tifs_result/direct_conv/multiinput
    for input_dim in ${input_dims[@]}
    do
        echo "Direct Convolution: (input, kernel) = ($input_dim, $kernel_dim)"
        output_path="tifs_result/direct_conv/multiinput/direct$input_dim.txt"
        env build/bin/direct_conv $input_dim $kernel_dim $poly_degree 2>$err  > $output_path
    done
}

function directconv_multipoly(){
    mkdir -p tifs_result/direct_conv/multipoly
    input_dim=$inputdim_multipoly
    for poly_degree in ${poly_degrees[@]}
    do
        echo "Direct Convolution: (input, kernel) = ($input_dim, $kernel_dim)"
        output_path="tifs_result/direct_conv/multipoly/direct$poly_degree.txt"
        env build/bin/direct_conv $input_dim $kernel_dim $poly_degree 2>$err  > $output_path
    done
}

function packedconv_multiinput(){
    mkdir -p tifs_result/packed_conv/multiinput
    for input_dim in ${input_dims[@]}
    do
        echo "Packed Convolution: (input, kernel, pack_num) = ($input_dim, $kernel_dim, $pack_num)"
        echo "Polynomial Length: $poly_degree"
        output_path="tifs_result/packed_conv/multiinput/pack$input_dim.txt"
        env build/bin/packed_conv $input_dim $kernel_dim $poly_degree $pack_num 2>$err > $output_path
    done
}

function packedconv_multipoly(){
    mkdir -p tifs_result/packed_conv/multipoly
    input_dim=$inputdim_multipoly
    for poly_degree in ${poly_degrees[@]}
    do
        echo "Packed Convolution: (input, kernel, pack_num) = ($input_dim, $kernel_dim, $pack_num)"
        echo "Polynomial Length: $poly_degree"
        output_path="tifs_result/packed_conv/multipoly/pack$poly_degree.txt"
        env build/bin/packed_conv $input_dim $kernel_dim $poly_degree $pack_num 2>$err > $output_path
    done
}

function general_lt_multiinput(){
    mkdir -p tifs_result/general_lt/multiinput
    # In general_lt, poly_degree 4096 or longer is too long to measure
    poly_degree=1024
    for input_dim in ${input_dims[@]}
    do
        echo "General Linear Transformation: ()"
        output_path="tifs_result/general_lt/multiinput/general_lt$input_dim.txt"
        env build/bin/general_lt $input_dim $kernel_dim $poly_degree 2>$err > $output_path
    done
}

function clean(){
    rm -rf tifs_result/*
}

directconv_multipoly
packedconv_multipoly
#general_lt
