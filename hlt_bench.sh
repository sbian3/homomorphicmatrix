#!/bin/bash

input_n=(4 8 16 32 64 128)
input_dims=()
kernel_dim=4
kernel_dims=(2 4 8 16 32 64 128)
poly_degrees=(1024 2048 4096 8192)
poly_degree=2048
poly_degree_multikernel=1024
pack_num=1
inputdim_multipoly=16
dirname_result_root=result
err=$dirname_result_root/err.log

# calculate input dimensions (n^2)
for n in ${input_n[@]}
do
    nn=$(($n ** 2))
    input_dims=("${input_dims[@]}" $nn)
done

# assume executable file is build/bin/direct_conv
function directconv_multiinput(){
    dirname=$dirname_result_root/direct_conv/multiinput/k$kernel_dim
    mkdir -p $dirname
    for input_dim in ${input_dims[@]}
    do
        echo "Direct Convolution: (input, kernel) = ($input_dim, $kernel_dim)"
        output_path="$dirname/direct_conv$input_dim.txt"
        env build/bin/direct_conv $input_dim $kernel_dim $poly_degree 2>$err  > $output_path
    done
}

function directconv_multipoly(){
    dirname=$dirname_result_root/direct_conv/multipoly/k$kernel_dim
    mkdir -p $dirname
    input_dim=$inputdim_multipoly
    for poly_degree in ${poly_degrees[@]}
    do
        echo "Direct Convolution: (input, kernel, poly_degree) = ($input_dim, $kernel_dim, $poly_degree)"
        output_path="$dirname/direct_conv$poly_degree.txt"
        env build/bin/direct_conv $input_dim $kernel_dim $poly_degree 2>$err  > $output_path
    done
}

function directconv_multikernel(){
    mkdir -p $dirname_result_root/direct_conv/multikernel
    input_dim=$inputdim_multipoly
    poly_degree=$poly_degree_multikernel
    for kernel_dim in ${kernel_dims[@]}
    do
        echo "Direct Convolution: (input, kernel, poly_degree) = ($input_dim, $kernel_dim, $poly_degree)"
        output_path="$dirname/direct_conv$kernel_dim.txt"
        env build/bin/direct_conv $input_dim $kernel_dim $poly_degree 2>$err  > $output_path
    done
}

function packedconv_multiinput(){
    dirname=$dirname_result_root/packed_conv/multiinput
    mkdir -p $dirname
    for input_dim in ${input_dims[@]}
    do
        echo "Packed Convolution: (input, kernel, pack_num) = ($input_dim, $kernel_dim, $pack_num)"
        echo "Polynomial Length: $poly_degree"
        output_path="$dirname/packed_conv$input_dim.txt"
        env build/bin/packed_conv $input_dim $kernel_dim $poly_degree $pack_num 2>$err > $output_path
    done
}

function packedconv_multipoly(){
    dirname=$dirname_result_root/packed_conv/multipoly/k$kernel_dim
    mkdir -p $dirname
    input_dim=$inputdim_multipoly
    for poly_degree in ${poly_degrees[@]}
    do
        echo "Packed Convolution: (input, kernel, poly_degree, pack_num) = ($input_dim, $kernel_dim, $poly_degree,  $pack_num)"
        output_path="$dirname/packed_conv$poly_degree.txt"
        env build/bin/packed_conv $input_dim $kernel_dim $poly_degree $pack_num 2>$err > $output_path
    done
}

function packedconv_multikernel(){
    mkdir -p $dirname_result_root/packed_conv/multikernel
    input_dim=$inputdim_multipoly
    poly_degree=$poly_degree_multikernel
    for kernel_dim in ${kernel_dims[@]}
    do
        echo "Packed Convolution: (input, kernel, poly_degree, pack_num) = ($input_dim, $kernel_dim, $poly_degree,  $pack_num)"
        output_path="$dirname_result_root/packed_conv/multikernel/packed_conv$kernel_dim.txt"
        env build/bin/packed_conv $input_dim $kernel_dim $poly_degree $pack_num 2>$err > $output_path
    done
}

function pc_toeplitz_multipoly(){
    dirname=$dirname_result_root/pc_toeplitz/multipoly/k$kernel_dim
    mkdir -p $dirname
    input_dim=$inputdim_multipoly
    for poly_degree in ${poly_degrees[@]}
    do
        echo "Packed Convolution(toeplitz multiplication): (input, kernel, poly_degree, pack_num) = ($input_dim, $kernel_dim, $poly_degree,  $pack_num)"
        output_path="$dirname/pc_toeplitz$poly_degree.txt"
        env build/bin/pc_toeplitz $input_dim $kernel_dim $poly_degree $pack_num 2>$err > $output_path
    done
}

function general_lt_multiinput(){
    mkdir -p $dirname_result_root/general_lt/multiinput
    # In general_lt, poly_degree 4096 or longer is too long to measure
    poly_degree=1024
    for input_dim in ${input_dims[@]}
    do
        echo "General Linear Transformation: ()"
        output_path="$dirname_result_root/general_lt/multiinput/general_lt$input_dim.txt"
        env build/bin/general_lt $input_dim $kernel_dim $poly_degree 2>$err > $output_path
    done
}

function general_lt_multipoly(){
    mkdir -p $dirname_result_root/general_lt/multipoly
    input_dim=$inputdim_multipoly
    for poly_degree in ${poly_degrees[@]}
    do
        echo "General Linear Transformation: (input, kernel, poly_degree) = ($input_dim, $kernel_dim, $poly_degree)"
        output_path="$dirname_result_root/general_lt/multipoly/general_lt$poly_degree.txt"
        env build/bin/general_lt $input_dim $kernel_dim $poly_degree 2>$err  > $output_path
    done
}

function clean(){
    rm -rf $dirname_result_root/*
}

#clean
directconv_multipoly
#packedconv_multipoly
pc_toeplitz_multipoly
#directconv_multiinput
#packedconv_multiinput
#directconv_multikernel
#packedconv_multikernel
#general_lt_multipoly
