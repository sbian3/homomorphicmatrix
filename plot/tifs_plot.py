import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum

poly_degrees = [ 1024, 2048, 4096, 8192, 16384]
#poly_degrees = [ 1024, 2048, 4096, 8192, 16384, 32768]
kernel_dims = [2, 4, 8, 16, 32, 64, 128]
input_dims = [16, 64, 256, 1024]

# for order reference
ref_nsq = []
ref_scale=10/(2048*2048)
ref_n15 = []
ref_n1 = []
scale_n15 = 10/pow(2048, 1.5)
scale_n1 = 10/2048
for poly_degree in poly_degrees:
    ref_nsq.append(poly_degree * poly_degree * ref_scale)
    ref_n15.append(pow(poly_degree, 1.5) * scale_n15)
    ref_n1.append(poly_degree * scale_n1)


class BenchType(Enum):
    multi_poly="multipoly"
    multi_kernel="multikernel"
    multi_input="multiinput"

class TimeUnit(Enum):
    ms="ms"
    us="us"
    ns="ns"

lt_prefix = "LinearTransformation:"
dec_prefix = "Decryption:"

def read_latency(filepath):
    for line in open(filepath):
        line = line.rstrip()
        line = line.split(" ")
        if(line[0] == lt_prefix):
            # current unit is [ms]
            if(line[2] == TimeUnit.ms.value):
                time_lt = float(line[1])
            elif(line[2] == TimeUnit.us.value):
                # if [us] in result txt, devide by 1000
                time_lt = float(line[1])/1000
        if(line[0] == dec_prefix):
            time_dec = float(line[1])
    return time_lt, time_dec

# read bench data from /home/work/kazumasita/relu/SEAL/tifs_result
def read_bench(prefix, benchtype):
    bench_lt = []
    if(benchtype == BenchType.multi_poly):
        for poly_degree in poly_degrees:
            filepath = "../tifs_result/" + prefix +   "/"+ benchtype.value+ "/" + prefix + str(poly_degree) + ".txt"
            time_lt, time_dec = read_latency(filepath)
            bench_lt.append(time_lt)
    elif(benchtype == BenchType.multi_kernel):
        for kernel_dim in kernel_dims:
            filepath = "../tifs_result/" + prefix +   "/"+ benchtype.value + "/" + prefix + str(kernel_dim) + ".txt"
            time_lt, time_dec = read_latency(filepath)
            bench_lt.append(time_lt)
    elif(benchtype == BenchType.multi_input):
        for input_dim in input_dims:
            filepath = "../tifs_result/" + prefix +   "/"+ benchtype.value + "/" + prefix + str(input_dim) + ".txt"
            time_lt, time_dec = read_latency(filepath)
            bench_lt.append(time_lt)
    return bench_lt

def read_data(benchtype):
    time_direct = read_bench("direct_conv", benchtype)
    time_direct_f = time_direct
    time_packed = read_bench("packed_conv", benchtype)
    # time_general is too long: put result directly
    time_general=[52570, 462696, 3234320]
    return time_direct_f,time_packed, time_general

def plot(x, x_label, time_direct, time_packed, time_general, benchtype):
    fig, ax = plt.subplots()
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=10)
    # direct conv
    ax.plot(x, time_direct, label="direct", color="red")
    # packed conv
    ax.plot(x, time_packed, label="packed conv", color="orange")
    # general_LT
    ax.plot(x[:3], time_general, label="general", color="blue")
    #ax.plot(x, ref_nsq, label="nsquare(ref)", color="green")
    #ax.plot(x, ref_n15, label="n1.5(ref)", color="purple")
    #ax.plot(x, ref_n1, label="n1(ref)", color="magenta")
    ax.legend()
    ax.set_xlabel(x_label)
    ax.set_xticks(x)
    ax.set_ylabel("latency [ms]")
    fig.savefig(benchtype.value + ".pdf", format="pdf")

def plot_tifs():
    time_direct,time_packed, time_general = read_data(BenchType.multi_poly)
    plot(poly_degrees, "polynomial degree", time_direct, time_packed, time_general, BenchType.multi_poly)
    #time_direct,time_packed, time_general = read_data(BenchType.multi_kernel)
    #plot(kernel_dims, "kernel dimensions", time_direct, time_packed, time_general, BenchType.multi_kernel)
    #time_direct,time_packed, time_general = read_data(BenchType.multi_input)
    #plot(input_dims, "input dimensions", time_direct, time_packed, time_general, BenchType.multi_input)

def main():
    plot_tifs()

main()
