import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum

poly_degrees = [ 1024, 2048, 4096, 8192]
#poly_degrees = [ 1024, 2048, 4096, 8192, 16384, 32768]
kernel_d=4
kernel_dims = [2, 4, 8, 16, 32, 64, 128]
#input_dims = [16, 64, 256, 1024]
input_dims = [16, 32, 64, 128, 256, 512, 1024]

# for order reference
def generate_order_reference(poly_degrees):
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
    return ref_nsq, ref_n15, ref_n1

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
            # current time unit is [ms]
            if(line[2] == TimeUnit.ms.value):
                time_lt = float(line[1])
                timeunit = TimeUnit.ms
            elif(line[2] == TimeUnit.us.value):
                # if time unit is [us] in result txt, devide by 1000
                time_lt = float(line[1])/1000
                timeunit = TimeUnit.us
        if(line[0] == dec_prefix):
            if(line[2] == TimeUnit.ms.value):
                time_dec = float(line[1])
                timeunit = TimeUnit.ms
            elif(line[2] == TimeUnit.us.value):
                time_dec = float(line[1])/1000
                timeunit = TimeUnit.us
    return time_lt, time_dec, timeunit

def join_path(result_dirname, algorithm_name, benchtype, param_bench):
    slash = "/"
    filepath = "../" + result_dirname + slash + algorithm_name + slash + benchtype.value + slash + str(param_bench) + ".txt"
    return filepath

class HltData():
    def __init__(self, lt_time, dec_time, timeunit):
        self.lt_time = lt_time
        self.dec_time = dec_time
        self.timeunit = timeunit

# read bench data from ../hlt_result
def read_bench(algorithm_name, benchtype, kernel_dim=0):
    lt_times = []
    dec_times = []
    result_dirname = "../hlt_result"
    if(benchtype == BenchType.multi_poly):
        for poly_degree in poly_degrees:
            filepath = result_dirname + "/" + algorithm_name +   "/"+ benchtype.value+ "/" + "k" + str(kernel_dim) + "/" +  algorithm_name + str(poly_degree) + ".txt"
            time_lt, time_dec, timeunit = read_latency(filepath)
            lt_times.append(time_lt)
            dec_times.append(time_dec)
    elif(benchtype == BenchType.multi_kernel):
        for kernel_dim in kernel_dims:
            filepath = result_dirname + "/" + algorithm_name +   "/"+ benchtype.value + "/" + algorithm_name + str(kernel_dim) + ".txt"
            time_lt, time_dec, timeunit = read_latency(filepath)
            lt_times.append(time_lt)
            dec_times.append(time_dec)
    elif(benchtype == BenchType.multi_input):
        for input_dim in input_dims:
            filepath = result_dirname + "/" + algorithm_name +   "/"+ benchtype.value + "/" + algorithm_name + str(input_dim) + ".txt"
            time_lt, time_dec, timeunit = read_latency(filepath)
            lt_times.append(time_lt)
            dec_times.append(time_dec)
    hltdata = HltData(lt_times, dec_times, timeunit)
    return hltdata

def read_data(benchtype, kernel_dim=0):
    data_direct = read_bench("direct_conv", benchtype, kernel_dim)
    data_packed = read_bench("packed_conv", benchtype, kernel_dim)
    data_pcdec = read_bench("pc_toeplitz", benchtype, kernel_dim)
    # time_general is too long: put result directly
    time_general=[52570, 462696, 3234320]
    return data_direct, data_packed, time_general, data_pcdec

def plot_lt(x, x_label, data_direct, data_packed, time_general, data_pcdec, benchtype, kernel_dim=0):
    fig, ax = plt.subplots()
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=10)
    # direct conv
    #ax.plot(x, data_direct.lt_time, label="direct", color="red")
    # packed conv
    ax.plot(x, data_packed.lt_time, label="optimized", color="orange")
    # general_LT
    ax.plot(x[:3], time_general, label="general_linear", color="blue")
    #ax.plot(x, ref_nsq, label="nsquare(ref)", color="green")
    #ax.plot(x, ref_n15, label="n1.5(ref)", color="purple")
    #ax.plot(x, ref_n1, label="n1(ref)", color="magenta")
    ax.legend(loc="upper left")
    ax.set_xlabel(x_label)
    ax.set_xticks(x)
    ax.set_ylabel("latency [ms]")
    fig_dirname = "./fig"
    fig.savefig(fig_dirname + "/" + "lt" + benchtype.value + str(kernel_dim) + ".pdf", format="pdf", bbox_inches="tight", pad_inches=0.1)

def plot_dec(x, x_label, data_direct, data_packed, time_general, data_pcdec, benchtype, kernel_dim=0):
    fig, ax = plt.subplots()
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=10)
    # direct conv
    #ax.plot(x, data_direct.lt_time, label="direct", color="red")
    # packed conv
    ax.plot(x, data_packed.dec_time, label="general_linear", color="blue")
    # toeplitz dec
    ax.plot(x, data_pcdec.dec_time, label="optimized", color="orange")
    ax.legend(loc="upper left")
    ax.set_xlabel(x_label)
    ax.set_xticks(x)
    ax.set_ylabel("latency [ms]")
    fig_dirname = "./fig"
    fig.savefig(fig_dirname + "/" + "dec" + benchtype.value + str(kernel_dim) + ".pdf", format="pdf", bbox_inches="tight", pad_inches=0.1)

def plot_tifs():
    data_direct ,data_packed, time_general, data_pcdec = read_data(BenchType.multi_poly, kernel_d)
    plot_lt(poly_degrees, "polynomial degree", data_direct, data_packed, time_general, data_pcdec, BenchType.multi_poly, kernel_d)
    plot_dec(poly_degrees, "polynomial degree", data_direct, data_packed, time_general, data_pcdec, BenchType.multi_poly, kernel_d)
    #time_direct,time_packed, time_general = read_data(BenchType.multi_kernel)
    #plot(kernel_dims, "kernel dimensions", time_direct, time_packed, time_general, BenchType.multi_kernel)
    #time_direct,time_packed, time_general = read_data(BenchType.multi_input)
    #plot(input_dims, "input dimensions", time_direct, time_packed, time_general, BenchType.multi_input)

def main():
    plot_tifs()

main()
