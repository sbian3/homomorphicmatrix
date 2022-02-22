import numpy as np
from enum import Enum

result_root_path="../result/compare"

class OpeType(Enum):
    Enc="encrypt:"
    LT="LinearTransformation:"
    Dec="Decryption:"

class BenchType(Enum):
    directconv="direct_conv"
    packedconv="packed_conv"

def read_from_file(opetype, filepath):
    time = 0
    iterates = 0
    for line in open(filepath):
        line = line.rstrip()
        line = line.split(" ")
        prefix = opetype.value
        if(line[0] == prefix):
            time += float(line[1])
            timeunit = line[2]
            iterates =iterates + 1
    # calc average
    time = time / iterates
    return time, timeunit

class ParamSet:

    def __init__(self, n, input_dim, kernel_dim, pack_num = 0):
        self.n = n
        self.input_dim = input_dim
        self.kernel_dim = kernel_dim
        if(pack_num == 0):
            self.benchtype = BenchType.directconv
        else:
            self.benchtype = BenchType.packedconv
            self.pack_num  = pack_num

slash = "/"
# rootpath do not end "/" ex:./hoge
def get_filepath(rootpath, parms):
    filepath = rootpath + slash + parms.benchtype.value + slash
    filepath += str(parms.input_dim) + "_" + str(parms.kernel_dim)+ ".txt"
    return filepath

class HomResult:

    def __init__(self, parms):
        self.parms = parms

    def read_bench(self, rootpath):
        filepath = get_filepath(rootpath, self.parms)
        self.enc_time, unit = read_from_file(OpeType.Enc, filepath)
        self.lt_time, _  = read_from_file(OpeType.LT, filepath)
        self.dec_time, _ = read_from_file(OpeType.Dec, filepath)
        self.timeunit = unit

    def multiply(self, scalar):
        self.enc_time *= scalar
        self.lt_time  *= scalar
        self.dec_time *= scalar

    def total(self):
        return self.enc_time + self.lt_time + self.dec_time

    def __str__(self):
        return "type: {}, enc: {} {}, lt: {} {}, dec: {} {}, total: {} {}".format(self.parms.benchtype.value, self.enc_time, self.timeunit, self.lt_time, self.timeunit, self.dec_time, self.timeunit, self.total(), self.timeunit)

def get_max_packing(n, input_dim, kernel_dim):
    return int(n/(input_dim + kernel_dim - 1))

def compare_direct_pack(n, input_dim, kernel_dim, rootpath):
    parms_dc = ParamSet(n, input_dim, kernel_dim)
    direct_result = HomResult(parms_dc)
    direct_result.read_bench(rootpath)
    pack_num = get_max_packing(n, input_dim, kernel_dim)
    parms_pack = ParamSet(n, input_dim, kernel_dim, pack_num)
    packed_result = HomResult(parms_pack)
    packed_result.read_bench(rootpath)
    print(direct_result)
    direct_result.multiply(pack_num)
    total_ratio = packed_result.total() / direct_result.total()
    print(direct_result)
    print(packed_result)
    print("total time ratio: {}".format(round(total_ratio, 3)))
    print("improvement     : {}%".format(round(float(1- total_ratio), 3)*100))
    
def compare_all_parms():
    inputs = [25, 64, 256, 784]
    kernels = [1, 9, 25]
    for input_dim in inputs:
        for kernel in kernels:
            print("compare input={}, kernel={}".format(input_dim, kernel))
            compare_direct_pack(1024, input_dim, kernel, result_root_path)

compare_all_parms()
