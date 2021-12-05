
import sys
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum

class GraphColor(Enum):
    BLUE = 'b'
    WHITE = 'w'
    GREEN = 'g'
    RED = 'r'
    ORANGE = 'orange'
    YELLOW = 'y'
    BLACK = 'k'
    CYAN = 'c'
    MAGENTA = 'm'

class GraphLineStyle(Enum):
    SOLID = 'solid'
    DASHED = 'dashed'

class PlotInfo:

    def __init__(self, label="hoge", color=GraphColor.RED, linewidth=2, linestyle=GraphLineStyle.SOLID):
        self.label = label
        self.color = color.value
        self.linewidth = linewidth
        self.linestyle = linestyle.value

    def setLabel(self, label):
        self.label = label

    def setColor(self, color):
        self.color = color.value

    def setLineStyle(self, linestyle):
        self.linestyle =  linestyle.value

class Graph_base2:
    savefig_path = "./fig"

    def __init__(self):
        #self.fig = plt.figure(figsize=(16, 9))
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(111)
        self.axes.set_xscale("log", base=2)
        self.axes.set_yscale("log", base=10)

    def set_title(self, title):
        self.axes.set_title(title)

    def set_xlabel(self, xlabel):
        self.axes.set_xlabel(xlabel)

    def set_ylabel(self, ylabel):
        self.axes.set_ylabel(ylabel)

    def set_xlim(self, xlim):
        self.axes.set_xlim(xlim)

    def set_ylim(self, ylim):
        self.axes.set_ylim(ylim)

    def plot(self, x, y):
        self.axes.plot(x, y)
        self.axes.legend()

    def addplot_withinfo(self, x, y, plotInfo):
        self.axes.plot(x, y, label=plotInfo.label, color=plotInfo.color, linewidth=plotInfo.linewidth, linestyle=plotInfo.linestyle)
        self.axes.legend(loc="upper left")

    def show(self):
        plt.show()

    def save(self, name):
        self.fig.savefig(self.savefig_path + "/" + name, bbox_inches="tight", pad_inches=0.1)

def print_numarr(arr):
    for item in arr:
        print(str(item))
