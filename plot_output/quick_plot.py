#!/usr/env python2

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    index = int(sys.argv[1])
except:
    index = 1
try:
    plot_all = bool(sys.argv[2])
except:
    plot_all = False

c = np.loadtxt("Chaste_control.txt")
d = np.loadtxt("Chaste_drug.txt")
#d = np.loadtxt("Chaste_drug_1000beat.txt")

if not plot_all:
    plt.plot(c[:,0], c[:,index], label='control')
    plt.plot(d[:,0], d[:,index], label='drug')
    plt.xlabel("time [ms]")
    plt.ylabel("state variable %d"%index)
    plt.legend()
    plt.savefig("quick_plot_variable_%d.png"%index)
else:
    fig, axes = plt.subplots(11, sharex=True, figsize=(6,15))
    axes[0].plot(c[:,0], c[:,1], label='control')
    axes[0].plot(d[:,0], d[:,1], label='drug')
    axes[0].set_ylabel("V [mV]")
    for i in xrange(10):
        axes[i+1].plot(c[:,0], c[:,35+i], label='control')
        axes[i+1].plot(d[:,0], d[:,35+i], label='drug')
        axes[i+1].set_ylabel("IKr variable %d"%i)
    axes[i+1].set_xlabel("time [ms]")
    plt.legend()
    plt.savefig("quick_plot_IKr.png", bbox_inch='tight')

