#!/usr/bin/env python

"""
This is a copy of the yields.py from runs-thc-ba
"""

import h5py
import numpy as np
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def normalize_yields(As, Ys):
    Anrm = np.arange(As.max()+1)
    Ynrm = np.zeros(int(As.max())+1)
    for i in range(Ynrm.shape[0]):
        Ynrm[i] = Ys[As == i].sum()
    return Anrm, Ynrm

if len(sys.argv) < 4:
    print("Usage: {} /path/to/tabulated_nucsyn.h5 " +\
        "/path/to/solar_r.dat /path/to/ejecta.h5".format(sys.argv[0]))
    exit(1)
tab_nucsyn = h5py.File(sys.argv[1], "r")
Asun, Ysun = np.loadtxt(sys.argv[2], unpack=True)
Asun, Ysun = normalize_yields(Asun, Ysun)
Ysun /= np.sum(Ysun)
tab_ejecta = h5py.File(sys.argv[3], "r")

mass = np.array(tab_ejecta["mass"])
Ye   = np.array(tab_ejecta["Ye"])
Ys   = np.array(tab_nucsyn["Y_final"])
As   = np.array(tab_nucsyn["A"])
Zs   = np.array(tab_nucsyn["Z"])

yields = []
for i in range(4):
    yields.append(np.zeros(Ys.shape[-1]))
for i in range(yields[0].shape[0]):
    yields[0][i] = np.sum(mass * Ys[:,:,:,i])
    for j, yeb in enumerate([(0.0, 0.25), (0.25, 0.35), (0.35, 0.5)]):
        idx = (Ye >= yeb[0]) & (Ye < yeb[1])
        yields[j+1][i] = np.sum(mass[idx,:,:] * Ys[idx,:,:,i])

outfile = h5py.File("yields.h5", "w")
outfile.create_dataset("/Y_final", data=yields[0])
outfile.create_dataset("/A", data=As)
outfile.create_dataset("/Z", data=Zs)

Ynrm = []
for i in range(4):
    Anrm, Ynrm_i = normalize_yields(As, yields[i])
    Ynrm.append(Ynrm_i)
norm = Ynrm[0].sum()
for i in range(4):
    Ynrm[i] /= norm

plt.plot(Asun, Ysun, 'k.', label='Solar')
plt.plot(Anrm, Ynrm[0], label='All')
plt.plot(Anrm, Ynrm[1], label='0 <= Ye < 0.25')
plt.plot(Anrm, Ynrm[2], label='0.25 <= Ye < 0.35')
plt.plot(Anrm, Ynrm[3], label='0.35 <= Ye < 0.5')
plt.yscale("log")
plt.ylim(ymin=1e-9, ymax=2e-1)
plt.xlim(xmin=0, xmax=275)
plt.ylabel("Relative final abundances")
plt.xlabel("A")
plt.legend(loc='best', numpoints=1)
plt.savefig("yields.png")
