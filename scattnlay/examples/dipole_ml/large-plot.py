#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

wl = 0.455
NL = 3


def read_data_from_files(file_start):
    from fs import open_fs
    from fs.walk import Walker
    home_fs = open_fs('./')
    walker = Walker(filter=[file_start+'*.txt'])
    data = []
    for path in walker.files(home_fs):
        with open('.'+path) as f:
            A=np.loadtxt((x.replace('[',' ').replace(']',' ').replace(',',' ') for x in f))
        # print(A)
        data.append(A)
        # else: data_list.append(np.loadtxt(dir_name+"/"+name, skiprows=9))
    #     pbar.update()
    # pbar.close()

    return np.array(data)


def set_pbar(total, desc):
    return tqdm(total=total, 
                desc=desc, 
                file=sys.stdout,
                bar_format='{l_bar}{bar}| elapsed: {elapsed} remaining: {remaining}')

def plot_data_all():
    plt.figure("all")
    for d in data:
        if len(d.shape)==1: continue
        max_iter = np.max(d[:,0])
        d = d[d[:,0]>max_iter-1]
        x = d[:,1]
        xm = np.max(d[:,4:4+NL], axis=1)/wl
        Dmax = d[:,2]#/xm
        print(max_iter)
        if max_iter==5000:
            plt.scatter(xm, Dmax, marker='x', s=34)
        else:
            plt.scatter(xm, Dmax, marker='o', s=24)
    plt.title(file_start)
    plt.xlabel("R_max/WL"   )
    plt.ylabel("Directivity")# * WL/R_max")

# file_start='out_index040-N50-iterations90'
#file_start='out_index040-N50-iterations50'
#file_start='out2_index040-N50-iterations6000'
#file_start='out2_index040-N50-NL2-iterations6000'
file_start='out2_index030-N50-NL3-iterations6000'
data = read_data_from_files(file_start)

plot_data_all()

d0 = data[0]
max_iter = np.max(d0[:,0])
d0 = d0[d0[:,0]>max_iter-1]
DDmax = d0[:,2]
for d in data:
    if len(d.shape)==1: continue
    max_iter = np.max(d[:,0])
    d = d[d[:,0]>max_iter-1]
    if len(d) != len(d0): continue
    Dmax = d[:,2]
    for i in range(len(d)):
        if DDmax[i] > Dmax[i]: continue
        DDmax[i] = Dmax[i]
        d0[i] = d[i]
        d0[i,1] = np.max(d0[i,4:4+NL])
        d0[i,3] = d0[i,3]/d0[i,1]
        d0[i,4:4+NL] =np.sort(d0[i,4:4+NL])/d0[i,1]
for i in range(len(d0)):
    if NL==3:
        print("%i, R_max=%5.5f, D_max=%3.1f \tdipole/R_max=%.5f \tR/R_max=[%8.7f, %8.7f, %1.1f]\tn=[% 3.1f, % 3.1f, % 3.1f]"%(tuple(d0[i])))
    if NL==2:
        print("%i, R_max=%5.5f, D_max=%3.1f \tdipole/R_max=%.5f \tR/R_max=[%8.7f, %1.1f]\tn=[% 3.1f, % 3.1f]"%(tuple(d0[i])))
#    print("%i, R_max=%16.15g, D_max=%3.1f \tdipole/R_max=%16.15g \tR/R_max=[%16.15g, %16.15g, %1.1f]\tn=[%16.15g, %16.15g, %16.15g]"%(tuple(d0[i])))
# xm = np.max(d0[:,4:7], axis=1)/wl

# plt.figure('max')
# plt.scatter(xm, DDmax)
plt.show()
