# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 11:57:24 2018

@author: colorbox
"""

import os
import numpy as np

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


directory=r'C:\Users\colorboxy\Documents\Github\episodic_sim_working_files\Current\multi_bump_multi_recall_mods\PRversions\out_for_cnnect'
fire_len=[]
file_list=os.listdir(directory)
out_f_list=[]
max_out=[]
for filen in file_list:
    if not('dots.csv' in filen):
        continue
    chunks=filen.split('-')
    gbarsyn = chunks[0]
    numba = chunks[1][:-8]
    out_f_list.append([int(numba),float(gbarsyn)])
    file_full=os.path.join(directory,filen)
    z=np.genfromtxt(file_full,delimiter=',')
    start=100
    try:
        z=z[(z[:,0]>100),:]
        fina=np.max(z[:,0])
        fire_len.append(len(z))  
    except:
        fina=0
        fire_len.append(0)  
       
    max_out.append(fina)
    print(filen)
out_arr=np.array(out_f_list)
fire_linarr=np.array(fire_len)
max_out=np.array(max_out)
            
    