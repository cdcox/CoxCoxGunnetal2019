# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 10:56:19 2018

@author: colorboxy
"""

import csv
import numpy as np
import os
in_dir = r'C:\Users\colorboxy\Documents\aliza_Training'
out_dir= r'C:\Users\colorboxy\Documents\aliza_Training\csvs'
f_list = os.listdir(in_dir)
for file_name in f_list:
    with open(os.path.join(in_dir,file_name),'r') as f:
        output=f.readlines()
    clean = []
    temp = []
    for i in output:
        i = i[:-1]
        if i==r'':
            clean.append(temp)
            temp=[]
        else:
            temp.append(float(i))
    file_name=file_name.replace(' ','_')
    file_name= file_name[:-4]
    for lnn,lists in enumerate(clean):
        lnn2=lnn+1
        out_name=os.path.join(out_dir,file_name+'_'+f"{lnn2:03}"+'.csv')
        lists=np.array(lists)
        np.savetxt(out_name,lists,delimiter=' ')
    print(file_name)