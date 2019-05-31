# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 09:57:42 2017

@author: colorbox
"""
import itertools
import os
import numpy as np
directory=r'C:\Users\colorboxy\Documents\Github\PRversions'
filename='recurrent_fire.py'
outdir=r'C:\Users\colorboxy\Documents\Github\episodic_sim_working_files\Current\multi_bump_multi_recall_mods\PRversions\out_for_cnnect'
with open(os.path.join(directory,filename),'r') as f:
    zed=f.readlines()

gbarsyn=np.arange(.2,1.1,.05)
targets_lists=gbarsyn
all_target_combo=list(itertools.product(gbarsyn))
names_list=['gbarsyn']
i=0
for zz in range(100):
    for combos in all_target_combo:
        nomen=r'outname="'
        for nn,name in enumerate(names_list):
            for ln,line in enumerate(zed):
                if not('=' in line):
                    continue
                parse=line.split('=')
                if name == parse[0]:
                    print('trig')
                    line=name+'='+str(np.round(combos[nn], 2))+'*nS\n'
                    zed[ln]=line
            nomen=nomen+str(np.round(combos[nn],2))
        zed[432]=nomen+'-'+str(zz)+r'"'
        i+=1
        with open(os.path.join(outdir,str(i)+'.py'),'w') as f:
            f.writelines(zed)
        
    
    
