# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 09:57:42 2017

@author: colorbox
"""
import itertools
import os
directory=r'C:\Users\colorboxy\Documents\Github\PRversions'
filename='recurrent_fire.py'
outdir=r'C:\Users\colorboxy\Documents\Github\PRversions\meta_temp'
with open(os.path.join(directory,filename),'r') as f:
    zed=f.readlines()

inhib_level=[0,.25,.5,.75,1]
on_off=[0,1]
targets_lists=[on_off]
all_target_combo=list(itertools.product(inhib_level,on_off))
names_list=['inhib_level','on_off']
i=0
for jjj in range(100):
    for combos in all_target_combo:
        nomen=r'outname="'
        for nn,name in enumerate(names_list):
            for ln,line in enumerate(zed):
                if not('=' in line):
                    continue
                parse=line.split('=')
                after=line.find('*')
                if name == parse[0]:
                    print('trig')
                    if after>0:
                        line=name+'='+str(combos[nn])+line[after:]
                    else:
                        line=line=name+'='+str(combos[nn])+'\n'
                    zed[ln]=line
            nomen=nomen+'-'+name+'-'+str(combos[nn])+'-'+str(jjj)
        zed[432]=nomen+str(i)+r'"'
        i+=1
        with open(os.path.join(outdir,str(i)+'.py'),'w') as f:
            f.writelines(zed)
        
                
    
    
