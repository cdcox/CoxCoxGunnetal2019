
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 15:59:46 2016

@author: colorbox
"""
from __future__ import division
from brian2 import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import csv
import random as rn
plt.switch_backend('agg')

def make_network(cores,number_of_neurons,netsize,connect_break): 
    out_list=[]
    test_net=nx.newman_watts_strogatz_graph(1000,netsize,0)
    factor=5
    
    for nn in test_net.nodes():
        if np.min(np.abs(cores-nn))<15:
            for i in range(factor):
                rand_core=cores[int(np.round(np.random.random()*6))]
                test_net.add_edge(nn,np.random.randint(0,1000))
    for ed in test_net.edges():
        if np.random.random()>connect_break:
            continue
        if np.abs(ed[0]-ed[1])>20 and np.abs(ed[0]-ed[1])<700:
            rng=np.random.choice([-1,0,0,1])
        elif np.abs(ed[0]-ed[1])>700:
            rng=1000
        else:
            rng=np.random.choice([0,0,0,0,0])
        if rng==0:
            out_list.append(ed)
            out_list.append(ed[::-1])
        elif rng==1:
            out_list.append(ed)
        elif rng==-1:
            out_list.append(ed[::-1])
        else:
            pass
    return out_list
    
def make_erc_timed_array(cores,times,indicies,cycles,count,stim_dt,active_counts):
    for corel in cores:
        forbidden=[]
        for i in range(10):
            for k in range(1):
                options=np.arange(corel-7,corel+7)
                terc=np.random.choice(options,50)
                for ex_count in range(active_counts):
                    for var in terc:
                        time_temp=((i*10+k+200+(count+ex_count)*200)*10+np.round((rand()-.5)*20))*.1*ms    
                        time_temp=np.round(time_temp,decimals=4)
                        if not((var,time_temp) in forbidden):
                            times.append(time_temp)
                            indicies.append(var)
                            forbidden.append((var,time_temp))
        count+=active_counts
    return indicies,times, count
    



def make_waving_stimulus(count,time_chunks,on_off,inhib_level):
    time_chunks=2500
    w_stimulus=np.zeros((time_chunks,number_of_neurons))
    for j in range(0,time_chunks,5):
        for i in range(number_of_neurons):
            if j<=110:
                if rn.random()<(np.abs(on_off-perc_reduce)):
                    w_stimulus[j,i]=(np.random.rand()-.5)*inhib_level+inhib_level 
            else:
                if rn.random()<perc_reduce:
                    w_stimulus[j,i]=(np.random.rand()-.5)*inhib_level+inhib_level
    for j in range(0,time_chunks,5):
        if j+5>=2500:
                continue
        for i in range(number_of_neurons):
            start=w_stimulus[j,i]
            stop=w_stimulus[j+5,i]
            slope=(stop-start)/5
            for nn in range(1,5):
                w_stimulus[j+nn,i]=start+slope*nn
    return w_stimulus
                
def make_easy_to_use_stimulus(sample_rate, current_w_sample,w_stimulus):
    upscale=current_w_sample/sample_rate
    start=np.array(np.shape(w_stimulus))
    start[0]=start[0]*upscale
    blank=np.zeros(start)
    for j in range(0,start[0]):
        blank[j,:]=w_stimulus[int(np.floor(j/upscale)),:]
    return blank


loop='0'
ltp_val='2'
excival='17'
latinhibyval='0'


stim_dt=100
number_of_neurons=1200
cores=np.array([np.arange(50,750,100)]).T
stimulus1=np.zeros((2500,number_of_neurons))
strength_of_wave=0
t_start=1600
perc_reduce= 0.7
time_chunks=100
connect_break=1
on_off=1
inhib_level=1

w_stimulus=np.zeros((350,number_of_neurons))

w_stimulus= make_waving_stimulus(number_of_neurons,time_chunks,on_off,inhib_level) #crucial
sample_rate=.5
current_w_sample=10
blank=make_easy_to_use_stimulus(sample_rate, current_w_sample,w_stimulus)

#w_stimulus=np.zeros((350,number_of_neurons))

stimulus1=stimulus1+w_stimulus #crucial line



stimulus1[5,:]=-200

stim=TimedArray(stimulus1*nA,dt=10*ms)  
cycles=2
count=4
times=[]
indicies=[]
active_counts=1
indicies,times,last_active=make_erc_timed_array(cores[3],times,indicies,cycles,count,stim_dt,active_counts)

erc_out=SpikeGeneratorGroup(1000,indicies,times)
defaultclock.dt = 0.01*ms

Cm = 1*uF # /cm**2
Iapp = 100*nA
gL = 0.3*msiemens
EL = -54.4*mV
ENa = 50*mV
EK = -77*mV
gNa = 35*msiemens
gK = 36*msiemens
tauCa=20*ms
tauinssynslo=1*ms

Isyn=0*uamp
areapyr=50000*um**2
test=areapyr
Ip0=-0.3*nA
gLs=0.1  *msiemens*cm**-2*test
gLd=0.1  *msiemens*cm**-2*test
gNa=30  *msiemens*cm**-2*test
gKdr=15  *msiemens*cm**-2*test
gCa=7  *msiemens*cm**-2*test
gKahp=0.8*msiemens*cm**-2*test  
gKC=15  *msiemens*cm**-2*test
VNa=60*mV
VCa=80*mV
VK=-75*mV
VL=-60*mV
Vsyn=0*mV 
vtsyn=0*mV 
gc=.5*msiemens*cm**-2*test
pp=0.5
Cm=3  *uF*cm**-2*test
gNMDA=0*msiemens*cm**-2*test
gampa=0*msiemens*cm**-2*test
gbarsyn=0.8*nS
gbaryisyn=0.6*nS
tausyn= 4*ms
taueyn= 2*ms
tauisyn=6*ms
scaler=(10**8*um**2/(test))

#interneuron stuff
interneuron_size=20000*um**2
Cmi = 1*uF/cm**2*interneuron_size
gL = 0.1*msiemens/cm**2*interneuron_size

VL2 = -65*mV
gNai = 35*msiemens/cm**2*interneuron_size
gKi = 9*msiemens/cm**2*interneuron_size
Esyn=-75*mV
taguisyn=2*ms
Iapp2=0*nA

Pinksy_rinzel_eqs='''
dVs/dt=(-gLs*(Vs-VL)-gNa*(Minfs**2)*hs*(Vs-VNa)-gKdr*ns*(Vs-VK)+(gc/pp)*(Vd-Vs)+(Ip0-Iinssyn)/pp)/Cm : volt
dVd/dt=(-gLd*(Vd-VL)-ICad-gKahp*qd*(Vd-VK)-gKC*cd*chid*(Vd-VK)+(gc*(Vs-Vd))/(1.0-pp)-Isyn/(1.0-pp))/Cm : volt
dCad/dt=  -0.13*ICad/uamp/ms*scaler-0.075*Cad/ms : 1
dhs/dt=  alphahs-(alphahs+betahs)*hs : 1
dns/dt=  alphans-(alphans+betans)*ns : 1
dsd/dt=  alphasd-(alphasd+betasd)*sd : 1
dcd/dt=  alphacd-(alphacd+betacd)*cd : 1
dqd/dt=  alphaqd-(alphaqd+betaqd)*qd : 1
ICad=     gCa*sd*sd*(Vd-VCa) : amp
alphams=  0.32*(-46.9-Vs/mV)/(exp((-46.9-Vs/mV)/4.0)-1.000001)/ms : Hz
betams=   0.28*(Vs/mV+19.9)/(exp((Vs/mV+19.9)/5.0)-1.000001)/ms : Hz
Minfs=    alphams/(alphams+betams) : 1
alphans=  0.016*(-24.9-Vs/mV)/(exp((-24.9-Vs/mV)/5.0)-1.000001)/ms : Hz
betans=   0.25*exp(-1.0-0.025*Vs/mV)/ms : Hz
alphahs=  0.128*exp((-43.0-Vs/mV)/18.0)/ms : Hz
betahs=   4.0/(1.00001+exp((-20.0-Vs/mV)/5.0))/ms : Hz
alphasd=  1.6/(1.0001+exp(-0.072*(Vd/mV-5.0)))/ms : Hz
betasd=   0.02*(Vd/mV+8.9)/(exp((Vd/mV+8.9)/5.0)-1.0)/ms : Hz
alphacd=((Vd/mV<=-10)*exp((Vd/mV+50.0)/11-(Vd/mV+53.5)/27)/18.975+(Vd/mV>-10)*2.0*exp((-53.5-Vd/mV)/27.0))/ms  : Hz
betacd=   ((Vd/mV<=-10)*(2.0*exp((-53.5-Vd/mV)/27.0)-alphacd*ms)+(Vd/mV>-10)*0)/ms : Hz
alphaqd=  clip(0.00002*Cad,0,0.01)/ms : Hz
betaqd=   0.001/ms : Hz
chid=     clip(Cad/250.0,0,1.0) : 1
INMDA=gNMDA*Si*(1+0.28*exp(-0.062*(Vs/mV)))**(-1)*(Vd-Vsyn) : amp
Hxs=0<=((Vs/mV+50)*1) : 1
Hxw=0<=((Vs/mV+40)*1) : 1
dSi/dt=(Hxs-Si/150)/second : 1
Isyn=Issyn+INMDA+Iapp+Iarray : amp
Issyn=gbarsyn*clip((ssyn+essyn),0,7000)*(Vd-Vsyn): amp
dssyn/dt=-ssyn/tausyn: 1
dessyn/dt=-essyn/taueyn :1
dinssyn/dt=-inssyn/tauisyn: 1
dinssynslo/dt=-inssynslo/tauinssynslo :1
Iinssyn=gbaryisyn*clip(inssyn+inssynslo,0,7000)*(Vs-VK): amp
dIapp/dt=0*amp/ms : amp
Iarray=stim(t,i) :amp
dCa/dt=-Ca/tauCa :1
counter:1
'''

buzsaki_eqs = '''
dv/dt = (-gNai*m**3*h*(v-VNa)-gKi*n**4*(v-VK)-gL*(v-VL2)-Issyn-Iinssyn-Iapp2)/Cmi : volt
m = alpha_m/(alpha_m+beta_m) : 1
alpha_m = -0.1/mV*(v+35*mV)/(exp(-0.1/mV*(v+35*mV))-1)/ms : Hz
beta_m = 4*exp(-(v+60*mV)/(18*mV))/ms : Hz
dh/dt = 5*(alpha_h*(1-h)-beta_h*h) : 1
alpha_h = 0.07*exp(-(v+58*mV)/(20*mV))/ms : Hz
beta_h = 1./(exp(-0.1/mV*(v+28*mV))+1)/ms : Hz
dn/dt = 5*(alpha_n*(1-n)-beta_n*n) : 1
alpha_n = -0.01/mV*(v+34*mV)/(exp(-0.1/mV*(v+34*mV))-1)/ms : Hz
beta_n = 0.125*exp(-(v+44*mV)/(80*mV))/ms : Hz
Issyn=gbarsyn*(clip(ssyn,0,7000))*(v-Vsyn): amp
dssyn/dt=-ssyn/tausyn: 1
dinssyn/dt=-inssyn/taguisyn: 1
Iinssyn=gbaryisyn*clip(inssyn,0,7000)*(v-VK):amp
'''
    
Pii=NeuronGroup(100,buzsaki_eqs,threshold='v>-20*mV',refractory='v>-60*mV',method='euler')
Pii.v='(randn()*60*.05-60)*mV'
Pii.h=.999
Pii.n=.0001
Pii.ssyn=0    
Pe = NeuronGroup(1000, Pinksy_rinzel_eqs,threshold='Vs>-20*mV',refractory='Vs>-60*mV',reset='Ca=50')
Pe.Vs='(randn()*60*.05-60)*mV'
Pe.Vd='(randn()*65*.05-65)*mV'
Pe.hs=.999
Pe.ns=.0001
Pe.sd=.009
Pe.cd=.007
Pe.qd=.01
Pe.Cad=.20
Pe.Si=0
Pe.ssyn=0
Pe.inssyn=0
Pe.essyn=0

erc_syn=Synapses(erc_out, Pe, on_pre='essyn +=12')
erc_syn.connect(j='i')
netsize=25
out_list=make_network(cores,number_of_neurons,netsize,connect_break)
out_list=np.array(out_list)


Pi = Pii

PI = PoissonInput(Pe,'ssyn',N=20,rate=10*Hz,weight='3')


erc_insyn=Synapses(erc_out,Pe,on_pre='inssyn+='+latinhibyval)
erc_insyn.connect(condition='(abs(i-j)>50)*(rand()<0.1)')
Ce = Synapses(Pe, Pe, ''' dw/dt=-0.001/ms : 1''',
              on_pre='''w=clip(w,1,'''+ltp_val+''')
              ssyn+='''+excival+'''*w
              w+=(Ca>20)*1
              w=clip(w,1,'''+ltp_val+''')
              ''',delay=1*rand()*ms)

Ce.connect(i=out_list[:,0],j=out_list[:,1])
Ce.w=1
Ci = Synapses(Pi, Pe, on_pre='inssyn +=35',delay=1*rand()*ms)
Ci.connect(p=0.3)
Ci2 = Synapses(Pi, Pi, on_pre='inssyn +=35',delay=1*rand()*ms)
Ci2.connect(p=0.3)
Ci.delay=2*ms
Cei=Synapses(Pe,Pi,on_pre='ssyn +=5')
Cei.connect(p=0.15)


M = SpikeMonitor(Pe)
SM = StateMonitor(Pe,'Vs',record=True,dt=.5*ms)

run(5*second,report='text')



plt.scatter(M.t,M.i,marker='.',s=.5)

outname='1'




zed=SM.Vs[200:400:30]/mV
np.savetxt(r'.//'+outname+'.csv',zed.T,delimiter=',')
b=np.vstack([M.t[:]/ms,M.i[:]])
np.savetxt(r'.//'+outname+'dots.csv',b.T,delimiter=',')



plt.savefig(r'.//'+outname+'.png')

