#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
"""

from numpy import loadtxt,average,std
import numpy as np
from matplotlib import pyplot as plt
from math import sqrt,log

def auto_corr(d):
    len_out = int(len(d)/2)
    out=np.zeros(len_out)
    avrg = average(d[:len_out])
    for t in range(len_out):
        for s in range(int(len_out)):
            out[t] = out[t] + d[t+s]*d[s] - avrg**2
        out[t] = out[t]/(len_out) #- avrg**2
    out = out/np.std(d)**2
    return(out)

def exp_time(correl_func):
    out=np.zeros(len(correl_func))
    for i in range(len(correl_func)):
        #print(i,log(abs(correl_func[i])))
        out[i] = -i/log(abs(correl_func[i]))
    return(out)
    
#root_dir = "D:\\git_repos\\SU3SimSuit"
root_dir = "/home/willian/git_repos/SU3SimSuit"
data=loadtxt(root_dir+'/output/avrg_plaquette.out')

#data=data[:5000]

plt.plot(range(len(data)),data)

#input("Press enter to continue.")

therm=50 #Thermalization time

c = auto_corr(data[therm:])
tau_exp = exp_time(c)
plt.figure()
plt.plot(range(len(c)),c)
#plt.figure()
#plt.plot(range(len(tau_exp)),tau_exp)


t=1
plt.figure()
plt.plot(range(0,len(data[therm::t])),data[therm::t])
print(average(data[therm::t]),'+/-',std(data[therm::t])/sqrt((len(data[therm::t]))))

def jacknife(data,nbins):
    bin_data=np.zeros(nbins,dtype=np.float64)
    bin_size=len(data)/nbins
    for i in range(len(data)):
        bin_index=int(i/bin_size)
        bin_data[bin_index] = bin_data[bin_index] + data[i]
    for i in range(nbins):
        bin_data[i] = bin_data[i]/bin_size
        
    avrg=average(data)
    error = 0.0
    for i in range(nbins):
        error=error+(bin_data[i]-avrg)**2
    
    error = error*(nbins-1.0)/nbins
    return(error)

print("Jacknife error:",jacknife(data[therm:],4))
plt.show()

