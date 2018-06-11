# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from numpy import loadtxt,average,std
import numpy as np
from matplotlib import pyplot as plt
from math import sqrt

def auto_corr(d):
    out=[]
    avrg = average(d)
    for i in range(int(len(d)/2)):
        out.append(0)
        for j in range(int(len(d)/2)):
            out[-1] =out[-1] + d[j]*d[i+j]
        out[-1] = out[-1]/int(len(d)/2)
        out[-1] = out[-1] - avrg**2
    return(out)            

root_dir = "D:\\git_repos\\SU3SimSuit"
#root_dir = "/home/willian/d/git_repos/SU3SimSuit"
data=loadtxt(root_dir+'/output/avrg_plaquette.out')


plt.plot(range(len(data)),data)

#input("Press enter to continue.")

therm=1 #Thermalization time

c = auto_corr(data[therm:])
plt.figure()
plt.plot(range(len(c)),c)

t=50
plt.figure()
plt.plot(range(len(data[therm::t])),data[therm::t])
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
        
print(average(data[therm:]),'+/-',std(data[therm:])/sqrt((len(data)-therm)))

print("Jacknife error:",jacknife(data[therm:],4))
#plt.plot(range(len(data)),data)

plt.show()

