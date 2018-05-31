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
data=loadtxt(root_dir+'/output/avrg_plaquette.out')


plt.plot(range(len(data)),data)

#input("Press enter to continue.")

therm=25 #Thermalization time

c = auto_corr(data[therm:])
plt.figure()
plt.plot(range(len(c)),c)


plt.figure()
plt.plot(range(len(data[therm:])),data[therm:])
print(average(data[therm::10]),'+/-',std(data[therm::10])/sqrt((len(data)-therm)/60))

#plt.plot(range(len(data)),data)

plt.show()

