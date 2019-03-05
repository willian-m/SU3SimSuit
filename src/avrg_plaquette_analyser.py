# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""

from matplotlib import pyplot as plt
import numpy as np


def correlation_func_computer(d):
    data_fft=np.fft.rfft(d)
    corr_fourier=[]
    for omega in data_fft:
        corr_fourier.append(omega*np.conjugate(omega))
    #end for
    return(np.fft.irfft(corr_fourier)/len(d)-float(data_fft[0])**2/len(d)**2)
    
def exp_corr_time(c):
    i = 0
    while (i < len(c)-1) & (c[i] > c[0]*np.exp(-1)):
        i =i + 1
    #end while
    return(i)

def int_corr_time(b_size,inp_data):
    nbins = len(inp_data)//b_size
    binned_data = []
    for i in range(nbins):
        binned_data.append(0)
        for j in range(bin_size):
            binned_data[-1] = binned_data[-1] + inp_data[j+i*bin_size]
        #end for
        binned_data[-1]=binned_data[-1]/b_size
    #end for
    
    return(np.std(binned_data)**(2)/np.std(inp_data)**(2))
    
def auto_widowing_int_corr_time(inp_corr,c):
    #c is the constant of the auto windowing algorithm
    #Typically c~6
    normalized_corr=inp_corr/inp_corr[0]
    m=1
    tau=1
    while (m < c*tau):
        tau=0
        for i in range(m):
            tau=tau+normalized_corr[i]
        m=m+1
    
    return(tau,m)
    
#Load action data
plaquette_data = np.loadtxt('output/avrg_plaquette.out',dtype=np.float64,delimiter=",")

corr_func=correlation_func_computer(plaquette_data)
exp_time=exp_corr_time(corr_func)

#Uses exp_time to estimate a thermalization time
therm_time=10*exp_time

#If resulting array has odd length, increase the therm_time by 1
if( (len(plaquette_data)-therm_time)%2 != 0):
    therm_time = therm_time + 1

delta_therm_time=10
#Recomputes correlation function, now without thermalization influence
while delta_therm_time > 0:
    corr_func=correlation_func_computer(plaquette_data[therm_time:])
    exp_time=exp_corr_time(corr_func)
    delta_therm_time=10*exp_time-therm_time
    therm_time=10*exp_time

int_time,m=auto_widowing_int_corr_time(corr_func,6)
#We bin the data on bins of size at least 2*exp_time
#bin_size = 2*exp_time
#while( len(plaquette_data[therm_time:])%bin_size != 0  & bin_size < len(plaquette_data[therm_time:])//2):
#    bin_size = bin_size + 1

#int_time = int_corr_time(bin_size,plaquette_data[therm_time:])

if (int_time < 1): #Decorrelated samples. Use standard techniques
    print("Thermalization time:",therm_time)
    print("Data can be treated as decorrelated")
    print("Average action:",np.average(plaquette_data[therm_time:]))
    print("Average error:",np.std(plaquette_data[therm_time:])/np.sqrt(len(plaquette_data[therm_time:])))
else: #Error is the naive times sqrt(int_time)
    print("Thermalization time:",therm_time)
    print("Correlation time:",int_time)
    print("Correlation time error:",np.sqrt(2*(2*m+1)*int_time**2/(len(plaquette_data[therm_time:]))))
    print("Average action:",np.average(plaquette_data[therm_time:]))
    print("Average error:",np.sqrt(2*int_time)*np.std(plaquette_data[therm_time:])/np.sqrt(len(plaquette_data[therm_time:])))
#plt.plot(range(len(corr_func)),corr_func)
plt.plot(range(len(plaquette_data)),plaquette_data)
plt.show()

plt.semilogy(range(len(corr_func)),corr_func/corr_func[0])
plt.show()