# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from numpy import loadtxt, zeros, sin
from numpy.fft import fftn, fftfreq
import xml.dom.minidom as md
from matplotlib import pyplot as plt

srcpath='./bin/'
input_xml = md.parse('./input.xml')
nx=int(input_xml.getElementsByTagName('nx')[0].firstChild.nodeValue)
ny=int(input_xml.getElementsByTagName('ny')[0].firstChild.nodeValue)
nz=int(input_xml.getElementsByTagName('nz')[0].firstChild.nodeValue)
nt=int(input_xml.getElementsByTagName('nt')[0].firstChild.nodeValue)

    
#Load data into memory
index1=loadtxt(srcpath+"fermion_prop.dat",usecols=(0,),dtype=int)
index2=loadtxt(srcpath+"fermion_prop.dat",usecols=(1,),dtype=int)
prop_re = loadtxt(srcpath+"fermion_prop.dat",usecols=(2,))
prop_im = loadtxt(srcpath+"fermion_prop.dat",usecols=(3,))

#takes the trace over dirac and color indexes
prop=zeros((nx,ny,nz,nt),dtype=complex)
for x in range(nx):
    for y in range(ny):
        for z in range(nz):
            for t in range(nt):
                for a in range(12):
                    for b in range(12):
                        i = b + 12*a + 144*x + 144*nx*y + 144*nx*ny*z + 144*nx*ny*nz*t 
                        prop[x,y,z,t] = prop[x,y,z,t] + complex( prop_re[i] , prop_im[i] )/12


prop_fftn = fftn(prop, s=(nx,ny,nz,nt), axes=(0,1,2,3))

freq_x = fftfreq(nx)
freq_y = fftfreq(ny)
freq_z = fftfreq(nz)
freq_t = fftfreq(nt)

orbits_dict = {}
for x in range(nx):
    for y in range(ny):
        for z in range(nz):
            for t in range(nt):
                p2 = freq_x[x]**2 + freq_y[y]**2 + freq_z[z]**2 + freq_t[t]**2
                p4 = freq_x[x]**4 + freq_y[y]**4 + freq_z[z]**4 + freq_t[t]**4
                if (p2,p4) not in orbits_dict.keys():
                    orbits_dict[p2,p4] = [[x,y,z,t]]
                else:
                    orbits_dict[p2,p4].append([x,y,z,t])
         
prop_orbit_averaged = []
for orbit in orbits_dict.keys():
    i = 0
    prop_orbit_averaged.append([orbit[0],orbit[1],0])
    for momentum in orbits_dict[orbit]:
        prop_orbit_averaged[-1][2] = prop_orbit_averaged[-1][2] + prop_fftn[momentum[0],momentum[1],momentum[2],momentum[3]]
        i = i+1
    #end for
    prop_orbit_averaged[-1][2]=prop_orbit_averaged[-1][2]/i
        

prop_theo_orbit_averaged = []
for orbit in orbits_dict.keys():
    i = 0
    prop_theo_orbit_averaged.append([orbit[0],orbit[1],0])
    for momentum in orbits_dict[orbit]:
        prop_theo_orbit_averaged[-1][2] = prop_theo_orbit_averaged[-1][2] + .4/(sin(momentum[0])**2 + sin(momentum[1])**2 + sin(momentum[2])**2 + sin(momentum[3])**2 + .5)
        i = i+1
    #end for
    prop_theo_orbit_averaged[-1][2]=prop_theo_orbit_averaged[-1][2]/i
                    
                    
                    
plt.plot([x[0] for x in prop_orbit_averaged[:-1]], [p[2] for p in prop_orbit_averaged[:-1]],"bo")
plt.plot([x[0] for x in prop_theo_orbit_averaged], [p[2] for p in prop_theo_orbit_averaged],"go")
plt.show()