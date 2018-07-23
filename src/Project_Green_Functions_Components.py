# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from numpy import loadtxt,sqrt
import xml.dom.minidom as md
from matplotlib import pyplot as plt

srcpath='./output/'
input_xml = md.parse('./input.xml')
nx=int(input_xml.getElementsByTagName('nx')[0].firstChild.nodeValue)
ny=int(input_xml.getElementsByTagName('ny')[0].firstChild.nodeValue)
nz=int(input_xml.getElementsByTagName('nz')[0].firstChild.nodeValue)
nt=int(input_xml.getElementsByTagName('nt')[0].firstChild.nodeValue)

    
#Load data into memory
green_function=[]
sigma=[]
int_time=[]
for i in [1,2,3]:
    green_function.append([])
    sigma.append([])
    int_time.append([])
    for j in [1,2,3]:
        momentum=loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(0,),dtype=int)
        green_function[i-1].append(loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(1,)))
        sigma[i-1].append(loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(2,)))
        #int_time[i-1].append(loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(3,)))
        for s in sigma[i-1][j-1]:
            if(s < 0):
                raise(ValueError,"Warning! Found negative error!")
            #end if
        #end for
    #end for
#end for
        
        
#Compute G_L and G_T as function of k

G_L = []
G_T = []
error_GL = []
error_GT = []
for k in range(len(momentum)):
    #Determines kx,ky,kz,omega
    p = momentum[k]
    omega = p//(nx*ny*nz)
    pz = (p - omega*nx*ny*nz)//(nx*ny)
    py = (p - omega*nx*ny*nz - pz*nx*ny)//nx
    px = p - omega*nx*ny*nz - pz*nx*ny - py*nx
    p_i = [px,py,pz,omega,px**2+py**2+pz**2,px**4+py**4+pz**4,px**6+py**6+pz**6]

    L = 0.0
    T = 0.0
    integrated_time=0
    err_L=0.0
    err_T=0.0
    #Computes longitudinal and transversal part
    if (omega == 0 and p_i[4] != 0):
        for i in range(3):
            for j in range(3):
                L = L + p_i[i]*p_i[j]*green_function[i][j][k]/p_i[4]
                err_L = err_L + p_i[i]*p_i[j]*sigma[i][j][k]/p_i[4]
                if(i == j):
                    T = T + green_function[i][j][k]
                    err_T = err_T + sigma[i][j][k]
                #end if
            #end for
        #end for
        T = (T - L)/2.0
        err_T = err_T + err_L
      
        G_L.append([p_i[4],p_i[5],p_i[6],L])
        G_T.append([p_i[4],p_i[5],p_i[6],T])
        error_GL.append([p_i[4],p_i[5],p_i[6],err_L])
        error_GT.append([p_i[4],p_i[5],p_i[6],err_T])
        
    #end if
#end for

plt.errorbar([point[0] for point in G_L[:]],[point[3] for point in G_L[:]],[point[3] for point in error_GL[:]],None,'bo',ms=3)
plt.errorbar([point[0] for point in G_T[:]],[point[3] for point in G_T[:]],[point[3] for point in error_GT[:]],None,'go',ms=3)
plt.show()

#Average over orbits (same p^2, p^[4] and p^[6])

#Create a vector with all points p^2, p^[4] and p^[6]

#sort the vectors
G_L.sort()
G_T.sort()
error_GL.sort()
error_GT.sort()

def avrg_over_orbit(data):
    p0=[data[0][0],data[0][1],data[0][2]]
    orbit_index=0
    orbit_count=0
    orb_avrgd=[list(p0)]
    orb_avrgd[0].append(0)
    for point in data:
        if ( [point[0],point[1],point[2]] == p0 ):
            orb_avrgd[orbit_index][3] = orb_avrgd[orbit_index][3] + point[3]
            orbit_count = orbit_count + 1
        else:
            orb_avrgd[orbit_index][3]=orb_avrgd[orbit_index][3]/orbit_count
            p0 = [point[0],point[1],point[2]]
            orbit_index = orbit_index + 1
            orbit_count = 1
            orb_avrgd.append(list(p0))
            orb_avrgd[orbit_index].append(point[3])
        #end if
    #end for
    return(orb_avrgd)
#end def
    
orb_avrgd_GL = avrg_over_orbit(G_L)
orb_avrgd_GT = avrg_over_orbit(G_T)
orb_avrgd_error_GL = avrg_over_orbit(error_GL)
orb_avrgd_error_GT = avrg_over_orbit(error_GT)

#Let's see how different they are comparing to before
plt.errorbar([point[0] for point in orb_avrgd_GL[:]],[point[3] for point in orb_avrgd_GL[:]],[point[3] for point in orb_avrgd_error_GL[:]],None,'bo',ms=3)
plt.errorbar([point[0] for point in orb_avrgd_GT[:]],[point[3] for point in orb_avrgd_GT[:]],[point[3] for point in orb_avrgd_error_GT[:]],None,'go',ms=3)
plt.show()

#Let us try an avergae over p^2. This should give us a first approximation of what to expect

def avrg_over_k2(data):
    p0=data[0][0]
    p_index=0
    p_count=0
    p_avrgd=[[p0]]
    p_avrgd[0].append(0)
    for point in data:
        if ( point[0] == p0 ):
            p_avrgd[p_index][1] = p_avrgd[p_index][1] + point[3]
            p_count = p_count + 1
        else:
            p_avrgd[p_index][1]=p_avrgd[p_index][1]/p_count
            p0 = point[0]
            p_index = p_index + 1
            p_count = 1
            p_avrgd.append([p0])
            p_avrgd[p_index].append(point[3])
        #end if
    #end for
    return(p_avrgd)
#end def
    
p2_avrgd_GL = avrg_over_k2(orb_avrgd_GL)
p2_avrgd_GT = avrg_over_k2(orb_avrgd_GT)
p2_avrgd_error_GL = avrg_over_k2(orb_avrgd_error_GL)
p2_avrgd_error_GT = avrg_over_k2(orb_avrgd_error_GT)

plt.errorbar([point[0] for point in p2_avrgd_GL[:]],[point[1] for point in p2_avrgd_GL[:]],[point[1] for point in p2_avrgd_error_GL[:]],None,'bo',ms=3)
plt.errorbar([point[0] for point in p2_avrgd_GT[:]],[point[1] for point in p2_avrgd_GT[:]],[point[1] for point in p2_avrgd_error_GT[:]],None,'go',ms=3)
plt.show()

#The right thing to do:
#   a) Select all points with the same p^2
#   b) Fit the values to F(p^[4]) = a + b*p^[4]
#   c) The value for that p^2 is the value of a

def extrapolator(data,error_data):