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
sigma2=[]
int_time=[]
for i in [1,2,3]:
    green_function.append([])
    sigma2.append([])
    int_time.append([])
    for j in [1,2,3]:
        momentum=loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(0,),dtype=int)
        green_function[i-1].append(loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(1,)))
        sigma2[i-1].append(loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(2,)))
        int_time[i-1].append(loadtxt(srcpath+'Corr4'+str(i)+'4'+str(j)+".dat",usecols=(3,)))
        for sigma in sigma2[i-1][j-1]:
            if(sigma < 0):
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
    p_i = [px,py,pz,omega,px**2+py**2+pz**2]

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
                if (int_time[i][j][k] != 0):
                    err_L = err_L + p_i[i]*p_i[j]*sqrt(sigma2[i][j][k]*2*int_time[i][j][k])/p_i[4]
                else:
                    err_L = err_L + p_i[i]*p_i[j]*sqrt(sigma2[i][j][k])/p_i[4]
                if(i == j):
                    T = T + green_function[i][j][k]
                    if (int_time[i][j][k] != 0):
                        err_T = err_T + sqrt(sigma2[i][j][k]*2*int_time[i][j][k])
                    else:
                        err_T = err_T + sqrt(sigma2[i][j][k])
                #end if
            #end for
        #end for
        T = (T - L)/2.0
        err_T = err_T + err_L
      
        G_L.append([p_i[4],L])
        G_T.append([p_i[4],T])
        error_GL.append([p_i[4],err_L])
        error_GT.append([p_i[4],err_T])
        
    #end if
#end for

plt.plot([point[0] for point in G_L[:]],[point[1] for point in G_L[:]],'bo',ms=3)
plt.plot([point[0] for point in G_T[:]],[point[1] for point in G_T[:]],'go',ms=3)
plt.show()

plt.errorbar([point[0] for point in G_L[:]],[point[1] for point in G_L[:]],[point[1] for point in error_GL[:]],None,'bo',ms=3)
plt.errorbar([point[0] for point in G_T[:]],[point[1] for point in G_T[:]],[point[1] for point in error_GT[:]],None,'go',ms=3)
plt.show()
plt.errorbar([point[0] for point in G_L[:]],[point[1] for point in G_L[:]],[point[1] for point in error_GL[:]],None,'bo',ms=3)
plt.show()
plt.errorbar([point[0] for point in G_T[:]],[point[1] for point in G_T[:]],[point[1] for point in error_GT[:]],None,'go',ms=3)
plt.show()