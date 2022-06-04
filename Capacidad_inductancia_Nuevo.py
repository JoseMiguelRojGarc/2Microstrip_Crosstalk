# -*- coding: utf-8 -*-
"""
Created on Wed May 11 19:18:44 2022

@author: josem
"""

#%%
import numpy as np
import scipy.constants as const
#%%

#%%
''' Constantes'''
W=0.002 #m
Wd=0.01
s=0.0013 #m
Hs=0.0016 #m
Ws=2*Wd+2*W+s
eps_r=1.000000000001
Cte=1/(2*np.pi*const.epsilon_0)


#%%

#%%
def en_linea1(n,dl):
    if (n>=Wd/dl-1)&(n<=(Wd+W)/dl-1):
        return 1
    else:
        return 0
def en_linea2(n,dl):
    if (n>(Wd+W+s)/dl-1)&(n<(Wd+2*W+s)/dl-1):
        return 1
    else:
        return 0

#%


#%%

'''Matriz S para los conductores '''
'''Se me queda estable hasta unos 2000'''
N=2000
xp=np.linspace(-Ws/2.0,Ws/2.0, num=N)
x=np.linspace(-Ws/2.0,Ws/2.0, num=N)
dl=Ws/N
x=x+dl/2.0
S=np.zeros((N,N))
n=0
while n<=N-1:
    m=0
    while m<=N-1:
        if (en_linea1(m,dl)==1) or (en_linea2(m,dl)==1):
            S[m][n]=Cte*(1/2)*dl*np.log((4*Hs**2+(x[m]-xp[n])**2)/(x[m]-xp[n])**2)
        else:
            if m==n:
                S[m][n]=Cte*dl*(-2*Hs/(4*Hs**2+(x[m]-xp[n])**2))+dl*((1+eps_r)/(2*const.epsilon_0*(1-eps_r)))
            else:
                S[m][n]=Cte*dl*(-2*Hs/(4*Hs**2+(x[m]-xp[n])**2))
        m=m+1
    n=n+1
#%%



#%%
'''Ecuaciones para el cálculo de las densidades de carga'''
'''definimos la matriz b de 1 y ceros'''
#OJOOOO, ME ESTÁ MACIENDO MAL EL INVERSO....................
#Si pongo una matriz demasiado grande me devielve que el inverso es cero 
b1=np.zeros(N)
b2=np.zeros(N)
i=0
while i<=N:
    if en_linea1(i, dl)==1:
        b1[i]=1
    elif en_linea2(i,dl)==1:
        b2[i]=1    
    i=i+1
    
inv_S=np.linalg.inv(S)
sigma1_n= np.dot(inv_S,b1)
sigma2_n=np.dot(inv_S,b2)

print(inv_S)
#%%

#%%
'''Calculo de la capacidad'''




#ARREGLAR ÉSTO, QUE ESTÁ MAL, LAS SIGMA_N SE INTEGRAN(SEGURO?????)
#ME ESTÁ DANDO CERO YA QUE EN CADA LÍNEA LA CARGA NETA ES NULA. 



'''Cálculo de C12'''

C12=0 
i=0
while i<=N-1:
    if en_linea1(i, dl)==1:
        C12=C12+((eps_r+1)/2)*sigma2_n[i]*dl
        j=0
        while j<=N-1:
            C12=C12+((1-eps_r)/(2*np.pi))*sigma2_n[j]*\
                    (-2*Hs/(4*Hs**2+(x[i]-xp[j])**2))*dl**2
            j=j+1
    i=i+1

print('C12=',C12)

C21=0 
i=0
while i<=N-1:
    if en_linea2(i, dl)==1:
        C21=C21+((eps_r+1)/2)*sigma1_n[i]*dl
        j=0
        while j<=N-1:
            C21=C21+((1-eps_r)/(2*np.pi))*sigma1_n[j]*\
                    (-2*Hs/(4*Hs**2+(x[i]-xp[j])**2))*dl**2
            j=j+1
    i=i+1


print('C21=',C21)

C11=0 
i=0
while i<=N-1:
    if en_linea1(i, dl)==1:
        C11=C11+((eps_r+1)/2)*sigma1_n[i]*dl
        j=0
        while j<=N-1:
            C11=C11+((1-eps_r)/(2*np.pi))*sigma1_n[j]*\
                    (-2*Hs/(4*Hs**2+(x[i]-xp[j])**2))*dl**2
            j=j+1
    i=i+1


print('C11=',C11)
C22=0  
i=0
while i<=N-1:
    if en_linea2(i, dl)==1:
        C22=C22+((eps_r+1)/2)*sigma2_n[i]*dl
        j=0
        while j<=N-1:
            C22=C22+((1-eps_r)/(2*np.pi))*sigma2_n[j]*\
                    (-2*Hs/(4*Hs**2+(x[i]-xp[j])**2))*dl**2
            j=j+1
    i=i+1

print('C22=',C22)
#%%
C0=np.zeros((2,2))
C0[0][0]=C11
C0[0][1]=C12
C0[1][0]=C21
C0[1][1]=C22
inv_C0=np.linalg.inv(C0)
print(inv_C0)
#%%

