# -*- coding: utf-8 -*-
"""
Created on Sun May 22 19:14:36 2022

@author: josem
"""

#%%
'''programa para graficas S13,S14'''
import numpy as np
import matplotlib.pyplot as plt



#%%


#%%
'''constantes'''
L=3.75e-7 #H/m
Lm=7.63e-8 #H/m
C=8.264e-10 #F/m
Cm=5.54e-11 #F/m
long=0.15 #m
#%%
w=np.linspace(10**8,3*10**9)
gmm= 1j*w*np.sqrt(L*Cm+Lm*C)
gm= 1j*w*np.sqrt(L*C)
V2= lambda x: (Lm/(2*L)-gmm**2/(4*gm**2)*(1+2*gm*x)*np.exp(-gm*x))+\
                    (gmm**2/(4*gm**2)-Lm/(2*L))*np.exp(-2*gm*long)*np.exp(gm*x)
                                                       
AbsV2=lambda x: abs(V2(x))
x=np.linspace(0,long)
plt.plot(x,AbsV2(x))
#%%
S13= 20*np.log10(abs(V2(0)))                
plt.plot(w, S13,label="S13")
plt.legend()
plt.show

S14=20*np.log10(abs(V2(long)))
plt.plot(w, S14, label="S14")
plt.legend()
plt.show
#%%
