#!/usr/bin/env python

import numpy as np
import f90nml as nr
import matplotlib.pyplot as plt

iparam=nr.read('int_params.nml')
mparam=nr.read('params.nml')
# print iparam['int_params']
# print mparam

f=open('mean_lyapunov.dat','r')
res=f.readlines()
f.close()

f=open('spectra-lyapexp-CO350-d1x10-8-llll.dat','r')
ref=f.readlines()
f.close()

mean=[]
var=[]
for x in res:
    y=x.split()
    if y[0]=='mean':
        z=mean
        for t in y[1:]:
            z.append(float(t))
    elif y[0]=='var':
        z=var
        for t in y[1:]:
            z.append(float(t))
    else:
        for t in y:
            z.append(float(t))

mean_ref=[]
var_ref=[]

for x in ref:
    y=x.split()
    mean_ref.append(float(y[1]))
    var_ref.append(float(y[2]))

mean=np.array(mean)
var=np.array(var)
mean_ref=np.array(mean_ref)
var_ref=np.array(var_ref)
n=len(mean)

facLE=mparam['aoscale']['f0']*24*3600

fig=plt.figure(num=10,figsize=(8,12))
ax=[]
ax.append(fig.add_subplot(2,1,1))
ax.append(fig.add_subplot(2,1,2))

ax[0].plot(range(1,n+1),mean*facLE,ls='',marker='x',color='b')
ax[1].plot(range(1,n+1),var*facLE**2,ls='',marker='x',color='b')

ax[0].plot(range(1,n+1),mean_ref,color='r')
ax[1].plot(range(1,n+1),var_ref,color='r')

ax[0].set_ylabel(r'[1/day]')
ax[1].set_ylabel(r'[1/day]')
ax[1].set_xlabel('Lyapunov exponent number')
plt.show()
