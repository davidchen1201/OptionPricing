#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm

##############################################################################
#Q1
S0 = 50
T = 0.25
r = 0.05 
sigma = 0.3

#Black-Scholes-Merton Put Option
def BSM_Put(S0,K,T,r,sigma):
    d1 = (np.log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*np.sqrt(T))
    d2 = d1-sigma*np.sqrt(T)
    put = K*np.exp(-r*T)*norm.cdf(-d2)-S0*norm.cdf(-d1)
    return put

M = np.power(10,5)
K = np.array([40, 50, 60])

BSM = np.zeros(3)

#naive MC
Y = np.zeros([3, M])
MC = np.zeros(3) 
Var = np.zeros(3)

#cv1
X1 = np.zeros([3, M])
Y1 = np.zeros([3, M])
MC1 = np.zeros(3)#estimator variance
Var1 = np.zeros(3)
rho1 = np.zeros(3)

#cv2
X2 = np.zeros([3, M])
Y2 = np.zeros([3, M])
MC2 = np.zeros(3)
Var2 = np.zeros(3)
rho2 = np.zeros(3)

#cv3
X3 = np.zeros([3, M])
Y3 = np.zeros([3, M])
MC3 = np.zeros(3)
Var3 = np.zeros(3)
rho3 = np.zeros(3)

for k in range(3):
    BSM[k] = BSM_Put(S0,K[k],T,r,sigma)
    
for k in range(3):
    for run in range(M):
        WT = np.random.normal(0,np.sqrt(T))   
        ST = S0*np.exp((r-0.5*sigma*sigma)*T+sigma*WT)
        X1[k, run] = ST
        X2[k, run] = WT
        X3[k, run] = WT*WT
        Y[k, run] = np.maximum(K[k]-ST,0)
    
    #naive MC
    MC[k] = np.mean(Y[k])*np.exp(-r*T)
    Var[k] = np.var(Y[k])/M*np.exp(-2*r*T)
    
    #cv1
    b1 = np.inner(X1[k]-np.mean(X1[k]),Y[k]-np.mean(Y[k]))/np.inner(X1[k]-np.mean(X1[k]),X1[k]-np.mean(X1[k]))   
    Y1[k] = Y[k]-b1*(X1[k]-np.mean(X1[k]))
    MC1[k] = np.mean(Y1[k])*np.exp(-r*T)
    Var1[k] = np.var(Y1[k])/M*np.exp(-2*r*T)
    rho1[k] = np.corrcoef(X1[k],Y[k])[0,1]
    
    #cv2
    b2 = np.inner(X2[k]-np.mean(X2[k]),Y[k]-np.mean(Y[k]))/np.inner(X2[k]-np.mean(X2[k]),X2[k]-np.mean(X2[k]))
    Y2[k] = Y[k]-b2*(X2[k]-np.mean(X2[k]))
    MC2[k] = np.mean(Y2[k])*np.exp(-r*T)
    Var2[k] = np.var(Y2[k])/M*np.exp(-2*r*T)
    rho2[k] = np.corrcoef(X2[k],Y[k])[0,1]
    
    #cv3
    b3 = np.inner(X3[k]-np.mean(X3[k]),Y[k]-np.mean(Y[k]))/np.inner(X3[k]-np.mean(X3[k]),X3[k]-np.mean(X3[k]))
    Y3[k] = Y[k]-b3*(X3[k]-np.mean(X3[k]))
    MC3[k] = np.mean(Y3[k])*np.exp(-r*T)
    Var3[k] = np.var(Y3[k])/M*np.exp(-2*r*T)
    rho3[k] = np.corrcoef(X3[k],Y[k])[0,1]   

    
##############################################################################
#Q2    
#Geometric brownian motion
S0 = 50
T = 0.25
r = 0.05 
sigma = 0.3

M = np.power(10,5)
K = np.array([40, 50, 60])


raw1 = np.zeros([3,M])
put1 = np.zeros(3)
var1 = np.zeros(3)
raw2 = np.zeros([3,M])
put2 = np.zeros(3)
var2 = np.zeros(3)
se = np.zeros(3)

A = np.linspace(0, 1, 11)  
for k in range(3): 
    #naive MC
    for m in range(M):
        Z= np.random.normal(0,1)
        ST= S0*np.exp((r-0.5*sigma*sigma)*T+sigma*np.sqrt(T)*Z)
        raw1[k,m] = np.maximum(K[k]-ST,0)
    put1[k] = np.mean(raw1[k])**np.exp(-r*T)
    var1[k] = np.inner(raw1[k]-np.mean(raw1[k]), raw1[k]-np.mean(raw1[k]))/(M-1)*np.exp(-2*r*T)    
    
    #stratifying sampling
    for i in range(A.size-1):
        for j in range(np.int(M/(A.size-1))):
            U = np.random.rand()
            V = A[i]+U*(A[i+1]-A[i])
            ST= S0*np.exp((r-0.5*sigma*sigma)*T+sigma*np.sqrt(T)*norm.ppf(V))
            raw2[k,i*np.int(M/(A.size-1))+j] = np.maximum(K[k]-ST,0)   
    put2[k] = np.mean(raw2[k])**np.exp(-r*T)
    var2[k] = np.inner(raw2[k]-np.mean(raw2[k]), raw2[k]-np.mean(raw2[k]))/(M-1)*np.exp(-2*r*T)  
    se[k] = np.sqrt(var2[k]/M)

##############################################################################
#Q3
S0 = 50
T = 0.25
r = 0.05 
sigma = 0.3

M = np.power(10,4)
K = np.array([10, 20, 60])
Y1 = np.zeros([3, M])
MC1 = np.zeros(3)
Var1 = np.zeros(3)
SE1 = np.zeros(3)
RE1 = np.zeros(3)
Y2 = np.zeros([3, M])
MC2 = np.zeros(3)
Var2 = np.zeros(3)
SE2 = np.zeros(3)
RE2 = np.zeros(3)

mu1 = np.log(S0)+(r-0.5*sigma*sigma)*T
#naive MC
for k in range(3): 
    for i in range(M):
        X = np.random.normal(mu1,sigma*np.sqrt(T))
        Y1[k,i] = 1*(X<=np.log(K[k]))
    MC1[k] = np.mean(Y1[k])
    Var1[k] =  np.inner(Y1[k]-np.mean(Y1[k]), Y1[k]-np.mean(Y1[k]))/(M-1) 
    SE1[k] = np.sqrt(Var1[k]/M)
    RE1[k] = SE1[k]/MC1[k]

#IS MC
for k in range(3):
    mu2 = np.log(K[k])
    for i in range(M):
        X = np.random.normal(mu2,sigma*np.sqrt(T))
        Y2[k,i] = 1*(X<=np.log(K[k]))*np.exp((X-mu1/2-mu2/2)*(mu1-mu2)/(sigma*sigma*T))
    MC2[k] = np.mean(Y2[k])
    Var2[k] =  np.inner(Y2[k]-np.mean(Y2[k]), Y2[k]-np.mean(Y2[k]))/(M-1) 
    SE2[k] = np.sqrt(Var2[k]/M)
    RE2[k] = SE2[k]/MC2[k]
