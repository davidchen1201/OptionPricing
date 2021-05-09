#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


##############################################################################
#Q1
hh1 = np.array([1,0.1,0.01])
hh2 = np.array([0.1,0.01,0.001])
n = 100000
S0 = 50 
r = 0.03
sigma = 0.2 
T = 1
K = 50

Delta = np.zeros(hh1.size)
SE1 = np.zeros(hh1.size)
Vega = np.zeros(hh2.size)
SE2 = np.zeros(hh1.size)

for num in range(hh1.size):
    h = hh1[num]
    X = np.zeros(n)
    Y = np.zeros(n)
    H = np.zeros(n)
    for i in range(n):
        Z= np.random.normal(0,1)
        X[i]= (S0+h)*np.exp((r-0.5*sigma**2)*T+sigma*np.sqrt(T)*Z)
        Y[i]= (S0-h)*np.exp((r-0.5*sigma**2)*T+sigma*np.sqrt(T)*Z)
        H[i] = 1/(2*h)*np.exp(-r*T)*(np.maximum(K-X[i],0)- np.maximum(K-Y[i],0))
    Delta[num] = np.mean(H)
    SE1[num] = np.sqrt(np.var(H)/(n-1))
print('Delta is', Delta)    
print('standard error is', SE1)

for num in range(hh2.size):
    h = hh2[num]
    X = np.zeros(n)
    Y = np.zeros(n)
    H = np.zeros(n)
    for i in range(n):
        Z= np.random.normal(0,1)
        X[i]= S0*np.exp((r-0.5*(sigma+h)**2)*T+(sigma+h)*np.sqrt(T)*Z)
        Y[i]= S0*np.exp((r-0.5*(sigma-h)**2)*T+(sigma-h)*np.sqrt(T)*Z)
        H[i] = 1/(2*h)*np.exp(-r*T)*(np.maximum(X[i]-K,0)- np.maximum(Y[i]-K,0))
    Vega[num] = np.mean(H)
    SE2[num] = np.sqrt(np.var(H)/(n-1))
print('Vega is', Vega)    
print('standard error is', SE2)
    
##############################################################################
#Q3
n = 100000
S0 = 50 
r = 0.03
Sigma = np.array([0.1,0.3,0.5])
T = 1
K = 50

Delta = np.zeros(Sigma.size)
SE = np.zeros(Sigma.size)

for num in range(Sigma.size):
    sigma = Sigma[num]
    print(sigma)
    D = np.zeros(n)
    for i in range(n):
        Z= np.random.normal(0,1)
        ST= S0*np.exp((r-0.5*sigma**2)*T+sigma*np.sqrt(T)*Z)
        D[i] = np.exp(-r*T)*(ST>=K)*Z/(S0*sigma*np.sqrt(T))
    Delta[num] = np.mean(D)
    SE[num] = np.sqrt(np.var(D)/(n-1))
print('Delta is', Delta)    
print('standard error is', SE)