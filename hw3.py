# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm
from scipy.special import factorial

#Q1
r = 0.01 
sigma1 = 0.3 
X0 = 50  
sigma2 = 0.2
Y0 = 30 
rho = 0.1
T = 0.5 
K = 15

mean = np.array([0, 0])
cov = np.array([[1, rho], [rho, 1]]) 

sample = 100000
XT = np.zeros(sample)
YT = np.zeros(sample)
ZT = np.zeros(sample)
for s in range(sample):   
    Z1, Z2 = np.random.multivariate_normal(mean, cov)
    XT[s] = X0*np.exp((r-0.5*sigma1*sigma1)*T+sigma1*np.sqrt(T)*Z1)
    YT[s] = Y0*np.exp((r-0.5*sigma2*sigma2)*T+sigma2*np.sqrt(T)*Z2)
    ZT[s] = np.maximum(XT[s]-YT[s]-K,0)
    
print('The price of spread option is', np.mean(ZT)*np.exp(-r*T))      


##############################################################################
#Q2
sample = 100000
steps = 100
T = 0.25
S0 = 52
K = 50
r = 0.05
sigma = 0.3
Lambda = 0.25
a = 0
b = 1

def VCS(S0,sigma,T,r,K): 
    dt = T/steps
    Xt = np.zeros(steps+1)
    Xt[0] = np.log(S0)
    vcs = np.zeros(sample)
    for s in range(sample):   
        for i in range(steps):
            Z = np.random.normal(0,1)
            N = np.random.poisson(Lambda*dt)
            M = 0
            for j in range(N):
                M = M+np.random.normal(a,b*b)
            Xt[i+1] = Xt[i]+(r-0.5*sigma*sigma)*dt+sigma*np.sqrt(dt)*Z+M
        vcs[s] = np.maximum(np.exp(Xt[-1])-K,0)
    MC= np.mean(vcs)*np.exp(-r*T)
    return(MC)   

def BSMC(S,sigma,T,r,K):
    d1 = (np.log(S/K)+(r+sigma*sigma/2)*T)/(sigma*np.sqrt(T))
    d2 = d1-sigma*np.sqrt(T)
    call = S*norm.cdf(d1)-K*np.exp(-r*T)*norm.cdf(d2)
    return(call)

m = np.exp(a+0.5*b*b) #value of E[Y]
def VCT(S0,sigma,T,r,K):
    vct = 0
    for n in range(1000):
        vct = vct + np.exp(-Lambda*m*T)*np.power(Lambda*m*T,n)/(factorial(n))*BSMC(S0,np.sqrt(sigma*sigma+b*b*n/T),T,r-Lambda*(m-1)+n*np.log(m)/T,K)
    return(vct)
    

print('simulated option value is', VCS(S0,sigma,T,r,K))
print('analytical option value is', VCT(S0,sigma,T,r,K))

##############################################################################
#Q3
a = 0.2
b = 0.05
sigma = 0.1
d = 4*b*a/sigma/sigma
X0 = 0.04
T = 0.5
steps = 100
dt = T/steps
M = 100000

#exact simulation
XX1 = np.zeros(M)
for s in range(M):    
    X1 = np.zeros(steps+1)
    X1[0] = X0

    for n in range(steps):
        c = sigma**2*(1-np.exp(-a*dt))/(4*a)
        Lambda = X1[n]*np.exp(-a*dt)/c
        Z = np.random.normal(0,1)
        X = np.random.chisquare(d-1)
        X1[n+1] = c*((Z+np.sqrt(Lambda))**2+X)       
    XX1[s] = X1[-1]
print('exact simulation is', np.mean(XX1))

#Euler scheme
XX2 = np.zeros(M)
for s in range(M):    
    X2 = np.zeros(steps+1)
    X2[0] = X0

    for n in range(steps):
        Z = np.random.normal(0,1)
        X2[n+1]=X2[n]+a*(b-X2[n])*dt+sigma*np.sqrt(X2[n]*dt)*Z
    XX2[s] = X2[-1]
print('Euler scheme is', np.mean(XX2))

#discretization error
error = np.mean(XX2-XX1)
print('discretization error is', error)
