# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import norm
import statsmodels.api as sm
import matplotlib.pyplot as plt

# Q2 Box-Muller method
sample = np.power(10,5)
X1 = np.zeros(sample)
X2 = np.zeros(sample)
i = 0
while i<sample:
    U1 = np.random.rand()
    U2 = np.random.rand()
    R = -2*np.log(U1)
    V = 2*np.pi*U2
    Z1 = np.sqrt(R)*np.cos(V)
    Z2 = np.sqrt(R)*np.sin(V)
    X1[i] = Z1
    X2[i] = Z2
    i = i+1

plt.hist(X1, density=True,bins = int(50))
plt.title('Histogram for Box-Muller Algorithm')
plt.xlabel(r'$Z_1$')
plt.ylabel('Frequency')
# plt.savefig('Histogram_BM1.png')
plt.savefig('Histogram_BM1.eps')
plt.show()
        
sm.qqplot(X1, line='s')
plt.title('Q-Q Plot for BM method')
# plt.savefig('QQ_BM1.png')
plt.savefig('QQ_BM1.eps')
plt.show()

plt.hist(X2, density=True,bins = int(50))
plt.title('Histogram for Box-Muller Algorithm')
plt.xlabel(r'$Z_2$')
plt.ylabel('Frequency')
# plt.savefig('Histogram_BM2.png')
plt.savefig('Histogram_BM2.eps')
plt.show()
        
sm.qqplot(X2, line='s')
plt.title('Q-Q Plot for BM method')
# plt.savefig('QQ_BM2.png')
plt.savefig('QQ_BM2.eps')
plt.show()

##############################################################################
#Q3
#Geometric brownian motion
def GBM(S0,T,r,sigma):
    Z= np.random.normal(0,1)
    ST = S0*np.exp((r-0.5*sigma*sigma)*T+sigma*np.sqrt(T)*Z)
    return ST

S0 = 50
T = 0.25
r = 0.05 
sigma = 0.3

#Q3(1)
Sample = 100
ST = np.zeros(Sample)
run = 0
for run in range(Sample):
    ST[run] = GBM(S0,T,r,sigma)   
    
plt.hist(ST)
plt.xlabel(r'$S_T$')
plt.ylabel('Count')
# plt.savefig('Histogram_ST.png')
plt.savefig('Histogram_ST.eps')
plt.show()

#Q3(2)
#Black-Scholes-Merton Put Option
def BSM_Put(S0,K,T,r,sigma):
    d1 = (np.log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*np.sqrt(T))
    d2 = d1-sigma*np.sqrt(T)
    put = K*np.exp(-r*T)*norm.cdf(-d2)-S0*norm.cdf(-d1)
    return put

M1 = np.power(10,5)
K = np.array([40, 50, 60])

Raw1 = np.zeros([3, M1])

for k in range(3):    
    for run in range(M1):
        Raw1[k, run] = np.maximum(K[k]-GBM(S0,T,r,sigma),0)

BSM = np.zeros(3)
MC1 = np.zeros(3)
Var = np.zeros(3)
cfl = np.zeros(3)
cfu = np.zeros(3)

for k in range(3):
    BSM[k] = BSM_Put(S0,K[k],T,r,sigma)
    MC1[k] = np.mean(Raw1[k,:])*np.exp(-r*T)
    Var[k] = np.var(Raw1[k,:])*np.exp(-2*r*T)
    cfl[k] = MC1[k]-1.96*np.sqrt(Var[k]/M1) #confidence interval lower bound
    cfu[k] = MC1[k]+1.96*np.sqrt(Var[k]/M1) #confidence interval upper bound
    
#Q3(3)
M2 = np.power(10,6)
Raw2 = np.zeros([3, M2])

for k in range(3):    
    for run in range(M2):
        Raw2[k, run] = np.maximum(K[k]-GBM(S0,T,r,sigma),0)

error = np.zeros([3,M2])
for k in range(3):    
    for m in range(M2):
        error[k,m] = np.mean(Raw2[k,:m+1])*np.exp(-r*T)-BSM[k]

plt.figure(figsize=(9,6))
plt.plot(range(M2), error[0,:], label='K=40')
plt.plot(range(M2), error[1,:], label='K=50')
plt.plot(range(M2), error[2,:], label='K=60')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('M',fontsize=15)
plt.ylabel('MC error',fontsize=15)
plt.legend(loc='upper right',fontsize=15)
plt.savefig('MCerrorlong.eps')
        
plt.figure(figsize=(9,6))
plt.plot(range(0,1000), error[0,:1000], label='K=40')
plt.plot(range(0,1000), error[1,:1000], label='K=50')
plt.plot(range(0,1000), error[2,:1000], label='K=60')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('M',fontsize=15)
plt.ylabel('MC error',fontsize=15)
plt.legend(loc='upper right',fontsize=15)
plt.savefig('MCerrorshort.eps')
