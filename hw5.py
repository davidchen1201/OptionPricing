import numpy as np
from sklearn.linear_model import LinearRegression

#Q2
S1 = np.array([9.9528, 7.8115, 9.1293, 8.3600, 8.7787, 10.9398])
X1 = np.zeros([S1.size, 2])
X1[:,0] = np.copy(S1)
X1[:,1] = S1**2

Y1 = np.array([5.1898, 1.3346, 4.4997, 2.6834, 2.8888, 3.2714])
reg1 = LinearRegression().fit(X1, Y1)
print(reg1.coef_)
print(reg1.intercept_)
Y1_pre = reg1.predict(X1)

S2 = np.array([8.3826, 11.9899, 6.8064, 7.0508, 11.2214, 8.9672, 11.5336])
X2 = np.zeros([S2.size, 2])
X2[:,0] = np.copy(S2)
X2[:,1] = S2**2
Y2 = np.array([5.1898,0,4.1885,4.4997,3.6400,2.8888,3.2714])*np.exp(-0.03*1/3)
reg2 = LinearRegression().fit(X2, Y2)
print(reg2.coef_)
print(reg2.intercept_)
Y2_pre = reg2.predict(X2)

p = (5.1382+0+0+5.1936+4.9492+3.6038+2.8601+3.2389)/8*np.exp(-0.03*1/3)
