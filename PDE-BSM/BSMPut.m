function price = BSMPut(S0,K,T,r,sigma)
d1 = (log(S0/K)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T);
price = K*exp(-r*T)*normcdf(-d2)-S0*normcdf(-d1);
