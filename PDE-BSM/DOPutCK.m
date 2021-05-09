function price = DOPutCK(S0,K,r,T,sigma,Smax,M,N)
% set up grid and adjust increments if necessary
dS = Smax/M;
dt = T/N;
matval = zeros(M+1,N+1);
vetS = linspace(0,Smax,M+1)';
veti = vetS / dS;
vetj = 0:N;
% set up boundary conditions
matval(: ,N+1) = max(K-vetS,0);
matval(1,:) = 0 ;
matval(M+1,:) = 0;
% set up the coefficients matrix
alpha = 0.25*dt*( sigma^2*(veti.^2) - r*veti);
beta = -dt*0.5*( sigma^2*(veti.^2) + r);
gamma = 0.25*dt*( sigma^2*(veti.^2) + r*veti);
M1 = -diag(alpha(3:M) ,-1) + diag(1-beta(2:M)) - diag(gamma(2:M-1) ,1) ;
[L,U] = lu(M1);
M2 = diag(alpha(3:M) ,-1) + diag(1+beta(2:M)) + diag(gamma(2:M-1) ,1);
% solve the sequence of linear systems
for j=N:-1:1
    matval(2:M,j) = U \ (L \ (M2*matval(2:M,j+1)));
end
% return price, possibly by linear interpolation outside the grid
price = interp1(vetS, matval(:,1), S0);

