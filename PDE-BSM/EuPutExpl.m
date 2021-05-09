function price = EuPutExpl(S0, K, r, T, sigma, Smax, M, N)
% set up grid and adjust increments if necessary
dS = Smax/M;
dt = T/N;
matval = zeros(M+1,N+1);
vetS = linspace(0,Smax,M+1)';
veti = 0:M;
vetj = 0:N;
% set up boundary conditions
matval(:,N+1) = max(K-vetS,0);
matval(1, :) = K*exp(-r*dt*(N-vetj)) ;
matval(M+1,:) = 0;
% set up coefficients
a = 0.5*dt*(sigma^2*veti - r).*veti;
b = 1- dt*(sigma^2*veti.^2 + r);
c = 0.5*dt*(sigma^2*veti + r).*veti;
% solve backward in time
for j = N:-1:1
    for i=2:M
        matval(i,j) = a(i)*matval(i-1,j+1) + b(i)*matval(i, j+1)+ c(i)*matval(i+1,j+1);
    end
end
% return price, possibly by linear interpolation outside the grid
price = interp1(vetS, matval(:,1), S0);
