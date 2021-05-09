function price = AmPutCK(S0,K,r,T,sigma,Smax,M,N,omega,tol)
% M = round(Smax/dS); 
dS = Smax/M; % set up grid
% N = round(T/dt); 
dt = T/N;
oldval = zeros(M-1,1); % vectors for Gauss-Seidel update
newval = zeros(M-1,1) ;
vetS = linspace(0,Smax,M+1)';
veti = 0:M; 
vetj = 0:N;
% set up boundary conditions
payoff = max (K-vetS(2:M),0);
pastval = payoff; % values for the last layer
boundval = K*exp(-r*dt*(N-vetj)); % boundary values
% set up the coefficients and the right hand side matrix
alpha = 0.25*dt*( sigma^2*(veti.^2) - r*veti);
beta = -dt*0.5*( sigma^2*(veti.^2) + r );
gamma = 0.25*dt*( sigma^2*(veti.^2) + r*veti);
M2 = diag(alpha(3:M) ,-1) + diag(1+beta(2:M)) + diag(gamma(2:M-1) ,1);
% solve the sequence of linear systems by SOR method
aux = zeros (M-1,1) ;
for j=N:-1:1    
    aux(1) = alpha(2) * (boundval(1,j) + boundval(1,j+1));
    % set up right hand side and initialize
    rhs = M2*pastval(:) + aux;
    oldval = pastval;
    error = realmax;
    while tol < error
        newval(1) = max(payoff(1),oldval(1)+omega/(1-beta(2))*(rhs(1)-(1-beta(2))*oldval(1)+gamma(2)*oldval(2)));
        for k=2:M-2
            newval(k) = max(payoff(k),oldval(k)+omega/(1-beta(k+1))*(rhs(k)+alpha(k+1)*newval(k-1)-(1-beta(k+1))*oldval(k)+gamma(k+1)*oldval(k+1)));
        end
        newval(M-1) = max(payoff(M-1),oldval(M-1)+omega/(1-beta(M))*(rhs(M-1)+alpha(M)*newval(M-2)-(1-beta(M))*oldval(M-1)));
        error = norm(newval - oldval);
        oldval = newval;
    end
    pastval = newval;
end          
newval = [boundval(1);newval;0]; % add missing values
% return price, possibly by linear interpolation outside the grid
price = interp1(vetS, newval, S0);
