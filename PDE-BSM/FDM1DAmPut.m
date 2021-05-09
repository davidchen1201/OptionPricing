function [P,Sf] = FDM1DAmPut(S,E,r,T,sigma,q,Smin,Smax,M,N,theta)
% Finite difference method for American put
f7 = 1;
[Am,Ap,B,Svec] = GridU1D(r,q,T,sigma,Smin,Smax,M,N,theta);
V0 = max(E - Svec,0); Payoff = V0(f7+(1:N-1));
Sf = zeros(M+1,1);
[L,U] = lu(Am); % Solve linear system for successive time steps
for j = M-1:-1:0
    V1 = V0;
    V0(f7+0) = E;
    b = Ap*V1(f7+(1:N-1)) + theta *B*V0(f7+[0 N]) + (1-theta)*B*V1(f7+[0 N]);
    solunc = U\(L\b);
    [V0(f7+(1:N-1)),Imax] = max([Payoff solunc],[],2);
    p = find(Imax == 2);
    i = p(1) - 1; % Grid line of last point below payoff
    Sf(f7+j) = interp1(solunc(i:i+2)-Payoff(i:i+2),Svec(f7+(i:i+2)),0);
end
P = interp1(Svec,V0,S);

function [Am,Ap,B,S] = GridU1D(r,q,T,sigma, Smin,Smax,M,N,theta)
% Uniform grid in original variables
f7 = 1;
dt = T/M;
S = linspace(Smin,Smax,N+1)';
dS = (Smax - Smin) / N;
a = (dt*sigma^2*S(f7+(1:N-1)).^2) ./ (2*dS^2);
b = (dt*(r-q)*S(f7+(1:N-1))) / (2*dS);
d = a - b;
m = - 2*a - dt*r;
u = a + b;
P = spdiags([d m u], 0:2, N-1, N+1);
Am = speye(N-1) - theta *P(:,f7+(1:N-1));
Ap = speye(N-1) + (1-theta)*P(:,f7+(1:N-1));
B = P(:,f7+[0 N]);
