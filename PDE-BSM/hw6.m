K = 50;
r = 0.05;
sigma1 = 0.2;
T = 0.5;
Smin  = 0;
Smax = 100;
S0 = linspace(5,15,30);

BSM1 = zeros(1,30);
for i = 1:30
    BSM1(i) = BSMPut(S0(i),K,T,r,sigma1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q1(a)
M1 = 100;
N1 = 30;
Exp1 = zeros(1,30);
Imp1 = zeros(1,30);
CK1 = zeros(1,30);
for i = 1:30
    Exp1(i) = EuPutExpl(S0(i), K, r, T, sigma1, Smax, M1, N1);
    Imp1(i) = EuPutImpl(S0(i),K,r,T,sigma1,Smax,M1,N1);
    CK1(i) = DOPutCK(S0(i),K,r,T,sigma1,Smax,M1,N1);
end
Errorexp1 =  abs(BSM1 - Exp1);
Errorimp1 = abs(BSM1 - Imp1);
Errorck1 = abs(BSM1 - CK1); 

figure(1)
plot(S0,Errorexp1,'g-^')
hold on;
plot(S0,Errorimp1,'b-o')
hold on;
plot(S0,Errorck1,'r-*')
hold on;
legend({'Explicit Scheme','Implicit Scheme','Crank-Nicolson Scheme'},'Location','NorthWest');
x=xlabel('$S_0$');
y=ylabel('$Error$');
z=title('$M=100, N=30, \sigma=0.2$');
set(x,'Interpreter','latex');
set(y,'Interpreter','latex');
set(z,'Interpreter','latex');
print(gcf,'-dpng','hw6Q1a.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q1(b)
M2 = 200;
N2 = 30;
Exp2 = zeros(1,30);
Imp2 = zeros(1,30);
CK2 = zeros(1,30);
for i = 1:30
    Exp2(i) = EuPutExpl(S0(i), K, r, T, sigma1, Smax, M2, N2);
    Imp2(i) = EuPutImpl(S0(i),K,r,T,sigma1,Smax,M2,N2);
    CK2(i) = DOPutCK(S0(i),K,r,T,sigma1,Smax,M2,N2);
end
Errorexp2 = abs(BSM1 - Exp2);
Errorimp2 = abs(BSM1 - Imp2);
Errorck2 = abs(BSM1 - CK2);

figure(2)
plot(S0,Errorexp2,'g-^')
hold on;
plot(S0,Errorimp2,'b-o')
hold on;
plot(S0,Errorck2,'r-*')
hold on;
legend({'Explicit Scheme','Implicit Scheme','Crank-Nicolson Scheme'},'Location','NorthWest');
x=xlabel('$S_0$');
y=ylabel('$Error$');
z=title('$M=200, N=30, \sigma=0.2$');
set(x,'Interpreter','latex');
set(y,'Interpreter','latex');
set(z,'Interpreter','latex');
print(gcf,'-dpng','hw6Q1b.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q1(c)
M3 = 700;
N3 = 30;
Exp3 = zeros(1,30);
Imp3 = zeros(1,30);
CK3 = zeros(1,30);
for i = 1:30
    Exp3(i) = EuPutExpl(S0(i), K, r, T, sigma1, Smax, M3, N3);
    Imp3(i) = EuPutImpl(S0(i),K,r,T,sigma1,Smax,M3,N3);
    CK3(i) = DOPutCK(S0(i),K,r,T,sigma1,Smax,M3,N3);
end
Errorexp3 = abs(BSM1 - Exp3);
Errorimp3 = abs(BSM1 - Imp3);
Errorck3 = abs(BSM1 - CK3);

figure(3)
plot(S0,Errorexp3,'g-^')
hold on;
plot(S0,Errorimp3,'b-o')
hold on;
plot(S0,Errorck3,'r-*')
hold on;
legend({'Explicit Scheme','Implicit Scheme','Crank-Nicolson Scheme'},'Location','NorthWest');
x=xlabel('$S_0$');
y=ylabel('$Error$');
z=title('$M=700, N=30, \sigma=0.2$');
set(x,'Interpreter','latex');
set(y,'Interpreter','latex');
set(z,'Interpreter','latex');
print(gcf,'-dpng','hw6Q1c.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q1(d)
sigma2 = 0.4;
BSM2 = zeros(1,30);
for i = 1:30
    BSM2(i) = BSMPut(S0(i),K,T,r,sigma2);
end

M4 = 200;
N4 = 30;
Exp4 = zeros(1,30);
Imp4 = zeros(1,30);
CK4 = zeros(1,30);
for i = 1:30
    Exp4(i) = EuPutExpl(S0(i), K, r, T, sigma2, Smax, M4, N4);
    Imp4(i) = EuPutImpl(S0(i),K,r,T,sigma2,Smax,M4,N4);
    CK4(i) = DOPutCK(S0(i),K,r,T,sigma2,Smax,M4,N4);
end
Errorexp4 = abs(BSM2 - Exp4);
Errorimp4 = abs(BSM2 - Imp4);
Errorck4 = abs(BSM2 - CK4);

figure(4)
plot(S0,Errorexp4,'g-^')
hold on;
plot(S0,Errorimp4,'b-o')
hold on;
plot(S0,Errorck4,'r-*')
hold on;
legend({'Explicit Scheme','Implicit Scheme','Crank-Nicolson Scheme'},'Location','NorthWest');
x=xlabel('$S_0$');
y=ylabel('$Error$');
z=title('$M=200, N=30, \sigma=0.4$');
set(x,'Interpreter','latex');
set(y,'Interpreter','latex');
set(z,'Interpreter','latex');
print(gcf,'-dpng','hw6Q1d.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q2
S = 50; E = 50; r = 0.10; T = 5/12; sigma = 0.40; q = 0;
Smin = 0; Smax = 100; N = 100; M = 100; theta = 1/2;
[P,Sf] = FDM1DAmPut(S,E,r,T,sigma,q,Smin,Smax,M,N,theta);

figure(5)
plot((1:N)*dt,Sf(1:N))
xlim([0,0.5])
ylim([30,50])
x=xlabel('$t$');
y=ylabel('$S_f(t)$');
set(x,'Interpreter','latex');
set(y,'Interpreter','latex');
print(gcf,'-dpng','hw6Q2sol.png')
