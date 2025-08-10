% soil column problem
% ODE approach
% Dankwertz b.c's
% B.C's: x=0, VC-DdC/dx=VCin; x=L, dc/dx=0;
function dankwerts_column_ode
D=100; Cin=100; v=1;
L=10; N=100; T=50;
dx=L/N;
x = (0:dx:L); 
n=N+1;
cinit = zeros(1,n);
odefun = @(t,c) column_ode(t,c,Cin,D,v,dx,n);
tspan = [0 T];
[~,C] = ode15s(odefun,tspan,cinit);
c = C(end,:);
plot(x,c,'bo');
legend('ODE approach');
hold off;
%%%%%%%
return

function [dcdt] = column_ode(t,c,Cin,D,v,dx,n)
dcdt = zeros(n,1);
c0=c(2)+2*v*dx*(Cin-c(1))/D;
dcdt(1) = D*(c(2) - 2*c(1)+c0)/dx^2 - v*(c(2) - c0)/(2*dx);
for i=2:n-1
dcdt(i) = D*(c(i+1) - 2*c(i)+c(i-1))/dx^2 - v*(c(i+1) - c(i-1))/(2*dx);
end
cn_plus_one=c(n-1);
dcdt(n) = D*(cn_plus_one - 2*c(n)+c(n-1))/dx^2 - v*(cn_plus_one - c(n-1))/(2*dx); % constant concentration at right end (dcdt=0)
return