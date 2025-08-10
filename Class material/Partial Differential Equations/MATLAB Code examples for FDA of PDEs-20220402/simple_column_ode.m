% soil column problem
% ODE approach
% B.C's: x=0, C=C0; x=L, C=0;
function simple_column_ode
v=0.36; alpha=5; c0=2; cL = 0;
L=100; N=100; T=50;
D=alpha*v;
dx=L/N;
x = (0:dx:L); 
n=N+1;
cinit = zeros(1,n);
cinit(1)=c0; cinit(n)=cL;
odefun = @(t,c) column_ode(t,c,D,v,dx,n);
tspan = [0 T];
[~,C] = ode45(odefun,tspan,cinit);
c = C(end,:);
plot(x,c,'-bo');
hold on;
%%%%%%%%%%%%%%%
% analytical solution
f1=(x-v*T)/(2*sqrt(D*T));
f2=(x+v*T)/(2*sqrt(D*T));
c=0.5*c0*(erfc(f1)+exp(v*x/D).*erfc(f2));
plot(x,c,'-r');
legend('ODE approach','analytical');
hold off;
%%%%%%%
return

function [dcdt] = column_ode(t,c,D,v,dx,n)
dcdt = zeros(n,1);
dcdt(1) = 0; % constant concentration at left end (dcdt=0)
for i=2:n-1
dcdt(i) = D*(c(i+1) - 2*c(i)+c(i-1))/dx^2 - v*(c(i+1) - c(i-1))/(2*dx);
end
dcdt(n) = 0; % constant concentration at right end (dcdt=0)
return