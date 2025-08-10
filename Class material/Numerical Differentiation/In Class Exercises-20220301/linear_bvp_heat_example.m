function [] = linear_bvp_heat_example
x=[0:0.5:10]';
T0=40; TL = 200;
[T] = heat_fda(x,T0,TL);
plot(x,T,'-go');
hold on;
x=[0:0.5:10];
solinit=bvpinit(x,[100 100]);
sol=bvp4c(@heat_ode,@heat_bc,solinit);
y=deval(sol,x);
T=y(1,:);
plot(x,T,'-rs');
hold on;
T=heat_shooting(@heat_ode,x,T0,TL);
plot(x,T,'-b^');
return

function [dydx] = heat_ode(x,y)
alpha=0.01; Ta=20;
dydx=zeros(2,1);
dydx(1)=y(2);
dydx(2)=alpha*(y(1)-Ta);
return

function [res] = heat_bc(ya,yb);
res = [ya(1)-40; yb(1)-200];
return

function [T] = heat_fda(x,T0,TL)
m=length(x); T=zeros(m,1);
alpha=0.01; Ta=20; dx=x(2)-x(1);
n=m-2;
e=ones(n,1);
ld=e;
d=e*(-2-alpha*dx^2);
ud=e;
b=e*(-alpha*Ta*dx^2);
b(1)=b(1)-T0;
b(n)=b(n)-TL;
A=spdiags([ld d ud], -1:1, n, n);
T(2:m-1)=A\b;
T(1)=T0;
T(m)=TL;
return

function [T] = heat_shooting(f,x,T0,TL)
m=length(x);
Ta=20; alpha=0.01;
% 1st attempt
q1=0;
[x,T1]=ode45(f,x,[T0 q1]);
p1=T1(m,1);
% 2nd attempt
q2=1;
[x,T2]=ode45(f,x,[T0 q2]);
p2=T2(m,1);
y2=q1+(TL-p1)*(q2-q1)/(p2-p1);
[x,Tx]=ode45(f,x,[T0 y2]);
T=Tx(:,1)';
return
