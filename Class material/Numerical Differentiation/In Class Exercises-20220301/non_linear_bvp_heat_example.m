function [] = non_linear_bvp_heat_example
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
legend('fda','bvp4c');
hold off;
return

function [dydx] = heat_ode(x,y)
alpha=5e-8; Ta=20;
dydx=zeros(2,1);
dydx(1)=y(2);
dydx(2)=alpha*(y(1)-Ta)^4;
return

function [res] = heat_bc(ya,yb);
res = [ya(1)-40; yb(1)-200];
return

function [T] = heat_fda(x,T0,TL)
% finite difference with Newton's method
m=length(x); T=zeros(m,1);
alpha=5e-8; Ta=20; dx=x(2)-x(1);
n=m-2;
Ttemp=zeros(n,1);
eps=1e-6; tol=inf; maxitr=100; k=0;
while (tol > eps & k <= maxitr)
    k=k+1;
    [J,b] = feval('heat_jacobian',Ttemp,alpha,dx,Ta,T0,TL);
    DT=J\b;
    Ttemp=Ttemp+DT;
    tol=norm(DT);
end
T(2:m-1)=Ttemp;
T(1)=T0;
T(m)=TL;
return

function [J,b] = heat_jacobian(T,alpha,dx,Ta,T0,TL)
% calculate jacobian and rhs for fda newton
n=length(T);
e=ones(n,1);
ld=e;
d=-2*e-4*alpha*dx^2*(T-Ta).^3;
ud=e;
J=spdiags([ld d ud],-1:1,n,n);
b=zeros(n,1);
for i=2:n-1
    b(i)=-T(i-1)+2*T(i)+alpha*dx^2*(T(i)-Ta)^4-T(i+1);
end
b(1)=-T0+2*T(1)+alpha*dx^2*(T(1)-Ta)^4-T(2);
b(n)=-T(n-1)+2*T(n)+alpha*dx^2*(T(n)-Ta)^4-TL;
return

