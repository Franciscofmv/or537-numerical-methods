function bvp_nonlinear_heat_example
x=0:0.5:10;
yinit=[0 100];
% matlab bvp solution
solinit = bvpinit(x,yinit);
sol = bvp4c(@heat, @heatbc, solinit);
T = deval(sol,x);
plot(x,T(1,:),'-bs');
hold on;
% fda newton
T0=40; TL=200;
T=heat_newton_fda(x,T0,TL); % analytical jacobian
plot(x,T,'-go');
hold on;
% T=heat_newton_numerical_jacobian_fda(@heat_df,x,T0,TL); %numerical jacobian
% shooting method
dT0=1;
y=bvp_shooting(@heat,x,T0,TL,dT0);
T=y(1,:);
plot(x,T,'-r+');
legend('matlab-bvp4c','fda-newton','shooting');
hold off;
return



function [y] = bvp_shooting(odefun,x,yi,yf,dy0)
n=length(x);
y0=zeros(2,1);
y0(1)=yi;
y0(2) = dy0; % initial guess for y0(2)
f=@(dy0)bvp_shoot(dy0,odefun,x,yi,yf);
y0(2)=fsolve(f,dy0); % root finding method
y0(2)=fsolve(@bvp_shoot,dy0,optimset('fsolve'),odefun,x,yi,yf);
y=runge_kutta(odefun,x,y0);
% comment the line above and uncomment 2 lines below to use ode45
% [X,Y]=ode45(odefun,x,y0);
% y = Y';
return


% need to write a func which has x, yi ->left bc, yf->right bc
function [res] = bvp_shoot(dy0,odefun,x,yi,yf) 
n=length(x);
x0=x(1);
xL=x(n);
y0=zeros(2,1);
y0(1)=yi;
y0(2)=dy0;
y=runge_kutta(odefun,x,y0); % can use ode45 instead
res=yf-y(1,n); % residua 
% comment the 2 lines above and uncomment 2 lines below to use ode45
% [X,Y]=ode45(odefun,x,y0);
% res=yf-Y(n,1);
return

function [res] = heatbc(ya,yb)
res = [ya(1)-40; yb(1)-200];
return

function [dydx] = heat(x,y)
Ta = 20;
alpha = 5e-8;
dydx=zeros(2,1);
dydx(1) = y(2);
dydx(2) = alpha*(y(1)-Ta)^4;
return

function [T] = heat_newton_fda(x,T0,TL)
maxitr=100; % maximum iterations
Ta=20; alpha=5e-8; dx=x(2)-x(1);
m=length(x);
n=m-2;
ld=ones(n,1);
ud=ones(n,1);
d=zeros(n,1);
b=zeros(n,1);
Tk=zeros(n,1);
T=zeros(m,1);
epsilon=inf;
tol=1e-3;
iter=0;
while (epsilon > tol && iter <= maxitr)
    d(1:n)=-2-4*alpha*dx^2*(Tk(1:n)-Ta).^3; % middle diagonal, diagonal of jacobian
    b(1)=Tk(2)-2*Tk(1)-alpha*dx^2*(Tk(1)-Ta)^4+T0;% boundary conditions
    b(n)=TL-2*Tk(n)-alpha*dx^2*(Tk(n)-Ta)^4+Tk(n-1); % boundary conditions
    for i=2:n-1
        b(i)=Tk(i+1)-2*Tk(i)-alpha*dx^2*(Tk(i)-Ta)^4+Tk(i-1); % interior point
    end
    J=my_spdiags([ld d ud],[-1,0,1]); % jacobian
    dT=-J\b;                           % correction vector
    epsilon=norm(dT,inf); % infinity norm of dt
    Tk=Tk+dT;
    iter=iter+1;
end
T(2:m-1)=Tk(1:n);
T(1)=T0;
T(m)=TL;
fprintf('number of iterations = %g\n',iter);
fprintf('residual error = %g\n',epsilon);
% plot(x,T);
return

function [T] = heat_semi_implicit_fda(x,T0,TL)
maxitr=100;
Ta=20; alpha=5e-8; dx=x(2)-x(1);
m=length(x);
n=m-2;
ld=ones(n,1);
ud=ones(n,1);
d=zeros(n,1);
b=zeros(n,1);
Tk=zeros(n,1);
T=zeros(m,1);
y=zeros(n,1);
epsilon=inf;
iter=0;
while (epsilon > 1e-3 && iter <= maxitr)
    d(1:n)=-2-alpha*dx^2*(Tk(1:n)-Ta).^3;
    b(1)=-T0-alpha*dx^2*Ta*(Tk(1)-Ta)^3;
    b(n)=-TL-alpha*dx^2*Ta*(Tk(n)-Ta)^3;
    for i=2:n-1
        b(i)=-alpha*dx^2*Ta*(Tk(i)-Ta)^3;
    end
    A=my_spdiags([ld d ud],[-1,0,1]); 
    dy=A\b-Tk;                           
    epsilon=norm(dy);
    Tk=Tk+dy;
    iter=iter+1;
end
T(2:m-1)=Tk(1:n);
T(1)=T0;
T(m)=TL;
% fprintf('number of iterations = %g\n',iter);
% fprintf('residual error = %g\n',epsilon);
return

function [T] = heat_newton_numerical_jacobian_fda(f,x,T0,TL)
maxitr=100;
Ta=20; alpha=5e-8; dx=x(2)-x(1);
m=length(x);
n=m-2;
b=zeros(n,1);
Tk=zeros(n,1);
T=zeros(m,1);
y=zeros(n,1);
epsilon=inf;
iter=0;
params(1)=dx; params(2)=alpha; params(3)=Ta;
while (epsilon > 1e-3 && iter <= maxitr)
    b=feval(f,Tk,params);
    b(1)=b(1)+T0;
    b(n)=b(n)+TL;
    J=my_jacobian(f,Tk,params); % jacobian
    dy=-J\b;                           % correction vector
    epsilon=norm(dy,inf);
    Tk=Tk+dy;
    iter=iter+1;
end
T(2:m-1)=Tk(1:n);
T(1)=T0;
T(m)=TL;
% fprintf('number of iterations = %g\n',iter);
% fprintf('residual error = %g\n',epsilon);
% plot(x,T);
return

function [J] = my_jacobian(f,x,params) 
m=length(x); h=1e-3;
for i = 1:m
    x1 = x;
    x2 = x;
    x1(i) = x(i)+h;
    x2(i) = x(i)-h;
    J(:,i) = (feval(f,x1,params)-feval(f,x2,params))/(2*h);
end
return

function [E] = heat_df(T,params)
n=length(T);
E=zeros(n,1);
dx=params(1); alpha=params(2); Ta=params(3);
for i=2:n-1
    E(i)=T(i-1)+(-2*T(i)-alpha*dx^2*(T(i)-Ta)^4)+T(i+1);
end
E(1)=(-2*T(1)-alpha*dx^2*(T(1)-Ta)^4)+T(2);
E(n)=T(n-1)+(-2*T(n)-alpha*dx^2*(T(n)-Ta)^4);
return


