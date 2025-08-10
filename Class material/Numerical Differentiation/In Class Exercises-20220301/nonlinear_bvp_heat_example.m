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
y0(2)=fsolve(f,dy0);
% y0(2)=fsolve(@bvp_shoot,dy0,optimset('fsolve'),odefun,x,yi,yf);
y=ode_solver('runge_kutta',odefun,x,y0);
return

function [res] = bvp_shoot(dy0,odefun,x,yi,yf)
n=length(x);
x0=x(1);
xL=x(n);
y0=zeros(2,1);
y0(1)=yi;
y0(2)=dy0;
y=ode_solver('runge_kutta',odefun,x,y0);
res=yf-y(1,n);
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
epsilon=inf;
iter=0;
while (epsilon > 1e-3 && iter <= maxitr)
    d(1:n)=-2-4*alpha*dx^2*(Tk(1:n)-Ta).^3;
    b(1)=Tk(2)-2*Tk(1)-alpha*dx^2*(Tk(1)-Ta)^4+T0;
    b(n)=TL-2*Tk(n)-alpha*dx^2*(Tk(n)-Ta)^4+Tk(n-1);
    for i=2:n-1
        b(i)=Tk(i+1)-2*Tk(i)-alpha*dx^2*(Tk(i)-Ta)^4+Tk(i-1);
    end
    J=spdiags([ld d ud],[-1,0,1],n,n); % jacobian
    dT=-J\b;                           % correction vector
    epsilon=norm(dT,inf);
    Tk=Tk+dT;
    iter=iter+1;
end
T(2:m-1)=Tk(1:n);
T(1)=T0;
T(m)=TL;
% fprintf('number of iterations = %g\n',iter);
% fprintf('residual error = %g\n',epsilon);
return


% numerical jacobian
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





