function [] = ode_examples
parachute_example;
generic_example;
pendulum_example;
return

function [] = parachute_example
hold off;
t = [0:12];
y0 = 0;
[y] = runge_kutta('parachute',t,y0);
figure(1);
plot(t,y,'b');
[T,Y]=ode45(@parachute,t,y0);
hold on;
plot(T,Y,'r');
hold on;
clear
syms g c m v vt
vt = dsolve('Dv=g-c*v/m','v(0)=0');
g=9.8; m=68.1; c=12.5;
% substitue g,m,c values
vt = subs(vt);
t = [0:12];
y = subs(vt);
%figure(1)
plot(t,y,'g');
legend('numerical','matlab','analytical');
hold off;
return

function [] = generic_example
hold off;
x = [0:0.1:4];
y0 = 1;
[y] = runge_kutta('generic_function',x,y0);
figure(2);
plot(x,y,'bo');
[x,y]=ode45(@generic_function,x,y0);
hold on;
plot(x,y,'r');
hold on;
clear
syms x y;
x = [0:0.1:4];
y = dsolve('Dy = -2*x^3+12*x^2-20*x+8.5','y(0)=1','x');
y = subs(y);
%figure(1)
plot(x,y,'g');
legend('numerical','matlab','analytical');
hold off;
return

function [] = pendulum_example
% pendulum problem
y0 = [pi/4 0]';
t0 = 0;
tn = 4;
n = 100;
t = linspace(t0,tn,n);
[y] = runge_kutta('pendulum',t,y0);
figure(3);
plot(t,y(1,:),'r',t,y(2,:),'b');
hold on;
% solution using ode45
[T,Y]=ode45(@pendulum,[t0 tn],y0);
plot(T,Y(:,1),'--g',T,Y(:,2),'--y');
xlim([0 4]);
legend('numerical-y1','numerical-y2','matlab-y1','matlab-y2');
hold off;
return

% and 'n' is the number of points at which the solution is evaluated
function [y] = runge_kutta(f,x,y0)
m = length(y0);
n = length(x);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h=x(i+1)-x(i);
    y(:,i+1)=y(:,i)+feval('rk4',f,h,x(i),y(:,i));
end
return

% 4th order Runge-Kutta function
% corrector is a column vector of size 'm'
% we have to use a column vector (ack..!)to be consistent
%               with matlab's ode45 function which requires
%               that f returns a column vector
function [corrector] = rk4(f,h,xi,yi)
k1=feval(f,xi,yi);
k2=feval(f,xi+0.5*h,yi+0.5*k1*h);
k3=feval(f,xi+0.5*h,yi+0.5*k2*h);
k4=feval(f,xi+h,yi+k3*h);
corrector=(k1+2*k2+2*k3+k4)*h/6;
return

% pendulum function
function [dydt] = pendulum(t,y)
g=32.2; l=2;
dydt=zeros(2,1);
dydt(1)= y(2);
dydt(2)= -g*sin(y(1))/l;
return

function [dv]=parachute(t,v)
g=9.8; c=12.5; m=68.1;
dv=g-c*v/m;
return

function [dydx] = generic_function(x,y)
dydx = -2*x^3+12*x^2-20*x+8.5;
return