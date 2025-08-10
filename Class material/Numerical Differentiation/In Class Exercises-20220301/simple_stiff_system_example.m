function [] = simple_stiff_system_example
t=[0:0.1:1]; A = [-5 3; 100 -301];
y0=[52.29; 83.82];
[y] = simple_backward_euler_system(t,A,y0);
plot(t,y(1,:),'-go');
hold on;
f=@stiff_system_function;
[t,y] = ode45(f,t,y0);
plot(t,y(:,1),'-b+');
clear;
syms t y1 y2;
s=dsolve('Dy1 = -5*y1+3*y2','Dy2 = 100*y1 - 301*y2','y1(0)=52.29','y2(0) = 83.82');
t=[0:0.1:1];
y1=subs(s.y1,t);
plot(t,y1,'-rs');
legend('backward-euler','runge-kutta','analytical');
hold off;
return

function [y] = simple_backward_euler_system(t,A,y0)
n=length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
I=eye(m);
for i=1:n-1
    dt=t(i+1)-t(i);
    y(:,i+1)=inv(I-dt*A)*y(:,i);
end
return

function [y] = runge_kutta(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    k1=feval(f,x(i),y(:,i));
    k2=feval(f,x(i)+dx/2,y(:,i)+k1*dx/2);
    k3=feval(f,x(i)+dx/2,y(:,i)+k2*dx/2);
    k4=feval(f,x(i)+dx,y(:,i)+k3*dx);
    y(:,i+1)=y(:,i)+(k1+2*k2+2*k3+k4)*dx/6;
end
return

function [dydt] = stiff_system_function(t,y)
dydt=zeros(2,1);
dydt(1) = -5*y(1)+3*y(2);
dydt(2) = 100*y(1)-301*y(2);
return