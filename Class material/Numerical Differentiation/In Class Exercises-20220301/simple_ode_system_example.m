% simple ode example
function [] = simple_ode_system_example
f=@test_func;
x = [0:0.2:4]; y0=[4; 6];
[y] = forward_euler(f,x,y0);
plot(x,y(1,:),'-go',x,y(2,:),'-go');
hold on;
[y] = heuns(f,x,y0);
plot(x,y(1,:),'-bd',x,y(2,:),'-bd');
hold on;
[y] = mid_point(f,x,y0);
plot(x,y(1,:),'-y^',x,y(2,:),'-y^');
hold on;
[y] = runge_kutta(f,x,y0);
plot(x,y(1,:),'-k+',x,y(2,:),'-k+');
hold on;
clear;
syms x y1 y2;
s=dsolve('Dy1 = -0.5*y1, Dy2 = 4-0.3*y2-0.1*y1','y1(0) = 4, y2(0) = 6','x');
x=[0:0.2:4];
y1=subs(s.y1);
y2=subs(s.y2);
% 
plot(x,y1,'-rs',x,y2,'-rs');
legend('euler','heuns','mid-point','runge_kutta','analytical');
hold off;
return

function [dydx] = test_func(x,y)
dydx=zeros(2,1);
dydx(1)=-0.5*y(1);
dydx(2)=4-0.3*y(2)-0.1*y(1);
return

function [y] = forward_euler(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    y(:,i+1)=y(:,i)+feval(f,x(i),y(:,i))*dx;
end
return

function [y] = heuns(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    ytemp=y(:,i)+feval(f,x(i),y(:,i))*dx;
    y(:,i+1)=y(:,i)+(feval(f,x(i),y(:,i))+feval(f,x(i+1),ytemp))*dx/2;
end
return

function [y] = mid_point(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    xtemp=(x(i+1)+x(i))/2;
    ytemp=y(:,i)+feval(f,x(i),y(:,i))*dx/2;
    y(:,i+1)=y(:,i)+feval(f,xtemp,ytemp)*dx;
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



