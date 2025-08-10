function [] = ode_system_example
hold off;
x = [0:0.1:2];
y0 = [4; 6]; % column vector
[y] = runge_kutta4('generic_function',x,y0);
y1=y(1,:);
y2=y(2,:);
plot(x,y1,'bo',x,y2,'ro');
hold on;
% analytical solution
clear;
s=dsolve('Dy1 = -0.5*y1','Dy2 = 4-0.3*y2-0.1*y1', 'y1(0) = 4', 'y2(0) = 6','x');
x=[0:0.1:2];
y1=subs(s.y1);
y2=subs(s.y2);
% y1=4*exp(-0.5*x);
% y2=2*exp(-0.5*x)-9.33*exp(-0.3*x)+13.33;
plot(x,y1,'b',x,y2,'r');
legend('numerical-y1','numerical-y2','analytical-y1','analytical-y2');
hold off;
return


function [dydx] = generic_function(x,y)
dydx=zeros(2,1);
dydx(1)= -0.5*y(1);
dydx(2) = 4-0.3*y(2)-0.1*y(1);
return

function [y] = runge_kutta4(df,x,y0)
m = length(y0);
n = length(x);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    k1 = feval(df,xi,yi);
    k2 = feval(df,xi+0.5*h,yi+0.5*k1*h);
    k3 = feval(df,xi+0.5*h,yi+0.5*k2*h);
    k4 = feval(df,xi+h,yi+k3*h);
    y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
end
return

function [y] = heuns(df,x,y0)
n = length(x);
m = length(y0);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    y_temp = yi+h*feval(df,xi,yi);
    y(:,i+1) = yi + h*(feval(df,xi,yi)+feval(df,x(i+1),y_temp))/2;
end
return

function [y] = mid_point(df,x,y0)
n = length(x);
m = length(y0);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    y_temp = yi+h*feval(df,xi,yi)/2;
    y(:,i+1) = yi + h*feval(df,xi+h/2,y_temp);
end
return

function [y] = euler(df,x,y0)
n = length(x);
m = length(y0);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    y(:,i+1) = yi+h*feval(df,xi,yi);
end
return
