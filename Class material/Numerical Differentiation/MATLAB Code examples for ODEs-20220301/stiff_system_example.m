function [] = stiff_system_linear_example
T=[0:0.05:2];
y0=[52.29; 83.82];
A=[-5 3; 100 -301];
y = backward_euler(T,A,y0);
% [y] = runge_kutta(@stiff_system_function,T,y0);
Y1 = y(1,:);
Y2 = y(2,:);
clear y;
syms t y;
s=dsolve('Dy1 = -5*y1+3*y2','Dy2 = 100*y1-301*y2','y1(0) = 52.29','y2(0) = 83.82','t');
t=[0:0.05:2];
y1=subs(s.y1,t);
y2=subs(s.y2,t);
plot(T,Y1,'bo',T,Y2,'ro',t,y1,'b',t,y2,'r');
legend('numerical-y1','numercial-y2','analytical-y1','analytical-y2');
hold off;
return

function [y] = backward_euler(t,A,y0)
n = length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
I=eye(m);
for i=1:n-1
    h = t(i+1)-t(i);
    y(:,i+1) = inv(I-h*A)*y(:,i);
end
return

function [y] = crank_nicolson(t,A,y0)
n = length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
I=eye(m);
for i=1:n-1
    h = t(i+1)-t(i);
    y(:,i+1) = (I-h*A)\y(:,i);
end
return

function [y] = runge_kutta(df,t,y0)
m = length(y0);
n = length(t);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = t(i+1)-t(i);
    ti=t(i);
    yi=y(:,i);
    k1 = feval(df,ti,yi);
    k2 = feval(df,ti+0.5*h,yi+0.5*k1*h);
    k3 = feval(df,ti+0.5*h,yi+0.5*k2*h);
    k4 = feval(df,ti+h,yi+k3*h);
    y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
end
return

% stiff system function
function [dydt] = stiff_system_function(t,y)
dydt=zeros(2,1);
dydt(1)=-5*y(1)+3*y(2);
dydt(2)=100*y(1)-301*y(2);
return