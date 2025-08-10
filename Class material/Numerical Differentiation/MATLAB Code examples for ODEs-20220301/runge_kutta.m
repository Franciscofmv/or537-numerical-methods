% extension of original runge_kutta for a system of ODE's

function [y] = runge_kutta(f,x,y0)
n=length(x);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    k1 = f(xi,yi);
    k2 = f(xi+0.5*h,yi+0.5*k1*h);
    k3 = f(xi+0.5*h,yi+0.5*k2*h);
    k4 = f(xi+h,yi+k3*h);
    y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6; % this is a matrix
end
end