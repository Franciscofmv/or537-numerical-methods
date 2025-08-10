function [y] = ode_solver(method,f,x,y0)
y = feval(method,f,x,y0);
return

function [y] = runge_kutta(df,x,y0)
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
