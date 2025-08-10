function reservoir_draining
f=@depth;
h=6:-0.1:0;
h0=0;
t=runge_kutta4(f,h,h0);
n=length(h);
tf=t(n);
plot(t,h);
% [H,T]=ode45(f,[6 0],h0);
% n=length(H);
% tf=T(n);
fprintf('Drain time = %g s\n',tf);
return

% dam reservoir function
function [dtdh] = depth(h,t)
g=9.81; d=0.25; e=1;
H=[6, 5, 4, 3, 2, 1, 0];
A=[1.17 0.97, 0.67, 0.45, 0.32, 0.18, 0];
Ah=interp1(H,A,h);
dtdh=-4*Ah/(pi*d^2*sqrt(2*g*(h+e)));
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
