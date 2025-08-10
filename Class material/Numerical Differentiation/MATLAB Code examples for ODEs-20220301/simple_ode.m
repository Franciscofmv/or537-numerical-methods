% simple ODE example
function simple_ode
f = @parachute;
v0=0;
t=[0:0.5:12];
v=rk4(f,t,v0);
plot(t,v,'-rs');
hold on;
g=9.8; c=12.5; m=68.1;
fplot(@(t) 9.8*68.1*(1-exp(-12.5*t/68.1))/12.5,[0,12]);
hold off;
return

function [y] = euler(f,x,y0)
n=length(x);
y=zeros(1,n);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    y(i+1)=y(i)+dx*f(x(i),y(i));
end
return

function [y] = heuns(f,x,y0)
n=length(x);
y=zeros(1,n);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    ytemp=y(i)+dx*f(x(i),y(i));
    y(i+1)=y(i)+(f(x(i),y(i))+f(x(i+1),ytemp))*dx/2;
end
return

function [y] = mid_point(f,x,y0)
n=length(x);
y=zeros(1,n);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    ytemp=y(i)+dx*f(x(i),y(i))/2;
    y(i+1)=y(i)+f(x(i)+dx/2,ytemp)*dx;
end
return

function [y] = rk4(f,x,y0)
n=length(x);
y=zeros(1,n);
y(1)=y0;
for i=1:n-1
    h=x(i+1)-x(i);
    k1=f(x(i),y(i));
    k2=f(x(i)+0.5*h,y(i)+0.5*k1*h);
    k3=f(x(i)+0.5*h,y(i)+0.5*k2*h);
    k4=f(x(i)+h,y(i)+k3*h);
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)*h/6;
end
return

function [dydx] = test_fnc(x,y)
dydx = -2*x^3+12*x^2-20*x+8.5;
return

function [dvdt] = parachute(t,v)
g=9.8; c=12.5; m=68.1;
dvdt = g-c*v/m;
return

function [dhdt] = reservoir(t,h)
H = [6 5 4 3 2 1];
A = [1.17 0.97 0.67 0.45 0.32 0.18 0];
Ah = interp1(H,A,h);
e=1; d =0.25;
dhdt = -pi*d^2*sqrt(2*g*h+e)/(4*Ah);
return
