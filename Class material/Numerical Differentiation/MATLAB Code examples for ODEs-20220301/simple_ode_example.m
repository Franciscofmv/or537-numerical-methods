function simple_ode_example
f = @simple_function;
y0 = 1;
x = (0:0.1:4);
figure;
y = euler(f,x,y0);
plot(x,y,'-ro');
hold on;
y = heuns(f,x,y0);
plot(x,y,'-go');
hold on;
y = mid_point(f,x,y0);
plot(x,y,'-yo');
hold on;
y = runge_kutta(f,x,y0);
plot(x,y,'-bo');
hold on;
[x,y] = ode45(f,x,y0);
plot(x,y,'-ms');
hold on;
f=@(x) -x.^4/2+4*x.^3-10*x.^2+8.5*x+1;
fplot(f,[0 4],'-k');
legend('euler','heuns','mid point','runge kutta','ode45','analytical');
hold off;
return

function [y] = euler(f,x,y0)
n = length(x);
y=zeros(1,n);
y(1) = y0;
for i=1:n-1
    xi=x(i);
    yi=y(i);
    h=x(i+1)-x(i);
    k1=f(xi,yi);
    y(i+1)=yi+h*k1;
end
return

function [y] = heuns(f,x,y0)
n = length(x);
y=zeros(1,n);
y(1) = y0;
for i=1:n-1
    xi=x(i);
    yi=y(i);
    h=x(i+1)-x(i);
    k1=f(xi,yi);
    k2=f(xi+h,yi+k1*h);
    y(i+1)=yi+h*(k1+k2)/2;
end
return

function [y] = mid_point(f,x,y0)
n = length(x);
y=zeros(1,n);
y(1) = y0;
for i=1:n-1
    xi=x(i);
    yi=y(i);
    h=x(i+1)-x(i);
    k1=f(xi,yi);
    k2=f(xi+0.5*h,yi+0.5*k1*h);
    y(i+1)=yi+h*k2;
end
return

function [y] = runge_kutta(f,x,y0)
n = length(x);
y=zeros(1,n);
y(1) = y0;
for i=1:n-1
    h = x(i+1)-x(i);
    k1 = f(x(i),y(i));
    k2 = f(x(i)+0.5*h,y(i)+0.5*k1*h);
    k3 = f(x(i)+0.5*h,y(i)+0.5*k2*h);
    k4 = f(x(i)+h,y(i)+k3*h);
    y(i+1) = y(i) + h*(k1+2*k2+2*k3+k4)/6;
end
return



function [dydx] = simple_function(x,y)
dydx = -2*x^3+12*x^2-20*x+8.5;
return