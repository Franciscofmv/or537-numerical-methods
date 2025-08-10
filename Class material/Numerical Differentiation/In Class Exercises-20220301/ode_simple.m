function ode_simple
f = @test_func;
x = 0:0.1:4;
y0 = 1;
y = euler(f,x,y0);
plot(x,y,'go');
hold on;
syms x y;
y = dsolve('Dy = -2*x^3+12*x^2-20*x+8.5','y(0) = 1','x');
x = 0:0.1:4;
y = double(subs(y));
plot(x,y,'-r');
hold off;
return

function [y] = euler(f,x,y0)
n = length(x);
y = zeros(1,n);
y(1) = y0;
for i=1:n-1
    dx = x(i+1)-x(i);
    y(i+1) = y(i) + dx*feval(f,x(i),y(i));
end
return

function [f] = test_func(x,y)
f = -2*x^3+12*x^2-20*x+8.5;
return