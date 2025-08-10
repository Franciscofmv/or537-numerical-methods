function simple_euler_example
f = @odefun
x = 0:0.1:4; % x vector
y0 = 1; % initial condition at x0 and x0 = 0 in this case
y = euler(f,x,y0)
plot(x,y)
end

function [dydx] = odefun(x,y)
dydx = -2*x^3+12*x^2-20*x+8.5;
end

function [y] = euler(f,x,y0) % y0 is initial condition
n = length(x); % for loop is equal to the length of x
y = zeros(1,n); % preallocating y in memory
y(1) = y0; % initial condition
    for i=1:n-1
        dx = x(i+1)-x(i); % in case of not having a linear vector
        y(i+1) = y(i) + f(x(i),y(i))*dx; 
    end
    
end