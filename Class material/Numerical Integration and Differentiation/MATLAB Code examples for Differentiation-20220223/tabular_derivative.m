function [dydx] = tabular_derivative(x,y) % x is IV y=f(x)
n = length(x); dydx = zeros(1,n);
% forward difference for first point
dydx(1) = (y(2)-y(1))/(x(2)-x(1));
% central difference for interior points
for i=2:n-1
    dydx(i) = (y(i+1)-y(i-1))/(x(i+1)-x(i-1));
end
% backward difference for last point
dydx(n) = (y(n)-y(n-1))/(x(n)-x(n-1));
return