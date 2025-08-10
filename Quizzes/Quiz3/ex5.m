
function exercise1

t = [0 0.1 0.2 0.3 0.5 0.7];
i = [0 0.16 0.32 0.56 0.84 2.0];
length(i)==length(t)
L=4;
format long;
a = L*tabular_derivative(t,i);
myout = matrix(a,t,i);
myout
end

function [dydx v] = tabular_derivative(x,y) % x is IV y=f(x)

n = length(x); dydx = zeros(1,n);

% forward difference for first point
dydx(1) = (y(2)-y(1))/(x(2)-x(1));

% central difference for interior points
for i=2:n-1
    dydx(i) = (y(i+1)-y(i-1))/((x(i+1)-x(i-1)));
end
% backward difference for last point

dydx(n) = (y(n)-y(n-1))/(x(n)-x(n-1));
end

function [v] = matrix(deriv,x,y)
n = length(x);
matr = zeros(n,3);
matr(1,1) = x(1);
matr(1,2)=y(1);
matr(1,3)=deriv(1);
for i=2:n-1
    matr(i,1) = x(i);
    matr(i,2) = y(i);
    matr(i,3) = deriv(i);
end
matr(n,1) = x(n);
matr(n,2) = y(n);
matr(n,3) = deriv(n);
v = matr;
end