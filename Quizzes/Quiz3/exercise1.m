
function exercise1

sigma = 10^(3)*[87.8 96.6 176 263 351 571 834 1229 1624 2107 2678 3380 4528];
eta = 10^(-3)*[153 204 255 306 357 408 459 510 561 612 663 714 765];
length(sigma)==length(eta)
format long;
a = tabular_derivative(eta,sigma);
myout = matrix(a,eta,sigma);
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