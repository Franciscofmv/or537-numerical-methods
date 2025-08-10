% chapra and canale 21.12 (p.541)
function jet_fighter
t = [0 0.52 1.04 1.75 2.37 3.25 3.83];
x = [153 185 210 249 261 271 273];
v = tabular_derivative(t,x);
figure(1);
plot(t,v);
a = tabular_derivative(t,v);
figure(2);
plot(t,a);
return

function [dydx] = tabular_derivative(x,y)
n = length(x); dydx = zeros(1,n);
dydx(1) = (y(2)-y(1))/(x(2)-x(1));
for i=2:n-1
    dydx(i) = (y(i+1)-y(i-1))/(x(i+1)-x(i-1));
end
dydx(n) = (y(n)-y(n-1))/(x(n)-x(n-1));
return