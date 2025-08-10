function exercise7
p = [0 10 20 30 40 50 60];

W = [122 130 135 160 175 190 200];
n = length(p);
tabular_trapezoidal(p)
f=@integrand;
I=0;
f=@(z)integrand(z,p,W)
integral(f,0,p(2))
end

function [y] = integrand(z,d,w)
rho = 10^3;
g = 9.8;

y = rho.*g.*w.*(d.-z);
end
function [I] = tabular_trapezoidal(x,y) % x and y vectors
if (length(x)==length(y))
I = 0;
n = length(x);
for i=1:n-1
h = x(i+1)-x(i); % width of segment i
I = I + h*(y(i+1)+y(i))/2;
end
else
disp("Vectors x and y must be same length.")
end