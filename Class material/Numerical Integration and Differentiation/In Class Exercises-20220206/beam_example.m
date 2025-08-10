% beam example (CC 24.20. p.685)
function beam_example
x = [0 0.375 0.75 1.125 1.5 1.875 2.25 2.625 3];
y = [0 -0.2571 -0.9484 -1.9689 -3.2262 -4.6414 -6.1503 -7.7051 -9.275]*1e-2;
E = 200e9; I = 0.0003;
slope = tabular_derivative(y,x);
moment = E*I*tabular_derivative(slope,x);
shear = tabular_derivative(moment,x);
distributed_load = -tabular_derivative(shear,x);
subplot(2,2,1);
plot(x,slope);
xlabel('x'); ylabel('slope'); title('slope vs x');
subplot(2,2,2);
plot(x,moment);
subplot(2,2,3);
plot(x,shear);
subplot(2,2,4);
plot(x,distributed_load);
return

% calculates derivatives at the same points supplied
function [dydx] = tabular_derivative(y,x)
n=length(x); dydx = zeros(1,n);
% use central for interior points
for i=2:n-1
    dydx(i)=(y(i+1)-y(i-1))/(x(i+1)-x(i-1));
end
% use forward for first point
dydx(1) = (y(2)-y(1))/(x(2)-x(1)); 
% use backward for last point
dydx(n)=(y(n)-y(n-1))/(x(n)-x(n-1));
return