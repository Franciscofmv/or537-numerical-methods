%% C & C 21.23 (p. 631)
function tabular_data_example1
% define tabular data
T=[0 1 2 3.25 4.5 6 7 8 8.5 9 10];
V=[0 5 6 5.5 7 8.5 8 6 7 7 5];
%% method 1: simple trapezoidal tabular data interpolation
d_trap_data = trap_data(T,V)
%% method 2: linearly interpolate and integrate using trap
% define integrand function with linear interpolation
f=@(t) vel_linear_interp(t,T,V);
% integrate to find distance 'd' from trap function
a=0; b=10; n=100;
d_linear_interp = trap(f,a,b,n)
return

% tabular data integration using trapezoidal rule
function [I] = trap_data(t,v)
n=length(t); I=0;
for i=1:n-1
    h=t(i+1)-t(i);
    I=I+h*(v(i+1)+v(i))/2;
end

% integrand function with linear interpolation
function [v] = vel_linear_interp(t,T,V)
v=interp1(T,V,t);
return


% trapezoidal integration of a function
function [I] = trap(f,a,b,n)
h = (b-a)/n; I=0;
for i=0:n
    x = a + i*h;
    if (i == 0 || i == n)  % end points
        I = I + feval(f,x);
    else                     % interior points
        I = I + 2*feval(f,x);
    end
end
I = I*h/2;
return

