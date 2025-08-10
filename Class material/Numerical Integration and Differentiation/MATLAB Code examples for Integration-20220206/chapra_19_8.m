% C & C 21.23 (p. 631)
function cc_21_23
% define tabular data
T=[0 1 2 3.25 4.5 6 7 8 8.5 9 10];
V=[0 5 6 5.5 7 8.5 8 6 7 7 5];
%% method 1: simple trapezoidal tabular data interpolation
d_trap_data = trap_data(T,V)
%% method 2: fit a cubic polynomial to data and integrate using trap
p=polyfit(T,V,3); % polynomial coefficients are in 'p'
% define uniformly spaced time data for evaluation of velocity
t=0:0.1:10;
% evaluate fitted velocity values V from T
v=polyval(p,t);
% plot data and fitted values
plot(T,V,'ro',t,v,'-g');
% define integrand function with polynomial fit
f1=@(t) vel_poly(t,p);
% define integrand function with linear interpolation
f2=@(t) vel_interp(t,T,V);
% integrate to find distance 'd' from trap function
a=0; b=10; n=100;
d_cubic_fit = trap(f1,a,b,n)
d_linear_interp = trap(f2,a,b,n)
return

% integrand function with polynomial fit
function [v] = vel_poly(t,p)
v = polyval(p,t);
return

% integrand function with linear interpolation
function [v] = vel_interp(t,T,V)
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

function [I] = trap_data(t,v)
n=length(t); I=0;
for i=1:n-1
    h=t(i+1)-t(i);
    I=I+h*(v(i+1)+v(i))/2;
end

