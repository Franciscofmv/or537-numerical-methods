function [] = simple_monte
% evaluate triple integral using Montecarlo integration
% function name to integrate
f='sinx';
% integration limits as vectors in x,y,z order
a=0;   % lower limit of integral
b=2;     % upper limit of integral
n=1000;  % number of montecarlo points
I3n=feval('monte',f,a,b,n);
disp(sprintf('The numerical value of single integral is %g\n',I3n));
% analytical evaluation of triple integral
clear f x;
syms x f;
f=sin(x);
I3a=subs(int(f,x,0,2));
fprintf('The analytical value of single integral is %g\n',I3a);
error=abs(I3n-I3a)/abs(I3a);
fprintf('Percentage error of numerical solution %g\n',error);
return

% montecarlo integration function
function [Im] = monte(f,a,b,n)
mu=b-a; % calculate total volume
y=rand(n,1); % random numbers within [0,1]
x=a+(b-a)*y; % convert to within [a,b]
Im=0;
for i=1:n
    Im=Im+feval(f,x(i));
end
Im=mu*Im/n;
return

function [f] = sinx(x)
f=sin(x);
return