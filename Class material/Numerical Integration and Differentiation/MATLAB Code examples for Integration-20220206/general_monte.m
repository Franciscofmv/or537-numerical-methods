function [] = general_monte
% evaluate triple integral using Montecarlo integration
% function name to integrate
f='x3_2yx_monte';
% integration limits as vectors in x,y,z order
a=[-1 0 -4]';   % lower limit of integral
b=[3 6 4]';     % upper limit of integral
n=100000;  % number of montecarlo points
m=3;           % number of dimensions = 3
I3n=feval('monte',f,a,b,m,n);
disp(sprintf('The numerical value of triple integral is %g\n',I3n));

% analytical evaluation of triple integral
clear f a b x n;
syms x y z f;
f=x^3-2*y*x;
I3a=subs(int(int(int(f,x,-1,3),y,0,6),z,-4,4));
disp(sprintf('The analytical value of triple integral is %g\n',I3a));
error=abs(I3n-I3a)/abs(I3a);
disp(sprintf('Percentage error of numerical solution %g\n',error));
%end first problem
return

% montecarlo integration function
function [Im] = monte(f,a,b,m,n)
x=zeros(m,n) % initialize sampling matrix
mu=prod(b-a); % calculate total volume
y=rand(m,n); % random numbers within [0,1]
% convert to within [a,b]
for i=1:m
    x(i,:)=a(i)+(b(i)-a(i))*y(i,:);
end
Im=mu*sum(feval(f,x))/n; % evaluate at all points at once
return

% function f=x^3-2yx for montecarlo integration
% input x,y,z as row vectors in x(m,n)
% use notation x = x(1,:); y=x(2,:); z=x(3,:) etc...
function [fx]=x3_2yx_monte(x)
fx=x(1,:).^3-2*x(2,:).*x(1,:);
return
