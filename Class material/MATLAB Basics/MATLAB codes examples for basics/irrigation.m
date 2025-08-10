% Palm Ex on p.137
% Optimization of an irrigation channel
% See explanation in class notes
function irrigation
A = 100;
% use anonymous function to pass in A in addition to default x
f = @(x) channel(x,A); % handle to irrigation function
x0=[20 1]; % initial guess for fminsearch
x=fminsearch(f,x0); % pass in f and x0; return solution in x
% output values
d = x(1)
theta = x(2)
b = A/d - d/tan(theta)
end

% this function will be used by fminsearch
% in an iterative fashion
function [L] = channel(x,A)
d=x(1);
theta=x(2);
b = A/d - d/tan(theta); % calculate b from A,d,theta
L = b +2*d/sin(theta);
end