% cc 2_7
function cc_2_7
a = 2;
a=input('Enter the number\n');
x = square_root(a);
fprintf('square root is %g\n',x);
end

% function to calculate square root of a positive number
% using an age old iterative method
% C & C Prob 2.7 (p.25)
% input 'a', output 'x'
function [x] = square_root(a)
% return zero for non-positive values of a
if (a <= 0)
    x = 0;
    return
end
% set initial values
tol=1e-6; % tolerance
e = inf; % set initial value of e to a large number  
x = a/2; % initialize x
while (e > tol)
    y = (x+a/x)/2;
    e = abs(y-x)/y;
    x = y;
end
end