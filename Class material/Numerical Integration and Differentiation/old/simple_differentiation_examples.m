function [] = simple_differentiation_examples
first_derivative_example;
second_derivative_example;
third_derivative_example;
tabular_data_derivative_example;
return

% first derivative example
function [] = first_derivative_example
f = 'xex';
x=1; h=1e-6;
method = 'forward';
diffx=diff1_fd(f,x,method,h);
fprintf('first derivative using simple %s differencing = %g\n',method, diffx);
method = 'central'; n = 5; order = 1;
diffx=diffn_fd(f,x,method,order,h,n);
fprintf('first derivative using %s differencing with %g points = %g\n',method, n, diffx);
clear;
syms x;
diffx = subs(diff('x*exp(x)'),1);
fprintf('first derivative analytical solution = %g\n\n',diffx);
return

function [] = second_derivative_example
f = 'xex';
x=1; h=1e-6;
method = 'forward';
diffx=diff2_fd(f,x,method,h);
fprintf('second derivative using simple %s differencing = %g\n',method, diffx);
method = 'central'; n = 5; order = 2;
diffx=diffn_fd(f,x,method,order,h,n);
fprintf('second derivative using %s differencing with %g points = %g\n',method, n, diffx);
diffx=diff_fd_recursive(f,x,method,h,0,2);
fprintf('second derivative using recursive %s differencing = %g\n',method, diffx);
clear;
syms x;
diffx = subs(diff(diff('x*exp(x)')),1);
fprintf('secod derivative analytical solution = %g\n\n',diffx);
return


function [] = third_derivative_example
f = 'xex';
x=1; h=1e-3; method='central';
diffx=diff_fd_recursive(f,x,method,h,0,3);
fprintf('third derivative using recursive %s differencing = %g\n',method,diffx);
clear;
syms x;
diffx = subs(diff(diff(diff('x*exp(x)'))),1);
fprintf('third derivative analytical solution = %g\n\n',diffx);
return

function [] = tabular_data_derivative_example;
n = 100;
x = rand(1,n);
y = x.*sin(x.^2);
diffx = tabular_derivative(x,y,0.23);
fprintf('tabular data derivative using central differencing = %g\n',diffx);
clear
syms x;
diffx = subs(diff('x*sin(x^2)'),0.23);
fprintf('tabular data derivative analytical solution = %g\n\n',diffx);
return

% first derivative function
function [diff1] = diff1_fd(f,x,method,h)
switch (method),
    case ('forward'),
        diff1=(feval(f,x+h)-feval(f,x))/h;
    case ('backward'),
        diff1=(feval(f,x)-feval(f,x-h))/h;
    otherwise,
        diff1=(feval(f,x+h)-feval(f,x-h))/(2*h);
end
return

% second derivative function
function [diff2] = diff2_fd(f,x,method,h)
switch (method)
    case ('forward'),
        diff2=(feval(f,x+2*h)-2*feval(f,x+h)+feval(f,x))/h^2;
    case ('backward'),
        diff2=(feval(f,x-2*h)-2*feval(f,x-h)+feval(f,x))/h^2;
    otherwise, % use central differencing
        diff2=(feval(f,x-h)-2*feval(f,x)+feval(f,x-h))/h^2;
end
return

% recursive derivative function
function [diff] = diff_fd_recursive(f,x,method,h,k,m)
fr='diff_fd_recursive';
k=k+1;
if (k == m)
    switch (method),
        case ('forward'),
            diff=(feval(f,x+h)-feval(f,x))/h;
        case ('backward'),
            diff=(feval(f,x)-feval(f,x-h))/h;
        otherwise,
            diff=(feval(f,x+h)-feval(f,x-h))/(2*h);
    end
else
    switch (method),
        case ('forward'),
            diff=(feval(fr,f,x+h,method,h,k,m)-feval(fr,f,x,method,h,k,m))/h;
        case ('backward'),
            diff=(feval(fr,f,x,method,h,k,m)-feval(fr,f,x-h,method,h,k,m))/h;
        otherwise,
            diff=(feval(fr,f,x+h,method,h,k,m)-feval(fr,f,x-h,method,h,k,m))/(2*h);
    end   
end
return


% general derivative function
function [diff] = diffn_fd(f,x,method,order,h,n)
switch (method)
    case ('forward'),
        if (order == 1)
            if (n <= 2)
                diff=(feval(f,x+h)-feval(f,x))/h;
            elseif (n <= 3)
                diff=(-1/2*feval(f,x+2*h)+2*feval(f,x+h)-3/2*feval(f,x))/h;
            else
                diff=(1/3*feval(f,x+3*h)-3/2*feval(f,x+2*h)+3*feval(f,x+h)...
                    -11/6*feval(x,h))/h;
            end
        elseif (order == 2)
            if (n <= 3)
                diff=(feval(f,x+2*h)-2*feval(f,x+h)+feval(f,x))/h^2;
            else
                diff=(-feval(f,x+3*h)+4*feval(f,x+2*h)-5*feval(f,x+h)...
                    +2*feval(f,x))/h^2;
            end
        end
    case ('backward'),
        if (order == 1)
            if (n <= 2)
                diff=(feval(f,x)-feval(f,x-h))/h;
            elseif (n <= 3)
                diff=(1/2*feval(f,x-2*h)-2*feval(f,x-h)+3/2*feval(f,x))/h;
            else
                diff=(-1/3*feval(f,x-3*h)+3/2*feval(f,x-2*h)-3*feval(f,x-h)...
                    +11/6*feval(x,h))/h;
            end
        elseif (order == 2)
            if (n <= 3)
                diff=(-feval(f,x-2*h)+2*feval(f,x-h)-feval(f,x))/h^2;
            else
                diff=(feval(f,x-3*h)-4*feval(f,x-2*h)+5*feval(f,x-h)...
                    -2*feval(f,x))/h^2;
            end
        end
    otherwise, % use central differencing
        if (order == 1)
            if (n <= 3)
                diff=(feval(f,x+h)-feval(f,x-h))/(2*h);
            else
                diff=(1/12*feval(f,x-2*h)-2/3*feval(f,x-h)...
                    +2/3*feval(f,x+h)-1/12*feval(f,x+2*h))/h;
            end
        elseif (order == 2)
            if (n <= 3)
                diff=(feval(f,x-h)-2*feval(f,x)+feval(f,x+h))/h^2;
            else
                diff=(-1/12*feval(f,x-2*h)+4/3*feval(f,x-h)-5/2*feval(f,x)...
                    +4/3*feval(f,x+h)-1/12*feval(f,x+2*h))/h^2;
            end
        end
end
return

function [f] = xex(x)
f = x*exp(x);
return

function [f] = tabular_derivative(x,y,x0)
[x, index]=sort(x);
y = y(index);
n=length(x); i=0;
for j=2:n-1
    if (x0 >= x(j) & x0 <= x(j+1))
        i=j;
    end
end
if (i == 0)
    if (x0 < x(2))
        i = 2;
    else
        i = n-1;
    end
end
% lagrange interpolation
% f = y(i-1)*(2*x0-x(i)-x(i+1))/(x(i-1)-x(i))*(x(i-1)-x(i+1))...
%     +y(i)*(2*x0-x(i-1)-x(i+1))/(x(i)-x(i-1))*(x(i)-x(i+1))...
%     +y(i+1)*(2*x0-x(i-1)-x(i))/(x(i+1)-x(i-1))*(x(i+1)-x(i));
y0=interp1(x,y,x0);
f = (y0-y(i-1))/(x0-x(i-1)); % backward difference
% f = (y(i)-y0)/(x(i)-x0); % forward difference
% f = (y(i+1)-y(i-1))/(x(i+1)-x(i-1)); % central difference
return