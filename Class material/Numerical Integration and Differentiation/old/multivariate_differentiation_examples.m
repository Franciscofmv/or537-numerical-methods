function [] = multivariate_differentiation_examples
gradient_example;
jacobian_example;
return

function [] = gradient_example
f='test_function1';
x = [2 1 -2]; h = 1e-3;
g=my_gradient(f,x,h);
fprintf('numerical gradient  = [ %g %g %g ]\n',g(1),g(2), g(3));
clear
syms x y z;
f = '10*log(((x-5)^2+(y-5)^2+(z-5)^2)/((x+5)^2+(y+5)^2+(z+5)^2))';
f1 = diff(f,x);
f2 = diff(f,y);
f3 = diff(f,z);
x = 2; y = 1; z = -2;
g(1) = subs(f1);
g(2) = subs(f2);
g(3) = subs(f3);
fprintf('analytical gradient  = [ %g %g %g ]\n',g(1),g(2), g(3));
return

function [] = jacobian_example
f='test_function2';
x = [2 1 -2]; h = 1e-3;
J=my_jacobian(f,x,h);
fprintf('numerical jacobian\n');
disp(J);
clear
syms x y z;
f = [x*y*z; y; x+z; x^2];
v = [x, y, z];
J = jacobian(f,v);
x = 2; y = 1; z = -2;
J =subs(J);
fprintf('analytical jacobian\n');
disp(J);
return

function [g] = my_gradient(f,x,h)
m=length(x);
for i = 1:m
    x1 = x;
    x2 = x;
    x1(i) = x(i)+h;
    x2(i) = x(i)-h;
    g(i) = (feval(f,x1)-feval(f,x2))/(2*h);
end
return

function [J] = my_jacobian(f,x,h)
m=length(x);
for i = 1:m
    x1 = x;
    x2 = x;
    x1(i) = x(i)+h;
    x2(i) = x(i)-h;
    J(:,i) = (feval(f,x1)-feval(f,x2))/(2*h);
end
return

function [f] = test_function1(x)
% x = x(1); y = x(2); z = x(3);
f = 10*log(((x(1)-5)^2+(x(2)-5)^2+(x(3)-5)^2)/((x(1)+5)^2+(x(2)+5)^2+(x(3)+5)^2));
return

function [f] = test_function2(x)
% make sure f is a column vector
f = zeros(3,1);
f(1) = x(1)*x(2)*x(3);
f(2) = x(2);
f(3) = x(1)+x(3);
f(4) = x(1)^2;
return