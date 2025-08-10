function [] = gradient_vector
f=@test_function;
x = [1 2 3]; h=1e-3;
g = my_gradient(f,x,h)
syms x y z;
f='sqrt(x^2+y^2+z^2)';
ga1=subs(diff(f,x),{x,y,z},{1,2,3});
ga2=subs(diff(f,y),{x,y,z},{1,2,3});
ga3=subs(diff(f,z),{x,y,z},{1,2,3});
ga = [ga1 ga2 ga3]
return

function [f] = test_function(x)
f = sqrt(x(1)^2+x(2)^2+x(3)^2);
return

function [g] = my_gradient(f,x,h)
n=length(x);
for i=1:n
    x1=x;
    x2=x;
    x1(i)=x(i)+h;
    x2(i)=x(i)-h;
    g(i)=(feval(f,x1)-feval(f,x2))/(2*h);
end
return