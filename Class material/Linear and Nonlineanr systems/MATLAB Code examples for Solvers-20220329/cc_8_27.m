% chapra and canale p. 218
function cc_8_27
x0=[6;];
f=@nl_func;
x_fsolve=fsolve(f,x0)
x_nr=newton_raphson(f,x0)
ezplot(f,[0 10]);
return

function [x] = newton_raphson(f,x0)
x=x0;
i=0; maxitr=100; epsilon=inf;
while (i < maxitr && epsilon > 1e-6)
    i=i+1;
    J=numerical_jacobian(f,x);
    b=feval(f,x);
    dx=-J\b;
    x=x+dx;
    epsilon=norm(dx);
end
return

function [J] = numerical_jacobian(f,x)
m=length(x); h=1e-6;
for i = 1:m
    xp = x;
    xn = x;
    xp(i) = x(i)+h;
    xn(i) = x(i)-h;
    J(:,i) = (feval(f,xp)-feval(f,xn))/(2*h);
end
return

function [f] = nl_func(x)
f=-(5/6)*(singularity(x,0)^4-singularity(x,5)^4)+...
    (15/6)*singularity(x,8)^3+...
    75*singularity(x,7)^2+...
    (57/6)*x^3-238.25*x;
return

function [f] = singularity(x,a)
if (x > a)
    f=x-a;
else
    f=0;
end
return
