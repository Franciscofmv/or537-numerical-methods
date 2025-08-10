function [] = nonlinear_solver_example
x0=[0; 0; 0;];
x_fsolve=fsolve(@nl_fun,x0,optimset('fsolve'))
x_nr=newton_raphson(x0)
return

function [f] = nl_fun(x)
f=zeros(3,1);
f(1)=exp(2*x(1))-x(2)-4;
f(2)=x(2)-x(3)^2-1;
f(3)=x(3)-sin(x(1));
return

function [x] = newton_raphson(x0)
x=x0;
i=0; maxitr=100; epsilon=inf;
while (i < maxitr & epsilon > 1e-6)
    i=i+1;
    J=jacobian(x);
    b=nl_fun(x);
    dx=-J\b;
    x=x+dx;
    epsilon=norm(dx);
end
return



function [J] = jacobian(x)
J=zeros(3);
J(1,1)=2*exp(2*x(1));
J(1,2)=-1;
J(2,2)=1;
J(2,3)=-2*x(3);
J(3,1)=-cos(x(1));
J(3,3)=1;
return
