function nonlinear_system_example
x0=[0; 0; 0;]; % initial guess
f=@nlfunc; % nonlinear function
% MATLAB's fsolve function
x_fsolve=fsolve(f,x0) 
% default newton raphson with numerical jacobian
x_nr=newton_raphson(f,x0) 
% newton raphson with analytical jacbian
jf=@jfunc; % analytical jacobian function
x_nra=newton_raphson_analytical_jacobian(f,jf,x0)
return

function [f] = nlfunc(x) % our system of equations
f=zeros(3,1);
f(1)=exp(2*x(1))-x(2)-4; % residuals?
f(2)=x(2)-x(3)^2-1;
f(3)=x(3)-sin(x(1));
return

function [J] = jfunc(x)
J=zeros(3);
J(1,1)=2*exp(2*x(1));
J(1,2)=-1;
J(2,2)=1;
J(2,3)=-2*x(3);
J(3,1)=-cos(x(1));
J(3,3)=1;
return
