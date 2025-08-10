% calculating J numerically
function [x] = newton_raphson(f,x0) %f and initial guess
x=x0;
i=0; maxitr=10; epsilon=inf;
while (i < maxitr && epsilon > 1e-6) 
    i=i+1;
    J=jacobian(f,x); % uses finite difference to calculate jacobian
    b=f(x);
    dx=-J\b; % calculating delta x
    x=x+dx; % x_k+1
    epsilon=norm(dx);
end
return


function [J] = jacobian(f,x)
n=length(x);
J=zeros(n);
h=1e-3;
for i=1:n
    x1=x; x2=x;
    x1(i)=x(i)+h;
    x2(i)=x(i)-h;
    J(:,i)=(feval(f,x1)-feval(f,x2))/(2*h);
end
return

