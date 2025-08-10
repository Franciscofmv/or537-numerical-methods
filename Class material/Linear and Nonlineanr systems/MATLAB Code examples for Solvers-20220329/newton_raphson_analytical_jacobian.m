function [x] = newton_raphson_analytical_jacobian(f,j,x0) % function, jacobian, and initial point
x=x0;
i=0; maxitr=100; epsilon=inf;
while (i < maxitr && epsilon > 1e-6)
    i=i+1;
    J=feval(j,x);
    b=feval(f,x);
    dx=-J\b;
    x=x+dx;
    epsilon=norm(dx);
end
return