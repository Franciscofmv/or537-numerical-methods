% line integration example
function line_integral_example
n=10;
a=0;
b=3;
f=@xex;
I_trap = trapezoidal(f,a,b,n);
fprintf('trapezoidal value of integral %g\n',I_trap);
I_simp = simpson(f,a,b,n);
fprintf('simpson value of integral %g\n',I_simp);
I_quad = quad(f,a,b);
fprintf('matlab quad function value of integral %g\n',I_quad);
clear;
syms x;
I_analytical=subs(int('x*exp(x)',x,0,3));
fprintf('analytical value of integral %g\n',I_analytical);
return

function [I_trap]= trapezoidal(f,a,b,n)
I_trap=0;
x=0;
h=(b-a)/n;
for i=0:n
    x=a+i*h;
    if (i == 0 || i == n) % end points
        I_trap=I_trap+feval(f,x);
    else % interior points
        I_trap=I_trap+2*feval(f,x);
    end
end
I_trap=I_trap*h/2;
return

function [I_simp] = simpson(f,a,b,n)
I_simp = 0;
x=0;
h=(b-a)/n;
for i=0:n
    x=a+i*h;
    if (i == 0 || i == n) % end points
        I_simp=I_simp+feval(f,x);
    elseif (mod(i,2) == 0) % even
        I_simp=I_simp+2*feval(f,x);
    else % odd
        I_simp=I_simp+4*feval(f,x);
    end
end
I_simp=I_simp*h/3;
return

function [f] = xex(x)
f=x.*exp(x);
return
