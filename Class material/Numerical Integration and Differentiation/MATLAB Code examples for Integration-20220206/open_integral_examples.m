% open integration examples
function open_integral_examples
example1;
example2;
example3;
return

function example1
n=10;
a=0;
b=1;
I = mid_point('inv_sqrt',a,b,n);
fprintf('mid point rule value of integral = %g\n',I);
I = gauss('inv_sqrt',a,b,5);
fprintf('gauss quadrature value of integral = %g\n',I);
% I = quad(@inv_sqrt,a,b);
% fprintf('matlab quad function value of integral = %g\n',I);
clear;
syms x;
I_analytical=subs(int('1/sqrt(x)',x,0,1));
fprintf('analytical value of integral = %g\n',I_analytical);
return

function example2
n=10;
a=0;
b=1;
I1 = simpson('fnc1',a,b,n);
I2 = simpson('fnc2',a,b,n);
I =I1+I2;
fprintf('simpson rule value of integral = %g\n',I);
I1 = mid_point('fnc1',a,b,n);
I2 = mid_point('fnc2',a,b,n);
I =I1+I2;
fprintf('mid point rule value of integral = %g\n',I);
clear;
syms x;
I_analytical=subs(int('4/(pi*(4+x^2))',x,0,inf));
fprintf('analytical value of integral = %g\n',I_analytical);
return

function example3
n=10;
f1=@fnc3;
I1 = mid_point('fnc3',-1,0,n);
I2 = simpson('fnc4',-1,1,n);
I =I1+I2;
fprintf('mid point + simpson value of integral = %g\n',I);
clear;
syms x;
I_analytical=subs(int('exp(-x^2/2)/(sqrt(2*pi))',x,-inf,1));
fprintf('analytical value of integral = %g\n',I_analytical);
return

function [I]= mid_point(f,a,b,n)
I=0;
x=0;
h=(b-a)/n;
for i=0:n-1
    x=a+(i+0.5)*h;
    I=I+feval(f,x);
end
I=I*h;
return

function [I] = gauss(f,a,b,n)
if (n == 1)
    zeta(1)=0;
    w(1)=2;
elseif (n == 2)
    zeta(1)=1/sqrt(3);
    zeta(2)=-zeta(1);
    w(1)=1;
    w(2)=1;
elseif (n == 3)
    zeta(1)=0;
    zeta(2)=sqrt(3/5);
    zeta(3)=-zeta(2);
    w(1)=8/9;
    w(2)=5/9;
    w(3)=w(2);
elseif (n == 4)
    zeta(1)=sqrt((3-2*sqrt(6/5))/7);
    zeta(2)=-zeta(1);
    zeta(3)=sqrt((3+2*sqrt(6/5))/7);
    zeta(4)=-zeta(3);
    w(1)=0.5+1/(6*sqrt(6/5));
    w(2)=w(1);
    w(3)=0.5-1/(6*sqrt(6/5));
    w(4)=w(3);
elseif (n == 5)
    zeta(1)=0;
    zeta(2)=sqrt(5-4*sqrt(5/14))/3;
    zeta(3)=-zeta(2);
    zeta(4)=sqrt(5+4*sqrt(5/14))/3;
    zeta(5)=-zeta(4);
    w(1)=128/225;
    w(2)=161/450+13/(180*sqrt(5/14));
    w(3)=w(2);
    w(4)=161/450-13/(180*sqrt(5/14));
    w(5)=w(4);
else
    disp('This value of n not implemented!');
    return
end
I=0;
for i=1:n
    x=0.5*((b-a)*zeta(i)+(b+a));
    I=I+w(i)*feval(f,x);
end
I=I*(b-a)/2;
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

function [f] = inv_sqrt(x)
f=1./sqrt(x);
return

function [f] = fnc1(x)
f=4./(pi*(4+x.^2));
return

function [f] = fnc2(z)
f=4./(pi*(4*z.^2+1));
return

function [f] = fnc3(z)
f=exp(-z^-2/2)/(sqrt(2*pi)*z^2);
return

function [f] = fnc4(x)
f=exp(-x^2/2)/sqrt(2*pi);
return

