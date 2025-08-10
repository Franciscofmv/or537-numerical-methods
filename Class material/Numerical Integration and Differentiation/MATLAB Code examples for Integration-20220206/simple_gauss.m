function simple_gauss
f=@sinsqx;
a=0;   % lower limit of integral
b=pi/2;     % upper limit of integral
n=2;  % number of gauss points
I3n=@gauss(f,a,b,n)
syms x;
I3a=double(int(sin(x)^2,x,0,pi/2))
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
alpha=(b-a)/2;
beta=(b+a)/2;
for i=1:n    
    x=alpha*zeta(i)+beta; % transformation
    I=I+w(i)*feval(f,x); % evaluate our function at x
end
I=I*alpha;
return

function [f] = sinsqx(x) % integrand function
f=sin(x)^2;
return