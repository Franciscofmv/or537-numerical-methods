% general gauss quadrature for an integral of m dimensions
function [I] = mgauss(f,a,b,n)
m = length(a);
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
elseif (n == 6)
    w(1)=0.3607615730481386; zeta(1)=0.6612093864662645;
    w(2)=0.3607615730481386; zeta(2)=-0.6612093864662645;
    w(3)=0.4679139345726910; zeta(3)=-0.2386191860831969;
    w(4)=0.4679139345726910; zeta(4)=0.2386191860831969;
    w(5)=0.1713244923791704; zeta(5)=-0.9324695142031521;
    w(6)=0.1713244923791704; zeta(6)=0.9324695142031521;
else
    disp('This value of n not implemented!')
    return
end
x(1:m)=0;
% in the beginning k=1
I=feval('kgauss',f,x,zeta,w,a,b,1,n);
return

% this function is recursively called to evaluate a multiple
% integral of dimension 'm' using gauss quadrature
% it is assumed that the outermost integral is x(1)
% for example in a triple integral, x=x(1), y=x(2), x=x(3)
% int(int(int(f(x,y,z)*dx)dy)dz
function [Ik] = kgauss(f,x,zeta,w,a,b,k,n)
Ik=0; m=length(a);
if isnumeric(a{k})
    ak=a{k};
else
    ak=feval(a{k},x);
end
if isnumeric(b{k})
    bk=b{k};
else
    bk=feval(b{k},x);
end
for i=1:n
    x(k)=0.5*((bk-ak)*zeta(i)+(bk+ak));
    if (k == m)
        Ik=Ik+w(i)*feval(f,x);
    else
        Ik=Ik+w(i)*feval(@kgauss,f,x,zeta,w,a,b,k+1,n);
    end
end
Ik=Ik*(bk-ak)/2;
return