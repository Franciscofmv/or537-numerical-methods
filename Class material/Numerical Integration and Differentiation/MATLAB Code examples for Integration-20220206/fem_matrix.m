% evaluate 4x4 fem element matrix
% for quadrilateral elements and poisson equation
% using Gauss quadrature
function fem_matrix
f='kme';
% need to input the limits as characters as we want to
% be able to handle function names as limits as well
a=char('-1','-1');   % lower limit of integral
b=char('1','1');     % upper limit of integral
n=2;  % number of gauss points
m=2;
KM = feval('m_integral_gauss',f,a,b,m,n)
return

% function for problem 3_3
function [fx] = kme(x)
eta=x(1);
eps=x(2);
f='gn_fnc';
GN=feval(f,eps,eta);
XY= [ 3 2; 7 1; 7 5; 3 4];
J = GN'*XY;
KE= [1 0; 0 2];
JINV = inv(J);
detj = det(J); %calculate determinant
fx = GN*JINV'*KE*JINV*GN'*detj; % matrix, integrand function
return

function [GN] = gn_fnc(eps,eta)
epsi=[-1 1 1 -1]; 
etai=[-1 -1 1 1];
gn_eps=epsi.*(1+etai*eta)/4;
gn_eta=etai.*(1+epsi*eps)/4;
GN = [gn_eps' gn_eta'];
return

% general gauss quadrature for an integral of m dimensions
function [Im] = m_integral_gauss(f,a,b,m,n)
x(1:m)=0;
[zeta, w] = gauss_table(n);
% in the beginning k=m
Im=feval('kth_integral_gauss',f,x,zeta,w,a,b,m,n);
return

function [zeta, w] = gauss_table(n)
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
    disp('This value of n not implemented!')
    return
end
return

% this function is recursively called to evaluate a multiple
% integral of dimension 'm' using gauss quadrature
% it is assumed that the outermost integral is x(1)
% for example in a triple integral, x=x(1), y=x(2), x=x(3)
% int(int(int(f(x,y,z)*dz)dy)dx
function [Ik] = kth_integral_gauss(f,x,zeta,w,a,b,k,n)
Ik=0;
ak=str2num(a(k,:));
bk=str2num(b(k,:));
if isempty(ak)
    fnc=deblank(a(k,:));
    ak=feval(fnc,x);
end
if isempty(bk)
    fnc=deblank(b(k,:));
    bk=feval(fnc,x);
end
for i=1:n
    x(k)=0.5*((bk-ak)*zeta(i)+(bk+ak));
    if (k == 1)
         Ik=Ik+w(i)*feval(f,x);
      else
         Ik=Ik+w(i)*feval('kth_integral_gauss',f,x,zeta,w,a,b,k-1,n);
    end
end
Ik=Ik*(bk-ak)/2;
return