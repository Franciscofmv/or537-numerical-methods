function general_gauss_example
f=@x4my2m2xyz;
% initialize upper and lower limits
a={-1,0,-4};   % lower limit of integral
b={3,6,4};     % upper limit of integral
n=3;  % number of gauss points
I3n=feval(@ggauss,f,a,b,n)
% analytical evaluation of triple integral
clear f a b x;
syms x y z f;
f=x^4-y^2+2*x*y*z;
I3a=double(int(int(int(f,x,-1,3),y,0,6),z,-4,4))
return

function [f]=x4my2m2xyz(x)
% function f=x^3-2yx
% use notation x = x(1); y=x(2); z=x(3) etc...
f=x(1)^4-x(2)^2+2*x(1)*x(2)*x(3);
return