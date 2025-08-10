%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate triple integral with variable limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function variable_limit_example
f=@fxyz;
a = {0,0,0}; % xmin, ymin, zmin
b = {8,@ymax,@zmax}; % xmax, ymax, zmax
n=10; % number of segments in each dimension
I=gtrap(f,a,b,n); % trapezoidal function
fprintf('The trapezoidal value of integral is %g\n',I);
I=gsimp(f,a,b,n); % simpson function
fprintf('The simpson value of integral is %g\n',I);
I=romberg(f,a,b,n); % romberg function
fprintf('The romberg value of integral is %g\n',I);
% define limits for integral3 (note: outer is x)
xmin=0; ymin=0; zmin=0;
xmax=8; ymax=@(x) x/2; zmax=@(x,y) x-2*y;
fxyz = @(x,y,z) exp(-(x.^2+y.^2+z.^2)/4);
I=integral3(fxyz,xmin,xmax,ymin,ymax,zmin,zmax);
fprintf('The integral3 value of integral is %g\n',I);
clear;
% analytical integraton
syms x y z;
f=exp(-(x^2+y^2+z^2)/4);
I=double(int(int(int(f,x,0,z-2*y),y,0,z/2),z,0,8));
fprintf('The analytical value of double is %g\n',I);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integrand function
function [f] = fxyz(X)
x=X(1); y=X(2); z=X(3);
f = exp(-(x^2+y^2+z^2)/4);
end

% ymax function x/2
function [f] = ymax(X)
x=X(1);
f = x/2;
end

% zmax function x-2y
function [f] = zmax(X)
x=X(1); y=X(2);
f = x-2*y;
end
