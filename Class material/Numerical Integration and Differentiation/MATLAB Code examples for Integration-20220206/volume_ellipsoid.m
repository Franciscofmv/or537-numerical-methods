% volume of an ellipsoid
function volume_ellipsoid
f=@ellipse;
A=4; B=2; C=1;
a = {-A,@(x)ymin(x,A,B,C),@(x)zmin(x,A,B,C)}; % xmin, ymin, zmin, lower limits
b = {A,@(x)ymax(x,A,B,C),@(x)zmax(x,A,B,C)}; % xmax, ymax, zmax. Upper limits
n=20;
I_simp = real(gsimp(f,a,b,n))% we are using simpson
% I_integral3 = integral3(f,-4,4,@ymin,@ymax,@zmin,@zmax) % can use this
% built in function

% simbolic integration:
syms x y z a b c; 
S=((int(int(int(1,x,...
    -a*sqrt(1-(z/c)^2-(y/b)^2),a*sqrt(1-(z/c)^2-(y/b)^2)),...
    y,-b*sqrt(1-(z/c)^2),b*sqrt(1-(z/c)^2)),...
    z,-c,c)));
a=4; b=2; c=1; % values for analytical
I_sym = double(subs(S))
I_anal = 4*pi*A*B*C/3
return

function [f] = ellipse(x)%integrand function
f = 1;
return

function [f] = ymin(x,A,B,C) 
X=x(1);
f = -B*sqrt(1-(X/A).^2);
return

function [f] = ymax(x,A,B,C)
X=x(1);
f = A*sqrt(1-(X/A).^2);
return

function [f] = zmin(x,A,B,C)
X=x(1); Y=x(2); 
f = -C*sqrt(1-(Y/B).^2-(X/A).^2);
return

function [f] = zmax(x,A,B,C)
X=x(1); Y=x(2);
f = C*sqrt(1-(Y/B).^2-(X/A).^2);
return

