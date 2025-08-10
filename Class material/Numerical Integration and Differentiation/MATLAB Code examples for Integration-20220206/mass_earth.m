% simple multi-dimensional integration
function mass_earth
% define limits
R=6380;
a = [0 0 0]; b = [R pi 2*pi];
n=10;
% define integrand function name for using our functions
f=@mass;
I_trap = mtrap(f,a,b,n)
I_simp = msimp(f,a,b,n)
% define integrand function name for using integral
f=@massm;
I_matlab = integral3(f,a(1),b(1),a(2),b(2),a(3),b(3))
return

% integrand function
function [f] = mass(x)
r=x(1); theta=x(2); beta=x(3);
rr=[0 1100 1500 2450 3400 3630 4500 5380 6060 6280 6380];
rhor=[13 12.4 12 11.2 9.7 5.7 5.2 4.7 3.6 3.4 3];
rho=interp1(rr,rhor,r);
f = rho*sin(theta)*r^2;
return

% integrand function for integral3
function [f] = massm(r,theta,beta)
rr=[0 1100 1500 2450 3400 3630 4500 5380 6060 6280 6380];
rhor=[13 12.4 12 11.2 9.7 5.7 5.2 4.7 3.6 3.4 3];
rho=interp1(rr,rhor,r);
f = rho.*sin(theta).*r.^2;
return