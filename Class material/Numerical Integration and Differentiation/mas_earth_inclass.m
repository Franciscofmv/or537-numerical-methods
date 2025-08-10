function mas_earth_inclass
R = 6380; % radius in km 
a = [0 0 0]; % lower limits
b = [pi*2 pi R]; % upper limits
f = @density% function that gives density at any given radius and agnle
n = 10;
f = @density;
mass = msimp(f,a,b,n)
end
% need to interpoalte

function [d] = density(x)
beta = x(1); alpha = x(2); r = x(3);
rr = [0 1100 1500 2450 3400 3630 4500 5380 6060 6280 6380];
rhor = [13 12.4 12 11.2 9.7 5.7 5.2 4.7 3.6 3.4 3];
rho = interp1(rr,rhor,r);
d = rho*r^2*sin(alpha);

end

