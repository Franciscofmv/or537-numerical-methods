function exercise5
x = [0 2 4 5 6 7 10]*10^(-2);
r = [1.35 1.34 1.6 1.58 1.42 2]*10^(-3);
dx = [2 2 1 1 1 3];
v=@pressure
v(dx,r)

end

function [p] = pressure(dx,r)
mu=0.005;%N · s/m^2 
rho = 1*10^3;%  kg/m3
Q = 10*10^(-6); % m^3/s 
dp = dx.*(-8*mu*Q)./(pi*r.^4);
sum(dp)
end