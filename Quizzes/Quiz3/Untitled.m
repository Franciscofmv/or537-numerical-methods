L = 3.0; % beam distance (m)
E=10^(9)*200;  % GPa
ag = 0.0003; % m^4
w0 = 10^(3)*250; %kN/m
x=1.5
format long
w0/(120*E*ag*L)*(-x^5 + 2*L^2*x^3-L^4*x)