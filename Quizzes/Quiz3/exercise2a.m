function exercise2a
syms theta(x) y(x) L E ag w0;
Y = diff(y(x),x);
eqs = diff(y(x),x) == (1/(120*E*ag*L))* w0*( -5*x^4 + 6*L^2*x^2 - L^4);
sol = dsolve(eqs, [y(0)==0,y(L)==0])
L = 3.0; % beam distance (m)
E=10^(9)*200;  % GPa
ag = 0.0003; % m^4
w0 = 10^(3)*250; %kN/m
x=1.5;
format long
double(subs(sol))
-0.07910156250000001

end