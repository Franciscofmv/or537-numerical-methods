%Velocity for person 1%%
m=70;t=9;c=12;
g=9.8; 
v1 = g*m - g*m*(exp(-(c*t)/m))/c
%%person2%%
syms v(t) t m c;
m = 80; c = 15;
subs(v)
v(t) = g*m - g*m*(exp(-(c*t)/m))/c-v1
double(solve(v,t))
%%what %%
d = 10*10^(-6);
rho_fluid = 1*(10^(-3))/(10^(-6));
rho_solid = 2.65*10^(-3)/(10^(-6));
mu = 0.014*(10^(-3))/(10^(-2)^2);
g = 9.8;

v = (g/18)*(rho_solid-rho_fluid)*d^2/mu