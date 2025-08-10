

f=@odefunc;
y0=0;

range=0:0.1:1.6;

[a b]=ode45(f,range,0)

function [g] = odefunc(x,y)
g =  (10^(3)*250/(120*3))*( -5*x^4 + 6*3^2*x^2 - 3^4);
end