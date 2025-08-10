function [] = stiff_example
t=[0:0.1:4];
y0=0;
[T,Y]=ode45(@stiff_function,t,y0);
syms t y;
y=dsolve('Dy = -1000*y+3000-2000*exp(-t)','y(0) = 0');
t=[0:0.1:4];
y=subs(y,t);
plot(t,y,'b',T,Y,'r');
legend('numerical','analytical');
return

function [dydt] = stiff_function(t,y)
dydt = -1000*y+3000-2000*exp(-t);
return