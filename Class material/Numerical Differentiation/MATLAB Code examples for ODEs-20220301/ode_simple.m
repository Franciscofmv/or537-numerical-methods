function ode_simple
f = @odefunc;
x = 0:0.2:4;
y0 = 1;
y = euler(f,x,y0);
plot(x,y,'g*');
hold on;
y = heuns(f,x,y0);
plot(x,y,'m^');
hold on;
y = mid_point(f,x,y0);
plot(x,y,'bo');
hold on;
y = runge_kutta(f,x,y0);
plot(x,y,'rs');
hold on;
[x,y] = ode45(f,x,y0);
plot(x,y,'kd');
hold on;
% solve with symbolic toolbox
syms y(x);
% define ODE
eqn = diff(y,x) == -2*x^3+12*x^2-20*x+8.5;
% define initial condition
cond = y(0) == 1;
% solve ODE
y = dsolve(eqn,cond);
x = 0:0.1:4; % define x vector
y = double(subs(y)); % conver to double
plot(x,y,'-k');
legend('euler','heuns','mid point','rk4','ode45','analytical');
hold off;
return


function [f] = odefunc(x,y)
f = -2*x^3+12*x^2-20*x+8.5;
return