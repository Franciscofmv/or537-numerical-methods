% simple ode example
function pendulum_ode_example
f=@pendulum;
t = (0:0.2:10); 
y0=[pi/4; 0]; % initial conditions vector
[y] = runge_kutta(f,t,y0);
plot(t,y(1,:),'-go',t,y(2,:),'-k+'); % solution runge kutta
hold on;
[t, y] = ode45(f,[0 10],y0'); 
plot(t,y(:,1),'-bs',t,y(:,2),'-r^');% solution ode45
legend('y1-rk4','y2-rk4','y1-ode45','y2-ode45');
hold off;
return

function [dydt] = pendulum(t,y)
dydt=zeros(2,1);
l=2; g=9.8;
dydt(1)=y(2);
dydt(2)=-g*sin(y(1))/l;
return

