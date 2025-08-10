% simple ode example
function pendulum_ode_example
f=@pendulum;
t = [0:0.2:10]; y0=[pi/4; 0];
[y] = runge_kutta(f,t,y0);
plot(t,y(1,:),'-go',t,y(2,:),'-k+');
hold on;
[t, y] = ode45(f,[0 10],y0');
plot(t,y(:,1),'-bs',t,y(:,2),'-r^');
legend('y1-rk','y2-rk','y1-matlab','y2-matlab');
hold off;
return

function [dydt] = pendulum(t,y)
dydt=zeros(2,1);
l=2; g=9.8;
dydt(1)=y(2);
dydt(2)=-g*sin(y(1))/l;
return


function [y] = runge_kutta(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    k1=feval(f,x(i),y(:,i));
    k2=feval(f,x(i)+dx/2,y(:,i)+k1*dx/2);
    k3=feval(f,x(i)+dx/2,y(:,i)+k2*dx/2);
    k4=feval(f,x(i)+dx,y(:,i)+k3*dx);
    y(:,i+1)=y(:,i)+(k1+2*k2+2*k3+k4)*dx/6;
end
return




