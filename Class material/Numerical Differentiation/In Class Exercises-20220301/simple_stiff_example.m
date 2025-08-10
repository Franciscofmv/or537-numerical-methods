function [] = simple_stiff_example
t=[0:0.1:1];
y0=0;
[y] = simple_backward_euler(t,y0);
plot(t,y,'-go');
hold on;
f=@stiff_function;
[t,y] = ode45(f,t,y0);
plot(t,y,'-b+');
clear;
syms t y;
y=dsolve('Dy = -1000*y+3000-2000*exp(-t)','y(0)=0');
t=[0:0.1:1];
y=subs(y,t);
plot(t,y,'-rs');
legend('backward-euler','runge-kutta','analytical');
hold off;
return

function [y] = simple_backward_euler(t,y0)
n=length(t);
y(1)=y0;
for i=1:n-1
    dt=t(i+1)-t(i);
    y(i+1)=(y(i)+3000*dt-2000*dt*exp(-t(i+1)))/(1+1000*dt);
end
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

function [dydt] = stiff_function(t,y)
dydt(1) = -1000*y(1)+3000-2000*exp(-t);
return