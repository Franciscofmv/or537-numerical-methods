% chemical reactor example
function reactor_ode_system_example
f=@reactor_func;
t = [0:0.2:10]; c0=[0; 0; 0; 0; 0];
[c] = runge_kutta(f,t,c0);
plot(t,c(1,:),'-go',t,c(2,:),'-k+');
hold on;
[t, c] = ode45(f,[0 10],c0);
plot(t,c(:,1),'-bs',t,c(:,2),'-r^');
legend('c1-rk','c2-rk','c1-matlab','c2-matlab');
hold off;
return

function [dcdt] = reactor_func(t,c)
dcdt=zeros(5,1); k1=0.1;
dcdt(1)=-0.12*c(1)+0.02*c(3)+1-k1*c(1);
dcdt(2)=0.15*c(1)-0.15*c(2);
dcdt(3)=0.02*c(2)-0.225*c(3)+4;
dcdt(4)=0.1*c(3)-0.137*c(4)+0.025*c(3);
dcdt(5)=0.03*c(1)+0.01*c(2)-0.04*c(5);
return

function [y] = forward_euler(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    y(:,i+1)=y(:,i)+feval(f,x(i),y(:,i))*dx;
end
return

function [y] = heuns(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    ytemp=y(:,i)+feval(f,x(i),y(:,i))*dx;
    y(:,i+1)=y(:,i)+(feval(f,x(i),y(:,i))+feval(f,x(i+1),ytemp))*dx/2;
end
return

function [y] = mid_point(f,x,y0)
n=length(x); m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    xtemp=(x(i+1)+x(i))/2;
    ytemp=y(:,i)+feval(f,x(i),y(:,i))*dx/2;
    y(:,i+1)=y(:,i)+feval(f,xtemp,ytemp)*dx;
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


