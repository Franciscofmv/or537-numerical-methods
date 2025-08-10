% CC 28.47
function forced_damping
f=@damping;
t = [0:0.001:4]; y0=[1; 0];
[y] = runge_kutta(f,t,y0);
plot(t,y(1,:),'-r',t,y(2,:),'-b');
hold on;
[t, y] = ode45(f,0:0.5:4,y0');
plot(t,y(:,1),'rs',t,y(:,2),'bo');
legend('y1-rk','y2-rk','y1-matlab','y2-matlab');
hold off;
return

function [dydt] = damping(t,y)
dydt=zeros(2,1);
m=2; a=5; k=6; F0=2.5; w=0.5;
y1=y(1); y2=y(2);
dydt(1)=y2;
dydt(2)=(-a*abs(y2)*y2-k*y1+F0*sin(w*t))/m;
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