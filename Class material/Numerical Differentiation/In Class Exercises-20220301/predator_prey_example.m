% chapra and canale example ex 28_2 p. 815
function predator_prey_example
t=0:0.1:30;
P0=[2; 1];
f=@predator_prey;
P=RK4(f,t,P0);
figure(1);
plot(t,P(1,:),'-g',t,P(2,:),'-k');
legend('Prey','Predator');
[T,P]=ode45(f,[0 30], P0);
figure(2);
plot(T,P(:,1),'-g',T,P(:,2),'-k');
return

function [dpdt] = predator_prey(t,p)
dpdt=zeros(2,1); 
x=p(1); y=p(2);
a=1.2; b=0.6; c=0.8; d=0.3;
dpdt(1)=a*x-b*x*y;
dpdt(2)=-c*y+d*x*y;
return

function [y] = Euler(f,x,y0)
n=length(x);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    y(:,i+1)=yi+feval(f,xi,yi)*dx;
end
return

function [y] = Heuns(f,x,y0)
n=length(x); 
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    f1=feval(f,xi,yi);
    ytemp=yi+f1*dx;
    f2=feval(f,x(i+1),ytemp);
    y(:,i+1)=yi+(f1+f2)*dx/2;
end
return


function [y] = RK4(f,x,y0)
n=length(x);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    k1=feval(f,xi,yi);
    k2=feval(f,xi+0.5*dx,yi+0.5*k1*dx);
    k3=feval(f,xi+0.5*dx,yi+0.5*k2*dx);
    k4=feval(f,xi+dx,yi+k3*dx);
    y(:,i+1)=yi+(k1+2*k2+2*k3+k4)*dx/6;
end
return


