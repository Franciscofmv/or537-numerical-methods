% chapra and canale p25.1 p. 750
function cc_25_1
t=0:0.1:1;
y0=1;
f=@func;
y=RK4(f,t,y0);
plot(t,y);
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

function [dydt] = func(t,y)
dydt=y*t^3-1.5*y;
return
