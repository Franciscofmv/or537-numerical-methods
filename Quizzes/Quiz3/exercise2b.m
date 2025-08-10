
f=@odefunc;
y0=0;

range=0:0.1:1.6;
t = Euler(f,range,0)
[a b]=ode45(f,range,0)



function [g] = odefunc(x,y)
g =  (10^(3)*250/(120*3))*( -20*x^3 + 2*6*3^2*x);
end



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
end

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
end


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
end