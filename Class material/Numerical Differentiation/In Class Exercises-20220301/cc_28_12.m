% chapra and canale p28.12 p. 828
function cc_28_12
C0=[1; 1; 0];
f=@reactor;
[T,C]=ode23s(f,[0 50],C0);
plot(T,C(:,1),'-g',T,C(:,2), '-r', T, C(:,3), '-b');
legend('c1','c2','c3');
return

function [dcdt] = reactor(t,c)
dcdt=zeros(3,1); 
c1=c(1); c2=c(2); c3=c(3);
dcdt(1)=-0.013*c1-1000*c1*c3;
dcdt(2)=-2500*c2*c3;
dcdt(3)=-0.013*c1-1000*c1*c3-2500*c2*c3;
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


