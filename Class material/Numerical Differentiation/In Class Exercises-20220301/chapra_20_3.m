% chapra 20_3
function chapra_20_3
f=@func;
t=0:0.25:3;
y0=1;
% runge kutta
y=rk4(f,t,y0);
plot(t,y,'-ro');
hold on;
%euler
y=euler(f,t,y0);
plot(t,y,'-gs');
hold on;
%heuns
y=heuns(f,t,y0);
plot(t,y,'-k^');
hold on;
% analytical
ezplot('t^2 - 1/exp(t) - 2*t + 2',[0,3]);
legend('rk4','euler','heuns','analytical');
hold off;
return

function [y] = euler(f,x,y0)
n=length(x);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    y(i+1)=y(i)+feval(f,x(i),y(i))*dx;
end
return

function [y] = rk4(f,x,y0)
n=length(x);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    k1=feval(f,x(i),y(i));
    k2=feval(f,x(i)+0.5*dx,y(i)+0.5*k1*dx);
    k3=feval(f,x(i)+0.5*dx,y(i)+0.5*k2*dx);
    k4=feval(f,x(i)+dx,y(i)+k3*dx);
    y(i+1)=y(i)+dx*(k1+2*k2+2*k3+k4)/6;
end
return

function [y] = ralston(f,x,y0)
n=length(x);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    k1=feval(f,x(i),y(i));
    k2=feval(f,x(i)+0.75*dx,y(i)+0.75*k1*dx);
    y(i+1)=y(i)+dx*(k1+2*k2)/3;
end
return

function [y] = heuns(f,x,y0)
n=length(x);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    % predictor
    ftemp=feval(f,x(i),y(i));
    ytemp=y(i)+ftemp*dx;
    % corrector
    y(i+1)=y(i)+(ftemp+feval(f,x(i+1),ytemp))*dx/2;
end
return

function [y] = mid_point(f,x,y0)
n=length(x);
y(1)=y0;
for i=1:n-1
    dx=x(i+1)-x(i);
    % predictor
    ftemp=feval(f,x(i),y(i));
    ytemp=y(i)+ftemp*dx/2;
    % corrector
    y(i+1)=y(i)+feval(f,x(i)+dx/2,ytemp)*dx;
end
return

function [dydt] = func(t,y)
dydt=-y+t^2;
return