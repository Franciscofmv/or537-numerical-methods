function nonlinear_stiff_system
y0=[52.29; 83.82]; t=0:0.005:1;
f=@stiff_system_function;
[y]=general_backward_euler_system(f,t,y0);
plot(t,y(1,:),'-g');
hold on;
[y]=multi_step_heuns(f,t,y0);
plot(t,y(1,:),'-r');
hold on;
[T,Y]=ode23s(f,[0 1],y0);
% y=runge_kutta(f,t,y0);
plot(T,Y(:,1),'-b');
hold on;
% hold on;
% clear
% syms t y1 y2;
% s=dsolve('Dy1 = -5*y1+3*y2','Dy2 = 100*y1-301*y2^3','y1(0) = 52.29','y2(0) = 83.82');
% t=[0:0.1:1];
% y1=subs(s.y1,t);
% plot(t,y1,'-rs');
legend('backward-euler','multi-step-heuns','ode23s');
% legend('backward-euler','ode23s');
hold off;
return

function [y] = general_backward_euler_system(df,t,y0)
n=length(t);
m=length(y0);
y=zeros(m,n);
I=eye(m);
y(:,1)=y0;
dy=zeros(1,m);
for i=1:n-1
    dt=t(i+1)-t(i);
    J=feval('numerical_jacobian',df,t(i+1),y(:,i));
    dy(:)=(I-dt*J)\(dt*feval(df,t(i+1),y(:,i)));
    y(:,i+1)=y(:,i)+dy(:);
end
return

function [y] = general_crank_nicolson_system(t,df,y0)
n=length(t);
m=length(y0);
y=zeros(m,n);
I=eye(m);
y(:,1)=y0;
dy=zeros(1,m);
for i=1:n-1
    dt=t(i+1)-t(i);
    t_half=(t(i+1)+t(i))/2;
    J=feval('analytical_jacobian',df,t_half,y(:,i));
    dy(:)=(I-0.5*dt*J)\(dt*feval(df,t_half,y(:,i)));
    y(:,i+1)=y(:,i)+dy(:);
end
return

function [dydt] = stiff_system_function(t,y)
dydt=zeros(2,1);
dydt(1) = -5*y(1)+3*y(2);
dydt(2) = 100*y(1)-301*y(2)^3;
return

function [y] = runge_kutta(df,x,y0)
m = length(y0);
n = length(x);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    k1 = feval(df,xi,yi);
    k2 = feval(df,xi+0.5*h,yi+0.5*k1*h);
    k3 = feval(df,xi+0.5*h,yi+0.5*k2*h);
    k4 = feval(df,xi+h,yi+k3*h);
    y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
end
return

function [J] = numerical_jacobian(f,t,y)
m=length(y);
dt=1e-4;
for i=1:m
    y1 = y;
    y2 = y;
    y1(i) = y(i)+dt;
    y2(i) = y(i)-dt;
    J(:,i) = (feval(f,t,y1)-feval(f,t,y2))/(2*dt);
end
return

function [y] = multi_step_heuns(df,x,y0)
n = length(x);
m = length(y0);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0; ms=5;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    dydx1=feval(df,xi,yi);
    ynext = yi+h*dydx1;
    for j=1:ms
        dydx2 = feval(df,x(i+1),ynext);
        ynext = yi + h*(dydx1+dydx2)/2;
    end
    y(:,i+1)=ynext;
end
return

function [J] = analytical_jacobian(f,t,y)
m=length(y);
J = [-5 3; 100 -903*y(2)^2];
return