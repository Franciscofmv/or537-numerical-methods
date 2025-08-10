init=[0 0]; % initial conditions
t=linspace(0,0.4);
sol = runge_kutta(@odesyst,t,init)
ode45(@odesyst,t,init)




function [dqdt] = odesyst(t,i)
dqdt=zeros(2,1);
L=1;c=0.25;w=sqrt(3.5);R=0.025;
dqdt(1) = i(2);
dqdt(2) = (1/L)*(EE(t)-R*i(2)-i(1)/c);
end
function [E] = EE(t)
w=sqrt(3.5);
if t==0
    E = 1;
else
    E = sin(w*t);
end
end


function [y] = runge_kutta(df,x,y0)
m = length(y0);
n = length(x);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
h = 0.1;
xi=x(i);
yi=y(:,i);
k1 = feval(df,xi,yi);
k2 = feval(df,xi+0.5*h,yi+0.5*k1*h);
k3 = feval(df,xi+0.5*h,yi+0.5*k2*h);
k4 = feval(df,xi+h,yi+k3*h);
y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
end
end