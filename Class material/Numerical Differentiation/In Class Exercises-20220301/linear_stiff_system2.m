function linear_stiff_system2
y0=[52.29; 83.82];
t=0:0.1:1;
A=@Afunc;
c=@cfunc;
y=linear_backward_euler(t,y0,A,c);
plot(t,y(1,:),'-go',t,y(2,:),'-rs');
hold on;
[t,y] = ode45(@odefunc, t,y0);
plot(t,y(:,1),'-b',t,y(:,2),'-k');
legend('y1','y2','y1-ode45','y2-ode45');
xlabel('t(s)');
ylabel('y(m)');
return

function [y] = linear_backward_euler(t,y0,Afunc,cfunc)
n=length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    dt=t(i+1)-t(i);
    I=eye(m);
    A=feval(Afunc,t(i+1));
    c=feval(cfunc,t(i+1));
    y(:,i+1)=inv(I-dt*A)*(y(:,i)+dt*c); % (I-dt*A)\(y(:,i)+dt*c)
end
return

function [A] = Afunc(t)
A=[-5, 3
    100 -301];
return

function [c] = cfunc(t)
c = [0
    0];
return

function [dydt] = odefunc(t,y)
dydt=zeros(2,1);
dydt(1)=-5*y(1)+3*y(2);
dydt(2)=100*y(1)-301*y(2);
return
