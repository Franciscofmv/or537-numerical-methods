function non_linear_stiff_system
y0=[52.29; 83.82];
t=0:0.001:1;
J=@Jfunc;
F=@odefunc;
y=non_linear_backward_euler(t,y0,J,F);
plot(t,y(1,:),'-g',t,y(2,:),'-r');
hold on;
t=0:0.1:1;
[t,y] = ode45(@odefunc, t,y0);
plot(t,y(:,1),'go',t,y(:,2),'rs');
legend('y1','y2','y1-ode45','y2-ode45');
xlabel('t(s)');
ylabel('y(m)');
hold off;
return

% basic nonlinear backward euler
function [y] = non_linear_backward_euler(t,y0,Jfunc,Ffunc)
n=length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
I=eye(m);
for i=1:n-1
    dt=t(i+1)-t(i);
    J=feval(Jfunc,t(i+1),y(:,i));
    F=feval(Ffunc,t(i+1),y(:,i));
    dy=inv(I-dt*J)*dt*F;
    y(:,i+1)=y(:,i)+dy;
end
return

% iterative version of backward euler
function [y] = non_linear_backward_euler_itr(t,y0,Jfunc,Ffunc)
n=length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
maxitr=2;
for i=1:n-1
    dt=t(i+1)-t(i);
    I=eye(m);
    y(:,i+1)=y(:,i);
    tol=inf; itr=0;
    while (tol > 1e-3 && itr < maxitr)
        J=feval(Jfunc,t(i+1),y(:,i+1));
        F=feval(Ffunc,t(i+1),y(:,i+1));
        dy=inv(I-dt*J)*dt*F;
        tol=norm(dy);
        y(:,i+1)=y(:,i)+dy;
        itr=itr+1;
    end
end
return

% function for analytical jacobian
function [J] = Jfunc(t,y)
J=[-5, 3
    100, -903*y(2)^2];
return


% ode function
function [dydt] = odefunc(t,y)
dydt=zeros(2,1);
dydt(1)=-5*y(1)+3*y(2);
dydt(2)=100*y(1)-301*y(2)^3;
return
