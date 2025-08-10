% basic nonlinear backward euler
function [y] = nonlinear_backward_euler(Ffunc,t,y0,Jfunc) % function, t vector, y0 initial cond., jacobian
n=length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
I=eye(m);
for i=1:n-1
    dt=t(i+1)-t(i);
    J=feval(Jfunc,t(i+1),y(:,i));
    F=feval(Ffunc,t(i+1),y(:,i)); % right hand side
    dy=(I-dt*J)\(dt*F);% backlash
    y(:,i+1)=y(:,i)+dy;
end
return