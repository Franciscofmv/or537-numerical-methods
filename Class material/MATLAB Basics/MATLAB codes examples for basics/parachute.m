% parachute example
function parachute
dt=0.1; tfinal=12;
t=0:dt:tfinal;
v=parachute_numerical(t);
hold off;
plot(t,v,'-g');
hold on;
t=0:tfinal;
v=parachute_analytical(t);
plot(t,v,'-rs');
legend('numerical','analytical');
xlabel('time (s)');
ylabel('velocity (m/s)');
return

function [v] = parachute_numerical(t)
n=length(t);
v=zeros(1,n);
% inital condtion
v(1)=0;
f=@parachute_func;
for i=1:n-1
    dt=t(i+1)-t(i);
    v(i+1)=v(i)+dt*feval(f,v(i));
end
return

function [v] = parachute_analytical(ta)
syms g m c v(t);
s=dsolve(diff(v) == g - c*v/m, v(0) == 0);
g=9.8; m=68.1; c=12.5; t=ta;
v=subs(s);
return

function [dvdt] = parachute_func(v)
g=9.8; m=68.1; c=12.5;
dvdt=g-c*v/m;
return
