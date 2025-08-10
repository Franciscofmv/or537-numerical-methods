function reservoir_draining
f=@depth;
h=6:-0.5:0;
t0=0;
t=runge_kutta(f,h,t0);
n=length(h);
tf=t(n);
plot(t,h,'-bo');
hold on;
[H,T]=ode45(f,[6 0],t0);
n=length(H);
plot(T,H,'r');
legend('rk4','ode45');
xlabel('time (s)');
ylabel('h (m)');
tf=T(n);
fprintf('Drain time = %g s\n',tf);
hold off;
return

% dam reservoir function
function [dtdh] = depth(h,t)
g=9.81; d=0.25; e=1;
H=[6, 5, 4, 3, 2, 1, 0];
A=[1.17 0.97, 0.67, 0.45, 0.32, 0.18, 0];
Ah=interp1(H,A,h); % doing interpolation inside to calculate ODE
dtdh=-4*Ah/(pi*d^2*sqrt(2*g*(h+e))); % ODE function, we are solving for t->dependant variable
return
