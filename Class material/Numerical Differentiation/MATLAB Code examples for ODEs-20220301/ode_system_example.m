% simple ode system example from lecture slides
function ode_system_example
x = (0:0.1:2); % x goes from 0 to 2, we can use any spacing we want 
f=@odefunc; % simple function
y0 = [4; 6]; % initial condition - column vector
[y] = runge_kutta(f,x,y0); % our rk4 solver
y1=y(1,:); y2=y(2,:); % extract values from solution
plot(x,y1,'-b',x,y2,'-r'); % plot rk4 results
hold on;
[x,y] = ode45(f,x,y0); % matlab's ode45 solver can use ode45(f,[0 2],y0)
y1=y(:,1); y2=y(:,2); % extract values from solution
plot(x,y1,'bo',x,y2,'ro'); %plot ode45 solution
% analytical solution, another way to solve an ODE
syms y1(x) y2(x);
eqns = [diff(y1,x)==-0.5*y1, diff(y2,x)==4-0.3*y2-0.1*y1];
conds = [y1(0)==4, y2(0)==6];
sol=dsolve(eqns,conds);
x=(0:0.1:2); % define numerical values for x
y1=double(subs(sol.y1,x)); % extract y1 solution at x
y2=double(subs(sol.y2,x)); % extract y2 solution at x
plot(x,y1,'-k',x,y2,'-g');  % plot analytical solution
legend('rk4-y1','rk4-y2','ode45-y1','ode45-y2','analytical-y1','analytical-y2');
hold off;
return


% simple function, each ODE is given in each row
function [dydx] = odefunc(x,y)
dydx=zeros(2,1);
dydx(1)= -0.5*y(1);
dydx(2) = 4-0.3*y(2)-0.1*y(1);
return

