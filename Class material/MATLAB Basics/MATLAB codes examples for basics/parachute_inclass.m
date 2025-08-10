%% Falling parachutist problem
%%
function parachute_inclass
% define constants
g = 9.81; c = 12.5; m = 68.1;
dt = 0.1;
% define t vector
t = 0:dt:12; % row vector
%% calculate numerical solution
[v] = feval(@numerical,g,m,c,t);
%% calculate analytical solution
[ta,va] = feval(@analytical);
%% plot results
plot(t,v,'-g');
hold on;
plot(ta,va,'ro');
xlabel('time (s)');
ylabel('velocity (m/s)');
legend('numerical','analytical');
hold off;
return

%% numerical function
function [v] = numerical(g,m,c,t)
n = length(t);
v = zeros(1,n);
v(1) = 0; % initial condition
% for loop to do time marching
for i = 1:n-1
    dt = t(i+1) - t(i);
    v(i+1) = v(i) + dt*(g - c*v(i)/m);
end
return

%% analytical function
function [t,v] = analytical
syms g m c t v;
v=dsolve('Dv = g - c*v/m','v(0) = 0'); % analytical solution by dsolve
g = 9.81; c = 12.5; m = 68.1; t =0:12;
v = subs(v);
%v = (g*m/c)*(1-exp(-c*t/m)); % analytical solution by hand
return