% solve robertson's DAE 
function dae_example
% ode function
f=@robertson;
% initial conditions
y0 = [1; 0; 0];
% time span
tspan = [0 4*logspace(-6,6)];
% Mass matrix
M=eye(3); % initialize to identity matrix
M(3,3)=0; % third equation is the algebraic equation so set to 0
% set ode options
options = odeset('Mass',M); % to know is an ode algebraic system
% solve using ode15s
[t,y] = ode15s(f,tspan,y0,options);
% plot solution after scaling y2
y(:,2)=1e4*y(:,2);
semilogx(t,y);
legend('y1','y2*1e4','y3');
xlabel('t');
ylabel('y');
end

function [dydt] = robertson(t,y)
dydt=zeros(3,1);
dydt(1)= -0.04*y(1) + 1e4*y(2).*y(3);
dydt(2)= 0.04*y(1) - 1e4*y(2).*y(3) - 3e7*y(2).^2;
dydt(3)= y(1) + y(2) + y(3) - 1;
end



