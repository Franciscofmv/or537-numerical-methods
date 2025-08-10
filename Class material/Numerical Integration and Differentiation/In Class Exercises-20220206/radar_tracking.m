% chapra 21_14 (p 541)
% or chapra and canale 23_22 (p 670)
function radar_tracking
t=[200, 202, 204, 206, 208, 210];
theta=[0.75, 0.72, 0.70, 0.68, 0.67, 0.66];
r=[5120, 5370, 5560, 5800, 6030, 6240];
r_dot=tabular_derivative(r,t);
theta_dot=tabular_derivative(theta,t);
r_doubledot=tabular_derivative(r_dot,t);
theta_doubledot=tabular_derivative(theta_dot,t);
% express v and a in terms of r and theta components
vr=r_dot; 
vtheta=r.*theta_dot;
ar=r_doubledot-r.*theta_dot.^2;
atheta=r.*theta_doubledot+2*r_dot.*theta_dot;
%determine v and a at t=206
v =[vr(4); vtheta(4)]
a = [ar(4); atheta(4)]
% plot values
x=r.*cos(theta);
y=r.*sin(theta);
vx=vr.*cos(theta)- vtheta.*sin(theta);
vy=vr.*sin(theta)+ vtheta.*cos(theta);
quiver(x,y,vx,vy);
return


% calculates derivatives at the same points supplied
function [dydx] = tabular_derivative(y,x)
n=length(x);
% use central for interior points
for i=2:n-1
    dydx(i)=(y(i+1)-y(i-1))/(x(i+1)-x(i-1));
end
% use forward for first point
dydx(1) = (y(2)-y(1))/(x(2)-x(1)); 
% use backward for last point
dydx(n)=(y(n)-y(n-1))/(x(n)-x(n-1));
return