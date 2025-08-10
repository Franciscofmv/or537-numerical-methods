% chapra 21_14 (p 541)
% or chapra and canale 23_22 (p 670)
function radar_tracking
t=[200, 202, 204, 206, 208, 210];
theta=[0.75, 0.72, 0.70, 0.68, 0.67, 0.66];
r=[5120, 5370, 5560, 5800, 6030, 6240];
% calculate the following:
r_dot=tabular_derivative(r,t); % first derivative of r
theta_dot=tabular_derivative(theta,t); % f' of theta
r_doubledot=tabular_derivative(r_dot,t); % tabular derivtive of theta_dot (f'' of theta
theta_doubledot=tabular_derivative(theta_dot,t); % second derivative of r
% express v and a in terms of r and theta components
vr=r_dot; 
vtheta=r.*theta_dot;
ar=r_doubledot-r.*theta_dot.^2;
atheta=r.*theta_doubledot+2*r_dot.*theta_dot;
%determine v and a at t=206
v =[vr(4); vtheta(4)] % velocity vector
a = [ar(4); atheta(4)] % acceleration vector
% plot values
x=r.*cos(theta);
y=r.*sin(theta);
vx=vr.*cos(theta)- vtheta.*sin(theta);
vy=vr.*sin(theta)+ vtheta.*cos(theta);
quiver(x,y,vx,vy);
return
