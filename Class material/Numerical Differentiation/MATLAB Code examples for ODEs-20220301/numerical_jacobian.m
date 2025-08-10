function [J] = numerical_jacobian(f,t,y)
m=length(y);
dt=1e-4;
for i=1:m
    y1 = y;
    y2 = y;
    y1(i) = y(i)+dt;
    y2(i) = y(i)-dt;
    J(:,i) = (feval(f,t,y1)-feval(f,t,y2))/(2*dt); % central differencin
end
return
