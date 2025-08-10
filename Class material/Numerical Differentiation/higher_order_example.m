function higher_order_example
x=0:0.5:5;
y0 = [2; 0]; % Initial Conditions
f = @odefun;
y = runge_kutta(f,x,y0);
y1 = y(1,:);
y2 = y(2,:);
plot(x,y1,'-go',x,y2,'-ro'); % what -go=- means line green circle and -ro means
hold on;
y1 = y(:,1);
y2 = y(:,2);
%plot(x,y1,'-bo',x,y2,'-yo'); % need to fix this line
end
function [dydx] = odefun(x,y)
dydx = zeros(2,1);
dydx(1) = y(2);
dydx(2) = -0.6*y(2)-8*y(1);
end


function [y] = runge_kutta(f,x,y0)
n=length(x);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    k1 = f(xi,yi);
    k2 = f(xi+0.5*h,yi+0.5*k1*h);
    k3 = f(xi+0.5*h,yi+0.5*k2*h);
    k4 = f(xi+h,yi+k3*h);
    y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6; % this is a matrix
end
end
