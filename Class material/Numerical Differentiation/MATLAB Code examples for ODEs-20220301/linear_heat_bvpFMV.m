function linear_heat_bvpFMV
odefun=@heatode;
bcfun=@heatbc
xinit = 0:10;
solinit=bvpinit(xinit,[100,0]); % bvp init is a build in function, initial temp is 100 degress,slope =0. 
x = 0:0.1:10;
sol = bvp4c(odefun,bcfun,solinit);
y=deval(sol,x);
plot(x,y(1,:),'ro')
xlabel('Distance (m)');
ylabel('temperature');

end

function [dydx] = heatode(x,y) % x IV, T DV; T=  y1, dT/dx = y2
dydx = zeros(2,1);
alpha = 0.01;
Ta = 20;
dydx(1) = y(2);
dydx(2) = alpha*(y(1)-Ta);
end

function [res] = heatbc(ya,yb) % boundary conditions
res=[ya(1)-40
    yb(1) - 200];
end