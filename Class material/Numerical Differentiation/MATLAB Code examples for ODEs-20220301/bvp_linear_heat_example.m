% heat example: linear two-point boundary value problem
function bvp_linear_heat_example
odefun = @heat;
bcfun = @heatbc;
xinit=[0:0.5:10];
yinit=[0 100]; % 
% matlab solution
solinit = bvpinit(xinit,yinit);% using bvpinit function
sol = bvp4c(odefun, bcfun, solinit);
xint=xinit;
yint = deval(sol,xinit); % function called deval, evaluating at xinit
% our solutions
T0=40;
TL=200;
x=[0:0.5:10];
T_fd=fda_heat(x,T0,TL); % fda
y=bvp_shooting(odefun,x,T0,TL); % shooting
T_sh=y(1,:);
% analytical solution
Ta=73.4523*exp(0.1*x)-53.4523*exp(-0.1*x)+20;
figure(1);
plot(xint,yint(1,:),'-ro',x,T_fd,'-g+',x,T_sh,'-k*',x,Ta,'-bs');
legend('bvp4c','fda','shooting','analytical');
hold off;
return

function [T] = fda_heat(x,T0,TL)
Ta=20; alpha=0.01; dx=x(2)-x(1);
m=length(x);
n=m-2; % only solving for the interior points
ld=ones(n,1); %lower diagonal is = 1
ud=ones(n,1);% upped diagonal is 1
d=zeros(n,1); % to initialize the vector

d(1:n)=-2-alpha*dx^2; % Ti

b=zeros(n,1);% initializing
b(1)=-dx^2*alpha*Ta-T0;

b(n)=-dx^2*alpha*Ta-TL; % compare to slide 79

b(2:n-1)=-dx^2*alpha*Ta;

%[ld d ud] marix nby3, n = size of matrix
A=my_spdiags([ld d ud],[-1,0,1]); % built in function that take the diagonals of matrix, use myspdiags instead
y=A\b;
T(2:m-1)=y;
T(1)=T0;
T(m)=TL;
return
 % see adjusted version in lecture

function [y] = bvp_shooting(odefun,x,yi,yf) %yi=T0 initial temp 40, yf=T_L 
n=length(x);
x0=x(1);
xL=x(n);
y0=zeros(2,1); % initial guess of IC
% first shot
q1=0;
y=runge_kutta(odefun,x,[yi; q1]);
p1=y(1,n); % end point
% second shot
q2=1; % assume, slope of 45 degrees
y=runge_kutta(odefun,x,[yi; q2]);
p2=y(1,n); % 
y2=q1+(yf-p1)*(q2-q1)/(p2-p1); % linear interpolation for slope
y=runge_kutta(odefun,x,[yi; y2]);
return

function [res] = heatbc(ya,yb)
res = [ya(1)-40; yb(1)-200];
return

function [dydx] = heat(x,y)
Ta = 20;
alpha = 0.01;
dydx=zeros(2,1);
dydx(1) = y(2);
dydx(2) = alpha*(y(1)-Ta);
return

