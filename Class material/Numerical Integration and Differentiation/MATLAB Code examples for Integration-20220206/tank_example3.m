%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% tank mass example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% varying height with decay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tank_example3
n=10;
T=3600; dt=10; nt=T/dt;
M=zeros(1,nt); t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
    f=@(x)tank_conc(x,t(i));
    zmax=@(x)tank_height(x,t(i));
    a={-2,-2,0};
    b={2,2,zmax};
    M(i)=gsimp(f,a,b,n);
end
plot(t,M);
fprintf('Mass in the tank at 1 hour is %g  Kg\n',M(nt));
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions to be evaluated follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f]=tank_conc(x,t)
% use notation x = x(1); y=x(2); z=x(3) etc...
C0=0.1; k=0.001;
f=C0*exp(-k*t)*(x(1)^2+x(2)^2)/sqrt(x(3)^2+1);
return

function [h] = tank_height(x,T)
A=pi*4;
h=gsimp(@qfunc,0,T,100)/A;
return

function [q] = qfunc(t)
Q0=0.1; omega=pi/6;
q=Q0*(1+sin(omega*t));
return
