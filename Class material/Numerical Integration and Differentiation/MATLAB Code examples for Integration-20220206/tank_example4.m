%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate average mass in tank over time T 
% with variable height and time dependendent concentration
% uses 'msimp' and 'gsimp' functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tank_example4
f=@tank_conc; T=3600;
a={-2,-2,0,0};
b={2,2,T,@tank_height};
n = 10;
MT=gsimp(f,a,b,n);
fprintf('Average mass in the tank through 1 hour is %g  Kg\n',MT/T);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions to be evaluated follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f]=tank_conc(X)
% use notation x = x(1); y=x(2); z=x(3) etc...
C0=0.1; k=1e-3; x=X(1); y=X(2); t=X(3); z=X(4);
f=C0*exp(-k*t)*(x^2+y^2)/sqrt(z^2+1);
return

function [h] = tank_height(x)
A=4*pi; T = x(3); n = 100;
h=msimp(@qfunc,0,T,n)/A;
return

function [q] = qfunc(t)
Q0=0.1; omega=pi/6;
q=Q0*(1+sin(omega*t));
return