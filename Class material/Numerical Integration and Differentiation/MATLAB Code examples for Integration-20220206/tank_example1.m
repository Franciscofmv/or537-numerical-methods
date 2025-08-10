%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% tank mass example with no decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% height is time dependent but obtained at a fixed time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tank_example1
n=100;
a=[-1,-1,0];
b=[1,1,3];
f=@tank_conc;
% evaluate mass with gsimp
M=msimp(f,a,b,n)
% evaluate with integral3
f=@tank_conc_integral3;
M=integral3(f,-1,1,-1,1,0,3)
fprintf('Mass in the tank at 1 hour is %g  Kg\n',M);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrand function for gsimp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f]=tank_conc(X)
% use notation x = X(1); y=X(2); z=X(3) etc...
C0=0.1; x=X(1); y=X(2); z=X(3);
f=C0*(x^2+y^2)/sqrt(z^2+1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrand function for integral3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f]=tank_conc_integral3(x,y,z)
C0=0.1;
f=C0*(x.^2+y.^2)./sqrt(z.^2+1);
return
