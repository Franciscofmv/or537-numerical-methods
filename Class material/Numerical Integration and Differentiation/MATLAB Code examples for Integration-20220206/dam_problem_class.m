function dam_problem_class
Z = [0 10 20 30 40 50 60];
W = [122 130 135 160 175 190 200];
rho = 1000; g = 9.81; D=60;
P = rho*g*W.*(D-Z);
M = rho*g*W.*Z.*(D-Z);
f1 = @(z) pressure_function(z,Z,P);
f2 = @(z) moment_function(z,Z,M);M
% limits
a=0; b=60; n=10;
force_trap_data = trap_data(Z,P)
depth_trap_data = trap_data(Z,M)/force_trap_data
force_simpson = simpson(f1,a,b,n)
depth_simpson = simpson(f2,a,b,n)/force_simpson
force_integral = integral(f1,a,b)
depth_integral = integral(f2,a,b)/force_integral
return

function [pz] = pressure_function(z,Z,P)
pz=interp1(Z,P,z);
return

function [mz] = moment_function(z,Z,M)
mz=interp1(Z,M,z);
return