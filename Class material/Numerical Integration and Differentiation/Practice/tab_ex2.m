function tab_ex2
T = [0 10 20 30 35 40 45 50]; % time(minutes)
Q = [4 4.8 5.2 5.0 4.6 4.3 4.3 5.0];% m^3/min
C = [10 35 55 52 40 37 32 34]; % mg/m^3
QC = Q.*C;
f1 = @trap_tabular;
mass_trapezoidal_tab = f1(T,QC)
f2 = @(t)mass(t,T,QC);
mass_interp1 = integral(f2,0,50)
mass_simpson = simp(f2, 0, 50, 10)
end

function [m] = mass(t,T,QC)
m = interp1(T,QC,t);
end