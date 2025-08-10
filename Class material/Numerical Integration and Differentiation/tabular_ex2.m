function tabular_ex2()
T = [0 10 20 30 35 40 45 50]; 
Q = [4 4.8 5.2 5 4.6 4.3 4.3 5]; 
C = [10 35 55 52 40 37 32 34];
d_tabular = trap_tab(T,Q.*C)
f =@(t)qc(t,T,Q,C);
interpolation = trap(f,0,50,200)
built_in_integral = integral(f,0,50)
d_simpson = simp(f,0,50,200)
end

function [v] = qc(t, T,Q,C)
    v = interp1(T,Q.*C,t)
end