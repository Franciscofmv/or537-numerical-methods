function tabular_ex1
T = [0 1 2 3.25 4.5 6 7 8 8.5 9 10]; % time vector in minutes
V = [0 5 6 5.5 7 8.5 8 6 7 7 5]; %  velocity vector
f = @trap_tabular; % arguments are x,y
distance_trap_tab = f(T,V)
f0 = @trapezoidal;
g = @(t)vel(t,T,V); % plug in T,v as defined above, Argument is t
distance_interpolated = f0(@(t)vel(t,T,V),0,10,200) % interpolating
distance_built = integral(@(t)vel(t,T,V),0,10)
distance_simpson = simp(g,0,10,100)
end

function [v] = vel(t,T,V)
v = interp1(T,V,t);
end

