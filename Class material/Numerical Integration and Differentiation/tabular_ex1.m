function tabular_ex1()
T = [0 1 2 3.25 4.5 6 7 8 8.5 9 10]; % time s
V = [0 5 6 5.5 7 8.5 8 6 7 7 5]; % velocity m/s
d_tabular = trap_tab(T,V)
f =@(t)vel(t,T,V);
interpolation = trap(f,0,10,200)
built_in_integral = integral(f,0,10)
d_simpson = simp(f,0,10,200)
end
function [I] = trap_tab(x,y) %x y are vectors, tabular data
n = length(x);
I = 0;
    for i = 1:n-1 % Tabular integration function
        h = x(i+1)-x(i);
        I = I+h*(y(i+1)+y(i))/2;
    end
end

function [v] = vel(t,T,V)
v = interp1(T,V,t);
end