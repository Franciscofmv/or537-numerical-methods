function wtf
quiz1(1,2,3)
end

function [v] = quiz1(ta,mass,dragc)
syms g m c v(t);
s=dsolve(diff(v) == g - c*v/m, v(0) == 0);
g = 9.8; m = mass; t = ta; c = dragc;
v = subs(s)
end