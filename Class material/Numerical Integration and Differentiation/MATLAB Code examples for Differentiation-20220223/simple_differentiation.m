function simple_differentiation
f = @xex;
x=1; h=1e-3;
diff_f=central_difference(f,x,h,3)
clear
syms x;
diff_anal= double(subs(diff(x*exp(x),3),x,1))

% we can see numerical and analytical are the same
return

function [f] = xex(x)
f = x*exp(x);
return
