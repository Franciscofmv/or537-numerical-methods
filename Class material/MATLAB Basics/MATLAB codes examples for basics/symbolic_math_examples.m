function [] = symbolic_math_examples
syms x;
f1=expand((x-1)^3) % expand (x-1)^3
f2=expand((2*x-1)^2) % expand (2x-1)^2
f = f1/f2;  % f = (x-1)^3/(2x-1)^2
% f = (x^3-3*x^2+3*x-1)/(1-4*x+4*x^2); % f = (x-1)^3/(2x-1)^2
simple_f=simple(f) % simplest form of f
diff_f=diff(f)     % differentiate f
int_f= int(f)      % integrate f
factor_f=factor(f) % factor f
taylor_f=taylor(f) % taylor series expansion about x=0
return
