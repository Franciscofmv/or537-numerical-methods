function poly_examples
% polynomial x^3-2x-5
p = [1 0 -2 -5];
% represent in symbolic notation
poly2sym(p)
% solve using roots
x = roots(p)
% solve using fzero
x0 = 1;
x = fzero(@poly_fun,x0)
% other polynomial functions
y = polyval(p,x) % evaluate polynomail p at a vector x
p = poly(x)      % construct polynomial with roots x
% multiply two polynomials
p = [1 0 -2 -5]; % p = x^3-2x-5
q = [1 -4 0 0 -2]; % q = x^4 - 4x^3 -2
r = conv(p,q) % result is a vector of length(p)+length(q)-1
% divide two polynomials
[u, v] = deconv(q,p)
% evaluate derivative of a polynomial
d = polyder(q)
return

function [p] = poly_fun(x)
p = x^3-2*x-5;
return
