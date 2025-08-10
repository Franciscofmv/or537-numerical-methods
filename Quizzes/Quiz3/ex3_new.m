
M(1.5)
shear(1.5)
function [M] = M(x)
w0=10^3*(250);L=3;
M = (w0/(120*L))*(-20*x^3 + 12*x*L^2);
end

function [V] = shear(x)
w0=10^3*(250);L=3;
V = (w0/(120*L))*(-60*x^2 + 12*L^2);
end