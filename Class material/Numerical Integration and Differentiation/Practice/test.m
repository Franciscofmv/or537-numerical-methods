
function test
g = @xex
f = @composite_simpson
f1 = @composite_simpson3_8
d1=f1(g,0,10,300)
d = integral(g,0,10)
end
function [y] = xex(x)
y = x+x.^2;
end