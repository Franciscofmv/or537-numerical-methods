function mdim_example
f=@fxy;
a=[-3 0]; b = [1 2]; n=10;
I=msimp(f,a,b,n)
I=integral2(@(x,y)y.^4.*(x.^2+x.*y),a(1),b(1),a(2),b(2))
end

function [f] = fxy(X)
x=X(1);
y=X(2);
f=y^4*(x^2+x*y);
end