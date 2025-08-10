function mid_point_example
f=@one_over_sqrtx;
a=0; b=1; n=10;
I=mid_point(f,a,b,n)
return

function [I] = mid_point(f,a,b,n)
I=0; h=(b-a)/n; a0=0.5*h;
for i=0:n-1
    x=a0+i*h;
    I=I+feval(f,x);
end
I=I*h;
return

function [f] = one_over_sqrtx(x)
f=1/sqrt(x);
return
