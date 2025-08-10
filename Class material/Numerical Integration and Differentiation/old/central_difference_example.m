function central_difference_example 
f=@test_function;
x=1; h=1e-3;
diffn=central_differencem(f,x,h,1,3)
syms x;
diffa=subs(diff(x*exp(x),3),1)
return

function [f] = test_function(x)
f = x*;
return

function [diff] = central_difference1(f,x,h)
diff = (feval(f,x+h)-feval(f,x-h))/(2*h);
return

function [diff] = central_difference2(f,x,h)
diff = (feval(f,x+h)-2*feval(f,x)+feval(f,x-h))/(h^2);
return

function [diff] = central_differencem(f,x,h,k,m)
cdf=@central_differencem;
if (k == m)
    diff = (feval(f,x+h)-feval(f,x-h))/(2*h);
else
    diff = (feval(cdf,f,x+h,h,k+1,m)...
        -feval(cdf,f,x-h,h,k+1,m))/(2*h);
end
return