% computes nth order central differencing
% x = independent variable
% h = spacing
% n = order of differenciation

% f=handle of function,x IV, h=spacing, n=order of differentiation
function [diff] = central_difference(f,x,h,n)
% starts from k=1
diff=central_differencen(f,x,h,1,n); 
return

% recursive central differencing for nth order
% when k is = to n it takes central diifference
function [diff] = central_differencek(f,x,h,k,n) 
cdf=@central_differencen;
if (k == n)
    diff = f((x+h)-f(,x-h))/(2*h);
else
    diff = (cdf(f,x+h,h,k+1,n)...
        -feval(cdf(f,x-h,h,k+1,n))/(2*h);
end
return