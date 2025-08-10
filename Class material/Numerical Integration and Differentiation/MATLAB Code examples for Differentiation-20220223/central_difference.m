% computes nth order central differencing
% x = independent variable
% h = spacing
% n = order of differenciation
function [diff] = central_difference(f,x,h,n)
% starts from k=1
diff=central_differencen(f,x,h,1,n);
return

% recursive central differencing for nth order
function [diff] = central_differencen(f,x,h,k,n)
cdf=@central_differencen;
if (k == n)
    diff = (feval(f,x+h)-feval(f,x-h))/(2*h);
else
    diff = (feval(cdf,f,x+h,h,k+1,n)...
        -feval(cdf,f,x-h,h,k+1,n))/(2*h);
end
return