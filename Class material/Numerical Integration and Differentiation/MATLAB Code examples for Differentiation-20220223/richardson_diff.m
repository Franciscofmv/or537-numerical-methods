function [D]=richardson_diff(f,x,h,n)
D1=central_difference(f,x,h,n);
D2=central_difference(f,x,2*h,n);
D=(4*D1-D2)/3; % weighted average
end

