function irrigation % main function
    x0 = [1 pi/4];% initial value of depth of water and initial angle thetha
    my_handle = @perimeter;
    x = fminsearch(my_handle,x0);
    d = x(1);
    theta=x(2);
    fprintf('x = %g\n',x)
    fprintf('d = %g\n',d)
    fprintf('theta = %g\n',theta)
end

function [L] = perimeter(x)
    d = x(1);
    theta=x(2);
    A = 100;
    b = A/d - d/tan(theta);
    L = b + 2*d/sin(theta);
end
