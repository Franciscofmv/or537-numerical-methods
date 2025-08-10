
function [I] = trap(f,a,b,n) % function,a,b,n is number of segments, trapezoidal rule
h = (b-a)/n; % segment size
I = 0;
    for i=0:n %summation for loop
        x = a + i*h; % defining x, i=0 then x =a, i = n then x = b
        if (i==0 || i==n)
            I = I + f(x); % end points
        else
            I = I + 2*f(x); % interior points
        end
    end
    I = I*h/2;
end