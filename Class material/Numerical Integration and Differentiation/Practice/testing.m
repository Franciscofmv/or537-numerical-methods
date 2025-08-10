function testing
f = @xex;
a = 0; b = 3; n=100; 
g = @composite_trapezoidal
 g(f,a,b,n)
 integral(f,a,b)
end


function [I] = composite_trapezoidal(f,a,b,n)
I = 0;
h = (b-a)/n;
    for i = 0:n
        x = a + i*h;
        if (i==0 || i==n)
            I = I + f(x);
        else
            I = I + 2*f(x);
        end 
    end
    I = I*h/2
end

function [f] = xex(x)
f = x.*exp(x);
end