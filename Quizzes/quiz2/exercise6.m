function exercise6
f=@simpson_one_third
built_in = integral(@integrand,0,30)
I4 = f(@integrand,0,30,4)
I8 = f(@integrand,0,30,8)
I16 = f(@integrand,0,30,16)

end

function [y] = integrand(x) 
y= (1.6*x - 0.045*x.^2).*cos(0.8 + 0.125*x-0.009*x.^2 + 0.0002*x.^3);
end


function [I] = simpson_one_third(f,a,b,n)
if (mod(n,2)==0)
    I=0;
    h = (b-a)/n;
    for i =0:n
        x = a + i*h;
        if (i==0||i==n)
            I = I + f(x);
        elseif (mod(i,2)==1)
            I = I + 4*f(x);
        else
            I = I + 2*f(x);
        end
    end
    I = I*h/3;
else
    disp("n has to be even")
end
end

