function exercise1
I1 = simpson_one_third(@integrand1,0,1,8)
I2 = mid_point(@integrand2,1,0,8)
I = I1+I2
I_built = integral(@integrand1,0,100)
end


function [v] = integrand1(y)
w= (1+y.^2).*(1+(y.^2/2));
v=1./(w);
end

function [v] = integrand2(z)
w=(1+(1./z.^2)).*(1+(1./(2*z.^2)))
v = (-1./z.^2).*1./(w);
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

function [I]= mid_point(f,a,b,n)
I=0;
x=0;
h=(b-a)/n;
for i=0:n-1
    x=a+(i+0.5)*h;
    I=I+f(x);
end
I=I*h;
end