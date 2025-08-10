function exercise4
f1=@integrand1;
f2=@integrand2;
f3=@integrand3;
g=@simpson_one_third
I_built = integral(f1,0,10)+integral(f2,10,20)+integral(f3,20,30)
I = g(f1,0,10,6)+g(f2,10,20,6)+g(f3,20,30,6)-26833.33
end

function [v] = integrand1(t)
v = 11*t.^2-5*t; % for 0<=t<=10
end
function [v] = integrand2(t)
v = 1100-5*t; % for 10<=t<=20
end
function [v] = integrand3(t)
v = 50*t+2*(t-20).^2; % for 20<=t<=30
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

