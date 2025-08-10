
function [I] = simp(f,a,b,n) % simpson
h = (b-a)/n; % segment size
I = 0;
    for i=0:n %summation for loop, end points have weight of 1
        x = a + i*h; % defining x, i=0 then x =a, i = n then x = b
        if (i==0 || i==n)
            I = I + f(x); % end points
        elseif(mod(i,2)==0) % check if i is even
            I = I + 2*f(x); % even interior point
        else
            I = I + 4*f(x); % odd interior points
        end
    end
    I = I*h/3;
end

function [I] = simpson(f,a,b,n)

    h=(b-a)/n;
    I=0;
    for i=0:n
        x = a + i*h;
        if(i==0||i==n)
            I = I + f(x);
        elseif (mod(n,2)==1 && i~=1)
            I = I + 4*f(x);
        else
            I = I + 2*f(x);
        end
    end
    I = I*h/3;
end