function [I] = composite_simpson3_8(f,a,b,n)
if (mod(n,3)==0)
    h = (b-a)/n;
    I = 0;
    for i=0:n
        x = a + i*h;
        if (i==0 || i==n)
            I = I +f(x);
        elseif (mod(i,3)==1)
            I = I + 3*f(x);
        elseif (mod(i,3)==2)
            I = I + 3*f(x);
        else
            I = I + 2*f(x);
        end
      
    end
    I = (3*h/8)*I
else
    disp("n must be multiple of 3")
end
end