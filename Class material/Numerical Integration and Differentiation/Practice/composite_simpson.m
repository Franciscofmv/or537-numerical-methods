function composite_simpson(f,a,b,n)
if (mod(n,2)==0)
    I =0;
    h = (b-a)/n;
    for i=0:n
        x = a + i*h;
        if (i==0||i==1)
            I = I + f(x);
        elseif (mod(i,2)==1)
            I = I + 4*f(x);
        else
            I = I+ 2*f(x);
        end
    end
else
    disp(" m must be even")
  
end
I = (h/3)*I
end