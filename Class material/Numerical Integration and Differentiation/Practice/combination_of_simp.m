function [I] = combination_of_simp(f,a,b,n)
m = n;
if n>4
f_1_3 = @composite_simpson;
f_1_8 = @composite_simpson3_8
    if (mod(n,2)==0)
        f_1_3(f,a,b,n)
    elseif (mod(n,2)==1)
        
        for i=1:(m-3)
            global I_1 = f_1_3(f,a,b,n);
        end
        
        for j=(m-3):m
           global I_2 = f_1_8(f,a,b,(m-3));
        end
    end
    I = I_1+I_2
else
    disp("n must be grater than 4")
end

end