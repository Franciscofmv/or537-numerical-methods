function romberg_test
n=10; a = 0; b =3; f = @xex; % since n=10, it will use 3 levels approxsqrt(10)
I_trap = trapezoidal(f,a,b,n)
I_romberg = romberg(f,a,b,n)
I_matlab = integral(f,a,b)
end

function [I] = romberg(f,a,b,N)
L=3; % set number of levels
IR=zeros(L,L); % initialize integral matrix
% loop over number of levels to calculate
% integral estimates at first level using trapezoidal rule
for i=1:L 
    n=N*2^(i-1);
    IR(i,1) = trapezoidal(f,a,b,n);
end
% calculate integral using Richardson's extrapolation
for j = 2:L
    for i = 1:L-j+1
        IR(i,j) = (4^(j-1)*IR(i+1,j-1) - IR(i,j-1))/(4^(j-1) - 1);
    end
end
I = IR(1,L); % final value in first row last column
end

function [I_trap] = trapezoidal(f,a,b,n)
I=0; h=(b-a)/n;
for i=0:n
    x=a+i*h;
    if (i == 0 || i == n)
        I=I+feval(f,x);
    else
        I=I+2*feval(f,x);
    end
end
I_trap=I*h/2;
end

function [f] = xex(x)
f = x.*exp(x);
end
