
function tabular_ex1
T = [0 1 2 3.25 4.5 6 7 8 8.5 9 10]; % time
V = [0 5 6 5.5 7 8.5 8 6 7 7 5]; %velocity
d= trap_tab(T,V);
a = 0; b = 10; n =10; f=@vel; % interpolation
f = @(t)vel(
d = trap(f,a,b,n)
d = integral(f,a,b) % matlan built in function
end

function [I]=trap_tab(x,y) % is IV y are the values
n = length(x);
I = 0;
for i = 1:n-1
    h = x(i+1)-x(i);
    I = I+h*(y(i+1)+y(i))/2;
end
end
function [v] = vel(T,)
v = interp1(T,V,t)
end

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