function simple_line_integral % dummy function
f = @xex;
a = 0; b = 3; n=100;  % integrating from 0 to 12
fplot(f,[a,b]); 
distance = trap(f,a,b,n)
distance_simp = simp(f,a,b,n) % simpson
distance_matlab = integral(f,a,b)
syms x;
s = int(x*exp(x),[0,3])
dist_anal=double(s); % distance analytical

% syms g m c t;
% s = int(sqrt(g*m/c)*tanh(sqrt(g*c/m)*t),t,[0,12]); % symbolic integration
% g = 9.81; m = 68.1;c = 0.25;% defining the values for constants
% dist_anal = double(subs(s)) % analytical value
end

% function [v] = bungee(t)% independent variable is t
% g = 9.81;%gravity
% m = 68.1;% mass
% c = 0.25;%
% v = sqrt(g*m/c)*tanh(sqrt(g*c/m)*t); % integrand function
% end

function [f] = xex(x)
f = x.*exp(x);
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