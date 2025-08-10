% line integral function
function bungee_jumper_distance
a=0; b = 10; n=11;
f = @vel;
distance_trap = trap(f,a,b,n)
distance_simpson = simpson(f,a,b,n)
distance_anal = distance(10)
distance_quad = quad(f,a,b)
return

function [I] = trap(f,a,b,n)
h = (b-a)/n; I=0;
for i=0:n
    x = a + i*h;
    if (i == 0 || i == n)  % end points
        I = I + feval(f,x);
    else                     % interior points
        I = I + 2*feval(f,x);
    end
end
I = I*h/2;
return

% simpson integration
function [I] = simpson(f,a,b,n)
h = (b-a)/n; I=0;
if (mod(n,2) == 0)
    for i=0:n
        x = a + i*h;
        if (i == 0 || i == n)  % end points
            I = I + feval(f,x);
        elseif (mod(i,2) == 1)                    % odd points
            I = I + 4*feval(f,x);
        else % even points
            I = I + 2*feval(f,x);
        end
    end
    I = I*h/3;
else
    for i=0:n-3
        x = a + i*h;
        if (i == 0 || i == n-3)  % end points
            I = I + feval(f,x);
        elseif (mod(i,2) == 1)         % odd points
            I = I + 4*feval(f,x);
        else % even points
            I = I + 2*feval(f,x);
        end
    end
     I1 = I*h/3;
     I=0;
    for i =n-3:n
        x = a + i*h;
        if (i == n-3 || i == n)  % end points
            I = I + feval(f,x);
        else % interior points
            I = I + 3*feval(f,x);
        end
    end
    I2 = I*3*h/8;
    I=I1+I2;
end
return

% integrand function
function [v] = vel(t)
g = 9.8; m = 68.1; c = 0.01;
v = sqrt(g*m/c)*tanh(sqrt(g*c/m)*t);
return

% analytical evaluation
function [dist] = distance(T)
% d = (m/c)*log(cosh(sqrt(g*c/m)*t));
syms g m c t v;
v = sqrt(g*m/c)*tanh(sqrt(g*c/m)*t);
d = int(v,t);
g = 9.8; m = 68.1; c = 0.01; t = T;
dist=subs(d);
return