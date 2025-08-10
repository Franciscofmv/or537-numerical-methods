function exercise3
format long
I_builtin=integral(@integrand,0,30)
I_numerical=segment_trapezoidal(@integrand,0,30,6)-10931.25
end

function [v] = integrand(t)
g = 9.8;%downward acceleration of gravity
q = 2500; % fuel compsumption kg
u = 1800; %velocity at which fuel is expelled
m0 = 160000; % initial mass of the rocket at time t
v = u.*log(m0./(m0-q*t))-g*t;
end

function [I] = segment_trapezoidal(f,a,b,n)
I=0; % initial value of integral
h=(b-a)/n;
for i=0:n
    x= a + i*h; % goin from a to b
    if (i==0 || i==n)
        I = I + f(x);
    else % if is not 0 or n
        I = I + 2*f(x);
    end
end
I = I*h/2;
end
