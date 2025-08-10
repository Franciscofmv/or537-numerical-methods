% chapra 19.16 (p. 494)
function mass_in_reactor
T=[0 10 20 30 35 40 45 50];
Q=[4 4.8 5.2 5.0 4.6 4.3 4.3 5.0];
c=[10 35 55 52 40 37 32 34];
mass_data = trap_data(T,Q,c)
f=@(t) mass_func(t,T,Q,c);
a=0; b=50; n=100;
mass_trap = trap(f,a,b,n)
mass_simpson = simpson(f,a,b,n)
mass_quad = quad(f,a,b)
return

function [I] = trap_data(t,Q,c)
I=0; n=length(t);
for i=1:n-1
    dt=t(i+1)-t(i);
    Qt=(Q(i)+Q(i+1))/2;
    ct=(c(i)+c(i+1))/2;
    I = I + Qt*ct*dt;
end
return

function [m] = mass_func(t,T,Q,c)
qt=interp1(T,Q,t);
ct=interp1(T,c,t);
m=qt.*ct;
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
return


    