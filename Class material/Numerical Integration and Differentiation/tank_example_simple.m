function tank_example_simple()
f = @conc; % concentration
a = [-1, -1, 0]; % lower limits
b = [1 1 3]; % upper limit
n = 10; % number of  segments
mass_trap = mtrap(f,a,b,n) % result of integral
mass_sim = msimp(f,a,b,n)

f2 = @concxyz;
mass_matlab  = integral3(f2,a(1),b(1),a(2),b(2),a(3),b(3)) % it doesn't take conc(X) with just one argumetn
syms f(x,y,z); % define symbolic variables
f(x,y,z) = 0.1*(x^2+y^2)/sqrt(z^2 + 1);
mass_ana = double(int(int(int(f,-1,1),-1,1),0,3))

end

function  [c] = concxyz(x,y,z) % for the matlab built in integral3
C0 = 0.1; % concentration in kg/m^3,
c = C0*(x.^2+y.^2)./sqrt(z.^2 + 1);

end

function [c] = conc(X); % X is a vector, X = [x,y,z]
x = X(1);y = X(2); z = X(3);
C0 = 0.1; % concentration in kg/m^3,
c = C0*(x^2+y^2)/sqrt(z^2 + 1);
end


% multi-dimensional trapezoidal main calling function
function [I] = mtrap(f,a,b,n)
m=length(a); x=zeros(1,m); k=1;
I=ktrap(f,x,a,b,n,k);
end
% recursive trapezoidal function
function [I] = ktrap(f,x,a,b,n,k)
I=0; m=length(a);
h=(b(k)-a(k))/n;
for i=0:n
    x(k)=a(k)+i*h;
    if (k == m)
        if (i == 0 || i == n)
            I=I+feval(f,x); % x is a vector
        else
            I=I+2*feval(f,x);
        end
    else
        if (i == 0 || i == n)
            I=I+feval(@ktrap,f,x,a,b,n,k+1);
        else
            I=I+2*feval(@ktrap,f,x,a,b,n,k+1);
        end
    end
end
I=I*h/2;
end
