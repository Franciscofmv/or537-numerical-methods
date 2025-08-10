%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate average mass in tank over time T 
% with variable height and time dependendent concentration
% uses 'msimp' and 'gsimp' functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tank_example4
f=@tank_conc; T=3600;
a={-2,-2,0,0};
b={2,2,T,@tank_height};
n = 10;
MT=gsimp(f,a,b,n);
fprintf('Average mass in the tank through 1 hour is %g  Kg\n',MT/T);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions to be evaluated follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f]=tank_conc(X)
C0=0.1; k=0.001; x=X(1); y=X(2); t=X(3); z=X(4);
f=C0*exp(-k*t)*(x^2+y^2)/sqrt(z^2+1);
return

function [h] = tank_height(X)
A=pi*4; T = X(3); n = 100;
h=msimp(@qfunc,0,T,n)/A;
return

function [q] = qfunc(t)
Q0=0.1; w=pi/6;
q=Q0*(1+sin(w*t));
return
% multi-dimensional simpson main calling function
% goes from outer to inner (i.e., x to z).
function [I] = gsimp(f,a,b,n)
m=length(a); x=zeros(1,m); k=1;
% convert limits to cell array if numbers are given
if (isnumeric(a))
    a=num2cell(a);
end
if (isnumeric(b))
    b=num2cell(b);
end
% call recursive integration function
I=ksimp(f,x,a,b,n,k);
return

% recursive simpson function
function [I] = ksimp(f,x,a,b,n,k)
I=0; m=length(a);
if isnumeric(a{k})
    ak=a{k};
else
    ak=feval(a{k},x);
end
if isnumeric(b{k})
    bk=b{k};
else
    bk=feval(b{k},x);
end
h =(bk-ak)/n;
for i=0:n
    x(k)=ak+i*h;
    if (k == m)
        if (i == 0 || i == n)
            I=I+feval(f,x);
        elseif (mod(i,2) == 1)
            I=I+4*feval(f,x);
        elseif (mod(i,2) == 0)
            I=I+2*feval(f,x);
        end
    else
        if (i == 0 || i == n)
            I=I+feval(@ksimp,f,x,a,b,n,k+1);
        elseif (mod(i,2) == 1)
            I=I+4*feval(@ksimp,f,x,a,b,n,k+1);
        elseif (mod(i,2) == 0)
            I=I+2*feval(@ksimp,f,x,a,b,n,k+1);
        end
    end
end
I=I*h/3;
return

