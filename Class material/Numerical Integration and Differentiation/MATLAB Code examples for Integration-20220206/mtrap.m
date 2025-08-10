% multi-dimensional trapezoidal main calling function
function [I] = mtrap(f,a,b,n)
m=length(a); x=zeros(1,m); k=1;
I=ktrap(f,x,a,b,n,k);
return

% recursive trapezoidal function
function [I] = ktrap(f,x,a,b,n,k)
I=0; m=length(a);
h=(b(k)-a(k))/n;
for i=0:n
    x(k)=a(k)+i*h;
    if (k == m)
        if (i == 0 || i == n)
            I=I+feval(f,x);
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
return
