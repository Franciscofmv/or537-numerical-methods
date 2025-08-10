% multi-dimensional simpson main calling function
function [I] = msimp(f,a,b,n)
m=length(a); x=zeros(1,m); k=1;
I=ksimp(f,x,a,b,n,k);
return

% recursive simpson function
function [I] = ksimp(f,x,a,b,n,k)
I=0; m=length(a);
h=(b(k)-a(k))/n;
for i=0:n
    x(k)=a(k)+i*h;
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