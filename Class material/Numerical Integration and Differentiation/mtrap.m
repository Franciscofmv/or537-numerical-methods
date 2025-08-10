% multi-dimensional trapezoidal main calling function
function [I] = mtrap(f,a,b,n)
m=length(a); x=zeros(1,m); k=1; % if two dimensions, m = 2, three m=3,...
I=ktrap(f,x,a,b,n,k);
return

% recursive trapezoidal function
function [I] = ktrap(f,x,a,b,n,k) % recursive function
I=0; m=length(a);
h=(b(k)-a(k))/n;
for i=0:n
    x(k)=a(k)+i*h;
    if (k == m) % it goes to the innermost ??
        if (i == 0 || i == n)
            I=I+feval(f,x); % we can use f
        else
            I=I+2*feval(f,x);
        end
    else
        if (i == 0 || i == n)
            I=I+feval(@ktrap,f,x,a,b,n,k+1);% k+1 making k bigger
        else
            I=I+2*feval(@ktrap,f,x,a,b,n,k+1);
        end
    end
end
I=I*h/2;
return
