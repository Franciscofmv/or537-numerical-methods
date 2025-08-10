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