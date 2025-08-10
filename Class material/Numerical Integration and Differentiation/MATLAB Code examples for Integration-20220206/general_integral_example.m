% m-dimensional integration
function general_integral_example
a=char('0','0','0');
b=char('8','xd2','xm2y');
n = [10 10 10];
f=@exyz;
I_trap = trap_g(f,a,b,n)
I_simp = simp_g(f,a,b,n)
clear
syms x y z;
I_anal=double(int(int(int(exp(-(x^2+y^2+z^2)/4),z,0,x-2*y),y,0,x/2),x,0,8))
return

function [f] = exyz(x)
% x(1) = x, x(2) = y, ...
f = exp(-(x(1)^2+x(2)^2+x(3)^2)/4);
return

function [f] = xm2y(x)
f = x(1)-2*x(2);
return

function [f] = xd2(x)
f = x(1)/2;
return


function [I] = trap_g(f,a,b,n)
m = length(n); x(1:m)=0; k=1;
I = trap_k(f,x,a,b,n,m,k);
return

function [I] = trap_k(f,x,a,b,n,m,k)
ak=str2num(a(k,:));
bk=str2num(b(k,:));
if (isempty(ak))
    fnc=deblank(a(k,:));
    ak=feval(fnc,x);
end
if (isempty(bk))
    fnc=deblank(b(k,:));
    bk=feval(fnc,x);
end
h(k)=(bk-ak)/n(k); 
I=0;
if (k == m)
    for i=0:n(k)
        x(k)=ak+i*h(k);
        if (i == 0 || i == n(k))
            I=I+feval(f,x);
        else
            I=I+2*feval(f,x);
        end
    end
else
    for i=0:n(k)
        x(k)=ak+i*h(k);
        if (i == 0 || i == n(k))
            I=I+feval(@trap_k,f,x,a,b,n,m,k+1);
        else
            I=I+2*feval(@trap_k,f,x,a,b,n,m,k+1);
        end
    end
end
I=I*h(k)/2;
return

function [I] = simp_g(f,a,b,n)
m = length(a); x(1:m)=0; k=1;
I = simp_k(f,x,a,b,n,m,k);
return

function [I] = simp_k(f,x,a,b,n,m,k)
ak=str2num(a(k,:));
bk=str2num(b(k,:));
if (isempty(ak))
    fnc=deblank(a(k,:));
    ak=feval(fnc,x);
end
if (isempty(bk))
    fnc=deblank(b(k,:));
    bk=feval(fnc,x);
end
h(k)=(bk-ak)/n(k); 
I=0;
if (k == m)
    for i=0:n(k)
        x(k)=ak+i*h(k);
        if (i == 0 || i == n(k))
            I=I+feval(f,x);
        elseif (mod(i,2) == 1)
            I=I+4*feval(f,x);
        else
            I=I+2*feval(f,x);
        end
    end
else
    for i=0:n(k)
        x(k)=ak+i*h(k);
        if (i == 0 || i == n(k))
            I=I+feval(@simp_k,f,x,a,b,n,m,k+1);
        elseif (mod(i,2) == 1)
            I=I+4*feval(@simp_k,f,x,a,b,n,m,k+1);
        else
            I=I+2*feval(@simp_k,f,x,a,b,n,m,k+1);
        end
    end
end
I=I*h(k)/3;
return

