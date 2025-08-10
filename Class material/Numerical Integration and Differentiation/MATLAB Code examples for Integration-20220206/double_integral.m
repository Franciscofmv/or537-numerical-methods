% simple double integration
function double_integral
ax=0; bx=3;
ay=-1; by=1;
n=10; f=@xexy;
I_trap = dbl_trap(f,ax,bx,ay,by,n)
I_simp = dbl_simp(f,ax,bx,ay,by,n)
I_matlab = integral2(f,ax,bx,ay,by)
syms x y;
I_anal = double(int(int(x*exp(x*y),x,0,3),y,-1,1))
return

function [f] = xexy(x,y)
f = x.*exp(x.*y);
return


function [I] = dbl_trap(f,ax,bx,ay,by,n)
I=0;
h=(by-ay)/n;
for i=0:n
    y=ay+i*h;
    if (i == 0 || i == n)
        I=I+trap(f,y,ax,bx,n);
    else
        I=I+2*trap(f,y,ax,bx,n);
    end
end
I=I*h/2;
return

function [I] = trap(f,y,a,b,n)
I=0;
h=(b-a)/n;
for i=0:n
    x=a+i*h;
    if (i == 0 || i == n)
        I=I+f(x,y);
    else
        I=I+2*f(x,y);
    end
end
I=I*h/2;
return


function [I] = dbl_simp(f,ax,bx,ay,by,n)
I=0;
h=(by-ay)/n;
for i=0:n
    y=ay+i*h;
    if (i == 0 || i == n)
        I=I+simp(f,y,ax,bx,n);
    elseif (mod(i,2) == 1)
        I=I+4*simp(f,y,ax,bx,n);
    elseif (mod(i,2) == 0)
        I=I+2*simp(f,y,ax,bx,n);
    end
end
I=I*h/3;
return

function [I] = simp(f,y,a,b,n)
I=0;
h=(b-a)/n;
for i=0:n
    x=a+i*h;
    if (i == 0 || i == n)
        I=I+f(x,y);
    elseif (mod(i,2) == 1)
        I=I+4*f(x,y);
    elseif (mod(i,2) == 0)
        I=I+2*f(x,y);
    end
end
I=I*h/3;
return
