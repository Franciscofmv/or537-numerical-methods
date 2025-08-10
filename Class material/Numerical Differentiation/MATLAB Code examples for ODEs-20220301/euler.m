function [y] = euler(f,x,y0)
n=length(x);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    xi=x(i);
    yi=y(:,i);
    h=x(i+1)-x(i);
    k1=f(xi,yi);
    y(:,i+1)=yi+h*k1;
end
end