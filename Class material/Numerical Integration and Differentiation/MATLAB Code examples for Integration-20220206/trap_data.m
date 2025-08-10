function [I] = trap_data(x,y)
I = 0; n = length(x);
for i=1:n-1
    h=x(i+1)-x(i);
    I=I+h*(y(i)+y(i+1))/2;
end
end