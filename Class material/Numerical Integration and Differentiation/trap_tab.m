function [I] = trap_tab(x,y) %x y are vectors, tabular data
n = length(x);
I = 0;
    for i = 1:n-1 % Tabular integration function
        h = x(i+1)-x(i);
        I = I+h*(y(i+1)+y(i))/2;
    end
end