function [I] = trap_tabular(x,y) % 

if length(x) == length(y)
    I = 0; % initial value for Integral I
    n = length(x);
    for i = 1:(n-1)
        h = x(i+1)-x(i); % width of segment 
        I = I + h*(y(i)+y(i+1))/2; % integration of tabular data formula    
    end
else
    disp("x and y vectors are of different length. To integrate tabular data, x and y must be of same length.")
end
end