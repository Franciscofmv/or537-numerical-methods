function [I] = romberg(f,a,b,N)
L=round(sqrt(N)); % set number of levels
IR=zeros(L,L); % initialize integral matrix
% loop over number of levels to calculate
% integral estimates at first level using trapezoidal rule
for i=1:L 
    n=N*2^(i-1);
    IR(i,1) = gtrap(f,a,b,n);
end
% calculate integral using Richardson's extrapolation
for j = 2:L
    for i = 1:L-j+1
        IR(i,j) = (4^(j-1)*IR(i+1,j-1) - IR(i,j-1))/(4^(j-1) - 1);
    end
end
I = IR(1,L); % final value in first row last column
end
