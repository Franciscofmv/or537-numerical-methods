% Palm 2.47 (p. 115)
function traffic
A=[ 1  0  1  0  0  0  0;
    1 -1  0 -1  0  0  0;
    0  1  0  0  1  0  0;
    0  0  1  0  0  1  0;
    0  0  0  1  0 -1  1;
    0  0  0  0  1  0  1;
    0  0  0  0  0  0  0];
n=size(A,1)
b = [300; -300; 600; 400; 200; 600; 0];
f=zeros(n,1);
% rank of A
rankA=rank(A)
if (rankA == n)
    % solve linear system
    f = A\b
else
    % construct augmented matrix AA
    AA=[A b];
    % reduce to row echelon form to get coefficients
    B=rref(AA)
    % measure f(6) and f(7) - set arbitrary values
    f(6)=300; f(7)=200;
    f(1:5) = B(1:5,8)-B(1:5,6)*f(6)-B(1:5,7)*f(7);
    f
end
% end example
end
