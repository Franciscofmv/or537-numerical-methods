% modification of spdiags to adopt linpack convention
% lower diagonals taken from bottom
% upper diagonals taken from top
% i.e., rows of diagonals directly correspond to rows of matrix
% matrix is assumed to be square
% input A = a tall nxm matrix containing the 'm' diagonals
% input id = a 1xm row vector containing the 'm' offsets for each diagonal
% output A = nxn sparse matrix
function [A] = my_spdiags(A,id)
[n m] = size(A);
B = zeros(n,m);
for i=1:m
    if (id(i) <= 0)
        B(1:n+id(i),i) = A(-id(i)+1:n,i);
    else
        B(1+id(i):n,i) = A(1:n-id(i),i);
    end
end
% construct matrix
A = spdiags(B,id,n,n);
return