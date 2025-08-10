function [x] = gauss_jordan(A,b)
n=size(A,1); % get size of matrix A
AA=[A b]; % add rhs to get augmented matrix AA
% loop over rows
for i=1:n 
    % if diagonal entry of row i is zero then swap it
    % with the row containing the largest entry in column i
    if (AA(i,i) == 0) 
        AA=pivot_row(AA,i); % swapping operation
    end
    % normalize by diagonal
    AA(i,:)=AA(i,:)/AA(i,i);
    % subtract a multiple of row j not equal to i
    % to zero out all entries in column i except diagonal entry
    for j=1:n
        if (j ~= i)
            AA(j,:)=AA(j,:)-AA(j,i)*AA(i,:);
        end
    end
end
x=AA(:,n+1);
return

% pivot function
function [A] = pivot_row(A,i)
% store row i in a temporary vector
row_temp=A(i,:);
% get the row index j containing the largest entry in column i
[~, j]=max(abs(A(:,i)));
% swap row j with i
A(i,:)=A(j,:);
A(j,:)=row_temp;
return
