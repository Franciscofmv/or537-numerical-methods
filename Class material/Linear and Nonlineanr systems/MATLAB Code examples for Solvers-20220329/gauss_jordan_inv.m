function [A_inv] = gauss_jordan_inv(A)
n=size(A,1);
I=eye(n);
AA=[A I];
for i=1:n
    if (AA(i,i) == 0)
        AA=pivot_row(AA,i);
    end
    AA(i,:)=AA(i,:)/AA(i,i);
    for j=1:n
        if (j ~= i)
            AA(j,:)=AA(j,:)-AA(j,i)*AA(i,:);
        end
    end
end
A_inv=AA(:,n+1:2*n);
return

function [A] = pivot_row(A,i)
[~,j]=max(abs(A(:,i)));
row_temp=A(i,:);
A(i,:)=A(j,:);
A(j,:)=row_temp;
return
