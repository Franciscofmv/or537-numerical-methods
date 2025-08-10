function gauss_jordan_test
A=rand(5);
A(1,1)=0;
b=rand(5,1);
x=A\b
x=gauss_jordan(A,b)
return


function [x] = gauss_jordan(A,b)
n=size(A,1);
AA=[A b];
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
x=AA(:,n+1);
return

function [A] = pivot_row(A,i)
[value,index]=max(abs(A(:,i)));
row_temp=A(i,:);
A(i,:)=A(index,:);
A(index,:)=row_temp;
return
