function linear_solver_example
% construct FDA matrix for a 2-D groundwater transport problem
[A,b] = gwflow2d_matrix;
n=size(b);
figure(1);
spy(A);
if (sum(diag(A)) < 0)
    A=-A;
    b=-b;
end
% solve using various methods
x=A\b;
% x=bicgstab(A,b);
% x=gmres(A,b);
% x=pcg(A,b);
ij=0;
nx=11; ny=7;
z=zeros(nx,ny);
for i=1:nx
    for j=1:ny
        ij=ij+1;
        z(i,j)=x(ij);
    end
end
figure(2);
[C, h]=contourf(z);
clabel(C,h);
[vx,vy] = gradient(z);
hold on;
h=streamslice(-vx,-vy);
set(h,'color','k');
axis tight;
hold off;
return

