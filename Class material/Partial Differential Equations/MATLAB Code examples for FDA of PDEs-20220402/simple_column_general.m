% soil column problem - FDA general implementation
% omega=0; explicit
% omega=1; fully implicit
% omega=0.5; crank-nicolson
% solves dc/dt=D*d^2c/dx^2-v*dc/dx
% dirichlet bc's: x=0, c=c0; x=L, c=0
% D=alpha*v
function simple_column_general
v=0.36; alpha=5; c0=2;
L=100; n=100; omega=1; T=50; nt=200;
dt=T/nt;
omemega=1-omega;
D=alpha*v;
dx=L/n;
x = (0:dx:L)';
b = zeros(n+1,1);
CN=v*dt/dx;
DN=D*dt/dx^2;
c =zeros(n+1,1);
P=v*dx/D;
if (P > 2)
    fprintf('Grid peclet number is greater than 2!\n');
    fprintf('Grid peclet number is %g\n',P);
    fprintf('Reduce grid spacing\n');
end
if (omega == 0)
    if (CN >= 1)
        fprintf('Courant number is greater than 1 for explicit scheme!\n');
        fprintf('Courant number is %g\n',CN);
        fprintf('Reduce time step size\n');
    end
    if (DN >= 0.5)
        fprintf('Diffusion number is greater than 0.5 for explicit scheme!\n');
        fprintf('Diffusion number is %g\n',DN);
        fprintf('Reduce time step size\n');
    end
end
for k=1:nt
    u=ones(n+1,1);
    d = -(DN+CN/2)*omega*u;
    e = (1+2*DN*omega)*u;
    f = (-DN+CN/2)*omega*u;
    % first node is a dirichlet b.c with c=c0;
    c(1)=c0; d(1)=0; e(1)=1; f(1)=0; b(1)=c0;
    % last node is a dirichlet b.c with c=0;
    c(n+1)=0; d(n+1)=0; e(n+1)=1; f(n+1)=0; b(n+1)=0;
    for i=2:n
        b(i)=c(i-1)*(DN+CN/2)*omemega+(1-2*DN*omemega)*c(i)+(DN-CN/2)*omemega*c(i+1);
    end
    if (omega ~=0) % implicit
        A = my_spdiags([d e f],-1:1);
        c = A\b;
    else       % explicit
        c=b;
    end
end
plot(x,c,'b');
hold on;
%%%%%%%%%%%%%%%
% analytical solution
f1=(x-v*T)/(2*sqrt(D*T));
f2=(x+v*T)/(2*sqrt(D*T));
c=0.5*c0*(erfc(f1)+exp(v*x/D).*erfc(f2));
plot(x,c,'--r');
legend('numerical','analytical');
hold off;
%%%%%%%
return

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

