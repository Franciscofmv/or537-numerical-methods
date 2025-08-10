% soil column problem - FDA crank nicolson implementation
% solves dc/dt=D*d^2c/dx^2-v*dc/dx
% dirichlet bc's: x=0, c=c0; x=L, c=0
% D=alpha*v
function simple_column_crank_nicolson
v=0.36; alpha=5; c0=2;
L=100; n=100; T=50; nt=200;
dt=T/nt;
D=alpha*v;
dx=L/n;
x = (0:dx:L)';
nx=length(x);
b = zeros(nx,1);
CN=v*dt/dx;
DN=D*dt/dx^2;
c =zeros(nx,1);
P=v*dx/D;
if (P > 2)
    fprintf('Grid peclet number is greater than 2!\n');
    fprintf('Grid peclet number is %g\n',P);
    fprintf('Reduce grid spacing\n');
end
for k=1:nt % foor loop
    % initialize matrix diagonals
    u=ones(nx,1); % ones vector
    % each diagonal is a vector:
    d = -(DN+CN/2)*0.5*u;% lower diagonal
    e = (1+DN)*u; % middle diagonal/main diagonal
    f = (-DN+CN/2)*0.5*u; % upper diagonal
    % first node is a dirichlet b.c with c=c0;
    c(1)=c0; d(1)=0; e(1)=1; f(1)=0; b(1)=c0; % for conditions
    % last node is a dirichlet b.c with c=0;
    c(nx)=0; d(nx)=0; e(nx)=1; f(nx)=0; b(nx)=0;
    % set rhs for interior nodes
    for i=2:nx-1
        b(i)=c(i-1)*(DN+CN/2)*0.5+(1-DN)*c(i)+(DN-CN/2)*0.5*c(i+1);
    end
    A = my_spdiags([d e f],-1:1); % modify sp_diags
    c = A\b;
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
return
%%%%%%%

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

