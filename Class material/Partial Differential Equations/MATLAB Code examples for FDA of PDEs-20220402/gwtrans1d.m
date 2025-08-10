% soil column problem
% solves dc/dt=D*d^2c/dx^2-v*dc/dx
% bc's: c=c0 at left end, dc/dx=0 at right end
% solution using crank-nicolson method
function gwtrans1d
v=0.36; alpha=2; c0=2;
L=100; n=100;
% case I:  T = 28 hrs
T=28; nt=28;
% numerical solution
[xn, cn] = column_fda(L, T, n, nt, v, alpha, c0);
% analytical solution
xa=[0:L]; D = alpha*v;
f1=(xa-v*T)/(2*sqrt(D*T));
f2=(xa+v*T)/(2*sqrt(D*T));
ca=0.5*c0*(erfc(f1)+exp(v*xa/D).*erfc(f2));
figure(1);
plot(xa,ca,'b', xn, cn, 'ro');
legend('numerical','analytical');
title('Concentration profile at T = 28 hrs');
hold off;
% case II: T = 167 hrs
T=167; nt=167;
% numerical solution
[xn, cn] = column_fda(L, T, n, nt, v, alpha, c0);
% analytical solution
xa=[0:L];
f1=(xa-v*T)/(2*sqrt(D*T));
f2=(xa+v*T)/(2*sqrt(D*T));
ca=0.5*c0*(erfc(f1)+exp(v*xa/D).*erfc(f2));
figure(2);
plot(xa,ca,'b', xn, cn, 'ro');
legend('numerical','analytical');
title('Concentration profile at T = 167 hrs');
hold off;
return
%%%%%%%


function [x, c] = column_fda(L, T, n, nt, v, alpha, c0)
dt=T/nt;
D=alpha*v;
dx=L/n;
x = (0:dx:L)';
nx=length(x); % number of segments
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
% initialize matrix diagonals
u=ones(nx,1);
d = -(DN+CN/2)*0.5*u;
e = (1+DN)*u;
f = (-DN+CN/2)*0.5*u;
for k=1:nt
    % first node is a dirichlet b.c with c=c0;
    c(1)=c0; d(1)=0; e(1)=1; f(1)=0; b(1)=c0;
    % last node is flux b.c with dc/dx=0;
    c(nx)=0; d(nx)=DN; e(nx)=1-DN; f(nx)=0;
    b(nx)=c(nx-1)*DN+(1-DN)*c(nx);
    % set rhs for interior nodes
    for i=2:nx-1
        b(i)=c(i-1)*(DN+CN/2)*0.5+(1-DN)*c(i)+(DN-CN/2)*0.5*c(i+1);
    end
    A = my_spdiags([d e f],-1:1);
    c = A\b;
end
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
