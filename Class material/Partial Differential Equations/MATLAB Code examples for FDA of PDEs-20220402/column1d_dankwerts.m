% 1-d column problem with Dankwerts bc's
% Crank Nicolson method
% bc's: vcin=vc-Ddc/dx at left end, dc/dx=0 at right end
function column1d_dankwerts
Cin=100; v = 1; D = 100;
L=100; N=100; nt=100;
% numerical solution for T = 5 min
T=50;
% FDA approach
[x1, C1] = column_fda(L, N, T, nt, Cin, D, v);
cout = C1(end);
fprintf('Exit concentration FDA = %g moles/m^3\n',cout);
% ODE approach
[x2, C2] = column_ode(L, N, T, Cin, D, v);
plot(x1,C1,'-r',x2,C2,'bo');
legend('FDA','ODE');
xlabel('distance (m)');
ylabel('concentration (mg/m^3)');
title('Concentration profile');
cout = C2(end);
fprintf('Exit concentration ODE = %g moles/m^3\n',cout);
return
%%%%%%%

% soil column problem
% ODE approach
% Dankwertz b.c's
% B.C's: x=0, VC-DdC/dx=VCin; x=L, dc/dx=0;
function [x, C] = column_ode(L, N, T, Cin, D, v)
dx=L/N;
x = (0:dx:L); 
n=N+1;
cinit = zeros(1,n);
odefun = @(t,c) ode_func(t,c,Cin,D,v,dx,n);
tspan = [0 T];
[~,c] = ode45(odefun,tspan,cinit);
C = c(end,:)'; %extract last time step results
%%%%%%%
return

function [dcdt] = ode_func(t,c,Cin,D,v,dx,n)
dcdt = zeros(n,1);
c0=c(2)+2*v*dx*(Cin-c(1))/D;
dcdt(1) = D*(c(2) - 2*c(1)+c0)/dx^2 - v*(c(2) - c0)/(2*dx);
for i=2:n-1
dcdt(i) = D*(c(i+1) - 2*c(i)+c(i-1))/dx^2 - v*(c(i+1) - c(i-1))/(2*dx);
end
cn_plus_one=c(n-1);
dcdt(n) = D*(cn_plus_one - 2*c(n)+c(n-1))/dx^2 - v*(cn_plus_one - c(n-1))/(2*dx); 
return


function [x, c] = column_fda(L, N, T, nt, cin, D, v)
dx=L/N; dt=T/nt;
x = (0:dx:L)';
t = (0:dt:T)';
Pe=v*dx/D;
CN=v*dt/dx;
DN=D*dt/dx^2;
if (Pe > 2)
    fprintf('Warning: grid peclet number is greater than 2!\n');
    fprintf('Results may be inaccurate!\n');
    m = L/(2*D/v);
    fprintf('Increase n to be greater than %g\n',m);
end
n=length(x);
u=ones(n,1);
c=zeros(n,1);
b=zeros(n,1);
for k=1:nt
    % set diagonals
    d = -0.5*(DN+CN/2)*u;
    e = (1+DN)*u;
    f = -0.5*(DN-CN/2)*u;
    % modify for left bc
    c0=c(2)+2*v*dx*(cin-c(1))/D;
    b(1)=0.5*(DN+CN/2)*c0+(1-DN)*c(1)+0.5*(DN-CN/2)*c(2);
    e(1)=e(1)-d(1)*2*v*dx/D;
    f(1)=f(1)+d(1);
    b(1)=b(1)-d(1)*2*v*dx*cin/D;
    for i=2:n-1
        b(i)=0.5*(DN+CN/2)*c(i-1)+(1-DN)*c(i)+0.5*(DN-CN/2)*c(i+1);
    end
    % modify for right bc
    d(n)=d(n)+f(n);
    c_n_plus_one=c(n-1);
    b(n)=0.5*(DN+CN/2)*c(n-1)+(1-DN)*c(n)+0.5*(DN-CN/2)*c_n_plus_one;
    % construct matrix
    A = my_spdiags([d e f],-1:1);
    % solve
    c = A\b;
end
return

% modification of spdiags to adopt linpack convention
% lower diagonals taken from bottom
% upper diagonals taken from top
% i.e., rows of diagonals directly correspond to rows of matrix
% matrix is assumed to be square
% input A = a tall nm matrix containing the 'm' diagonals
% input id = a 1xm row vector containing the 'm' offsets for each diagonal
% output A = nn sparse matrix
function [A] = my_spdiags(A,id)
[n,m] = size(A);
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

