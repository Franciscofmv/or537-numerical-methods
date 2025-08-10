% cooling fin: linear two-point boundary value problem
% this is a singular problem
% chapra and canale 28_48
function cooling_fin
% for bvp4c singular option does not work
% so start x from 0.01 and assume u0 = 0.01 in @coolbc
xinit= 0.01:0.01:1;
yinit=[0 0];
% matlab bvp4c solution
solinit = bvpinit(xinit,yinit);
p = 10;
% options = bvpset('SingularTerm',[0 0; 2*p -2]);
% sol = bvp4c(@cool, @coolbc, solinit, options);
% options = bvpset('SingularTerm',[0 0; 2*p -2]);
sol = bvp4c(@coolode, @coolbc, solinit);
X = 0.1:0.05:1;
Y = deval(sol,X);
% FDA solution
x = 0:0.01:1;
u0=0;
u1=1;
u=fda_cool(x,u0,u1); % fda
plot(X,Y(1,:),'-ro',x,u,'-g');
legend('bvp4c solution','FDA solution');
xlabel('x');
ylabel('u');
% syms y(x);
% s=dsolve(diff(y,2) == (2/x)*y-(20/x)*diff(y), y(0) == 0, y(1) == 1);
% x = 0:0.01:1;
% y = subs(s);
hold off;
return

function [u] = fda_cool(X,u0,u1)
p=10; dx=X(2)-X(1);
m=length(X);
n=m-2;
x = X(2:m-1)';
e = ones(n,1);
ld=x.*e-dx;
ud=x.*e+dx;
d=zeros(n,1);
d(1:n)=-2*x-2*p*dx^2;
b=zeros(n,1);
b(1)=-u0*(x(1)-dx);
b(n)=-u1*(x(n)+dx);
A=my_spdiags([ld d ud],[-1,0,1]);
y=A\b;
u(2:m-1)=y;
u(1)=u0;
u(m)=u1;
return


function [res] = coolbc(ya,yb)
res = [ya(1) - 0.01; yb(1)-1];
return

function [dydx] = cool(x,y)
dydx=zeros(2,1);
dydx(1) = y(2);
dydx(2) = 0;
return


function [y] = bvp_shooting(odefun,x,yi,yf)
n=length(x);
x0=x(1);
xL=x(n);
y0=zeros(2,1);
% first shot
q1=0;
[X,Y]=ode45(odefun,x,[yi; q1]);
p1=Y(n,1);
% second shot
q2=1; % assume
[X,Y]=ode45(odefun,x,[yi; q2]);
p2=Y(n,1);
y2=q1+(yf-p1)*(q2-q1)/(p2-p1);
[X,Y]=ode45(odefun,x,[yi; y2]);
y = Y';
return

function [dydx] = coolode(x,y)
dydx=zeros(2,1);
p = 10;
dydx(1) = y(2);
dydx(2) = -2*(y(2)-p*y(1))/x;
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

