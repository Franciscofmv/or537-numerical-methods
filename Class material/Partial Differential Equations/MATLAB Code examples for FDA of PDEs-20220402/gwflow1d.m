% solves steady state 1d gwflow equation
% of the form d/dx(Khdh/dx)+N=0
% K and N are assumed constant
function gwflow1d
K = 1; N=1e-5; dx=10; L=1000;
x=(0:dx:L)';
h0=10; hL=0; % dirichlet conditions at both ends
m=length(x);
h=zeros(m,1);
n=m-2;
ld=ones(n,1);
ud=ones(n,1);
d=-2*ones(n,1);
b=-2*N*ones(n,1)*dx^2/K;
b(1)=b(1)-h0;
b(n)=b(n)-hL;
A=spdiags([ld d ud],[-1,0,1],n,n);
y=A\b;
h(2:m-1)=sqrt(y);
h(1)=sqrt(h0);
h(m)=sqrt(hL);
plot(x,h,'-r');
return