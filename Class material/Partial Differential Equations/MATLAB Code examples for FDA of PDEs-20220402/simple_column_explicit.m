% soil column problem - FDA explicit implementation
% solves dc/dt=D*d^2c/dx^2-v*dc/dx
% dirichlet bc's: x=0, c=c0; x=L, c=0
% D=alpha*v
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
% use notation c^k+1=c and c^k=cold
for k=1:nt
    c(1)=c0;
    c(nx)=0;
    cold=c;
    for i=2:nx-1
        c(i)=cold(i)+cold(i-1)*(DN+CN/2)+(-2*DN)*cold(i)+(DN-CN/2)*cold(i+1);
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

