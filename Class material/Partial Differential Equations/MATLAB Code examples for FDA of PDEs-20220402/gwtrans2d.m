function [] = gwtrans2d
alphaL=2; alphaT=0.5; Dm=0.01; vx=1;, vy=0;
Lx=100; Ly=60; T=50; nx=51; ny=31; nt=25;
R=1; lambda=0; isol = [1 nt];
[csol] = gwtrans2d_solve(Lx, Ly, T, nx, ny, nt, isol, vx, vy, ...
    alphaL, alphaT, Dm, R, lambda);
ij=0;
for i=1:nx
    for j=1:ny
        ij=ij+1;
        c(i,j)=csol(ij,2);
    end
end
[C, h] = contour(c');
clabel(C,h);
colormap cool;
return
%%%%%%%

function [csol] = gwtrans2d_solve(Lx, Ly, T, nx, ny, nt, isol, vx, vy, ...
    alphaL, alphaT, Dm, R, lambda)
% hardcode initial condition
x1=14; x2=30; y1=16; y2=28; c0=70;
t=zeros(nt,1);
V=sqrt(vx^2+vy^2);
Dxx=alphaL*V+Dm;
Dyy=alphaT*V+Dm;
dt=T/nt;
dx=Lx/(nx-1);
dy=Ly/(ny-1);
n = nx*ny;
b = zeros(n,1);
c = zeros(n,1);
nsol = length(isol);
m = 1;
csol = zeros(n,nsol);
% intialize source
c=initialize(x1,y1,x2,y2,dx,dy,nx,ny,c0);
% Convection or Courant number
CNx=vx*dt/dx;
CNy=vy*dt/dy;
CN=max(CNx,CNy);
% Diffusion number
DNx=Dxx*dt/dx^2;
DNy=Dyy*dt/dy^2;
DN=max(DNx,DNy);
% grid peclet number
P=max(vx*dx/Dxx,vy*dy/Dyy);
if (P > 2)
    fprintf('Grid peclet number is greater than 2!\n');
    fprintf('Grid peclet number is %g\n',P);
    fprintf('Reduce grid spacing\n');
end
u=ones(n,1);
for k=1:nt
    t(k)=k*dt;
    % internal nodes
    for i=2:nx-1
        for j=2:ny-1
            ij1=feval(@ind,ny,i-1,j);
            ij2=feval(@ind,ny,i,j-1);
            ij3=feval(@ind,ny,i,j);
            ij4=feval(@ind,ny,i,j+1);
            ij5=feval(@ind,ny,i+1,j);
            % rhs
            b(ij3)=0.5*(DNx+0.5*CNx)*c(ij1)+...
                 0.5*(DNy+0.5*CNy)*c(ij2)+...
                (R-DNx-DNy-lambda*0.5)*c(ij3)+...
                 0.5*(DNy-0.5*CNy)*c(ij4)+...              
                 0.5*(DNx-0.5*CNx)*c(ij5);
            
        end
    end
    % matrix diagonals
    if (k == 1)
        d = -0.5*(DNx+0.5*CNx)*u;
        e = -0.5*(DNy+0.5*CNy)*u;
        f = (R+DNx+DNy+lambda*0.5)*u;
        g = -0.5*(DNy-0.5*CNy)*u;
        h = -0.5*(DNx-0.5*CNx)*u;
    end
    % special treatment for boundary nodes
    % left end (x = 0)
    i=1;
    for j=2:ny-1
         ij1=feval(@ind,ny,i-1,j);
         ij2=feval(@ind,ny,i,j-1);
         ij3=feval(@ind,ny,i,j);
         ij4=feval(@ind,ny,i,j+1);
         ij5=feval(@ind,ny,i+1,j);   
         b(ij3)= 0.5*(DNy+0.5*CNy)*c(ij2)+...
                (R-DNx-DNy-lambda*0.5)*c(ij3)+...
                 0.5*(DNy-0.5*CNy)*c(ij4)+...              
                DNx*c(ij5);
         d(ij3)=0;
         h(ij3)=-DNx;
    end
    % right end (x = Lx)
    i=nx;
    for j=2:ny-1
         ij1=feval(@ind,ny,i-1,j);
         ij2=feval(@ind,ny,i,j-1);
         ij3=feval(@ind,ny,i,j);
         ij4=feval(@ind,ny,i,j+1);
         ij5=feval(@ind,ny,i+1,j);   
         b(ij3)= DNx*c(ij1)+...
                 0.5*(DNy+0.5*CNy)*c(ij2)+...
                (R-DNx-DNy-lambda*0.5)*c(ij3)+...
                 0.5*(DNy-0.5*CNy)*c(ij4);
         d(ij3)=-DNx;
         h(ij3)=0;
    end
    % bottom end (y = 0)
    j=1;
    for i=2:nx-1
         ij1=feval(@ind,ny,i-1,j);
         ij2=feval(@ind,ny,i,j-1);
         ij3=feval(@ind,ny,i,j);
         ij4=feval(@ind,ny,i,j+1);
         ij5=feval(@ind,ny,i+1,j);   
         b(ij3)= 0.5*(DNx+0.5*CNx)*c(ij1)+...
                (R-DNx-DNy-lambda*0.5)*c(ij3)+...
                 DNy*c(ij4)+...
                 0.5*(DNx-0.5*CNx)*c(ij5);
         e(ij3)=0;
         g(ij3)=-DNy;
    end
    % top end (y = Ly)
    j=ny;
    for i=2:nx-1
         ij1=feval(@ind,ny,i-1,j);
         ij2=feval(@ind,ny,i,j-1);
         ij3=feval(@ind,ny,i,j);
         ij4=feval(@ind,ny,i,j+1);
         ij5=feval(@ind,ny,i+1,j);   
         b(ij3)=0.5*(DNx+0.5*CNx)*c(ij1)+...
                 DNy*c(ij2)+...
                (R-DNx-DNy-lambda*0.5)*c(ij3)+...              
                 0.5*(DNx-0.5*CNx)*c(ij5);
         e(ij3)=-DNy;
         g(ij3)=0;
    end
    % bottom left corner (x=0, y=0)
    i=1;j=1;
    ij1=feval(@ind,ny,i-1,j);
    ij2=feval(@ind,ny,i,j-1);
    ij3=feval(@ind,ny,i,j);
    ij4=feval(@ind,ny,i,j+1);
    ij5=feval(@ind,ny,i+1,j);
    b(ij3)=(R-DNx-DNy-lambda*0.5)*c(ij3)+...
                 DNy*c(ij4)+...              
                 DNx*c(ij5);
    if (k == 1)
        d(ij3)=0;
        e(ij3)=0;
        g(ij3)=-DNy;
        h(ij3)=-DNx;
	end
    % bottom right corner (x=Lx, y=0)
    i=nx;j=1;
    ij1=feval(@ind,ny,i-1,j);
    ij2=feval(@ind,ny,i,j-1);
    ij3=feval(@ind,ny,i,j);
    ij4=feval(@ind,ny,i,j+1);
    ij5=feval(@ind,ny,i+1,j);
    b(ij3)=DNx*c(ij1)+...
           (R-DNx-DNy-lambda*0.5)*c(ij3)+...
           DNy*c(ij4);
    if (k == 1)
        d(ij3)=-DNx;
        e(ij3)=0;
        g(ij3)=-DNy;
        h(ij3)=0;
    end
    % top left corner (x=0, y=Ly)
    i=1;j=ny;
    ij1=feval(@ind,ny,i-1,j);
    ij2=feval(@ind,ny,i,j-1);
    ij3=feval(@ind,ny,i,j);
    ij4=feval(@ind,ny,i,j+1);
    ij5=feval(@ind,ny,i+1,j);
    b(ij3)=DNy*c(ij2)+...
           (R-DNx-DNy-lambda*0.5)*c(ij3)+... 
            DNx*c(ij5);
    if (k == 1)
        d(ij3)=0;
        e(ij3)=-DNy;
        g(ij3)=0;
        h(ij3)=-DNx;
    end
    % top right corner (x=Lx, y=Ly)
    i=nx;j=ny;
    ij1=feval(@ind,ny,i-1,j);
    ij2=feval(@ind,ny,i,j-1);
    ij3=feval(@ind,ny,i,j);
    ij4=feval(@ind,ny,i,j+1);
    ij5=feval(@ind,ny,i+1,j);
    b(ij3)=DNx*c(ij1)+...
           DNy*c(ij2)+...
           (R-DNx-DNy-lambda*0.5)*c(ij3);
    if (k == 1)
        d(ij3)=-DNx;
        e(ij3)=-DNy;
        g(ij3)=0;
        h(ij3)=0;
        A = my_spdiags([d e f g h],[-ny,-1,0,1,ny]);
%         spy(A);
    end
    % solve matrix
    c = A\b;
    if (k == isol(m))
        csol(:,m) = c;
        m=m+1;
    end
end
return

function [ij] = ind(ny,i,j)
ij=(i-1)*ny+j;
return

function [c] = initialize(x1,y1,x2,y2,dx,dy,nx,ny,c0)
% initial conditions
n=nx*ny;
c=zeros(n,1);
continuous=1;
if (continuous == 1)
 	ix1=2+floor(x1/dx);
	ix2=max(ceil(x2/dx),1);
	iy1=2+floor(y1/dy);
	iy2=max(ceil(y2/dy),1);
	% initialize interior source nodes
	for i=ix1:ix2
        for j=iy1:iy2
            c(feval(@ind,ny,i,j))=c0;
        end
	end
	%initialize left and right edges of source nodes
	cx1=c0*((ix1-1)*dx-x1)/dx;
	cx2=c0*(x2-(ix2-1)*dx)/dx;
	for j=iy1:iy2
        i=ix1-1;
        c(feval(@ind,ny,i,j))=cx1;
        i=ix2+1;
        c(feval(@ind,ny,i,j))=cx2;
	end
	cy1=c0*((iy1-1)*dy-y1)/dy;
	cy2=c0*(y2-(iy2-1)*dy)/dy;
	% initialize top and bottom edges of source nodes
	for i=ix1:ix2
        j=iy1-1;
        c(feval(@ind,ny,i,j))=cy1;
        j=iy2+1;
        c(feval(@ind,ny,i,j))=cy2;
	end
	% initialize corner nodes
	i=ix1-1; j=iy1-1;
	c(feval(@ind,ny,i,j))=(cx1+cy1)/2;
	i=ix1-1; j=iy2+1;
	c(feval(@ind,ny,i,j))=(cx1+cy2)/2;
	i=ix2+1; j=iy1-1;
	c(feval(@ind,ny,i,j))=(cx2+cy1)/2;
	i=ix2+1; j=iy2+1;
	c(feval(@ind,ny,i,j))=(cx2+cy2)/2;
else
    ix1=1+round(x1/dx);
	ix2=1+round(x2/dx);
	iy1=1+round(y1/dy);
	iy2=1+round(y2/dy);
	for i=ix1:ix2
        for j=iy1:iy2
            c(feval(@ind,ny,i,j))=c0;
        end
	end
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