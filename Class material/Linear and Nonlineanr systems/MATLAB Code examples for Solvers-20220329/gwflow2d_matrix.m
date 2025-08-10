function [A,b] = gwflow2d_matrix
global d e f g h b nx ny bctype_face bcval_face
global dx dy kx ky q xw yw nw;
Lx=20; Ly=12; nx=11; ny=7;
kx=1; ky=1; % hydraulic conductivity
nw=2; % number of wells
bctype_face=zeros(4,1);
bcval_face=zeros(4,1);
% specify boundary conditions
% left
bctype_face(1)=1;
bcval_face(1)=1;
% bottom
bctype_face(2)=2;
bcval_face(2)=0;
% right
bctype_face(3)=1;
bcval_face(3)=0;
% top
bctype_face(4)=2;
bcval_face(4)=0;
% wells
q=zeros(nw,1);
xw=zeros(nw,1);
yw=zeros(nw,1);
q(1)=-1; q(2)=-1;
xw(1)=2; xw(2)=8;
yw(1)=2; yw(2)=4;
dx=Lx/(nx-1);
dy=Ly/(ny-1);
n = nx*ny;
b = zeros(n,1);
u = ones(n,1);
% set default values for matrix diagonals
d = kx*u*dy^2;
e = ky*u*dx^2;
f = -2*kx*u*dy^2-2*ky*u*dx^2;
g = kx*u*dy^2;
h = ky*u*dx^2;
% adjustments for boundary conditions
for i=1:nx
    for j=1:ny
        ij=feval(@ind,ny,i,j);
        % assign well fluxes to rhs
        wells(i,j,ij);
        dirichlet(i,j,ij);
        neumann(i,j,ij);
    end
end
% solve matrix
A = my_spdiags([d e f g h],[-ny -1 0 1 ny]);
return
%%%%%%%

function [ij] = ind(ny,i,j)
ij=(i-1)*ny+j;
return

function wells(i,j,ij)
global d e f g h b nx ny bctype_face bcval_face
global dx dy kx ky q xw yw nw;
for k=1:nw
    if (i == xw(k) & j == yw(k))
        % rhs
        b(ij)=-q(k)*dx*dy;
    end
end
return

function neumann(i,j,ij)
global d e f g h b nx ny bctype_face bcval_face
global dx dy kx ky q xw yw nw;
% adjust for neumann boundary conditions        
if (i == 1 & bctype_face(1) == 2)
    b(ij)=b(ij)+2*kx*bcval_face(1)*dx*dy^2;
    f(ij)=f(ij)+kx*dy^2;
    d(ij)=0;
elseif (j == 1 & bctype_face(2) == 2)
    b(ij)=b(ij)+2*ky*bcval_face(2)*dy*dx^2;
    f(ij)=f(ij)+ky*dx^2;
    e(ij)=0;
elseif (i == nx & bctype_face(3) == 2)
    b(ij)=b(ij)+2*kx*bcval_face(1)*dx*dy^2;
    f(ij)=f(ij)+kx*dy^2;
    h(ij)=0;
elseif (j == ny & bctype_face(4) == 2)
    b(ij)=b(ij)+2*ky*bcval_face(2)*dy*dx^2;
    f(ij)=f(ij)+ky*dx^2;
    g(ij)=0;
end
return

function dirichlet(i,j,ij)
global d e f g h b nx ny bctype_face bcval_face
global dx dy kx ky q xw yw nw;
% adjust for dirichlet boundary conditions        
if (i == 2 & bctype_face(1) == 1)
    b(ij)= b(ij)-dy^2*kx*bcval_face(1);
    d(ij)=0;
end
if (j == 2 & bctype_face(2) == 1)
    b(ij)= b(ij)-dx^2*ky*bcval_face(2);
    e(ij)=0;
end
if (i == nx-1 & bctype_face(3) == 1)
    b(ij)= b(ij)-dy^2*kx*bcval_face(3);
    h(ij)=0;
end
if (j == ny-1 & bctype_face(4) == 1)
    b(ij)=b(ij)-dx^2*ky*bcval_face(4);
    g(ij)=0;
end
if (i == 1 & bctype_face(1) == 1)
    b(ij)=bcval_face(1);
    d(ij)=0;
    e(ij)=0;
    f(ij)=1;
    g(ij)=0;
    h(ij)=0;
elseif (j == 1 & bctype_face(2) == 1)
    b(ij)=bcval_face(2);
    d(ij)=0;
    e(ij)=0;
    f(ij)=1;
    g(ij)=0;
    h(ij)=0;
elseif (i == nx & bctype_face(3) == 1)
    b(ij)=bcval_face(3);
    d(ij)=0;
    e(ij)=0;
    f(ij)=1;
    g(ij)=0;
    h(ij)=0;
elseif (j == ny & bctype_face(4) == 1)
    b(ij)=bcval_face(4);
    d(ij)=0;
    e(ij)=0;
    f(ij)=1;
    g(ij)=0;
    h(ij)=0;
end
return

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
