% chapra and canale p. 223 (8.44)
function pipe_example
x0=1e-3*ones(6,1);
f=@nl_func;
x=fsolve(f,x0);
Q1=x(1)
Q3=x(2)
z2=x(3) % elevation of tank b

% newton raphson needs a good guess to work for this problem
x0=[0.1; 0.1; 100; 0.01; 0.01; 0.01]; 
x_nr=newton_raphson(f,x0);
Q1=x(1)
Q3=x(2)
z2=x(3)

return

function [fx] = nl_func(x)
fx=zeros(6,1); zv=zeros(3,1);
Q=zv; D=zv; f=zv; z=zv; L=zv; v=zv;
g=9.82; eps=0.0012; nu=1e-6;
Q(1)=x(1); Q(2)= 0.1; Q(3)=x(2);
z(1)=200; z(2)=x(3); z(3)=172.5;
f(1)=x(4); f(2)=x(5); f(3)=x(6);
D(1)=0.4; D(2)=0.25; D(3)=0.2;
L(1)=1800; L(2)=500; L(3)=1400;
v=Q./(pi*D.^2/4);
hL=f.*(L./D).*(v.^2/(2*g));
Re=v.*D/nu;

% first three equations that we wrote
fx(1)=Q(1)-Q(2)-Q(3);
fx(2)=z(1)-z(2)-hL(1)-hL(2);
fx(3)=z(1)-z(3)-hL(1)-hL(3);

fx(4:6)=(1./sqrt(f))+2*log10((eps./(3.7*D))+(2.51./(Re.*sqrt(f))));% nonlinear equation
return