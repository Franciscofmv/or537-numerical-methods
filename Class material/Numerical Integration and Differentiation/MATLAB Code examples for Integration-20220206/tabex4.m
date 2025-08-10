% reactor problem C&C 21.20 (p. 631)
function tabex4
T=[0 1 5.5 10 12 14 16 18 20 24];
C=[1 1.5 2.3 2.1 4 5 5.5 5 3 1.2];
Q = 20 + 10*sin(pi*(T-10)/12);
QC = Q.*C;
I1=trap_data(T,QC);
I2=trap_data(T,Q);
I_trap=I1/I2 % average concentration
f1=@(t)qc(t,T,C); % mass
f2=@Q; % volume
a=0; b=24; n=20;
I1=integral(f1,a,b);
I2=integral(f2,a,b);
I_matlab=I1/I2
I1=simpson(f1,a,b,n);
I2=simpson(f2,a,b,n);
I_simp=I1/I2
end

function [f] = Q(t)
f = 20 + 10*sin(pi*(t-10)/12);
end

function [f] = qc(t,T,C)
q=Q(t); 
c=interp1(T,C,t);% to find small c
f=q.*c;
end

function [I] = trap_data(x,y)
I = 0; n = length(x);
for i=1:n-1
    h=x(i+1)-x(i);
    I=I+h*(y(i)+y(i+1))/2;
end
end
