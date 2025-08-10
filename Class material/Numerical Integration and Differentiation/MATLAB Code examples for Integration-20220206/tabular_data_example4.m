% reactor problem C&C 21.20 (p. 631)
function tabular_data_example4
T=[0 1 5.5 10 12 14 16 18 20 24];
C=[1 1.5 2.3 2.1 4 5 5.5 5 3 1.2];
I=tab_data_trap(T,C)
f1=@(t)qc(t,T,C);
f2=@Q;
a=0; b=24; n=20;
I1=simpson(f1,a,b,n);
I2=simpson(f2,a,b,n);
I=I1/I2
end

function [I] = tab_data_trap(t,c)
n=length(t);
I1=0; I2=0;
for i=1:n-1
    dt=t(i+1)-t(i);
    I1=I1+dt*(c(i)*Q(t(i))+c(i+1)*Q(t(i+1)))/2;
    I2=I2+dt*(Q(t(i))+Q(t(i+1)))/2;
end
I=I1/I2;
end

function [f] = Q(t)
f = 20 + 10*sin(pi*(t-10)/12);
end

function [f] = qc(t,T,C)
q=Q(t);
c=interp1(T,C,t);
f=q*c;
end


