function tabular_data_example2
T=[0 10 20 30 35 40 45 50];
Q=[4 4.8 5.2 5.0 4.6 4.3 4.3 5.0];
c=[10 35 55 52 40 37 32 34];
I_trap=trap_data(T,Q,c)
f=@(t)qc_interp(t,T,Q,c);
I_interp=integral(f,0,50)
end

function [I] = trap_data(T,Q,c)
n=length(T);
I=0;
for i=1:n-1
    h=T(i+1)-T(i);
    qt=(Q(i+1)+Q(i))/2;
    ct=(c(i+1)+c(i))/2;
    I=I+h*qt*ct;
end
end

function [qc] = qc_interp(t,T,Q,c)
q=interp1(T,Q,t);
c=interp1(T,c,t);
qc=q.*c;
end