function tabex2
T=[0 10 20 30 35 40 45 50];
Q=[4 4.8 5.2 5.0 4.6 4.3 4.3 5.0];
c=[10 35 55 52 40 37 32 34];
Qc = Q.*c;
I_trap=trap_data(T,Qc)
f=@(t)qc_interp(t,T,Qc);
I_interp=integral(f,0,50)
end

function [qc] = qc_interp(t,T,Qc)
qc=interp1(T,Qc,t);
end