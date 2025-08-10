function tabex1
T=[0 1 2 3.25 4.5 6 7 8 8.5 9 10];
V=[0 5 6 5.5 7 8.5 8 6 7 7 5];
f=@(t)vel(t,T,V);
a=0; b=10; n=10;
I=trap(f,a,b,n)
I=trap_data(T,V)
end

function [vt] = vel(t,T,V)
vt = interp1(T,V,t);
end