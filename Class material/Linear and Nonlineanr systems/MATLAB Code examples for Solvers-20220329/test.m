function test
x0 = [300 300 3000 5000];
f=@homew;
fsolve(f,x0)

end



function [fx] = homew(x)
fx=zeros(4,1);
Tc=x(1);Th=x(2);Jc=x(3);Jh=x(4);
fx(1)=5.67*10^(-8)*Tc^4+17.41*Tc-Jc-5188.8;
fx(2)=Jc-0.71*Jh+7.46*Tc-2357.71;
fx(3)=5.67*10^(-8)*Th^4+1.865*Th-Th-2250;
fx(4)=Jh-0.71*Jc+7.46*Th-11093; 
end